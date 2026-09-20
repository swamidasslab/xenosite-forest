#!/usr/bin/env python3
"""Three-way path-search H2H: archive BFS, archive DFS, live find_path.

Archive (read-only ``xenosite._archive_forest``):
  classic ``bfs`` / ``dfs`` metabolite enumeration (PhaseOneQF), single-reactant
  mode. Cost bound only at the harness: stop after ``MAX_MOLS`` yielded
  metabolites (``max_paths=MAX_MOLS`` + consumer break). DFS passes
  ``all_paths=True`` so an empty endpoint set does not trip the archive's
  early-return quirk. No archive source edits.

Live (``xenosite.forest``):
  ``find_path`` with PhaseOne and its native ``max_nodes`` budget.

Cases emphasize mid-size / multi-edit targets where ordering blow-up makes
BFS/DFS miss under ``MAX_MOLS`` while live find_path returns a long plan fast.

Re-run:
  uv run python tests/forest/bench_find_path_h2h.py
"""

from __future__ import annotations

import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from xenosite._archive_forest.bfs import bfs as archive_bfs
from xenosite._archive_forest.bfs import dfs as archive_dfs
from xenosite._archive_forest.utils import canon_smi
from xenosite.forest.find_path import PathCounters
from xenosite.forest.find_path import find_path as live_find_path
from xenosite.forest.rdkitutil import canon_smiles
from xenosite.forest.rulesets import PhaseOne

warnings.filterwarnings("ignore", category=UserWarning, module="xenosite._archive_forest")

ROOT = Path(__file__).resolve().parents[2]
ART_OUT = ROOT / "artifacts" / "bench_find_path_h2h_3way.out"
ART_LIVE = ROOT / "artifacts" / "bench_find_path_h2h_3way.live.log"

# Shared archive yield cap. Large enough that a short path can appear; small
# enough that multi-edit ordering blow-up hits the cap instead of hanging.
MAX_MOLS = 200
ARCHIVE_DEPTH = 4
LIVE_MAX_NODES = 800
MAX_PATHS_LIVE = 1

# Mid-size / multi-edit stories (not toy 2–4 atom cases).
CASES: list[tuple[str, str, str]] = [
    # eugenol → allyl-quinone (4-step plan; BFS/DFS cap)
    (
        "eugenol→allyl-quinone",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
    # 3,4-dimethoxyphenethylamine → catechol (2 unordered dealkylations)
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    # 4-methoxyphenol → hydroxyquinone (4 concurrent edits)
    (
        "MeOPhOH→hydroxyquinone",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
    # TBA (22 heavy atoms): frontier fills MAX_MOLS before the cleavage product
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        "CC(C)(C)C#CC=CC=O",
    ),
    # 2-methoxynaphthalene → 1,2-NQ (reachable 3-step plan).
    # Contrast: 2-MeO → 1,4-NQ has no PhaseOne path (see docs/forest/PERFORMANCE.md).
    (
        "2-MeO-naph→1,2-NQ",
        "COc1ccc2ccccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
    ),
]


@dataclass
class ArchiveResult:
    hit: bool
    valid: bool
    product: str | None
    seconds: float
    path_len: int | None
    mols_yielded: int
    capped: bool


@dataclass
class LiveResult:
    hit: bool
    valid: bool
    product: str | None
    seconds: float
    path_len: int | None
    plan_str: str | None
    nodes: int
    mol_edits: int
    rule_expansions: int
    billed: int
    budget: int
    budget_exhausted: bool


def _git_sha() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=ROOT,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        return "unknown"


def run_archive(search_fn, reactant: str, target_c: str) -> ArchiveResult:
    """Enumerate metabolites; stop at target hit or MAX_MOLS yields."""

    t0 = time.perf_counter()
    n = 0
    product: str | None = None
    path_len: int | None = None
    valid = False

    kwargs: dict = dict(
        ruleset="PhaseOneQF",
        depth=ARCHIVE_DEPTH,
        max_paths=MAX_MOLS,
        phase1=True,
    )
    # Archive DFS early-returns when the endpoint set is empty unless all_paths.
    if search_fn is archive_dfs:
        kwargs["all_paths"] = True

    for _smi, sites, mols in search_fn([reactant], **kwargs):
        n += 1
        prod = canon_smi(mols[-1]) if mols else None
        product = prod
        if prod == target_c:
            valid = True
            path_len = len(sites) if sites is not None else None
            break
        if n >= MAX_MOLS:
            break

    capped = (not valid) and n >= MAX_MOLS
    return ArchiveResult(
        hit=valid,
        valid=valid,
        product=product if valid else None,
        seconds=time.perf_counter() - t0,
        path_len=path_len,
        mols_yielded=n,
        capped=capped,
    )


def run_live(reactant: str, target: str, target_c: str) -> LiveResult:
    counters = PathCounters()
    t0 = time.perf_counter()
    hits = list(
        live_find_path(
            reactant,
            target,
            ruleset=PhaseOne,
            counters=counters,
            max_paths=MAX_PATHS_LIVE,
            max_nodes=LIVE_MAX_NODES,
        )
    )
    elapsed = time.perf_counter() - t0
    product = None
    path_len = None
    plan_str = None
    if hits:
        product = canon_smiles(hits[0].smiles)
        plan = hits[0].plan
        plan_str = str(plan)
        steps = getattr(plan, "steps", None)
        path_len = len(steps) if steps is not None else None
    valid = product == target_c if product is not None else False
    nodes = int(counters.nodes)
    exhausted = (not valid) and nodes >= LIVE_MAX_NODES
    return LiveResult(
        hit=bool(hits),
        valid=valid,
        product=product,
        seconds=elapsed,
        path_len=path_len if valid else None,
        plan_str=plan_str if valid else None,
        nodes=nodes,
        mol_edits=int(counters.mol_edits),
        rule_expansions=int(counters.rule_expansions),
        billed=int(counters.billed),
        budget=LIVE_MAX_NODES,
        budget_exhausted=exhausted,
    )


def _fmt_archive(label: str, r: ArchiveResult) -> str:
    if r.valid:
        status = f"ok hops={r.path_len}"
    elif r.capped:
        status = "CAP"
    else:
        status = "miss"
    return f"{label} {r.seconds:.3f}s  {status}  mols={r.mols_yielded}/{MAX_MOLS}"


def _fmt_live(r: LiveResult) -> str:
    if r.valid:
        status = f"ok steps={r.path_len}"
    elif r.budget_exhausted:
        status = "EXH"
    elif r.hit:
        status = "INVALID"
    else:
        status = "miss"
    return (
        f"live {r.seconds:.3f}s  {status}  "
        f"nodes={r.nodes}/{r.budget}  ed={r.mol_edits}  "
        f"re={r.rule_expansions}  bill={r.billed}"
    )


def main() -> int:
    ART_OUT.parent.mkdir(parents=True, exist_ok=True)
    live_f = ART_LIVE.open("w")
    lines: list[str] = []
    sha = _git_sha()

    def log(msg: str = "") -> None:
        print(msg, flush=True)
        live_f.write(msg + "\n")
        live_f.flush()
        lines.append(msg)

    log(f"find_path H2H 3-way  sha={sha}")
    log(
        f"archive: PhaseOneQF bfs/dfs enum  depth={ARCHIVE_DEPTH}  "
        f"MAX_MOLS={MAX_MOLS} (harness yield cap)"
    )
    log(
        f"live:    PhaseOne find_path  max_nodes={LIVE_MAX_NODES}  "
        f"max_paths={MAX_PATHS_LIVE}"
    )
    log("archive path length = reaction hops on the metabolite walk")
    log("live path length    = len(plan.steps) on PathOutcome")
    log(f"live log: {ART_LIVE}")
    log()

    rows: list[tuple[str, str, ArchiveResult, ArchiveResult, LiveResult]] = []
    for i, (label, reactant, target) in enumerate(CASES, 1):
        target_c = canon_smiles(Chem.MolFromSmiles(target))
        ha = Chem.MolFromSmiles(reactant).GetNumHeavyAtoms()
        log(f"[{i}/{len(CASES)}] {label}  (reactant heavy atoms={ha}) ...")
        log(f"  reactant={reactant}")
        log(f"  product ={target}")
        log("  archive bfs...")
        br = run_archive(archive_bfs, reactant, target_c)
        log(f"    {_fmt_archive('bfs', br)}")
        log("  archive dfs...")
        dr = run_archive(archive_dfs, reactant, target_c)
        log(f"    {_fmt_archive('dfs', dr)}")
        log("  live find_path...")
        lr = run_live(reactant, target, target_c)
        log(f"    {_fmt_live(lr)}")
        if lr.plan_str:
            log(f"    plan={lr.plan_str}")
        log(f"  target={target_c}")
        log()
        rows.append((label, target_c, br, dr, lr))

    log("=" * 110)
    log(
        f"{'case':26}  {'bfs_s':>7}  {'dfs_s':>7}  {'live_s':>7}  "
        f"{'bfs':^16}  {'dfs':^16}  {'live':^24}"
    )
    log("-" * 110)

    tot_b = tot_d = tot_l = 0.0
    for label, _t, br, dr, lr in rows:
        tot_b += br.seconds
        tot_d += dr.seconds
        tot_l += lr.seconds

        def ac(r: ArchiveResult) -> str:
            if r.valid:
                return f"ok hops={r.path_len} n={r.mols_yielded}"
            if r.capped:
                return f"CAP n={r.mols_yielded}"
            return "miss"

        def lc(r: LiveResult) -> str:
            if r.valid:
                return f"ok steps={r.path_len} nd={r.nodes}"
            if r.budget_exhausted:
                return f"EXH nd={r.nodes}"
            return f"miss nd={r.nodes}"

        log(
            f"{label:26}  {br.seconds:7.3f}  {dr.seconds:7.3f}  {lr.seconds:7.3f}  "
            f"{ac(br):^16}  {ac(dr):^16}  {lc(lr):^24}"
        )

    log("=" * 110)
    log()
    log(f"Totals wall: bfs={tot_b:.3f}s  dfs={tot_d:.3f}s  live={tot_l:.3f}s")
    n_b = sum(1 for _, _, br, _, _ in rows if br.valid)
    n_d = sum(1 for _, _, _, dr, _ in rows if dr.valid)
    n_l = sum(1 for _, _, _, _, lr in rows if lr.valid)
    log(f"Valid hits: bfs={n_b}/{len(rows)}  dfs={n_d}/{len(rows)}  live={n_l}/{len(rows)}")
    log()
    log("Metric definitions:")
    log(f"  MAX_MOLS={MAX_MOLS}        archive metabolite yield cap (harness)")
    log(f"  ARCHIVE_DEPTH={ARCHIVE_DEPTH}     classic hop ceiling")
    log(f"  LIVE_MAX_NODES={LIVE_MAX_NODES}   live find_path node budget")
    log("  CAP                  archive hit MAX_MOLS without seeing the target")
    log("  EXH                  live nodes reached max_nodes without a valid hit")
    log("  Conjugation is not in PhaseOneQF / PhaseOne.")

    body = "\n".join(lines) + "\n"
    ART_OUT.write_text(body)
    log(f"wrote {ART_OUT}")
    live_f.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
