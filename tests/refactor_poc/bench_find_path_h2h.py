#!/usr/bin/env python3
"""Head-to-head find_path: forest (PhaseOneQF) vs refactor_poc (PhaseOne).

Same reactants/targets. Reports work counters and wall time. Does not paste
guided_path rankers; each side uses its own search defaults under a shared
budget ceiling large enough that both sides can finish when a path exists.

Rulesets (closest Phase I + QF pair; neither default includes conjugation):
  forest: PhaseOneQF  (find_path default)
  poc:    PhaseOne    (rulesets.PhaseOne; includes QuinoneFormation)

Re-run:
  uv run python tests/refactor_poc/bench_find_path_h2h.py
"""

from __future__ import annotations

import sys
import time
import warnings
from dataclasses import dataclass

from rdkit import Chem

from xenosite.forest import PathSearchCounters
from xenosite.forest import find_path as forest_find_path
from xenosite.forest.utils import canon_smi
from xenosite.refactor_poc.find_path import PathCounters, find_path as poc_find_path
from xenosite.refactor_poc.rdkitutil import canon_smiles
from xenosite.refactor_poc.rulesets import PhaseOne

warnings.filterwarnings("ignore", category=UserWarning, module="xenosite.forest")

# Shared ceilings: high enough for hard cases; not a contest target.
FOREST_MAX_EXPANSIONS = 200
POC_MAX_NODES = 800
MAX_PATHS = 1

CASES: list[tuple[str, str, str]] = [
    # label, reactant, target
    ("anisole→phenol", "COc1ccccc1", "Oc1ccccc1"),
    ("benzene→quinone", "c1ccccc1", "O=C1C=CC(=O)C=C1"),
    ("phenol→quinone", "Oc1ccccc1", "O=C1C=CC(=O)C=C1"),
    ("butylbenzene→ω-OH", "c1ccc(CCCC)cc1", "OCCCCc1ccccc1"),
    # Prior hard: terbinafine → enyne aldehyde (many expansions / cleavages)
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        "CC(C)(C)C#CC=CC=O",
    ),
    # Prior hard: acetate + O-dealkylation to catechol
    ("acetate→catechol", "CC(=O)Oc1ccc(OC)cc1", "Oc1ccc(O)cc1"),
    # Multi-hop oxygenation + dearomatization under PhaseOne
    ("PhCH2OH→quinone", "OCc1ccccc1", "O=C1C=CC(=O)C=C1"),
]


@dataclass
class SideResult:
    hit: bool
    product: str | None
    seconds: float
    work: dict


def _forest_product(outcome) -> str | None:
    if not outcome.smiles:
        return None
    last = outcome.smiles[-1]
    return canon_smi(last) if last else None


def run_forest(reactant: str, target: str) -> SideResult:
    counters = PathSearchCounters()
    t0 = time.perf_counter()
    hits = list(
        forest_find_path(
            reactant,
            target,
            ruleset="PhaseOneQF",
            max_paths=MAX_PATHS,
            max_expansions=FOREST_MAX_EXPANSIONS,
            counters=counters,
        )
    )
    elapsed = time.perf_counter() - t0
    product = _forest_product(hits[0]) if hits else None
    return SideResult(
        hit=bool(hits),
        product=product,
        seconds=elapsed,
        work={
            "site_applies": counters.site_applies,
            "billed": counters.billed(),
            "linearizations": counters.linearizations_applied,
            "rule_expansions": counters.rule_expansions,
            "mol_edits": counters.mol_edits,
            "budget_exhausted": counters.budget_exhausted,
        },
    )


def run_poc(reactant: str, target: str) -> SideResult:
    counters = PathCounters()
    t0 = time.perf_counter()
    hits = list(
        poc_find_path(
            reactant,
            target,
            ruleset=PhaseOne,
            counters=counters,
            max_paths=MAX_PATHS,
            max_nodes=POC_MAX_NODES,
        )
    )
    elapsed = time.perf_counter() - t0
    product = hits[0].smiles if hits else None
    return SideResult(
        hit=bool(hits),
        product=product,
        seconds=elapsed,
        work={
            "mol_edits": counters.mol_edits,
            "billed": counters.billed,
            "nodes": counters.nodes,
            "sites_considered": counters.sites_considered,
            "sites_skipped": counters.sites_skipped,
            "rule_expansions": counters.rule_expansions,
        },
    )


def _fmt_forest(r: SideResult) -> str:
    # Guided PhaseOneQF bills plan linearizations; site_applies often stays 0.
    return "bill=%s lin=%s ed=%s sa=%s  %.3fs%s" % (
        r.work["billed"],
        r.work["linearizations"],
        r.work["mol_edits"],
        r.work["site_applies"],
        r.seconds,
        " EXH" if r.work.get("budget_exhausted") else "",
    )


def _fmt_poc(r: SideResult) -> str:
    return "ed=%s bill=%s nd=%s  %.3fs" % (
        r.work["mol_edits"],
        r.work["billed"],
        r.work["nodes"],
        r.seconds,
    )


def _forest_idle(r: SideResult) -> bool:
    """True only when forest did no billed work (not merely site_applies==0)."""
    return int(r.work.get("billed", 0) or 0) == 0 and int(
        r.work.get("mol_edits", 0) or 0
    ) == 0


def main() -> int:
    want = canon_smiles
    print("find_path H2H  forest=PhaseOneQF  poc=PhaseOne")
    print(
        "ceilings: forest max_expansions=%s  poc max_nodes=%s  max_paths=%s"
        % (FOREST_MAX_EXPANSIONS, POC_MAX_NODES, MAX_PATHS)
    )
    print(flush=True)
    rows = []
    for i, (label, reactant, target) in enumerate(CASES, 1):
        target_c = want(Chem.MolFromSmiles(target))
        print("[%s/%s] %s ..." % (i, len(CASES), label), flush=True)
        print("  forest...", flush=True)
        fr = run_forest(reactant, target)
        print(
            "    %s hit=%s product=%s" % (_fmt_forest(fr), fr.hit, fr.product),
            flush=True,
        )
        if _forest_idle(fr) and not fr.hit:
            print("    WARN: forest idle (billed=0, no hit)", flush=True)
        print("  poc...", flush=True)
        pr = run_poc(reactant, target)
        print(
            "    %s hit=%s product=%s" % (_fmt_poc(pr), pr.hit, pr.product),
            flush=True,
        )

        match = False
        if fr.hit and pr.hit and fr.product and pr.product:
            match = fr.product == pr.product == target_c
        rows.append((label, reactant, target_c, fr, pr, match))
        print("  match=%s  target=%s" % (match, target_c), flush=True)
        print(flush=True)

    print("=" * 110)
    print(
        "%-22s  %-36s  %-28s  %s"
        % ("case", "forest (bill/lin/ed / s)", "poc (ed/bill / s)", "match")
    )
    print("-" * 110)
    for label, _r, _target_c, fr, pr, match in rows:
        fcell = "b=%s L=%s ed=%s %.2fs" % (
            fr.work["billed"],
            fr.work["linearizations"],
            fr.work["mol_edits"],
            fr.seconds,
        )
        if fr.work.get("budget_exhausted"):
            fcell += " EXH"
        if not fr.hit:
            fcell += " miss"
        pcell = "ed=%s b=%s %.2fs" % (
            pr.work["mol_edits"],
            pr.work["billed"],
            pr.seconds,
        )
        if not pr.hit:
            pcell += " miss"
        print("%-22s  %-36s  %-28s  %s" % (label, fcell, pcell, match))
    print("=" * 110)
    print()
    print("Notes:")
    print(
        "  forest: billed=linearizations+site_applies; guided PhaseOneQF often"
        " has site_applies=0 with linearizations>0 (not idle)."
    )
    print("  poc:    mol_edits; billed=mol_edits+nodes")
    print("  Units differ; compare within-side and wall time across sides.")
    print("  Conjugation is not in PhaseOneQF / PhaseOne (closest Phase I+QF pair).")
    n_match = sum(1 for *_, m in rows if m)
    n_both = sum(1 for *_, fr, pr, _m in rows if fr.hit and pr.hit)
    print(
        "  both hit: %s/%s   products match: %s/%s"
        % (n_both, len(rows), n_match, len(rows))
    )
    idle = [label for label, _r, _t, fr, _pr, _m in rows if _forest_idle(fr)]
    if idle:
        print("  IDLE forest cases (invalid for H2H): %s" % ", ".join(idle))
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
