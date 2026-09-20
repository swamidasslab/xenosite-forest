#!/usr/bin/env python3
"""Head-to-head find_path: archived forest (PhaseOneQF) vs live forest (PhaseOne).

Same reactants/targets. Reports work counters that exist on both sides plus
side-specific billed units, with labels that say what each number means.

Comparable columns (both sides expose these):
  mol_edits   — accepted reaction applies / kekulé overlays
  rule_exp    — frontier rule expansions (metabolize / enumerate calls)
  nodes       — poc: queue pops; forest: nodes_enqueued
  wall_s      — wall time

Side-specific (not the same unit — do not equate):
  forest billed = linearizations_applied + site_applies
    (guided PhaseOneQF often has site_applies=0; linearizations are the work)
  poc billed    = mol_edits + nodes

Rulesets (closest Phase I + QF pair; neither default includes conjugation):
  forest: PhaseOneQF  (find_path default)
  poc:    PhaseOne    (rulesets.PhaseOne; includes QuinoneFormation)

Re-run:
  uv run python tests/forest/bench_find_path_h2h.py
"""

from __future__ import annotations

import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from xenosite._archive_forest import PathSearchCounters
from xenosite._archive_forest.guided_path import find_path as forest_find_path
from xenosite._archive_forest.utils import canon_smi
from xenosite.forest.find_path import PathCounters
from xenosite.forest.find_path import find_path as poc_find_path
from xenosite.forest.rdkitutil import canon_smiles
from xenosite.forest.rulesets import PhaseOne

warnings.filterwarnings("ignore", category=UserWarning, module="xenosite._archive_forest")

ROOT = Path(__file__).resolve().parents[2]
ART_OUT = ROOT / "artifacts" / "bench_find_path_h2h_post_swap.out"
ART_LIVE = ROOT / "artifacts" / "bench_find_path_h2h_post_swap.live.log"

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
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        "CC(C)(C)C#CC=CC=O",
    ),
    ("acetate→catechol", "CC(=O)Oc1ccc(OC)cc1", "Oc1ccc(O)cc1"),
    ("PhCH2OH→quinone", "OCc1ccccc1", "O=C1C=CC(=O)C=C1"),
    # Forest flaky / high-budget under PhaseOneQF; poc should solve.
    (
        "MeOPhOH→hydroxyquinone",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
    (
        "MeOPhOH→orthocarbonate Q",
        "COc1ccc(O)cc1",
        "O=C1C=CC(OC(O)O)=CC1=O",
    ),
    ("naphthalene→1,4-NQ", "c1ccc2ccccc2c1", "O=C1C=CC(=O)c2ccccc12"),
]


@dataclass
class SideResult:
    hit: bool
    product: str | None
    seconds: float
    work: dict
    valid: bool  # hit and product canonically equals target


def _forest_product(outcome) -> str | None:
    if not outcome.smiles:
        return None
    last = outcome.smiles[-1]
    return canon_smi(last) if last else None


def _valid_product(product: str | None, target_c: str) -> bool:
    """Correctness gate: a hit counts only when the product is the target."""

    return product is not None and product == target_c


def run_forest(reactant: str, target: str, target_c: str) -> SideResult:
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
    hit = bool(hits)
    return SideResult(
        hit=hit,
        product=product,
        seconds=elapsed,
        work={
            "mol_edits": counters.mol_edits,
            "rule_expansions": counters.rule_expansions,
            "nodes": counters.nodes_enqueued,
            "linearizations": counters.linearizations_applied,
            "site_applies": counters.site_applies,
            "billed": counters.billed(),
            "budget_exhausted": counters.budget_exhausted,
        },
        valid=_valid_product(product, target_c) if hit else False,
    )


def run_poc(reactant: str, target: str, target_c: str) -> SideResult:
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
    if product is not None:
        product = canon_smiles(product)
    hit = bool(hits)
    return SideResult(
        hit=hit,
        product=product,
        seconds=elapsed,
        work={
            "mol_edits": counters.mol_edits,
            "rule_expansions": counters.rule_expansions,
            "nodes": counters.nodes,
            "sites_considered": counters.sites_considered,
            "sites_skipped": counters.sites_skipped,
            "billed": counters.billed,
        },
        valid=_valid_product(product, target_c) if hit else False,
    )


def _gate(fr: SideResult, pr: SideResult) -> str:
    """Classify correctness for speed celebration.

    both_ok       — both hit with product == target (speed comparable)
    poc_only_ok   — poc valid; forest miss/EXH/invalid (poc win, note separately)
    forest_only_ok — forest valid; poc miss/invalid (poc regression)
    invalid_hit   — a side claimed hit but product != target
    both_miss     — neither found a valid path
    """

    if fr.valid and pr.valid:
        return "both_ok"
    if pr.valid and not fr.valid:
        return "poc_only_ok"
    if fr.valid and not pr.valid:
        return "forest_only_ok"
    if (fr.hit and not fr.valid) or (pr.hit and not pr.valid):
        return "invalid_hit"
    return "both_miss"


def _fmt_side(label: str, r: SideResult) -> str:
    """Print comparable counters first; side billed after."""

    if label == "forest":
        return (
            "ed=%s re=%s nd=%s  bill(lin+sa)=%s (L=%s sa=%s)  %.3fs%s"
            % (
                r.work["mol_edits"],
                r.work["rule_expansions"],
                r.work["nodes"],
                r.work["billed"],
                r.work["linearizations"],
                r.work["site_applies"],
                r.seconds,
                " EXH" if r.work.get("budget_exhausted") else "",
            )
        )
    return "ed=%s re=%s nd=%s  bill(ed+nd)=%s  %.3fs" % (
        r.work["mol_edits"],
        r.work["rule_expansions"],
        r.work["nodes"],
        r.work["billed"],
        r.seconds,
    )


def _forest_idle(r: SideResult) -> bool:
    """True only when forest did no billed work (not merely site_applies==0)."""
    return int(r.work.get("billed", 0) or 0) == 0 and int(
        r.work.get("mol_edits", 0) or 0
    ) == 0


def main() -> int:
    ART_OUT.parent.mkdir(parents=True, exist_ok=True)
    live = ART_LIVE.open("w")
    lines: list[str] = []

    def log(msg: str = "") -> None:
        print(msg, flush=True)
        live.write(msg + "\n")
        live.flush()
        lines.append(msg)

    want = canon_smiles
    log("find_path H2H  forest=PhaseOneQF  poc=PhaseOne  (after xf)")
    log(
        "ceilings: forest max_expansions=%s  poc max_nodes=%s  max_paths=%s"
        % (FOREST_MAX_EXPANSIONS, POC_MAX_NODES, MAX_PATHS)
    )
    log(
        "comparable: mol_edits (ed) / rule_expansions (re) / nodes (nd) / wall_s"
    )
    log(
        "billed units differ: forest=lin+site_applies; poc=mol_edits+nodes"
    )
    log(
        "correctness gate: valid = hit AND product canonically equals target"
    )
    log(f"live log: {ART_LIVE}")
    log()
    rows = []
    for i, (label, reactant, target) in enumerate(CASES, 1):
        target_c = want(Chem.MolFromSmiles(target))
        log("[%s/%s] %s ..." % (i, len(CASES), label))
        log("  forest...")
        fr = run_forest(reactant, target, target_c)
        log(
            "    %s hit=%s valid=%s product=%s"
            % (_fmt_side("forest", fr), fr.hit, fr.valid, fr.product)
        )
        if _forest_idle(fr) and not fr.hit:
            log("    WARN: forest idle (billed=0, no hit)")
        log("  poc...")
        pr = run_poc(reactant, target, target_c)
        log(
            "    %s hit=%s valid=%s product=%s"
            % (_fmt_side("poc", pr), pr.hit, pr.valid, pr.product)
        )

        gate = _gate(fr, pr)
        match = fr.valid and pr.valid and fr.product == pr.product == target_c
        rows.append((label, reactant, target_c, fr, pr, match, gate))
        log("  gate=%s  match=%s  target=%s" % (gate, match, target_c))
        log()

    log("=" * 140)
    log(
        "%-24s  %-12s  %8s  %8s  %-55s  %-40s"
        % ("case", "gate", "forest_s", "poc_s", "forest ed/re/nd", "poc ed/re/nd")
    )
    log("-" * 140)
    forest_total = poc_total = 0.0
    both_ok_f = both_ok_p = 0.0
    n_both_ok = n_poc_only = n_forest_only = n_invalid = n_miss = 0
    for label, _r, _target_c, fr, pr, match, gate in rows:
        forest_total += fr.seconds
        poc_total += pr.seconds
        if gate == "both_ok":
            n_both_ok += 1
            both_ok_f += fr.seconds
            both_ok_p += pr.seconds
        elif gate == "poc_only_ok":
            n_poc_only += 1
        elif gate == "forest_only_ok":
            n_forest_only += 1
        elif gate == "invalid_hit":
            n_invalid += 1
        else:
            n_miss += 1
        fcell = "ed=%s re=%s nd=%s b=%s(L=%s sa=%s)%s%s" % (
            fr.work["mol_edits"],
            fr.work["rule_expansions"],
            fr.work["nodes"],
            fr.work["billed"],
            fr.work["linearizations"],
            fr.work["site_applies"],
            " EXH" if fr.work.get("budget_exhausted") else "",
            " miss" if not fr.hit else ("" if fr.valid else " INVALID"),
        )
        pcell = "ed=%s re=%s nd=%s b=%s%s" % (
            pr.work["mol_edits"],
            pr.work["rule_expansions"],
            pr.work["nodes"],
            pr.work["billed"],
            " miss" if not pr.hit else ("" if pr.valid else " INVALID"),
        )
        log(
            "%-24s  %-12s  %8.3f  %8.3f  %-55s  %-40s"
            % (label, gate, fr.seconds, pr.seconds, fcell, pcell)
        )
    log("=" * 140)
    log()
    log("Totals: forest=%.3fs  poc=%.3fs  (all cases)" % (forest_total, poc_total))
    if n_both_ok:
        speedup = both_ok_f / both_ok_p if both_ok_p > 0 else float("inf")
        log(
            "Speed (both_ok only, n=%s): forest=%.3fs  poc=%.3fs  ratio=%.2fx"
            % (n_both_ok, both_ok_f, both_ok_p, speedup)
        )
    log(
        "Gates: both_ok=%s  poc_only_ok=%s  forest_only_ok=%s  "
        "invalid_hit=%s  both_miss=%s"
        % (n_both_ok, n_poc_only, n_forest_only, n_invalid, n_miss)
    )
    log()
    log("Metric definitions:")
    log("  mol_edits (ed)     both: accepted reaction apply / overlay")
    log("  rule_expansions    both: frontier rule metabolize/enumerate")
    log("  nodes (nd)         poc queue pops; forest nodes_enqueued")
    log("  forest billed      linearizations_applied + site_applies")
    log("                     (sa often 0 on guided; lin is the real work)")
    log("  poc billed         mol_edits + nodes")
    log("  Conjugation is not in PhaseOneQF / PhaseOne.")
    log(
        "  Speed celebration only on both_ok; poc_only_ok noted separately "
        "(forest EXH/miss)."
    )
    idle = [label for label, _r, _t, fr, _pr, _m, _g in rows if _forest_idle(fr)]
    if idle:
        log("  IDLE forest cases (invalid for H2H): %s" % ", ".join(idle))
    body = "\n".join(lines) + "\n"
    ART_OUT.write_text(body)
    log(f"wrote {ART_OUT}")
    live.close()
    if n_forest_only or n_invalid:
        return 1
    if idle:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
