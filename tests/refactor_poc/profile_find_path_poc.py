#!/usr/bin/env python3
"""cProfile harness for refactor_poc find_path (poc only).

Reuses substrate/target pairs and PhaseOne wiring from bench_find_path_h2h.
Writes text dump + pstats under artifacts/.

  uv run python tests/refactor_poc/profile_find_path_poc.py
"""

from __future__ import annotations

import cProfile
import io
import pstats
import sys
import time
from pathlib import Path

from xenosite.refactor_poc.find_path import PathCounters, find_path
from xenosite.refactor_poc.rulesets import PhaseOne

# Subset of bench_find_path_h2h.CASES: one cheap, one mid, two hard hits.
CASES: list[tuple[str, str, str]] = [
    ("anisole→phenol", "COc1ccccc1", "Oc1ccccc1"),
    ("PhCH2OH→quinone", "OCc1ccccc1", "O=C1C=CC(=O)C=C1"),
    ("acetate→catechol", "CC(=O)Oc1ccc(OC)cc1", "Oc1ccc(O)cc1"),
    (
        "MeOPhOH→hydroxyquinone",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
]

POC_MAX_NODES = 800
MAX_PATHS = 1
TOP_N = 30

ROOT = Path(__file__).resolve().parents[2]
# Post-xf / no-canonicalize / lazy-csmi snapshot (prior baseline: poc_find_path_profile.out).
ART_OUT = ROOT / "artifacts" / "poc_find_path_profile_after_xf.out"
ART_PSTATS = ROOT / "artifacts" / "poc_find_path_profile_after_xf.pstats"
ART_LIVE = ROOT / "artifacts" / "poc_find_path_profile_after_xf.live.log"


def run_case(label: str, reactant: str, target: str) -> tuple[bool, float, PathCounters]:
    counters = PathCounters()
    t0 = time.perf_counter()
    hits = list(
        find_path(
            reactant,
            target,
            ruleset=PhaseOne,
            counters=counters,
            max_paths=MAX_PATHS,
            max_nodes=POC_MAX_NODES,
        )
    )
    elapsed = time.perf_counter() - t0
    return bool(hits), elapsed, counters


def dump_stats(pr: cProfile.Profile, stream: io.StringIO) -> None:
    ps = pstats.Stats(pr, stream=stream)
    ps.strip_dirs()
    stream.write("\n=== sort by cumulative (top %d) ===\n" % TOP_N)
    ps.sort_stats("cumulative")
    ps.print_stats(TOP_N)
    stream.write("\n=== sort by tottime (top %d) ===\n" % TOP_N)
    ps.sort_stats("tottime")
    ps.print_stats(TOP_N)
    stream.write("\n=== callers for top cumulative (top 15) ===\n")
    ps.sort_stats("cumulative")
    ps.print_callers(15)


def deepcopy_verdict(pr: cProfile.Profile) -> str:
    """Summarize copy.deepcopy share vs prior ~59% / 1.17s cum baseline."""

    buf = io.StringIO()
    ps = pstats.Stats(pr, stream=buf)
    ps.strip_dirs()
    wall = float(ps.total_tt)
    deep = 0.0
    deep_ncalls = 0
    canon = 0.0
    for (filename, _line, func), (_cc, nc, _tt, ct, _callers) in ps.stats.items():
        if filename.endswith("copy.py") and func == "deepcopy":
            deep = ct
            deep_ncalls = nc
        if func == "cannonicalize_order":
            canon = ct
    pct = (100.0 * deep / wall) if wall > 0 else 0.0
    return (
        f"deepcopy verdict: cum={deep:.3f}s  ncalls={deep_ncalls}  "
        f"~{pct:.1f}% of profile wall={wall:.3f}s  "
        f"cannonicalize_order cum={canon:.3f}s  "
        f"(prior baseline ~1.17s / ~59% via canonicalize+trace)"
    )


def main() -> int:
    ART_OUT.parent.mkdir(parents=True, exist_ok=True)
    live = ART_LIVE.open("w")
    lines: list[str] = []

    def log(msg: str) -> None:
        print(msg, flush=True)
        live.write(msg + "\n")
        live.flush()
        lines.append(msg)

    log("poc find_path cProfile (after xf / no-canonicalize / lazy csmi)")
    log(f"cases={len(CASES)}  max_nodes={POC_MAX_NODES}  max_paths={MAX_PATHS}")
    log(f"live log: {ART_LIVE}")
    log("")

    # Warm import / RDKit once outside the profile window for cleaner signal.
    log("[warmup] anisole→phenol ...")
    hit, sec, c = run_case(*CASES[0])
    log(f"  hit={hit}  {sec:.3f}s  ed={c.mol_edits} re={c.rule_expansions} nd={c.nodes}")

    pr = cProfile.Profile()
    log("")
    log("[profile] enabling cProfile")
    pr.enable()
    case_rows: list[str] = []
    for i, (label, reactant, target) in enumerate(CASES, 1):
        log(f"[{i}/{len(CASES)}] {label} ...")
        hit, sec, c = run_case(label, reactant, target)
        row = (
            f"  hit={hit}  {sec:.3f}s  "
            f"ed={c.mol_edits} re={c.rule_expansions} nd={c.nodes} "
            f"sites={c.sites_considered}/{c.sites_skipped}skip"
        )
        log(row)
        case_rows.append(
            f"{label}: hit={hit} {sec:.3f}s ed={c.mol_edits} "
            f"re={c.rule_expansions} nd={c.nodes}"
        )
    pr.disable()
    log("[profile] disabled")

    verdict = deepcopy_verdict(pr)
    log("")
    log(verdict)

    buf = io.StringIO()
    dump_stats(pr, buf)
    stats_text = buf.getvalue()
    print(stats_text, flush=True)
    live.write(stats_text)
    live.flush()

    ART_PSTATS.parent.mkdir(parents=True, exist_ok=True)
    pr.dump_stats(str(ART_PSTATS))

    body = "\n".join(lines) + "\n" + stats_text
    body += f"\npstats written to {ART_PSTATS}\n"
    ART_OUT.write_text(body)
    log(f"wrote {ART_OUT}")
    log(f"wrote {ART_PSTATS}")
    live.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
