#!/usr/bin/env python3
"""cProfile off vs on for opt-in canonical lex-orbit site emission.

Same find_path case shape as ``profile_forest_copy.py``.

  uv run python tests/forest/profile_canonical_emitted_sites.py

Writes ``artifacts/canonical_emitted_sites_profile.{out,pstats,live.log}``.
"""

from __future__ import annotations

import cProfile
import io
import pstats
import time
from pathlib import Path

from xenosite.forest.find_path import PathCounters, find_path
from xenosite.forest.rulesets import PhaseOne

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

MAX_NODES = 800
MAX_PATHS = 1
TOP_N = 30
ROOT = Path(__file__).resolve().parents[2]
ART_OUT = ROOT / "artifacts" / "canonical_emitted_sites_profile.out"
ART_PSTATS = ROOT / "artifacts" / "canonical_emitted_sites_profile.pstats"
ART_LIVE = ROOT / "artifacts" / "canonical_emitted_sites_profile.live.log"


def _log(msg: str) -> None:
    line = msg if msg.endswith("\n") else msg + "\n"
    print(msg, flush=True)
    ART_LIVE.parent.mkdir(parents=True, exist_ok=True)
    with ART_LIVE.open("a", encoding="utf-8") as fh:
        fh.write(line)


def run_batch(*, canonical: bool) -> tuple[float, list[tuple[str, bool, float, PathCounters]]]:
    """Run cases with explicit ``canonical_emitted_sites`` on ``find_path``."""

    rows: list[tuple[str, bool, float, PathCounters]] = []
    t0 = time.perf_counter()
    for label, reactant, target in CASES:
        counters = PathCounters()
        c0 = time.perf_counter()
        hits = list(
            find_path(
                reactant,
                target,
                ruleset=PhaseOne,
                counters=counters,
                max_paths=MAX_PATHS,
                max_nodes=MAX_NODES,
                canonical_emitted_sites=canonical,
            )
        )
        elapsed = time.perf_counter() - c0
        rows.append((label, bool(hits), elapsed, counters))
        _log(
            f"  [{'ON' if canonical else 'OFF'}] {label}: "
            f"hit={bool(hits)} wall={elapsed:.3f}s "
            f"expansions={counters.rule_expansions} edits={counters.mol_edits}"
        )
    return time.perf_counter() - t0, rows


def dump_stats(pr: cProfile.Profile, stream: io.StringIO) -> None:
    ps = pstats.Stats(pr, stream=stream)
    ps.strip_dirs()
    stream.write("\n=== sort by cumulative (top %d) ===\n" % TOP_N)
    ps.sort_stats("cumulative")
    ps.print_stats(TOP_N)
    stream.write("\n=== sort by tottime (top %d) ===\n" % TOP_N)
    ps.sort_stats("tottime")
    ps.print_stats(TOP_N)


def _canon_cost(pr: cProfile.Profile) -> str:
    buf = io.StringIO()
    ps = pstats.Stats(pr, stream=buf)
    ps.strip_dirs()
    wall = float(ps.total_tt)
    interesting = (
        "ensure_lexical_orbit_representatives",
        "canonicalize_smarts_match",
        "canonicalize_pair_match",
        "automorphism_to_representative",
        "lexical_orbit_representatives",
        "restamp_product_forest_last_layer",
    )
    lines = [f"profile total_tt={wall:.4f}s"]
    for (filename, _line, func), (_cc, nc, _tt, ct, _callers) in ps.stats.items():
        if func in interesting:
            pct = (100.0 * ct / wall) if wall > 0 else 0.0
            lines.append(f"  {func}: cum={ct:.4f}s ({pct:.2f}%) ncalls={nc}")
    return "\n".join(lines)


def main() -> None:
    ART_OUT.parent.mkdir(parents=True, exist_ok=True)
    if ART_LIVE.exists():
        ART_LIVE.unlink()
    _log("=== canonical emitted sites profile (off vs on) ===")

    pr_off = cProfile.Profile()
    pr_off.enable()
    wall_off, rows_off = run_batch(canonical=False)
    pr_off.disable()

    pr_on = cProfile.Profile()
    pr_on.enable()
    wall_on, rows_on = run_batch(canonical=True)
    pr_on.disable()

    stream = io.StringIO()
    stream.write("canonical_emitted_sites profile (find_path PhaseOne)\n")
    stream.write(f"cases={len(CASES)} max_nodes={MAX_NODES} max_paths={MAX_PATHS}\n\n")
    stream.write(f"OFF wall={wall_off:.3f}s\n")
    for label, hit, elapsed, counters in rows_off:
        stream.write(
            f"  {label}: hit={hit} {elapsed:.3f}s "
            f"exp={counters.rule_expansions} edits={counters.mol_edits}\n"
        )
    stream.write(f"\nON wall={wall_on:.3f}s\n")
    for label, hit, elapsed, counters in rows_on:
        stream.write(
            f"  {label}: hit={hit} {elapsed:.3f}s "
            f"exp={counters.rule_expansions} edits={counters.mol_edits}\n"
        )
    delta = wall_on - wall_off
    pct = (100.0 * delta / wall_off) if wall_off > 0 else 0.0
    stream.write(f"\nDELTA on-off = {delta:+.3f}s ({pct:+.1f}%)\n\n")
    stream.write("=== OFF cProfile highlights ===\n")
    stream.write(_canon_cost(pr_off) + "\n")
    dump_stats(pr_off, stream)
    stream.write("\n=== ON cProfile highlights ===\n")
    stream.write(_canon_cost(pr_on) + "\n")
    dump_stats(pr_on, stream)

    ART_OUT.write_text(stream.getvalue(), encoding="utf-8")
    pr_on.dump_stats(str(ART_PSTATS))
    _log(f"wrote {ART_OUT}")
    _log(f"DELTA on-off = {delta:+.3f}s ({pct:+.1f}%)")


if __name__ == "__main__":
    main()
