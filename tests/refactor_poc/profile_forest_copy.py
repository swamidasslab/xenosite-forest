#!/usr/bin/env python3
"""Before/after cProfile harness for ``forest_copy`` vs ``copy.deepcopy`` on ``_forest``.

Same find_path case shape as ``profile_find_path_poc.py`` — stresses forest
copy on every work-copy / product path.

  uv run python tests/refactor_poc/profile_forest_copy.py before
  uv run python tests/refactor_poc/profile_forest_copy.py after
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


def _paths(tag: str) -> tuple[Path, Path, Path]:
    stem = f"forest_copy_{tag}"
    art = ROOT / "artifacts"
    return (
        art / f"{stem}.out",
        art / f"{stem}.pstats",
        art / f"{stem}.live.log",
    )


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


def copy_verdict(pr: cProfile.Profile) -> str:
    buf = io.StringIO()
    ps = pstats.Stats(pr, stream=buf)
    ps.strip_dirs()
    wall = float(ps.total_tt)
    deep = 0.0
    deep_ncalls = 0
    forest_cp = 0.0
    forest_ncalls = 0
    for (filename, _line, func), (_cc, nc, _tt, ct, _callers) in ps.stats.items():
        if filename.endswith("copy.py") and func == "deepcopy":
            deep = ct
            deep_ncalls = nc
        if func == "forest_copy":
            forest_cp = ct
            forest_ncalls = nc
    deep_pct = (100.0 * deep / wall) if wall > 0 else 0.0
    fc_pct = (100.0 * forest_cp / wall) if wall > 0 else 0.0
    return (
        f"copy verdict: deepcopy cum={deep:.3f}s ncalls={deep_ncalls} "
        f"~{deep_pct:.1f}% of profile wall={wall:.3f}s; "
        f"forest_copy cum={forest_cp:.3f}s ncalls={forest_ncalls} "
        f"~{fc_pct:.1f}% of wall"
    )


def main(argv: list[str]) -> int:
    if len(argv) != 2 or argv[1] not in {"before", "after"}:
        print("usage: profile_forest_copy.py {before|after}", file=sys.stderr)
        return 2
    tag = argv[1]
    art_out, art_pstats, art_live = _paths(tag)
    art_out.parent.mkdir(parents=True, exist_ok=True)
    live = art_live.open("w")
    lines: list[str] = []
    wall0 = time.perf_counter()

    def log(msg: str) -> None:
        stamp = time.perf_counter() - wall0
        line = f"[{stamp:8.3f}s] {msg}"
        print(line, flush=True)
        live.write(line + "\n")
        live.flush()
        lines.append(line)

    log(f"forest_copy cProfile tag={tag}")
    log(f"cases={len(CASES)}  max_nodes={POC_MAX_NODES}  max_paths={MAX_PATHS}")
    log(f"live log: {art_live}")
    log("")

    log("[warmup] anisole→phenol ...")
    hit, sec, c = run_case(*CASES[0])
    log(f"  hit={hit}  {sec:.3f}s  ed={c.mol_edits} re={c.rule_expansions} nd={c.nodes}")

    pr = cProfile.Profile()
    log("")
    log("[profile] enabling cProfile")
    pr.enable()
    case_wall0 = time.perf_counter()
    for i, (label, reactant, target) in enumerate(CASES, 1):
        log(f"[{i}/{len(CASES)}] {label} ...")
        hit, sec, c = run_case(label, reactant, target)
        log(
            f"  hit={hit}  {sec:.3f}s  "
            f"ed={c.mol_edits} re={c.rule_expansions} nd={c.nodes} "
            f"sites={c.sites_considered}/{c.sites_skipped}skip"
        )
    pr.disable()
    case_wall = time.perf_counter() - case_wall0
    log("[profile] disabled")
    log(f"perf_counter wall (profiled cases): {case_wall:.3f}s")

    verdict = copy_verdict(pr)
    log("")
    log(verdict)

    buf = io.StringIO()
    dump_stats(pr, buf)
    stats_text = buf.getvalue()
    print(stats_text, flush=True)
    live.write(stats_text)
    live.flush()

    pr.dump_stats(str(art_pstats))
    body = "\n".join(lines) + "\n" + stats_text
    body += f"\npstats written to {art_pstats}\n"
    art_out.write_text(body)
    log(f"wrote {art_out}")
    log(f"wrote {art_pstats}")
    live.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
