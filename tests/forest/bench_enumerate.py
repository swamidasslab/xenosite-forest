#!/usr/bin/env python3
"""Live Python ``bfs`` / ``dfs`` wall times vs Rust ``enumerate_bench``.

Pairs with::

  cargo run -p xenosite-forest --example enumerate_bench --release

  uv run python tests/forest/bench_enumerate.py
  uv run python tests/forest/bench_enumerate.py --phase-one
  uv run python tests/forest/bench_enumerate.py --compare   # run Rust too, print side-by-side
"""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

from xenosite.forest.find_path import bfs, dfs
from xenosite.forest.rules import Hydroxylation
from xenosite.forest.rulesets import PhaseOne, RuleSet

ROOT = Path(__file__).resolve().parents[2]
REPEATS = 5

OH = RuleSet((Hydroxylation,), name="OH")

OH_CASES: list[tuple[str, str, int, str]] = [
    ("OH ethane d1 bfs", "CC", 1, "bfs"),
    ("OH ethane d1 dfs", "CC", 1, "dfs"),
    ("OH ethane d2 bfs", "CC", 2, "bfs"),
    ("OH ethane d2 dfs", "CC", 2, "dfs"),
    ("OH ethane d3 bfs", "CC", 3, "bfs"),
    ("OH benzene d2 bfs", "c1ccccc1", 2, "bfs"),
    ("OH benzene d3 bfs", "c1ccccc1", 3, "bfs"),
    ("OH butylbenzene d2 bfs", "c1ccc(CCCC)cc1", 2, "bfs"),
    ("OH butylbenzene d2 dfs", "c1ccc(CCCC)cc1", 2, "dfs"),
    ("OH butylbenzene d3 bfs", "c1ccc(CCCC)cc1", 3, "bfs"),
    ("OH toluene d2 bfs", "Cc1ccccc1", 2, "bfs"),
    ("OH toluene d3 bfs", "Cc1ccccc1", 3, "bfs"),
]

PHASE_ONE_CASES: list[tuple[str, str, int, str]] = [
    ("P1 ethane d2 bfs", "CC", 2, "bfs"),
    ("P1 ethane d3 bfs", "CC", 3, "bfs"),
    ("P1 anisole d1 bfs", "COc1ccccc1", 1, "bfs"),
    ("P1 anisole d2 bfs", "COc1ccccc1", 2, "bfs"),
    ("P1 anisole d2 dfs", "COc1ccccc1", 2, "dfs"),
    ("P1 anisole d3 bfs", "COc1ccccc1", 3, "bfs"),
    ("P1 eugenol d1 bfs", "COc1ccc(CC=C)cc1O", 1, "bfs"),
    ("P1 eugenol d2 bfs", "COc1ccc(CC=C)cc1O", 2, "bfs"),
    ("P1 veratrole d2 bfs", "COc1ccc(OC)cc1", 2, "bfs"),
    ("P1 veratrole d3 bfs", "COc1ccc(OC)cc1", 3, "bfs"),
    ("P1 phenacetin d2 bfs", "CCOc1ccc(NC(C)=O)cc1", 2, "bfs"),
]


def _git_sha() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
            )
            .strip()
        )
    except Exception:
        return "?"


def run_one(smiles: str, depth: int, order: str, ruleset) -> tuple[int, float, set[str]]:
    enum = bfs if order == "bfs" else dfs
    # Warmup
    list(enum(smiles, ruleset, depth=depth))
    best = float("inf")
    n = 0
    smiles_set: set[str] = set()
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        hits = list(enum(smiles, ruleset, depth=depth))
        elapsed = time.perf_counter() - t0
        n = len(hits)
        smiles_set = {mol.xf.csmi for mol, _ in hits}
        if elapsed < best:
            best = elapsed
    return n, best, smiles_set


def bench_suite(title: str, cases: list[tuple[str, str, int, str]], ruleset) -> dict[str, tuple[int, float]]:
    print(f"\n=== {title} (best-of-{REPEATS}) ===")
    print(f"{'case':<28} {'n':>6} {'seconds':>12} {'µs/hit':>8}")
    out: dict[str, tuple[int, float]] = {}
    total_n = 0
    total_t = 0.0
    for name, smiles, depth, order in cases:
        n, secs, _ = run_one(smiles, depth, order, ruleset)
        us_per = (secs * 1e6 / n) if n else 0.0
        print(f"{name:<28} {n:>6} {secs:>12.6f} {us_per:>8.1f}")
        out[name] = (n, secs)
        total_n += n
        total_t += secs
    print(f"{'TOTAL':<28} {total_n:>6} {total_t:>12.6f}")
    return out


def run_rust(extra: list[str]) -> dict[str, tuple[int, float]]:
    cmd = [
        "cargo",
        "run",
        "-p",
        "xenosite-forest",
        "--example",
        "enumerate_bench",
        "--release",
        "--",
        *extra,
    ]
    print("+", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)
    print(proc.stdout)
    parsed: dict[str, tuple[int, float]] = {}
    for line in proc.stdout.splitlines():
        if len(line) < 28:
            continue
        name = line[:28].rstrip()
        rest = line[28:].split()
        if name in {"case", "TOTAL"} or name.startswith("=") or name.startswith("Rust"):
            continue
        if len(rest) >= 2:
            try:
                n = int(rest[0])
                secs = float(rest[1])
            except ValueError:
                continue
            parsed[name] = (n, secs)
    return parsed


def compare(py: dict[str, tuple[int, float]], rs: dict[str, tuple[int, float]]) -> None:
    print("\n=== Python vs Rust ===")
    print(
        f"{'case':<28} {'py_n':>6} {'rs_n':>6} {'n_ok':>5} "
        f"{'py_s':>10} {'rs_s':>10} {'speedup':>8}"
    )
    names = sorted(set(py) | set(rs), key=lambda k: (k.startswith("P1"), k))
    for name in names:
        if name not in py or name not in rs:
            print(f"{name:<28}  (missing one side)")
            continue
        pn, pt = py[name]
        rn, rt = rs[name]
        n_ok = "yes" if pn == rn else "NO*"
        speed = pt / rt if rt > 0 else float("inf")
        speed_s = f"{speed:>7.2f}x" if speed < 9999 else "  >9999x"
        print(
            f"{name:<28} {pn:>6} {rn:>6} {n_ok:>5} "
            f"{pt:>10.6f} {rt:>10.6f} {speed_s}"
        )
    print(
        "\n* PhaseOne product counts can differ (Python RDKit PhaseOne vs Rust "
        "chematic phase_one catalog/doors). Hydroxylation-only rows should match."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--phase-one", action="store_true")
    group.add_argument("--oh", action="store_true")
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Also run Rust enumerate_bench and print a speedup table",
    )
    args = parser.parse_args()

    print(f"Python enumerate bfs/dfs  best-of-{REPEATS}  sha={_git_sha()}")
    py: dict[str, tuple[int, float]] = {}
    rust_extra: list[str] = []
    if args.phase_one:
        py.update(bench_suite("PhaseOne", PHASE_ONE_CASES, PhaseOne))
        rust_extra = ["--phase-one"]
    elif args.oh:
        py.update(bench_suite("Hydroxylation-only", OH_CASES, OH))
        rust_extra = ["--oh"]
    else:
        py.update(bench_suite("Hydroxylation-only", OH_CASES, OH))
        py.update(bench_suite("PhaseOne", PHASE_ONE_CASES, PhaseOne))

    if args.compare:
        rs = run_rust(rust_extra)
        compare(py, rs)


if __name__ == "__main__":
    main()
