#!/usr/bin/env python3
"""Live Python ``find_path`` wall times for Rust H2H (same cases / knobs).

Pairs with ``cargo run -p xenosite-forest --example find_path_bench --release``.

  uv run python tests/forest/bench_find_path_rust_h2h.py
  uv run python tests/forest/bench_find_path_rust_h2h.py --larger
  uv run python tests/forest/bench_find_path_rust_h2h.py --hard
"""

from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

from xenosite.forest.find_path import PathCounters, find_path
from xenosite.forest.rdkitutil import canon_smiles
from xenosite.forest.rulesets import PhaseOne

ROOT = Path(__file__).resolve().parents[2]
MAX_NODES = 800
MAX_PATHS = 1
REPEATS = 5

# Same rows as crates/xenosite-forest/examples/find_path_bench.rs
CASES: list[tuple[str, str, str]] = [
    ("eugenol→allyl-quinone", "COc1ccc(CC=C)cc1O", "O=C1C=CC(=O)C(CC=C)=C1"),
    ("dimethoxy-PEA→catechol", "COc1ccc(CCN)cc1OC", "NCCc1ccc(O)c(O)c1"),
    ("MeOPhOH→hydroxyquinone", "COc1ccc(O)cc1", "O=C1C=C(O)C(=O)C(O)=C1"),
    # Stereo form matches chematic product; RDKit canon of this spelling hits.
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        r"C(#C/C=C/C=O)C(C)(C)C",
    ),
    ("2-MeO-naph→1,2-NQ", "COc1ccc2ccccc2c1", "O=C1C(=O)c2ccccc2C=C1"),
]

LARGER: list[tuple[str, str, str]] = [
    (
        "tBu-bis-ND→dialdehyde",
        "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
        "O=Cc1ccc(C=O)cc1",
    ),
    (
        "macrocycle-ND→aminoK",
        "C1CCCCCCNC2CCCC(CC2)NCCCC1",
        "NC1CCCC(=O)CC1",
    ),
    (
        "tribenzyl→PhCHO",
        "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
        "O=Cc1ccccc1",
    ),
    (
        "triPh-butyl→OH",
        "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
        "Oc1ccccc1CCCCc2ccccc2CCCCc3ccccc3",
    ),
    (
        "MeO-diphenyl→catechol",
        "COc1ccc(Cc2ccc(OC)cc2)cc1",
        "Oc1ccc(Cc2ccc(O)cc2)cc1",
    ),
]

# Larger scaffolds with ≥3–8 PhaseOne hops (paired with find_path_bench --hard).
HARD: list[tuple[str, str, str]] = [
    (
        "trimethoxy-PEA→catechol",
        "COc1cc(OC)c(OC)c(CCN)c1",
        "NCCc1cc(O)c(O)c(O)c1",
    ),
    (
        "eugenol-MeO→allylQ",
        "COc1cc(CC=C)cc(OC)c1O",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "bisMeO-naph→1,2NQ",
        "COc1ccc2c(OC)cccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
    ),
    (
        "tetraMeO-biphenyl→tetraOH",
        "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC",
        "Oc1ccc(-c2ccc(O)c(O)c2)cc1O",
    ),
    (
        "veratrole-allyl→allylQ",
        "COc1ccc(CC=C)c(OC)c1OC",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "tetraMeO-naph→polyOH-NQ",
        "COc1cc(OC)c2c(OC)cc(OC)cc2c1",
        "O=C1C=C(O)C(=O)c2c(O)cc(O)cc12",
    ),
]


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


def run_one(reactant: str, target: str) -> tuple[bool, float, int, int, int, int]:
    want = canon_smiles(target)
    # Warmup
    find_path(
        reactant,
        target,
        ruleset=PhaseOne,
        counters=PathCounters(),
        max_paths=MAX_PATHS,
        max_nodes=MAX_NODES,
    )
    best: tuple[bool, float, int, int, int, int] | None = None
    for _ in range(REPEATS):
        counters = PathCounters()
        t0 = time.perf_counter()
        hits = list(
            find_path(
                reactant,
                target,
                ruleset=PhaseOne,
                counters=counters,
                max_paths=MAX_PATHS,
                max_nodes=MAX_NODES,
            )
        )
        elapsed = time.perf_counter() - t0
        product = canon_smiles(hits[0].smiles) if hits else None
        hit = product == want
        steps = len(hits[0].plan.steps) if hits and hit else 0
        row = (
            hit,
            elapsed,
            steps,
            int(counters.nodes),
            int(counters.mol_edits),
            int(counters.billed),
        )
        if best is None or elapsed < best[1]:
            best = row
    assert best is not None
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--larger", action="store_true")
    group.add_argument("--hard", action="store_true")
    args = parser.parse_args()
    if args.hard:
        cases = HARD
        title = "hard HA≈14–20 · multi-step (≥3–8 hops)"
    elif args.larger:
        cases = LARGER
        title = "larger HA≈17–26"
    else:
        cases = CASES
        title = "mid-size / multi-edit"

    print(
        f"Python live find_path PhaseOne  max_nodes={MAX_NODES}  "
        f"max_paths={MAX_PATHS}  best-of-{REPEATS}  sha={_git_sha()}"
    )
    print(f"\n=== {title} (use_filters=True / atom_diff) ===")
    print(
        f"{'case':<32} {'hit':>4} {'seconds':>9} {'steps':>5} "
        f"{'nodes':>6} {'edits':>7} {'bill':>6}"
    )
    total = 0.0
    for name, reactant, target in cases:
        hit, seconds, steps, nodes, edits, billed = run_one(reactant, target)
        total += seconds
        print(
            f"{name:<32} {'ok' if hit else 'MISS':>4} {seconds:9.3f} "
            f"{steps:5d} {nodes:6d} {edits:7d} {billed:6d}"
        )
    print(f"{'TOTAL':<32} {'':>4} {total:9.3f}")


if __name__ == "__main__":
    main()
