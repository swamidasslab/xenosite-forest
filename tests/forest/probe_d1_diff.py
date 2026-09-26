#!/usr/bin/env python3
"""Dump / compare PhaseOne depth-1 products vs Rust (same chematic keys).

Python emits ``rule\\tpattern\\trdkit_csmi``. Rekey through the Rust
``probe_d1_diff rekey`` helper so both sides share Chematic
``canonical_smiles_stable_key`` before set-diff.

Usage::

  uv run python tests/forest/probe_d1_diff.py dump ibuprofen 'CC(C)Cc1ccc(C(C)C(=O)O)cc1' > /tmp/py.tsv
  cargo run -p xenosite-forest --example probe_d1_diff --release -- dump ibuprofen '…' > /tmp/rs.tsv
  cargo run -p xenosite-forest --example probe_d1_diff --release -- rekey < /tmp/py.tsv > /tmp/py_keyed.tsv
  cargo run -p xenosite-forest --example probe_d1_diff --release -- compare ibuprofen /tmp/py_keyed.tsv /tmp/rs.tsv
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from xenosite.forest.find_path import bfs
from xenosite.forest.rulesets import PhaseOne

ROOT = Path(__file__).resolve().parents[2]


def _rule_name(info: dict) -> str:
    leaf = info["rule"][0]
    return getattr(leaf, "name", None) or type(leaf).__name__


def _pattern_name(info: dict) -> str:
    pat = info.get("pattern") or {}
    if isinstance(pat, dict):
        return str(pat.get("name") or "?")
    return "?"


def dump(name: str, smiles: str) -> None:
    # depth-1 bfs; unique by display CSMI already. Emit every path's rule.
    # Same product from two rules → two rows (compare joins on key).
    by_csmi: dict[str, dict[str, str]] = defaultdict(dict)
    for mol, info in bfs(smiles, PhaseOne, depth=1):
        csmi = mol.xf.csmi
        rule = _rule_name(info)
        pat = _pattern_name(info)
        by_csmi[csmi].setdefault(rule, pat)
    print(f"# python {name}: bfs_d1_unique={len(by_csmi)}", file=sys.stderr)
    for csmi, rules in sorted(by_csmi.items()):
        for rule, pat in sorted(rules.items()):
            print(f"{rule}\t{pat}\t{csmi}")


def run_compare(name: str, smiles: str) -> None:
    py_raw = ROOT / "artifacts" / f"d1_py_{name}.tsv"
    py_key = ROOT / "artifacts" / f"d1_py_{name}_keyed.tsv"
    rs_tsv = ROOT / "artifacts" / f"d1_rs_{name}.tsv"
    py_raw.parent.mkdir(exist_ok=True)

    # Python dump
    proc = subprocess.run(
        [sys.executable, __file__, "dump", name, smiles],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    sys.stderr.write(proc.stderr)
    py_raw.write_text(proc.stdout)

    rust = [
        "cargo",
        "run",
        "-p",
        "xenosite-forest",
        "--example",
        "probe_d1_diff",
        "--release",
        "--",
    ]
    # Rust dump
    proc = subprocess.run(
        [*rust, "dump", name, smiles],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    sys.stderr.write(proc.stderr)
    rs_tsv.write_text(proc.stdout)

    # Rekey Python through chematic
    proc = subprocess.run(
        [*rust, "rekey"],
        cwd=ROOT,
        input=py_raw.read_text(),
        capture_output=True,
        text=True,
        check=True,
    )
    py_key.write_text(proc.stdout)

    # Compare
    proc = subprocess.run(
        [*rust, "compare", name, str(py_key), str(rs_tsv)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)


def main() -> None:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dump")
    d.add_argument("name")
    d.add_argument("smiles")
    c = sub.add_parser("run")
    c.add_argument("name")
    c.add_argument("smiles")
    args = p.parse_args()
    if args.cmd == "dump":
        dump(args.name, args.smiles)
    else:
        run_compare(args.name, args.smiles)


if __name__ == "__main__":
    main()
