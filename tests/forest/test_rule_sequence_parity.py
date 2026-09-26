"""Multi-rule sequence parity: catch Kekulé / structure cache corruption.

Apply two or three distinct paired leaf rules in order on a small mol.
After each hop, compare RDKit CSMI product bags (Python vs Rust). Failures
are appended to :mod:`tests.forest.sequence_parity_replay` for replay
(commit that file when a new case appears).
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest
from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.forest.find_path_rust import native_available
from xenosite.forest.rdkit_api import MolFromSmiles, MolToSmiles

from .pattern_info_inventory import instantiate_rule
from .rule_parity_pairs import paired_rule_names, python_leaf_classes
from .sequence_parity_replay import SEQUENCE_PARITY_REPLAY

pytestmark = pytest.mark.skipif(
    not native_available(),
    reason="xenosite-forest-native not installed",
)

_REPLAY_PATH = Path(__file__).resolve().parent / "sequence_parity_replay.py"

_STARTS = (
    "CCO",
    "CCN",
    "CC",
    "c1ccccc1",
    "Oc1ccccc1",
    "COc1ccccc1",
    "CC(=O)Nc1ccc(O)cc1",
    "c1ccc2[nH]ccc2c1",
    "c1ccsc1",
)

# Leaves that are cheap and commonly chained (avoid conjugates / rare rules).
_SEQUENCE_RULES = (
    "Hydroxylation",
    "Dehydrogenation",
    "Dealkylation",
    "NDealkylation",
    "Epoxidation",
    "Dehydration",
)


def _rdkit_csmi(smiles: str) -> str | None:
    mol = MolFromSmiles(smiles)
    if mol is None:
        return None
    RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _python_bag(rule_name: str, smiles: str) -> set[str]:
    cls = python_leaf_classes()[rule_name]
    rule = instantiate_rule(cls)
    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    out: set[str] = set()
    for products, _info in rule.metabolize(mol):
        for product in products:
            csmi = _rdkit_csmi(
                MolToSmiles(product, canonical=True, isomericSmiles=False)
            )
            if csmi:
                out.add(csmi)
    return out


def _rust_bag(rule_name: str, smiles: str) -> set[str]:
    import xenosite_forest as native  # type: ignore[import-not-found]

    rs = native.RuleSet.leaf(rule_name)
    mol = native.ForestMol(smiles)
    out: set[str] = set()
    for row in rs.metabolize(mol):
        for product in row[4]:
            for piece in product.split("."):
                csmi = _rdkit_csmi(piece)
                if csmi:
                    out.add(csmi)
    return out


def _record_failure(start: str, rules: tuple[str, ...]) -> None:
    """Append a replay case if missing (caller should commit the file)."""

    case = (start, rules)
    if case in SEQUENCE_PARITY_REPLAY:
        return
    existing = list(SEQUENCE_PARITY_REPLAY)
    existing.append(case)
    lines = [
        '"""Committed replay cases for multi-rule sequence parity.',
        "",
        "Appended automatically when :mod:`test_rule_sequence_parity` discovers a new",
        "failure (write the tuple, commit this file). Each entry is",
        "``(start_smiles, (rule_name, ...))`` — apply the named leaves in order and",
        "assert Rust↔Python product bags match after every hop.",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "SEQUENCE_PARITY_REPLAY: tuple[tuple[str, tuple[str, ...]], ...] = (",
    ]
    for smi, rs in existing:
        rs_repr = ", ".join(repr(r) for r in rs)
        lines.append(f"    ({smi!r}, ({rs_repr},)),")
    lines.append(")")
    lines.append("")
    _REPLAY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _assert_sequence(start: str, rules: tuple[str, ...]) -> None:
    paired = set(paired_rule_names())
    for name in rules:
        if name not in paired:
            pytest.skip(f"{name} not paired")
    current_py = start
    current_rs = start
    for hop, rule_name in enumerate(rules):
        py = _python_bag(rule_name, current_py)
        rs = _rust_bag(rule_name, current_rs)
        if py != rs:
            _record_failure(start, rules)
            only_py = sorted(py - rs)
            only_rs = sorted(rs - py)
            raise AssertionError(
                f"sequence parity fail at hop {hop} ({rule_name}) "
                f"start={start!r} rules={rules}:\n"
                f"  input_py={current_py!r} input_rs={current_rs!r}\n"
                f"  only Python ({len(only_py)}): {only_py[:12]}\n"
                f"  only Rust   ({len(only_rs)}): {only_rs[:12]}\n"
                f"  recorded into {_REPLAY_PATH.name} — commit for replay"
            )
        # Advance along a shared product when both sides agree and non-empty.
        shared = sorted(py & rs)
        if not shared:
            return
        # Prefer a mono-component organic product for the next hop.
        nxt = next((s for s in shared if "." not in s and "[" not in s), shared[0])
        current_py = nxt
        current_rs = nxt


@pytest.mark.parametrize("start,rules", SEQUENCE_PARITY_REPLAY)
def test_sequence_parity_replay(start: str, rules: tuple[str, ...]) -> None:
    _assert_sequence(start, rules)


def test_sequence_parity_random_hops() -> None:
    """Random 2–3 distinct paired leaves on a handful of starts."""

    paired = [n for n in _SEQUENCE_RULES if n in set(paired_rule_names())]
    if len(paired) < 2:
        pytest.skip("need ≥2 paired sequence rules")
    rng = random.Random(int(os.environ.get("XENOSITE_SEQ_SEED", "42")))
    n_trials = int(os.environ.get("XENOSITE_SEQ_TRIALS", "24"))
    for _ in range(n_trials):
        start = rng.choice(_STARTS)
        k = rng.choice((2, 3))
        rules = tuple(rng.sample(paired, k=min(k, len(paired))))
        _assert_sequence(start, rules)
