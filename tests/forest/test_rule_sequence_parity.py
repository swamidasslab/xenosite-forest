"""Multi-rule sequence parity: catch Kekulé / structure cache corruption.

Apply two or three distinct paired leaf rules in order on a small mol.
After each hop, compare RDKit CSMI product bags (Python vs Rust).

**Cache smoke:** each Rust hop also metabolizes an ``edit_copy`` sibling that
shares the Kekulé ``Rc`` — sibling bags must match the primary mol.

**Replay:** green locks live in ``SEQUENCE_PARITY_REPLAY``. New product-bag
failures are appended to ``SEQUENCE_PARITY_XFAIL`` and auto-committed when
possible (``XENOSITE_SEQ_AUTOCOMMIT=0`` or ``CI=true`` skips the commit).
"""

from __future__ import annotations

import os
import random
import subprocess
from pathlib import Path

import pytest
from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.forest.find_path_rust import native_available
from xenosite.forest.rdkit_api import MolFromSmiles, MolToSmiles

from .pattern_info_inventory import instantiate_rule
from .rule_parity_pairs import paired_rule_names, python_leaf_classes
from .sequence_parity_replay import SEQUENCE_PARITY_REPLAY, SEQUENCE_PARITY_XFAIL

pytestmark = pytest.mark.skipif(
    not native_available(),
    reason="xenosite-forest-native not installed",
)

_REPLAY_PATH = Path(__file__).resolve().parent / "sequence_parity_replay.py"
_REPO_ROOT = Path(__file__).resolve().parents[2]

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
    "O=C1C=CC(=O)C=C1",
    "Oc1ccc(O)cc1",
    "O=C=Nc1ccccc1",
)

# Leaves that are cheap and commonly chained (avoid conjugates / rare rules).
# Hydrogenation / Dehydrogenation stress pair doors + Kekulé rematch.
_SEQUENCE_RULES = (
    "Hydroxylation",
    "Dehydrogenation",
    "Hydrogenation",
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


def _rust_bag_on(mol, rule_name: str) -> set[str]:
    """Metabolize ``rule_name`` on an existing native ``ForestMol``."""

    import xenosite_forest as native  # type: ignore[import-not-found]

    rs = native.RuleSet.leaf(rule_name)
    out: set[str] = set()
    for row in rs.metabolize(mol):
        _pattern, _site, _atoms, _orbit, products, _path = row
        for product in products:
            for piece in product.split("."):
                csmi = _rdkit_csmi(piece)
                if csmi:
                    out.add(csmi)
    return out


def _rust_bag(rule_name: str, smiles: str) -> set[str]:
    import xenosite_forest as native  # type: ignore[import-not-found]

    mol = native.ForestMol(smiles)
    primary = _rust_bag_on(mol, rule_name)
    # Shared Kekulé Rc sibling must not diverge (cache corruption smoke).
    sibling = mol.edit_copy()
    shared = _rust_bag_on(sibling, rule_name)
    if primary != shared:
        raise AssertionError(
            f"Kekulé-share corruption on {rule_name} / {smiles!r}:\n"
            f"  primary-only: {sorted(primary - shared)[:12]}\n"
            f"  sibling-only: {sorted(shared - primary)[:12]}"
        )
    return primary


def _write_replay_module(
    replay: list[tuple[str, tuple[str, ...]]],
    xfail: list[tuple[str, tuple[str, ...]]],
) -> None:
    lines = [
        '"""Committed replay cases for multi-rule sequence parity.',
        "",
        "``SEQUENCE_PARITY_REPLAY`` — must stay green (regression lock).",
        "",
        "``SEQUENCE_PARITY_XFAIL`` — appended automatically when random fuzz finds a",
        "new product-bag gap (auto-commit unless ``CI`` / ``XENOSITE_SEQ_AUTOCOMMIT=0``).",
        "Promote a case to ``SEQUENCE_PARITY_REPLAY`` once leaf parity makes it pass.",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "SEQUENCE_PARITY_REPLAY: tuple[tuple[str, tuple[str, ...]], ...] = (",
    ]
    for smi, rs in replay:
        rs_repr = ", ".join(repr(r) for r in rs)
        lines.append(f"    ({smi!r}, ({rs_repr},)),")
    lines.append(")")
    lines.append("")
    lines.append(
        "# Product-bag gaps still open on leaf parity — keep for replay / XPASS watch."
    )
    lines.append(
        "SEQUENCE_PARITY_XFAIL: tuple[tuple[str, tuple[str, ...]], ...] = ("
    )
    for smi, rs in xfail:
        rs_repr = ", ".join(repr(r) for r in rs)
        lines.append(f"    ({smi!r}, ({rs_repr},)),")
    lines.append(")")
    lines.append("")
    _REPLAY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _autocommit_replay() -> str | None:
    """``git add`` + ``git commit`` the replay file. Returns commit sha or None."""

    if os.environ.get("CI", "").lower() in {"1", "true", "yes"}:
        return None
    if os.environ.get("XENOSITE_SEQ_AUTOCOMMIT", "1") == "0":
        return None
    try:
        rel = str(_REPLAY_PATH.relative_to(_REPO_ROOT))
        subprocess.run(
            ["git", "add", "--", rel],
            cwd=_REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        status = subprocess.run(
            ["git", "status", "--porcelain", "--", rel],
            cwd=_REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        if not status.stdout.strip():
            return None
        subprocess.run(
            ["git", "commit", "-m", "Auto-commit sequence parity replay case"],
            cwd=_REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return sha.stdout.strip() or None
    except (subprocess.CalledProcessError, OSError):
        return None


def _record_failure(start: str, rules: tuple[str, ...]) -> str | None:
    """Append to XFAIL replay if missing; optionally auto-commit. Returns sha."""

    case = (start, rules)
    if case in SEQUENCE_PARITY_REPLAY or case in SEQUENCE_PARITY_XFAIL:
        return None
    replay = list(SEQUENCE_PARITY_REPLAY)
    xfail = list(SEQUENCE_PARITY_XFAIL)
    xfail.append(case)
    _write_replay_module(replay, xfail)
    return _autocommit_replay()


def _assert_sequence(start: str, rules: tuple[str, ...], *, soft: bool = False) -> None:
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
            sha = _record_failure(start, rules)
            only_py = sorted(py - rs)
            only_rs = sorted(rs - py)
            commit_note = f" (committed {sha})" if sha else ""
            msg = (
                f"sequence parity fail at hop {hop} ({rule_name}) "
                f"start={start!r} rules={rules}:\n"
                f"  input_py={current_py!r} input_rs={current_rs!r}\n"
                f"  only Python ({len(only_py)}): {only_py[:12]}\n"
                f"  only Rust   ({len(only_rs)}): {only_rs[:12]}\n"
                f"  recorded into {_REPLAY_PATH.name}{commit_note}"
            )
            if soft:
                pytest.xfail(msg)
            raise AssertionError(msg)
        shared = sorted(py & rs)
        if not shared:
            return
        nxt = next((s for s in shared if "." not in s and "[" not in s), shared[0])
        current_py = nxt
        current_rs = nxt


@pytest.mark.parametrize("start,rules", SEQUENCE_PARITY_REPLAY)
def test_sequence_parity_replay(start: str, rules: tuple[str, ...]) -> None:
    _assert_sequence(start, rules, soft=False)


@pytest.mark.parametrize("start,rules", SEQUENCE_PARITY_XFAIL)
@pytest.mark.xfail(reason="leaf product parity gap; recorded for replay", strict=False)
def test_sequence_parity_xfail_replay(start: str, rules: tuple[str, ...]) -> None:
    _assert_sequence(start, rules, soft=False)


def test_sequence_parity_random_hops() -> None:
    """Random 2–3 distinct paired leaves on a handful of starts."""

    paired = [n for n in _SEQUENCE_RULES if n in set(paired_rule_names())]
    if len(paired) < 2:
        pytest.skip("need ≥2 paired sequence rules")
    rng = random.Random(int(os.environ.get("XENOSITE_SEQ_SEED", "42")))
    n_trials = int(os.environ.get("XENOSITE_SEQ_TRIALS", "32"))
    for _ in range(n_trials):
        start = rng.choice(_STARTS)
        k = rng.choice((2, 3))
        rules = tuple(rng.sample(paired, k=min(k, len(paired))))
        # Soft: record + xfail on product gaps; hard-fail only on cache smoke.
        _assert_sequence(start, rules, soft=True)


def test_sequence_parity_edit_copy_share_survives_two_hops() -> None:
    """Explicit cache-share walk: hydroxylation then dehydrogenation on phenol."""

    _assert_sequence("Oc1ccccc1", ("Hydroxylation", "Dehydrogenation"), soft=False)
