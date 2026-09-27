"""Rust↔RDKit leaf-rule product parity (parametric over cover intents).

Default: full paired×corpus cartesian via :func:`parity_param_cases`.
Focused CoverIntent-only with ``XENOSITE_PARITY_FULL=0`` /
``pytest --parity-focused`` (CI may disable for speed).

For each ``(rule, mol)``:

1. Collect emissions on both engines.
2. Collapse sites by **that engine's** topological ranks (orbits).
3. Assert equal unique topological site→product bags.
4. Assert equal product sets after **RDKit** CSMI.

Known open gaps (C18 chemistry-first plan) are ``xfail`` via
:mod:`parity_xfail_cases` — remove entries as Phase 1–3 closes them against
**approved chemistry**, not literal Python match.

Pairing / exceptions: :mod:`test_rule_parity_pairs` (attribute data on rules).
Corpus coverage: :mod:`test_rule_parity_corpus` (every pattern / every when).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from typing import Any

import pytest
from rdkit import Chem
from rdkit.Chem.rdmolops import RemoveStereochemistry

from xenosite.forest.find_path_rust import native_available
from xenosite.forest.rdkit_api import MolFromSmiles, MolToSmiles
from xenosite.forest.rules import ReactionRule

from .parity_xfail_cases import (
    C18_OPEN_PRODUCT_XFAIL_COUNT,
    PARITY_PRODUCT_XFAIL,
    parity_xfail_reason,
)
from .pattern_info_inventory import instantiate_rule
from .rule_parity_corpus import parity_full_enabled, parity_param_cases
from .rule_parity_pairs import paired_rule_names, python_leaf_classes

pytestmark = pytest.mark.skipif(
    not native_available(),
    reason="xenosite-forest-native not installed",
)


def _rdkit_csmi(smiles_or_mol: Any) -> str | None:
    """Canonical non-isomeric RDKit SMILES; ``None`` if unparseable."""

    if smiles_or_mol is None:
        return None
    if isinstance(smiles_or_mol, str):
        mol = MolFromSmiles(smiles_or_mol)
    else:
        mol = Chem.Mol(smiles_or_mol)
    if mol is None:
        return None
    RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _site_kind(rule: ReactionRule) -> str:
    kind = getattr(rule, "site_kind", "atom")
    assert kind in {"atom", "bond", "directed_bond", "atom_pair"}, kind
    return kind


def _topo_key(ranks: dict[int, int] | list[int], atoms: Iterable[int], kind: str) -> tuple:
    """Engine-local topological site key (ranks are not cross-engine)."""

    def rank(i: int) -> int:
        if isinstance(ranks, dict):
            return int(ranks[i])
        return int(ranks[i])

    vals = [rank(i) for i in atoms]
    if kind in {"bond", "atom_pair", "directed_bond"}:
        return tuple(sorted(vals))
    if len(vals) == 1:
        return (vals[0],)
    return tuple(sorted(vals))


def _python_emissions(
    rule: ReactionRule, smiles: str, kind: str
) -> tuple[Counter, set[str]]:
    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    ranks = mol.xf.topol_equiv
    by_topo: dict[tuple, set[str]] = {}
    all_products: set[str] = set()
    for products, info in rule.metabolize(mol):
        site = info.get("discovered_site", info["site"])
        if isinstance(site, (set, frozenset, tuple)):
            atoms = list(site)
        else:
            atoms = [int(site)]
        key = _topo_key(ranks, atoms, kind)
        bucket = by_topo.setdefault(key, set())
        for product in products:
            csmi = _rdkit_csmi(product)
            if csmi is None:
                continue
            bucket.add(csmi)
            all_products.add(csmi)
    site_product_bags = Counter(
        frozenset(prods) for prods in by_topo.values() if prods
    )
    return site_product_bags, all_products


def _rust_emissions(
    rule_name: str, smiles: str, kind: str
) -> tuple[Counter, set[str]]:
    import xenosite_forest as native  # type: ignore[import-not-found]

    rs = native.RuleSet.leaf(rule_name)
    mol = native.ForestMol(smiles)
    ranks = mol.ranks()
    by_topo: dict[tuple, set[str]] = {}
    all_products: set[str] = set()
    for _pattern, _site, site_atoms, _orbit, products, _path in rs.metabolize(mol):
        atoms = list(site_atoms) if site_atoms else [_site]
        key = _topo_key(ranks, atoms, kind)
        bucket = by_topo.setdefault(key, set())
        for product in products:
            for piece in product.split("."):
                csmi = _rdkit_csmi(piece)
                if csmi is None:
                    continue
                bucket.add(csmi)
                all_products.add(csmi)
    site_product_bags = Counter(
        frozenset(prods) for prods in by_topo.values() if prods
    )
    return site_product_bags, all_products


def _assert_parity(rule_name: str, smiles: str) -> None:
    cls = python_leaf_classes()[rule_name]
    rule = instantiate_rule(cls)
    kind = _site_kind(rule)
    py_sites, py_prods = _python_emissions(rule, smiles, kind)
    rs_sites, rs_prods = _rust_emissions(rule_name, smiles, kind)

    if py_prods != rs_prods:
        only_py = sorted(py_prods - rs_prods)
        only_rs = sorted(rs_prods - py_prods)
        raise AssertionError(
            f"product CSMI mismatch for {rule_name} on {smiles!r}:\n"
            f"  only Python ({len(only_py)}): {only_py[:12]}\n"
            f"  only Rust   ({len(only_rs)}): {only_rs[:12]}"
        )

    if sum(py_sites.values()) != sum(rs_sites.values()):
        raise AssertionError(
            f"topological site count mismatch for {rule_name} on {smiles!r}: "
            f"python={sum(py_sites.values())} rust={sum(rs_sites.values())} "
            f"(product sets already match)"
        )

    if py_sites != rs_sites:
        raise AssertionError(
            f"site→product bag mismatch for {rule_name} on {smiles!r}: "
            f"python_bags={len(py_sites)} rust_bags={len(rs_sites)} "
            f"(counts matched, bag multiset differed)"
        )


def _paired_names() -> list[str]:
    if not native_available():
        return []
    return paired_rule_names()


def _parity_param_ids() -> list:
    """Parametrize with C18 xfail marks on known open gaps."""

    raw = (
        parity_param_cases(_PAIRED)
        if _PAIRED
        else (("Hydroxylation", "CCO"),)
    )
    out = []
    for rule_name, smiles in raw:
        reason = parity_xfail_reason(rule_name, smiles)
        if reason:
            out.append(
                pytest.param(
                    rule_name,
                    smiles,
                    marks=pytest.mark.xfail(reason=reason, strict=False),
                    id=f"{rule_name}-{smiles}",
                )
            )
        else:
            out.append(pytest.param(rule_name, smiles, id=f"{rule_name}-{smiles}"))
    return out


_PAIRED = _paired_names()
_CASES = _parity_param_ids()


@pytest.mark.parametrize("rule_name,smiles", _CASES)
def test_leaf_rule_product_parity(rule_name: str, smiles: str) -> None:
    """Paired leaf × cover mol (or full cartesian when parity_full enabled)."""

    if not _PAIRED:
        pytest.skip("no paired leaf rules")
    _assert_parity(rule_name, smiles)


@pytest.mark.parametrize("rule_name", _PAIRED or ["Hydroxylation"])
def test_leaf_parity_on_example_substrate(rule_name: str) -> None:
    """One guaranteed hit per paired leaf (``_example_substrates``)."""

    if not _PAIRED:
        pytest.skip("no paired leaf rules")
    cls = python_leaf_classes()[rule_name]
    examples = getattr(cls, "_example_substrates", ()) or ()
    assert examples, f"{rule_name} lacks _example_substrates"
    reason = parity_xfail_reason(rule_name, examples[0])
    if reason:
        pytest.xfail(reason)
    _assert_parity(rule_name, examples[0])


def test_parity_param_mode_defaults_to_full() -> None:
    """Default collection uses full cartesian unless focused is requested."""

    if not _PAIRED:
        pytest.skip("no paired leaf rules")
    from .rule_parity_corpus import (
        parity_rule_mol_cases,
        parity_rule_mol_cases_full,
    )

    focused = parity_rule_mol_cases(_PAIRED)
    full = parity_rule_mol_cases_full(_PAIRED)
    # Count underlying (rule, smiles) pairs — marks wrap pytest.param.
    n_cases = len(_CASES)
    if parity_full_enabled():
        assert n_cases == len(full)
        assert len(full) > len(focused)
    else:
        assert n_cases == len(focused)
        assert len(focused) < len(full)


def test_parity_product_xfail_table_covers_live_gaps_only() -> None:
    """Xfail table keys are unique; no stale PASS entries (strict cleanup)."""

    assert len(PARITY_PRODUCT_XFAIL) == len(set(PARITY_PRODUCT_XFAIL))
    assert len(PARITY_PRODUCT_XFAIL) == C18_OPEN_PRODUCT_XFAIL_COUNT, (
        f"update C18_OPEN_PRODUCT_XFAIL_COUNT when editing PARITY_PRODUCT_XFAIL "
        f"(baseline tag py-rust-parity-c18-baseline was {C18_OPEN_PRODUCT_XFAIL_COUNT})"
    )
    # Every xfail key must still fail; otherwise remove it (Phase close).
    stale: list[str] = []
    for (rule_name, smiles), reason in sorted(PARITY_PRODUCT_XFAIL.items()):
        try:
            _assert_parity(rule_name, smiles)
        except AssertionError:
            continue
        stale.append(f"{rule_name} {smiles!r} ({reason})")
    assert not stale, (
        "parity xfail entries now PASS — remove from PARITY_PRODUCT_XFAIL "
        "and decrement C18_OPEN_PRODUCT_XFAIL_COUNT:\n  "
        + "\n  ".join(stale)
    )


def test_parity_product_xfail_marks_every_param_gap() -> None:
    """No unmarked product-parity failure on the parametric grid."""

    if not _PAIRED:
        pytest.skip("no paired leaf rules")
    unmarked: list[str] = []
    for rule_name, smiles in parity_param_cases(_PAIRED):
        if parity_xfail_reason(rule_name, smiles):
            continue
        try:
            _assert_parity(rule_name, smiles)
        except AssertionError as err:
            unmarked.append(f"{rule_name} {smiles!r}: {str(err).splitlines()[0]}")
    assert not unmarked, (
        "product parity gaps missing from PARITY_PRODUCT_XFAIL "
        "(add with C18 reason, or fix):\n  " + "\n  ".join(unmarked)
    )
