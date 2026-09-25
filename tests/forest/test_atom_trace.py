"""Atom trace on ``mol.xf.tracing``. Indexes are 0-based. Depth 0 is the reactant.

The old suite spoke 1-based atom-map numbers. This API does not stamp those
maps. ``atom_root`` is the depth-0 index, or None when the atom was created
later. ``atom_origin`` stays the earliest recorded index, including a birth
index. Metabolize always installs a trace; there is no ``do_not_tag_atoms``.
"""

from __future__ import annotations

import inspect

import pytest
from rdkit import Chem
from test_rules import examples

from xenosite.forest.find_path import bfs, dfs
from xenosite.forest.rules import (
    Hydroxylation,
    ReactionRule,
    SmirksReactionRule,
    describe,
)
from xenosite.forest.rulesets import PhaseOne

_SKIP = frozenset(
    {
        "ReactionRule",
        "SmirksReactionRule",
        "ResonanceRule",
        "ResonancePairRule",
        "ConjugationRule",
        "TautomerRule",
    }
)

_FALLBACK = {
    "Acetylation": "C1=CC(=C(C=C1N)C(=O)O)O",
    "SulfurReduction": "c1ccccc1SSc1ccccc1",
    "NitroaromaticReduction": "[O-][N+](C1=CC2=C(C=C1)NC(CN=C2C3=CC=CC=C3Cl)=O)=O",
    "ThiopheneSulfurOxidation": "O=C(c1ccc(cc1)C(C(=O)O)C)c2sccc2",
    "NDealkylation": "CN(C)Cc1ccccc1",
}


def _rule_classes():
    from xenosite.forest import rules as forest_rules

    return [
        cls
        for _, cls in inspect.getmembers(forest_rules, inspect.isclass)
        if issubclass(cls, ReactionRule)
        and cls.__module__ == forest_rules.__name__
        and cls.__name__ not in _SKIP
    ]


def _reactant_for(name: str) -> str:
    if name in examples:
        return examples[name][0][1]
    return _FALLBACK[name]


def _heavy(mol):
    return [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]


def test_untraced_queries_are_none():
    mol = Chem.MolFromSmiles("CCO")
    tracing = mol.xf.tracing
    assert not tracing.active
    assert tracing.depth is None
    assert tracing.atom_indices(0) is None
    assert tracing.atom_depths(0) is None
    assert tracing.atom_root(0) is None
    assert tracing.atom_origin(0) is None
    assert tracing.atom_added_by(0) is None
    assert tracing.removed_roots() == frozenset()


def test_stamp_is_identity_and_skips_hydrogen():
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    stamped = mol.xf.tracing._stamp()
    assert stamped.xf.tracing.active
    assert stamped.xf.tracing.depth == 0
    heavies = _heavy(stamped)
    assert len(heavies) == 1
    carbon = heavies[0].GetIdx()
    assert stamped.xf.tracing.atom_root(carbon) == carbon
    assert stamped.xf.tracing.atom_indices(carbon) == (carbon,)
    assert stamped.xf.tracing.atom_depths(carbon) == (0,)
    for atom in stamped.GetAtoms():
        if atom.GetAtomicNum() == 1:
            assert stamped.xf.tracing.atom_indices(atom.GetIdx()) is None


def test_hydroxylation_cco_roots_and_added_oxygen():
    # Carbon 0 of CCO gains OH. The new oxygen has no depth-0 root.
    mol = Chem.MolFromSmiles("CCO")
    products, info = next(Hydroxylation().metabolize(mol))
    product = products[0]
    assert info["site"] == frozenset({0})
    tracing = product.xf.tracing
    assert tracing.active
    assert tracing.depth == 1
    assert tracing.removed_roots() == frozenset()

    roots = {}
    added = []
    for atom in _heavy(product):
        idx = atom.GetIdx()
        assert tracing.atom_indices(idx) is not None
        assert tracing.atom_depths(idx) is not None
        assert len(tracing.atom_indices(idx)) == len(tracing.atom_depths(idx))
        root = tracing.atom_root(idx)
        if root is None:
            added.append(atom)
            assert tracing.atom_added_by(idx) == ("Hydroxylation", frozenset({0}))
            assert tracing.atom_depths(idx) == (1,)
        else:
            roots[root] = idx
            assert tracing.atom_added_by(idx) is None
            assert tracing.atom_depths(idx) == (0, 1)
    assert set(roots) == {0, 1, 2}
    assert len(added) == 1
    assert added[0].GetAtomicNum() == 8
    # The carbon that was index 0 is still that root, now possibly moved.
    assert product.GetAtomWithIdx(roots[0]).GetAtomicNum() == 6


def test_restamp_does_not_rewrite_history():
    mol = Chem.MolFromSmiles("CCO")
    products, _info = next(Hydroxylation().metabolize(mol))
    product = products[0]
    before = {
        atom.GetIdx(): product.xf.tracing.atom_indices(atom.GetIdx())
        for atom in _heavy(product)
    }
    product.xf.tracing._stamp()
    after = {
        atom.GetIdx(): product.xf.tracing.atom_indices(atom.GetIdx())
        for atom in _heavy(product)
    }
    assert before == after


class _CarbonToOxygen(SmirksReactionRule):
    smirks = (
        (
            "[C:1]Cl>>[O:1].[Cl:2]",
            describe(
                site_map=1,
                adds="O",
                removes="C",
                cleaves=True,
                leave_count=1,
                leave_formula={"Cl": 1},
            ),
        ),
    )

    def __init__(self):
        super().__init__(name="CarbonToOxygen")


def test_element_change_keeps_the_root():
    mol = Chem.MolFromSmiles("CCl")
    pieces = []
    for products, _info in _CarbonToOxygen().metabolize(mol):
        pieces.extend(products)
    oxygens = [p for p in pieces if p.GetAtomWithIdx(0).GetAtomicNum() == 8]
    assert oxygens
    oxygen = next(
        atom for atom in oxygens[0].GetAtoms() if atom.GetAtomicNum() == 8
    )
    tracing = oxygens[0].xf.tracing
    assert tracing.atom_root(oxygen.GetIdx()) == 0
    assert tracing.atom_added_by(oxygen.GetIdx()) is None
    assert 0 not in {
        atom.GetIdx()
        for atom in _heavy(oxygens[0])
        if tracing.atom_root(atom.GetIdx()) is None
    }


def _first_product(rule_cls, smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    try:
        products, _info = next(rule_cls().metabolize(mol))
        product = products[0]
    except StopIteration:
        return None
    return product


@pytest.mark.parametrize("rule_cls", _rule_classes(), ids=lambda c: c.__name__)
def test_trace_every_rule(rule_cls):
    try:
        smiles = _reactant_for(rule_cls.__name__)
    except KeyError:
        pytest.skip(f"{rule_cls.__name__}: no example reactant")
    product = _first_product(rule_cls, smiles)
    if product is None:
        pytest.skip(f"{rule_cls.__name__}: no product; chemistry gap, not a trace gap")
    tracing = product.xf.tracing
    assert tracing.active
    assert tracing.depth == 1
    roots = set()
    for atom in _heavy(product):
        idx = atom.GetIdx()
        indices = tracing.atom_indices(idx)
        depths = tracing.atom_depths(idx)
        assert indices is not None and depths is not None
        assert len(indices) == len(depths)
        assert indices[-1] == idx
        root = tracing.atom_root(idx)
        if root is None:
            assert tracing.atom_added_by(idx) is not None
            assert 0 not in depths
        else:
            assert tracing.atom_added_by(idx) is None
            assert depths[0] == 0
            roots.add(root)
    assert tracing.removed_roots().isdisjoint(roots)


def test_two_hydroxylations_keep_original_carbons():
    ethane = Chem.MolFromSmiles("CC")
    ethanol = next(
        product
        for products, _info in Hydroxylation().metabolize(ethane)
        for product in products
        if product.xf.csmi == "CCO"
    )
    glycol = next(
        product
        for products, _info in Hydroxylation().metabolize(ethanol)
        for product in products
        if product.xf.csmi == "OCCO"
    )
    tracing = glycol.xf.tracing
    assert tracing.depth == 2
    carbon_roots = set()
    for atom in _heavy(glycol):
        idx = atom.GetIdx()
        root = tracing.atom_root(idx)
        if atom.GetAtomicNum() == 6:
            assert root in {0, 1}
            assert tracing.atom_depths(idx) == (0, 1, 2)
            assert len(tracing.atom_indices(idx)) == 3
            carbon_roots.add(root)
        else:
            assert root is None
            assert tracing.atom_added_by(idx)[0] == "Hydroxylation"
    assert carbon_roots == {0, 1}


def test_bfs_and_dfs_products_are_traced():
    for enumerator in (bfs, dfs):
        products = list(enumerator("CC", PhaseOne, depth=1))
        assert products
        for product, _info in products:
            assert "." not in product.xf.csmi
            assert product.xf.tracing.active
            assert product.xf.tracing.depth == 1
            present = {
                atom.GetIdx()
                for atom in _heavy(product)
                if product.xf.tracing.atom_indices(atom.GetIdx())
            }
            assert {atom.GetIdx() for atom in _heavy(product)} <= present
