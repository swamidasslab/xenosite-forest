"""Partial kekulé parents. Helpers take a dict. The rule stores it."""

import copy

from rdkit import Chem

from xenosite.refactor_poc.rdkitutil import (
    ensure_forest,
    ensure_kekule_parents,
    get_forest,
    parents_for_ends,
)
from xenosite.refactor_poc.rules import Epoxidation, OxygenReduction

_ANTHRACENE = "c1ccc2cc3ccccc3cc2c1"
_POLYPHENYL = "c1ccc(-c2ccc(-c3ccc(-c4ccc(-c5ccc(-c6ccccc6)cc5)cc4)cc3)cc2)cc1"
_CC = Chem.MolFromSmarts("[#6:1]=,:[#6:2]")


def _canon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _products(rule, smiles):
    found = set()
    for mol, _info in rule.metabolize(Chem.MolFromSmiles(smiles)):
        text = Chem.MolToSmiles(mol, isomericSmiles=False)
        parsed = Chem.MolFromSmiles(text)
        assert parsed is not None, text
        found.add(Chem.MolToSmiles(parsed, isomericSmiles=False))
    return found


def _fill(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    cache = {}
    for left, right in mol.GetSubstructMatches(_CC):
        ensure_kekule_parents(mol, left, right, cache)
    return mol, cache


def test_helpers_do_not_touch_forest():
    mol, cache = _fill("c1ccccc1")
    assert getattr(mol, "_forest", None) is None
    assert len(cache["parents"]) == 2


def test_anthracene_is_four_parents_not_sixteen():
    mol, cache = _fill(_ANTHRACENE)
    assert getattr(mol, "_forest", None) is None
    assert len(cache["parents"]) == 4
    assert all(
        not atom.GetIsAromatic()
        for parent in cache["parents"]
        for atom in parent.GetAtoms()
    )


def test_polyphenyl_is_twelve_parents_not_sixty_four():
    mol, cache = _fill(_POLYPHENYL)
    assert getattr(mol, "_forest", None) is None
    assert len(cache["parents"]) == 12
    for parent in cache["parents"]:
        aromatic = sum(atom.GetIsAromatic() for atom in parent.GetAtoms())
        assert aromatic == 30


def test_rule_stores_the_dict_and_deepcopy_keeps_the_mols():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert mol is not None
    list(Epoxidation().metabolites(mol))
    cache = get_forest(ensure_forest(mol))["structure"]["kekule_parents"]
    assert len(cache["parents"]) == 2
    copied = copy.deepcopy(get_forest(ensure_forest(mol)))
    assert len(copied["structure"]["kekule_parents"]["parents"]) == 2


def test_epoxidation_keeps_the_measured_products():
    anthracene = _products(Epoxidation(), _ANTHRACENE)
    assert _canon("C1=c2cc3ccccc3cc2=CC2OC12") in anthracene
    assert _canon("C1=CC2OC2c2cc3ccccc3cc21") in anthracene
    poly = _products(Epoxidation(), _POLYPHENYL)
    kept = _canon(
        "C1=CC=C(C2=CC=C(C3=CC=C(C4=CC=C(C5=CC6OC6(C6=CC=CC=C6)C=C5)C=C4)C=C3)C=C2)C=C1"
    )
    assert kept in poly
    assert _products(OxygenReduction(), "C([O-])=O") == {"[O-]CO"}


def test_pair_ends_in_different_systems_are_not_a_product():
    mol, cache = _fill(_POLYPHENYL)
    systems = list(cache["systems"])
    assert len(systems) == 6
    start = next(iter(systems[0]))
    end = next(iter(systems[1]))
    ends = parents_for_ends(mol, start, end, cache)
    assert ends.same_system is False
    assert len(ends.parents) == 4
    assert len(cache["parents"]) == 12
