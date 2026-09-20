"""Focused tests for ``forest_copy`` layering and rule identity."""

from __future__ import annotations

from types import MappingProxyType

from rdkit import Chem

from xenosite.refactor_poc.forest_copy import forest_copy, start_labels_of
from xenosite.refactor_poc.rdkitutil import copy_mol, ensure_forest
from xenosite.refactor_poc.rules import Hydroxylation


def test_forest_copy_shares_rules_by_id():
    mol = Chem.MolFromSmiles("CCO")
    product, _info = next(Hydroxylation().metabolize(mol))
    forest = product._forest
    addition = next(iter(forest["atom_trace"]["additions"].values()))
    rule = addition["rules"][0]
    copied = forest_copy(forest, same_structure=False)
    copied_addition = next(iter(copied["atom_trace"]["additions"].values()))
    assert copied_addition["rules"][0] is rule
    assert id(copied_addition["rules"][0]) == id(rule)
    # Mutable containers are new.
    assert copied["atom_trace"] is not forest["atom_trace"]
    assert copied["atom_trace"]["additions"] is not forest["atom_trace"]["additions"]


def test_forest_copy_shallow_shares_immutable():
    mol = Chem.MolFromSmiles("CCO")
    mol.GetAtomWithIdx(0).SetProp("atomLabel", "alpha")
    stamped = mol.xf.tracing._stamp()
    forest = stamped._forest
    labels = start_labels_of(forest)
    assert labels is not None
    assert dict(labels) == {0: "alpha"}
    assert isinstance(forest["immutable"], MappingProxyType)

    same = forest_copy(forest, same_structure=True)
    diff = forest_copy(forest, same_structure=False)
    # New immutable shell; nested start_labels shared by identity.
    assert same["immutable"] is not forest["immutable"]
    assert diff["immutable"] is not forest["immutable"]
    assert same["immutable"]["start_labels"] is labels
    assert diff["immutable"]["start_labels"] is labels


def test_forest_copy_cache_dropped_vs_kept():
    parent = ensure_forest(Chem.MolFromSmiles("CCO"))
    _ = parent.xf.csmi
    forest = parent._forest
    cache = forest["cache"]
    assert cache.get("csmi") == "CCO"

    kept = forest_copy(forest, same_structure=True)
    dropped = forest_copy(forest, same_structure=False)
    assert kept["cache"] is cache
    assert "cache" not in dropped

    via_copy_mol = copy_mol(parent)._forest
    assert via_copy_mol["cache"] is cache


def test_forest_copy_shares_pattern_info_identity():
    mol = Chem.MolFromSmiles("CCO")
    product, info = next(Hydroxylation().metabolize(mol))
    pattern = info["pattern"]
    addition = next(iter(product._forest["atom_trace"]["additions"].values()))
    assert addition["pattern"] is pattern
    copied = forest_copy(product._forest, same_structure=True)
    copied_addition = next(iter(copied["atom_trace"]["additions"].values()))
    assert copied_addition["pattern"] is pattern
