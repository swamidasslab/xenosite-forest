"""Parity: forest AtomTracker vs thin POC AtomTracker (xf proxy)."""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.base import AtomTracker as ForestAtomTracker
from xenosite.forest.rules import Hydroxylation as ForestHydroxylation
from xenosite.refactor_poc.atom_tracker import AtomTracker as PocAtomTracker
from xenosite.refactor_poc.rules import Hydroxylation as PocHydroxylation


def test_topol_equiv_parity_on_same_smiles():
    smi = "COc1ccc(O)cc1"
    forest_mol = Chem.MolFromSmiles(smi)
    poc_mol = Chem.MolFromSmiles(smi)
    assert forest_mol is not None and poc_mol is not None
    assert ForestAtomTracker.topol_equiv(forest_mol) == PocAtomTracker.topol_equiv(
        poc_mol
    )


def test_initialize_tags_depths_parity():
    forest_mol = Chem.MolFromSmiles("CCO")
    poc_mol = Chem.MolFromSmiles("CCO")
    assert forest_mol is not None and poc_mol is not None
    ForestAtomTracker().initialize_tags(forest_mol)
    PocAtomTracker().initialize_tags(poc_mol)
    assert ForestAtomTracker.depths(forest_mol) == PocAtomTracker.depths(poc_mol) == [0]
    f_tags = ForestAtomTracker.tags(forest_mol)
    p_tags = PocAtomTracker.tags(poc_mol)
    assert len(f_tags) == len(p_tags)
    # Same number of heavy atoms at depth 0.
    assert {tuple(v["depth"]) for v in f_tags.values()} == {(0,)}
    assert {tuple(v["depth"]) for v in p_tags.values()} == {(0,)}


def test_hydroxylation_product_depths_and_added_atom_parity():
    """Side-by-side: both trackers see depth 0+1 and an atom without depth 0."""

    forest_parent = Chem.MolFromSmiles("CCO")
    poc_parent = Chem.MolFromSmiles("CCO")
    assert forest_parent is not None and poc_parent is not None

    forest_rows = list(ForestHydroxylation().metabolize(forest_parent))
    poc_product, _info = next(PocHydroxylation().metabolize(poc_parent))

    # Forest metabolize yields ((rule, site), [mols], ...).
    forest_product = None
    for item in forest_rows:
        if not isinstance(item, tuple) or len(item) < 2:
            continue
        mols = item[1]
        if isinstance(mols, list) and mols:
            forest_product = mols[0]
            break
    assert forest_product is not None

    f_depths = ForestAtomTracker.depths(forest_product)
    p_depths = PocAtomTracker.depths(poc_product)
    assert f_depths == p_depths == [0, 1]

    f_tags = ForestAtomTracker.tags(forest_product)
    p_tags = PocAtomTracker.tags(poc_product)
    f_new = [v for v in f_tags.values() if 0 not in v["depth"]]
    p_new = [v for v in p_tags.values() if 0 not in v["depth"]]
    assert f_new and p_new  # added oxygen has no depth-0 root


def test_poc_tags_strict_untraced_raises():
    mol = Chem.MolFromSmiles("CC")
    assert mol is not None
    with pytest.raises(KeyError):
        PocAtomTracker.tags(mol, strict=True)
    assert PocAtomTracker.tags(mol, strict=False) == {}


def test_site_to_topol_site_parity():
    mol = Chem.MolFromSmiles("CCC")
    assert mol is not None
    te = PocAtomTracker.topol_equiv(mol)
    site = ("Hydroxylation_SmartsReactionRuleRxn0", (0, 2))
    assert ForestAtomTracker.site_to_topol_site(
        site, te
    ) == PocAtomTracker.site_to_topol_site(site, te)
