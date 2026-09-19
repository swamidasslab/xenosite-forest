"""Cache identity for the RDKit door. Callers read NamedTuple attributes."""

import copy

from rdkit import Chem

from xenosite.refactor_poc.rdkitutil import (
    copy_mol,
    get_csmi,
    get_forest,
    mcs_matches,
    molecule_formula,
    resonance_bond_maps,
    sanitized_fragments,
    smarts_matches,
    split_fragments,
    topol_equiv,
)


def test_cached_answers_are_the_same_object():
    mol = Chem.MolFromSmiles("CCC")
    assert topol_equiv(mol) is topol_equiv(mol)
    assert get_csmi(mol) is get_csmi(mol)
    smarts = "[C:1][C:2]"
    assert smarts_matches(mol, smarts) is smarts_matches(mol, smarts)


def test_mcs_matches_cache_on_the_reactant():
    reactant = Chem.MolFromSmiles("c1ccccc1")
    target = Chem.MolFromSmiles("c1ccccc1")
    first = mcs_matches(reactant, target)
    second = mcs_matches(reactant, target)
    assert first is second
    assert len(first.embeddings) > 1
    held = get_forest(reactant)["structure"]["mcs_matches"][get_csmi(target)]
    assert held is first
    assert isinstance(held.embeddings, tuple)


def test_new_structure_does_not_reuse_parent_caches():
    parent = Chem.MolFromSmiles("c1ccccc1O")
    parent_csmi = get_csmi(parent)
    parent_maps = resonance_bond_maps(parent)

    product = Chem.Mol(parent)
    built = get_forest(product, new_structure=True)
    assert "structure" not in built
    structure = get_forest(product)["structure"]
    assert "csmi" not in structure
    assert "resonance_bonds" not in structure
    assert structure is not get_forest(parent)["structure"]
    assert get_csmi(product) == parent_csmi
    assert resonance_bond_maps(product) is not parent_maps

    child = copy_mol(parent)
    stripped = get_forest(child, new_structure=True)
    assert "csmi" not in stripped.get("structure", {})
    assert "resonance_bonds" not in stripped.get("structure", {})


def test_forest_stays_a_dict():
    mol = Chem.MolFromSmiles("CC")
    forest = get_forest(mol)
    forest["structure"]["csmi"] = "CC"
    forest["atom_trace"] = {"depth": 0}
    forest["not_a_schema_key"] = 1
    assert forest["structure"]["csmi"] == "CC"
    assert forest["atom_trace"]["depth"] == 0
    assert forest["not_a_schema_key"] == 1
    assert isinstance(forest, dict)

    formula = molecule_formula(mol)
    assert isinstance(formula, dict)
    assert formula["counts"]["C"] == 2
    assert formula["charge"] == 0
    assert formula is molecule_formula(mol)

    cloned = copy.deepcopy(forest)
    assert cloned["structure"]["formula"]["counts"]["C"] == 2
    assert cloned["atom_trace"]["depth"] == 0


def test_fragment_split_uses_pieces():
    mol = Chem.MolFromSmiles("CCO")
    split = split_fragments(mol)
    assert split.pieces == (mol,)
    fragments = sanitized_fragments(mol)
    assert len(fragments.pieces) == 1
