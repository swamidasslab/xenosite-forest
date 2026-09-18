"""Closed mols reject structural edits; edit_mol clears resonance after mutation."""

from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolfiles import MolFromSmiles

import pytest

from xenosite.forest.base import _resonance_cache, copy_mol
from xenosite.forest.edit_guard import MolClosedError, edit_mol


def test_closed_mol_rejects_bond_edit():
    mol = MolFromSmiles("C=C")
    with pytest.raises(MolClosedError):
        mol.GetBondWithIdx(0).SetBondType(BondType.SINGLE)


def test_closed_mol_allows_props():
    mol = MolFromSmiles("C")
    mol.GetAtomWithIdx(0).SetProp("note", "ok")
    mol.GetAtomWithIdx(0).SetAtomMapNum(3)
    assert mol.GetAtomWithIdx(0).GetProp("note") == "ok"
    assert mol.GetAtomWithIdx(0).GetAtomMapNum() == 3


def test_edit_mol_allows_mutation_and_clears_resonance():
    mol = MolFromSmiles("c1ccccc1")
    cache = _resonance_cache(mol)
    assert mol._forest["resonance"] is cache
    with edit_mol(mol):
        mol.GetBondWithIdx(0).SetBondType(BondType.SINGLE)
    assert "resonance" not in mol._forest


def test_edit_mol_without_mutation_keeps_resonance():
    mol = MolFromSmiles("C")
    cache = _resonance_cache(mol)
    with edit_mol(mol):
        mol.SetProp("still", "closed-chem")
    assert mol._forest["resonance"] is cache


def test_closed_mol_rejects_kekulize():
    from rdkit.Chem import Kekulize

    mol = MolFromSmiles("c1ccccc1")
    with pytest.raises(MolClosedError):
        Kekulize(mol)


def test_copy_while_open_does_not_stay_open():
    mol = MolFromSmiles("C=C")
    with edit_mol(mol):
        twin = copy_mol(mol)
        with pytest.raises(MolClosedError):
            twin.GetBondWithIdx(0).SetBondType(BondType.SINGLE)
        mol.GetBondWithIdx(0).SetBondType(BondType.SINGLE)
