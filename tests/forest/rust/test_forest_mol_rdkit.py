"""Lazy RDKit bridge on product ForestMol (optional extra)."""

from __future__ import annotations

import pytest

from xenosite.forest import available

pytest.importorskip("rdkit")

assert available()


def test_forest_mol_rdkit_roundtrip_smiles():
    from rdkit import Chem

    from xenosite.forest import ForestMol

    rd = Chem.MolFromSmiles("CCO")
    mol = ForestMol(rd)
    assert mol.csmi
    back = mol.to_rdkit()
    assert Chem.MolToSmiles(back) == Chem.MolToSmiles(rd)
