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


def test_forest_mol_from_rdkit_preserves_atom_index_frame():
    """RDKit ingest must not rewrite through canonical SMILES.

    Anisole ``c1ccccc1OC`` is C0…C5, O6, C7 in RDKit. Canonical
    ``MolToSmiles`` yields ``COc1ccccc1`` (methyl at 0). Demethylation sites
    must stay on RDKit index 7 so predict SOM scores attach correctly.
    """
    from rdkit import Chem

    from xenosite.forest import ForestMol, dealkylation

    rd = Chem.MolFromSmiles("c1ccccc1OC")
    assert [a.GetSymbol() for a in rd.GetAtoms()] == list("CCCCCCOC")
    assert Chem.MolToSmiles(rd) == "COc1ccccc1"  # canonical reorder exists

    mol = ForestMol(rd)
    methyl_sites = {
        em.site
        for em in dealkylation().metabolize(mol)
        if "methyl" in (em.pattern_name or "")
    }
    assert 7 in methyl_sites
    assert 0 not in methyl_sites


def test_require_identity_smiles_atom_order_rejects_canonical_reorder():
    from rdkit import Chem

    from xenosite.forest.mol import _require_identity_smiles_atom_order

    rd = Chem.MolFromSmiles("c1ccccc1OC")
    Chem.MolToSmiles(rd, canonical=True)
    assert Chem.MolToSmiles(rd) == "COc1ccccc1"
    with pytest.raises(ValueError, match="reordered heavy atoms"):
        _require_identity_smiles_atom_order(rd)
