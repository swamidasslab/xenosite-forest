"""Conjugation star adducts, thiol filtering, and CXSMILES labels."""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest import load_ruleset
from xenosite.forest.rules import (
    STAR_ONLY_LABELS,
    Glucuronidation,
    Glutathionation,
)
from xenosite.forest.utils import mol_to_cxsmiles

PHENOL = "c1ccccc1O"
STYRENE_OXIDE = "c1ccccc1C1OC1"
BENZYL_CHLORIDE = "ClCc1ccccc1"
TERMINAL_ALKENE = "C=CC"
ALIPHATIC_THIOL = "CCS"
THIOPHENOL = "Sc1ccccc1"
CYSTEINE_LIKE = "SC[C@H](N)C(=O)O"


def _smiles(mol):
    return Chem.MolToSmiles(mol)


def _has_products(rule, smi):
    mol = Chem.MolFromSmiles(smi)
    return any(True for _ in rule.metabolites(mol))


def test_glucuronidation_defaults_to_bare_star():
    mol = Chem.MolFromSmiles(PHENOL)
    _, products = next(Glucuronidation().metabolites(mol))
    smi = _smiles(products[0])
    assert "*" in smi
    assert "O=C(O)C1OC" not in smi
    assert smi == "*Oc1ccccc1"


def test_glucuronidation_full_structure_opt_in():
    mol = Chem.MolFromSmiles(PHENOL)
    _, products = next(Glucuronidation(as_star=False).metabolites(mol))
    smi = _smiles(products[0])
    assert "*" not in smi
    assert "O=C(O)" in smi


def test_glucuronidation_star_label_glca():
    mol = Chem.MolFromSmiles(PHENOL)
    _, products = next(Glucuronidation(star_label="GlcA").metabolites(mol))
    cx = mol_to_cxsmiles(products[0])
    assert cx is not None
    assert cx.split()[0] == "*Oc1ccccc1"
    assert "GlcA" in cx
    dummy = next(a for a in products[0].GetAtoms() if a.GetAtomicNum() == 0)
    assert dummy.GetProp("atomLabel") == "GlcA"


@pytest.mark.parametrize("smi", [ALIPHATIC_THIOL, THIOPHENOL, CYSTEINE_LIKE])
def test_thiols_match_gsh_not_no_thiol(smi):
    assert _has_products(Glutathionation(), smi)
    assert not _has_products(Glutathionation(include_thiol=False), smi)
    assert not _has_products(
        Glutathionation(include_thiol=False, star_label="DNA"), smi
    )


@pytest.mark.parametrize("smi", [STYRENE_OXIDE, BENZYL_CHLORIDE, TERMINAL_ALKENE])
def test_electrophiles_match_both_thiol_modes(smi):
    assert _has_products(Glutathionation(), smi)
    assert _has_products(Glutathionation(include_thiol=False), smi)


def test_glutathionation_defaults_to_star_not_peptide():
    mol = Chem.MolFromSmiles(STYRENE_OXIDE)
    _, products = next(Glutathionation().metabolites(mol))
    smi = _smiles(products[0])
    assert "*" in smi
    assert "NC(CCC(=O)N" not in smi


@pytest.mark.parametrize(
    "label",
    ["GSH", "Protein", "DNA", "Cyanide"],
)
def test_glutathionation_star_labels(label):
    include_thiol = label not in ("DNA", "Cyanide")
    rule = Glutathionation(include_thiol=include_thiol, star_label=label)
    mol = Chem.MolFromSmiles(STYRENE_OXIDE)
    _, products = next(rule.metabolites(mol))
    cx = mol_to_cxsmiles(products[0])
    assert cx is not None
    assert "*" in cx.split()[0]
    assert label in cx
    dummy = next(a for a in products[0].GetAtoms() if a.GetAtomicNum() == 0)
    assert dummy.GetProp("atomLabel") == label


@pytest.mark.parametrize("label", sorted(STAR_ONLY_LABELS))
def test_star_only_labels_reject_full_structure(label):
    with pytest.raises(ValueError, match="full conjugate"):
        Glutathionation(star_label=label, as_star=False)


def test_load_glutathionation_no_thiol_ruleset():
    rs = load_ruleset("GlutathionationNoThiol")
    assert rs.name == "GlutathionationNoThiol"
    assert _has_products(rs.rules[0], STYRENE_OXIDE)
    assert not _has_products(rs.rules[0], ALIPHATIC_THIOL)
