"""Conjugation star adducts, thiol filtering, and CXSMILES labels."""

from __future__ import annotations

import re

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
BENZOIC_ACID = "c1ccccc1C(=O)O"
STYRENE_OXIDE = "c1ccccc1C1OC1"
BENZYL_CHLORIDE = "ClCc1ccccc1"
TERMINAL_ALKENE = "C=CC"
ALIPHATIC_THIOL = "CCS"
THIOPHENOL = "Sc1ccccc1"
CYSTEINE_LIKE = "SC[C@H](N)C(=O)O"

# One probe substrate per conjugation SMARTS (UGT / GSH / Protein / DNA / Cyanide).
# GlutathionationNoThiol drops the thiol SMARTS, so DNA/Cyanide alkene is rxn index 2.
_CONJUGATION_SMARTS_PROBES = [
    # (id, rule factory, smiles, smarts rxn index)
    ("ugt_acid", lambda: Glucuronidation(star_label="GlcA"), BENZOIC_ACID, 0),
    ("ugt_phenol", lambda: Glucuronidation(star_label="GlcA"), PHENOL, 1),
    ("gsh_epoxide", lambda: Glutathionation(star_label="GSH"), STYRENE_OXIDE, 0),
    ("gsh_chloride", lambda: Glutathionation(star_label="GSH"), BENZYL_CHLORIDE, 1),
    ("gsh_thiol", lambda: Glutathionation(star_label="GSH"), ALIPHATIC_THIOL, 2),
    ("gsh_alkene", lambda: Glutathionation(star_label="GSH"), TERMINAL_ALKENE, 3),
    ("protein_epoxide", lambda: Glutathionation(star_label="Protein"), STYRENE_OXIDE, 0),
    ("protein_chloride", lambda: Glutathionation(star_label="Protein"), BENZYL_CHLORIDE, 1),
    ("protein_thiol", lambda: Glutathionation(star_label="Protein"), ALIPHATIC_THIOL, 2),
    ("protein_alkene", lambda: Glutathionation(star_label="Protein"), TERMINAL_ALKENE, 3),
    (
        "dna_epoxide",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        STYRENE_OXIDE,
        0,
    ),
    (
        "dna_chloride",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        BENZYL_CHLORIDE,
        1,
    ),
    (
        "dna_alkene",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        TERMINAL_ALKENE,
        2,
    ),
    (
        "cyanide_epoxide",
        lambda: Glutathionation(include_thiol=False, star_label="Cyanide"),
        STYRENE_OXIDE,
        0,
    ),
    (
        "cyanide_chloride",
        lambda: Glutathionation(include_thiol=False, star_label="Cyanide"),
        BENZYL_CHLORIDE,
        1,
    ),
    (
        "cyanide_alkene",
        lambda: Glutathionation(include_thiol=False, star_label="Cyanide"),
        TERMINAL_ALKENE,
        2,
    ),
]


def _smiles(mol):
    return Chem.MolToSmiles(mol)


def _rxn_index(site_name: str) -> int:
    m = re.search(r"Rxn(\d+)$", site_name)
    if not m:
        raise AssertionError(f"unexpected site name {site_name!r}")
    return int(m.group(1))


def _sites_for_rxn(rule, smi, rxn_idx):
    mol = Chem.MolFromSmiles(smi)
    out = []
    for (name, site), products in rule.metabolites(mol):
        if _rxn_index(name) == rxn_idx:
            out.append((site, products))
    return out


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


@pytest.mark.parametrize(
    "probe_id, rule_factory, smi, rxn_idx",
    _CONJUGATION_SMARTS_PROBES,
    ids=[p[0] for p in _CONJUGATION_SMARTS_PROBES],
)
def test_conjugation_smarts_som_is_single_atom(probe_id, rule_factory, smi, rxn_idx):
    """Every conjugation SMARTS must report a single-atom site of metabolism."""
    hits = _sites_for_rxn(rule_factory(), smi, rxn_idx)
    assert hits, f"{probe_id}: no products for SMARTS rxn {rxn_idx} on {smi}"
    for site, _products in hits:
        assert len(site) == 1, (
            f"{probe_id}: expected single-atom SOM, got {sorted(site)} "
            f"(len={len(site)}) for {smi}"
        )


@pytest.mark.parametrize(
    "probe_id, smi, rxn_idx, expected_o",
    [
        # Canonical SMILES atom order: phenolic / acid OH oxygen.
        ("ugt_phenol", "Oc1ccccc1", 1, 0),
        ("ugt_acid", "O=C(O)c1ccccc1", 0, 2),
        ("ugt_apap", "CC(=O)Nc1ccc(O)cc1", 1, 8),
        ("ugt_ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 0, 12),
    ],
)
def test_ugt_som_is_oxygen(probe_id, smi, rxn_idx, expected_o):
    """UGT SOM must be the oxygen that receives GlcA (matches high ugt scores)."""
    mol = Chem.MolFromSmiles(smi)
    assert mol.GetAtomWithIdx(expected_o).GetAtomicNum() == 8
    hits = _sites_for_rxn(Glucuronidation(star_label="GlcA"), smi, rxn_idx)
    assert hits, f"{probe_id}: no UGT products for rxn {rxn_idx}"
    for site, _ in hits:
        assert site == frozenset({expected_o}), (
            f"{probe_id}: expected SOM {{{expected_o}}}, got {sorted(site)}"
        )
