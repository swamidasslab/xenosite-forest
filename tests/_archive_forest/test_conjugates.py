"""Conjugation star adducts, thiol filtering, and CXSMILES labels."""

from __future__ import annotations

import re

import pytest
from rdkit import Chem

from xenosite._archive_forest import load_ruleset
from xenosite._archive_forest.rules import (
    STAR_ONLY_LABELS,
    Glucuronidation,
    Glutathionation,
)
from xenosite._archive_forest.utils import mol_to_cxsmiles

PHENOL = "c1ccccc1O"
BENZOIC_ACID = "c1ccccc1C(=O)O"
STYRENE_OXIDE = "c1ccccc1C1OC1"
BENZYL_CHLORIDE = "ClCc1ccccc1"
BENZYL_BROMIDE = "BrCc1ccccc1"
BENZYL_IODIDE = "ICc1ccccc1"
TERMINAL_ALKENE = "C=CC"
ALIPHATIC_THIOL = "CCS"
THIOPHENOL = "Sc1ccccc1"
CYSTEINE_LIKE = "SC[C@H](N)C(=O)O"
BENZOQUINONE = "O=C1C=CC(=O)C=C1"
NAPQI = "CC(=O)N=C1C=CC(=O)C=C1"
MENADIONE = "CC1=CC(=O)c2ccccc2C1=O"
BETA_ENONE = "CC=CC(=O)CO"
CINNAMALDEHYDE = "O=CC=Cc1ccccc1"
FORMALDEHYDE = "C=O"
ACETALDEHYDE = "CC=O"
BENZALDEHYDE = "O=Cc1ccccc1"
METHYLGLYOXAL = "CC(=O)C=O"
FURFURAL = "O=Cc1ccco1"
N_PHENYLAZIRIDINE = "C1CN1c1ccccc1"
METHYL_MESYLATE = "COS(=O)(=O)C"
BENZYL_MESYLATE = "CS(=O)(=O)OCc1ccccc1"
PHENYL_ISOCYANATE = "O=C=Nc1ccccc1"
PHENYL_ISOTHIOCYANATE = "S=C=Nc1ccccc1"
# Lookalikes that must not match the dedicated SMARTS (other rules may still fire).
AZETIDINE_PHENYL = "C1CCN1c1ccccc1"
PYRROLIDINE = "C1CCNC1"
SULFONAMIDE = "CS(=O)(=O)Nc1ccccc1"
SULFONE = "CS(=O)(=O)C"
SULFONIC_ACID = "CS(=O)(=O)O"
BENZONITRILE = "N#Cc1ccccc1"
ACETANILIDE = "CC(=O)Nc1ccccc1"
CARBODIIMIDE = "N=C=Nc1ccccc1"

# One probe substrate per conjugation SMARTS (UGT / GSH / Protein / DNA / Cyanide).
# Glutathionation indices: 0 epoxide, 1 halide, 2 thiol, 3 alkene, 4 Michael,
# 5 aldehyde, 6 aziridine, 7 sulfonate, 8 isocyanate.
# NoThiol drops thiol → aldehyde=4, aziridine=5, sulfonate=6, isocyanate=7.
_CONJUGATION_SMARTS_PROBES = [
    # (id, rule factory, smiles, smarts rxn index)
    ("ugt_acid", lambda: Glucuronidation(star_label="GlcA"), BENZOIC_ACID, 0),
    ("ugt_phenol", lambda: Glucuronidation(star_label="GlcA"), PHENOL, 1),
    ("gsh_epoxide", lambda: Glutathionation(star_label="GSH"), STYRENE_OXIDE, 0),
    ("gsh_chloride", lambda: Glutathionation(star_label="GSH"), BENZYL_CHLORIDE, 1),
    ("gsh_bromide", lambda: Glutathionation(star_label="GSH"), BENZYL_BROMIDE, 1),
    ("gsh_thiol", lambda: Glutathionation(star_label="GSH"), ALIPHATIC_THIOL, 2),
    ("gsh_alkene", lambda: Glutathionation(star_label="GSH"), TERMINAL_ALKENE, 3),
    ("gsh_michael", lambda: Glutathionation(star_label="GSH"), BENZOQUINONE, 4),
    ("protein_epoxide", lambda: Glutathionation(star_label="Protein"), STYRENE_OXIDE, 0),
    ("protein_chloride", lambda: Glutathionation(star_label="Protein"), BENZYL_CHLORIDE, 1),
    ("protein_thiol", lambda: Glutathionation(star_label="Protein"), ALIPHATIC_THIOL, 2),
    ("protein_alkene", lambda: Glutathionation(star_label="Protein"), TERMINAL_ALKENE, 3),
    ("protein_michael", lambda: Glutathionation(star_label="Protein"), BENZOQUINONE, 4),
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
        "dna_michael",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        BENZOQUINONE,
        3,
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
    (
        "cyanide_michael",
        lambda: Glutathionation(include_thiol=False, star_label="Cyanide"),
        BENZOQUINONE,
        3,
    ),
    ("gsh_aldehyde", lambda: Glutathionation(star_label="GSH"), ACETALDEHYDE, 5),
    (
        "dna_aldehyde",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        ACETALDEHYDE,
        4,
    ),
    ("gsh_aziridine", lambda: Glutathionation(star_label="GSH"), N_PHENYLAZIRIDINE, 6),
    ("gsh_sulfonate", lambda: Glutathionation(star_label="GSH"), METHYL_MESYLATE, 7),
    ("gsh_isocyanate", lambda: Glutathionation(star_label="GSH"), PHENYL_ISOCYANATE, 8),
    (
        "dna_aziridine",
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        N_PHENYLAZIRIDINE,
        5,
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


@pytest.mark.parametrize(
    "smi, expected_soms, product_substr",
    [
        # Canonical SMILES; Michael β-carbons match high reactivity.gsh atoms.
        (BENZOQUINONE, {2, 3, 6, 7}, "*C1C=C(O)C=CC1=O"),
        (NAPQI, {5, 6, 9, 10}, None),
        (MENADIONE, {2}, "*C1C(=O)c2ccccc2C(O)=C1C"),
        (BETA_ENONE, {1}, None),
        (CINNAMALDEHYDE, {3}, None),
    ],
)
def test_glutathionation_michael_acceptors(smi, expected_soms, product_substr):
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 4)
    assert hits, f"no Michael products for {can}"
    soms = {next(iter(site)) for site, _ in hits}
    assert expected_soms <= soms, f"{can}: expected {expected_soms} ⊆ {soms}"
    for site, products in hits:
        assert len(site) == 1
        assert "*" in _smiles(products[0])
    if product_substr:
        assert any(product_substr == _smiles(p[0]) for _, p in hits)


@pytest.mark.parametrize("smi", [BENZYL_BROMIDE, BENZYL_IODIDE])
def test_glutathionation_benzyl_bromide_iodide(smi):
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 1)
    assert hits
    for site, products in hits:
        assert site == frozenset({1})
        assert _smiles(products[0]) == "*Cc1ccccc1"


@pytest.mark.parametrize(
    "rule_factory",
    [
        lambda: Glutathionation(star_label="GSH"),
        lambda: Glutathionation(star_label="Protein"),
        lambda: Glutathionation(include_thiol=False, star_label="DNA"),
        lambda: Glutathionation(include_thiol=False, star_label="Cyanide"),
    ],
)
def test_reactivity_heads_enumerate_benzoquinone(rule_factory):
    assert _has_products(rule_factory(), BENZOQUINONE)


@pytest.mark.parametrize(
    "smi, carbonyl_idx, product_smi",
    [
        # Canonical SMILES; SOM is the aldehyde carbon (high reactivity.gsh).
        (FORMALDEHYDE, 0, "*CO"),
        (ACETALDEHYDE, 1, "*C(C)O"),
        (BENZALDEHYDE, 1, "*C(O)c1ccccc1"),
        (METHYLGLYOXAL, 3, "*C(O)C(C)=O"),
        (FURFURAL, 1, "*C(O)c1ccco1"),
    ],
)
def test_glutathionation_aldehyde_thiohemiacetal(smi, carbonyl_idx, product_smi):
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    mol = Chem.MolFromSmiles(can)
    assert mol.GetAtomWithIdx(carbonyl_idx).GetSymbol() == "C"
    assert mol.GetAtomWithIdx(carbonyl_idx).GetTotalNumHs() >= 1
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 5)
    assert hits, f"no aldehyde products for {can}"
    soms = {next(iter(site)) for site, _ in hits}
    assert carbonyl_idx in soms, f"{can}: expected carbonyl {carbonyl_idx} in {soms}"
    for site, products in hits:
        assert len(site) == 1
        assert site == frozenset({carbonyl_idx})
        assert _smiles(products[0]) == product_smi


def test_benzaldehyde_includes_carbonyl_thiohemiacetal():
    """Aromatic aldehydes must include carbonyl C adduct (ring hits OK; low specificity)."""
    can = Chem.MolToSmiles(Chem.MolFromSmiles(BENZALDEHYDE))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 5)
    assert hits
    soms = {next(iter(site)) for site, _ in hits}
    assert 1 in soms
    assert any(_smiles(p[0]) == "*C(O)c1ccccc1" for _, p in hits)


@pytest.mark.parametrize(
    "smi, rxn_idx, expected_soms, product_smi",
    [
        (N_PHENYLAZIRIDINE, 6, {5, 6}, "*CCNc1ccccc1"),
        (METHYL_MESYLATE, 7, {0}, "*C"),
        (BENZYL_MESYLATE, 7, {5}, "*Cc1ccccc1"),
        (PHENYL_ISOCYANATE, 8, {2}, "*C(=O)Nc1ccccc1"),
        (PHENYL_ISOTHIOCYANATE, 8, {2}, "*C(=S)Nc1ccccc1"),
    ],
)
def test_glutathionation_aziridine_sulfonate_isocyanate(
    smi, rxn_idx, expected_soms, product_smi
):
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, rxn_idx)
    assert hits, f"no products for rxn {rxn_idx} on {can}"
    soms = {next(iter(site)) for site, _ in hits}
    assert expected_soms <= soms
    for site, products in hits:
        assert len(site) == 1
        assert _smiles(products[0]) == product_smi


@pytest.mark.parametrize(
    "smi, rxn_idx, reason",
    [
        # Aziridine lookalikes (4+/acyclic N, pyridine, epoxide).
        (AZETIDINE_PHENYL, 6, "azetidine"),
        (PYRROLIDINE, 6, "pyrrolidine"),
        ("CCNc1ccccc1", 6, "aniline"),
        ("c1ccncc1", 6, "pyridine"),
        (STYRENE_OXIDE, 6, "epoxide"),
        # Sulfonate lookalikes (no C–O–SO2 alkylator).
        (SULFONAMIDE, 7, "sulfonamide"),
        (SULFONE, 7, "sulfone"),
        (SULFONIC_ACID, 7, "sulfonic acid"),
        ("CS(C)=O", 7, "sulfoxide"),
        # Isocyanate lookalikes.
        (BENZONITRILE, 8, "nitrile"),
        (ACETANILIDE, 8, "amide"),
        ("O=C=O", 8, "CO2"),
        ("CC#N", 8, "acetonitrile"),
        (CARBODIIMIDE, 8, "carbodiimide"),
    ],
)
def test_glutathionation_lookalikes_miss_dedicated_smarts(smi, rxn_idx, reason):
    """Similar motifs must not accidentally match aziridine/sulfonate/isocyanate SMARTS."""
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, rxn_idx)
    assert hits == [], f"{reason} ({can}) unexpectedly matched rxn {rxn_idx}: {hits}"

