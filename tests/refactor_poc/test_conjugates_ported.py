"""Port of ``tests/test_conjugates.py`` for poc conjugation rules.

Probe SMILES are imported from ``test_conjugates`` (not copied). Forest
``include_thiol`` / ``load_ruleset`` / ``STAR_ONLY_LABELS`` are unfinished on
the poc surface — those cases are skipped. GSH SMARTS indices differ because
the poc expands epoxide and aziridine into three patterns each; helpers remap
forest indices to the matching poc ``rxn_num`` set.
"""

from __future__ import annotations

import pytest

from test_conjugates import (
    ACETALDEHYDE,
    ACETANILIDE,
    AZETIDINE_PHENYL,
    BENZALDEHYDE,
    BENZOIC_ACID,
    BENZOQUINONE,
    BENZYL_BROMIDE,
    BENZYL_CHLORIDE,
    BENZYL_IODIDE,
    BENZYL_MESYLATE,
    BENZONITRILE,
    BETA_ENONE,
    CARBODIIMIDE,
    CINNAMALDEHYDE,
    CYSTEINE_LIKE,
    FORMALDEHYDE,
    FURFURAL,
    MENADIONE,
    METHYL_MESYLATE,
    METHYLGLYOXAL,
    NAPQI,
    N_PHENYLAZIRIDINE,
    PHENOL,
    PHENYL_ISOCYANATE,
    PHENYL_ISOTHIOCYANATE,
    PYRROLIDINE,
    STYRENE_OXIDE,
    SULFONAMIDE,
    SULFONE,
    SULFONIC_ACID,
    TERMINAL_ALKENE,
    THIOPHENOL,
    ALIPHATIC_THIOL,
)
from xenosite.forest.utils import mol_to_cxsmiles
from xenosite.refactor_poc.rdkit_api import Mol, MolFromSmiles, MolToSmiles
from xenosite.refactor_poc.rules import Glucuronidation, Glutathionation

# Forest GSH index → poc rxn_num values (epoxide/aziridine are three patterns).
_GSH_FOREST_TO_POC: dict[int, frozenset[int]] = {
    0: frozenset({0, 1, 2}),  # epoxide
    1: frozenset({3}),  # halide
    2: frozenset({4}),  # thiol
    3: frozenset({5}),  # alkene
    4: frozenset({6}),  # Michael
    5: frozenset({7}),  # aldehyde
    6: frozenset({8, 9, 10}),  # aziridine
    7: frozenset({11}),  # sulfonate
    8: frozenset({12}),  # isocyanate
}

# Forest UGT: acid=0 phenol=1. Poc: phenol=0 acid=1.
_UGT_FOREST_TO_POC: dict[int, frozenset[int]] = {
    0: frozenset({1}),
    1: frozenset({0}),
}


def _smiles(mol: Mol) -> str:
    return MolToSmiles(mol)


def _sites_for_rxn(rule, smi: str, forest_rxn_idx: int, *, family: str):
    mol = MolFromSmiles(smi)
    assert mol is not None, smi
    if family == "gsh":
        wanted = _GSH_FOREST_TO_POC[forest_rxn_idx]
    elif family == "ugt":
        wanted = _UGT_FOREST_TO_POC[forest_rxn_idx]
    else:
        raise AssertionError(family)
    out = []
    for por in rule.metabolites(mol):
        if por.info["rxn_num"] in wanted:
            site = por.info["site"]
            site_fs = frozenset((site,) if isinstance(site, int) else site)
            out.append((site_fs, por.products))
    return out


def _has_products(rule, smi: str) -> bool:
    mol = MolFromSmiles(smi)
    assert mol is not None, smi
    return any(True for _ in rule.metabolites(mol))


def test_glucuronidation_defaults_to_bare_star():
    mol = MolFromSmiles(PHENOL)
    product, _info = next(Glucuronidation().metabolize(mol))
    smi = _smiles(product)
    assert "*" in smi
    assert "O=C(O)C1OC" not in smi
    assert smi == "*Oc1ccccc1"


def test_glucuronidation_full_structure_opt_in():
    mol = MolFromSmiles(PHENOL)
    product, _info = next(Glucuronidation(as_star=False).metabolize(mol))
    smi = _smiles(product)
    assert "*" not in smi
    assert "O=C(O)" in smi


def test_glucuronidation_star_label_glca():
    mol = MolFromSmiles(PHENOL)
    product, _info = next(Glucuronidation(star_label="GlcA").metabolize(mol))
    cx = mol_to_cxsmiles(product)
    assert cx is not None
    assert cx.split()[0] == "*Oc1ccccc1"
    assert "GlcA" in cx
    dummy = next(a for a in product.GetAtoms() if a.GetAtomicNum() == 0)
    assert dummy.GetProp("atomLabel") == "GlcA"


@pytest.mark.parametrize("smi", [ALIPHATIC_THIOL, THIOPHENOL, CYSTEINE_LIKE])
def test_thiols_match_gsh(smi):
    """Default Glutathionation hits thiols. ``include_thiol=False`` is unfinished."""

    assert _has_products(Glutathionation(), smi)


@pytest.mark.parametrize("smi", [STYRENE_OXIDE, BENZYL_CHLORIDE, TERMINAL_ALKENE])
def test_electrophiles_match_gsh(smi):
    assert _has_products(Glutathionation(), smi)


def test_glutathionation_defaults_to_star_not_peptide():
    mol = MolFromSmiles(STYRENE_OXIDE)
    product, _info = next(Glutathionation().metabolize(mol))
    smi = _smiles(product)
    assert "*" in smi
    assert "NC(CCC(=O)N" not in smi


@pytest.mark.parametrize("label", ["GSH", "Protein", "DNA", "Cyanide"])
def test_glutathionation_star_labels(label):
    rule = Glutathionation(star_label=label)
    mol = MolFromSmiles(STYRENE_OXIDE)
    product, _info = next(rule.metabolize(mol))
    cx = mol_to_cxsmiles(product)
    assert cx is not None
    assert "*" in cx.split()[0]
    assert label in cx
    dummy = next(a for a in product.GetAtoms() if a.GetAtomicNum() == 0)
    assert dummy.GetProp("atomLabel") == label


@pytest.mark.parametrize("label", ["Protein", "DNA", "Cyanide"])
def test_star_only_labels_reject_full_structure(label):
    with pytest.raises(ValueError, match="full conjugate"):
        Glutathionation(star_label=label, as_star=False)


# Skipped: include_thiol=False probes, GlutathionationNoThiol load_ruleset,
# and STAR_ONLY_LABELS import — unfinished poc public surface.


_SMARTS_PROBES = [
    ("ugt_acid", lambda: Glucuronidation(star_label="GlcA"), BENZOIC_ACID, 0, "ugt"),
    ("ugt_phenol", lambda: Glucuronidation(star_label="GlcA"), PHENOL, 1, "ugt"),
    ("gsh_epoxide", lambda: Glutathionation(star_label="GSH"), STYRENE_OXIDE, 0, "gsh"),
    ("gsh_chloride", lambda: Glutathionation(star_label="GSH"), BENZYL_CHLORIDE, 1, "gsh"),
    ("gsh_bromide", lambda: Glutathionation(star_label="GSH"), BENZYL_BROMIDE, 1, "gsh"),
    ("gsh_thiol", lambda: Glutathionation(star_label="GSH"), ALIPHATIC_THIOL, 2, "gsh"),
    ("gsh_alkene", lambda: Glutathionation(star_label="GSH"), TERMINAL_ALKENE, 3, "gsh"),
    ("gsh_michael", lambda: Glutathionation(star_label="GSH"), BENZOQUINONE, 4, "gsh"),
    ("protein_epoxide", lambda: Glutathionation(star_label="Protein"), STYRENE_OXIDE, 0, "gsh"),
    ("protein_chloride", lambda: Glutathionation(star_label="Protein"), BENZYL_CHLORIDE, 1, "gsh"),
    ("protein_thiol", lambda: Glutathionation(star_label="Protein"), ALIPHATIC_THIOL, 2, "gsh"),
    ("protein_alkene", lambda: Glutathionation(star_label="Protein"), TERMINAL_ALKENE, 3, "gsh"),
    ("protein_michael", lambda: Glutathionation(star_label="Protein"), BENZOQUINONE, 4, "gsh"),
    ("gsh_aldehyde", lambda: Glutathionation(star_label="GSH"), ACETALDEHYDE, 5, "gsh"),
    ("gsh_aziridine", lambda: Glutathionation(star_label="GSH"), N_PHENYLAZIRIDINE, 6, "gsh"),
    ("gsh_sulfonate", lambda: Glutathionation(star_label="GSH"), METHYL_MESYLATE, 7, "gsh"),
    ("gsh_isocyanate", lambda: Glutathionation(star_label="GSH"), PHENYL_ISOCYANATE, 8, "gsh"),
]


@pytest.mark.parametrize(
    "probe_id, rule_factory, smi, forest_rxn, family",
    _SMARTS_PROBES,
    ids=[p[0] for p in _SMARTS_PROBES],
)
def test_conjugation_smarts_som_is_single_atom(probe_id, rule_factory, smi, forest_rxn, family):
    hits = _sites_for_rxn(rule_factory(), smi, forest_rxn, family=family)
    assert hits, f"{probe_id}: no products for forest SMARTS {forest_rxn} on {smi}"
    for site, _products in hits:
        assert len(site) == 1, (
            f"{probe_id}: expected single-atom SOM, got {sorted(site)} "
            f"(len={len(site)}) for {smi}"
        )


@pytest.mark.parametrize(
    "probe_id, smi, forest_rxn, expected_o",
    [
        ("ugt_phenol", "Oc1ccccc1", 1, 0),
        ("ugt_acid", "O=C(O)c1ccccc1", 0, 2),
        ("ugt_apap", "CC(=O)Nc1ccc(O)cc1", 1, 8),
        ("ugt_ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 0, 12),
    ],
)
def test_ugt_som_is_oxygen(probe_id, smi, forest_rxn, expected_o):
    mol = MolFromSmiles(smi)
    assert mol is not None
    assert mol.GetAtomWithIdx(expected_o).GetAtomicNum() == 8
    hits = _sites_for_rxn(
        Glucuronidation(star_label="GlcA"), smi, forest_rxn, family="ugt"
    )
    assert hits, f"{probe_id}: no UGT products for forest rxn {forest_rxn}"
    for site, _ in hits:
        assert site == frozenset({expected_o}), (
            f"{probe_id}: expected SOM {{{expected_o}}}, got {sorted(site)}"
        )


@pytest.mark.parametrize(
    "smi, expected_soms, product_substr",
    [
        # Symmetry collapse (DIVERGENCES.md): one beta carbon stands for the class.
        # Forest lists {2, 3, 6, 7}.
        (BENZOQUINONE, {2}, "*C1C=C(O)C=CC1=O"),
        # Forest lists {5, 6, 9, 10}. Ranks keep one carbon from each class.
        (NAPQI, {5, 6}, None),
        (MENADIONE, {2}, "*C1C(=O)c2ccccc2C(O)=C1C"),
        (BETA_ENONE, {1}, None),
        (CINNAMALDEHYDE, {3}, None),
    ],
)
def test_glutathionation_michael_acceptors(smi, expected_soms, product_substr):
    can = MolToSmiles(MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 4, family="gsh")
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
    can = MolToSmiles(MolFromSmiles(smi))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 1, family="gsh")
    assert hits
    for site, products in hits:
        assert site == frozenset({1})
        assert _smiles(products[0]) == "*Cc1ccccc1"


@pytest.mark.parametrize(
    "rule_factory",
    [
        lambda: Glutathionation(star_label="GSH"),
        lambda: Glutathionation(star_label="Protein"),
        lambda: Glutathionation(star_label="DNA"),
        lambda: Glutathionation(star_label="Cyanide"),
    ],
)
def test_reactivity_heads_enumerate_benzoquinone(rule_factory):
    assert _has_products(rule_factory(), BENZOQUINONE)


@pytest.mark.parametrize(
    "smi, carbonyl_idx, product_smi",
    [
        (FORMALDEHYDE, 0, "*CO"),
        (ACETALDEHYDE, 1, "*C(C)O"),
        (BENZALDEHYDE, 1, "*C(O)c1ccccc1"),
        (METHYLGLYOXAL, 3, "*C(O)C(C)=O"),
        # progression: furfural carbonyl thiohemiacetal is emitted.
        (FURFURAL, 1, "*C(O)c1ccco1"),
    ],
)
def test_glutathionation_aldehyde_thiohemiacetal(smi, carbonyl_idx, product_smi):
    can = MolToSmiles(MolFromSmiles(smi))
    mol = MolFromSmiles(can)
    assert mol is not None
    assert mol.GetAtomWithIdx(carbonyl_idx).GetSymbol() == "C"
    assert mol.GetAtomWithIdx(carbonyl_idx).GetTotalNumHs() >= 1
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 5, family="gsh")
    assert hits, f"no aldehyde products for {can}"
    soms = {next(iter(site)) for site, _ in hits}
    assert carbonyl_idx in soms, f"{can}: expected carbonyl {carbonyl_idx} in {soms}"
    for site, products in hits:
        assert len(site) == 1
        assert site == frozenset({carbonyl_idx})
        assert _smiles(products[0]) == product_smi


def test_benzaldehyde_includes_carbonyl_thiohemiacetal():
    can = MolToSmiles(MolFromSmiles(BENZALDEHYDE))
    hits = _sites_for_rxn(Glutathionation(star_label="GSH"), can, 5, family="gsh")
    assert hits
    soms = {next(iter(site)) for site, _ in hits}
    assert 1 in soms
    assert any(_smiles(p[0]) == "*C(O)c1ccccc1" for _, p in hits)


@pytest.mark.parametrize(
    "smi, forest_rxn, expected_soms, product_smi",
    [
        # Both aziridine carbons are one rank class. Forest lists {5, 6}.
        (N_PHENYLAZIRIDINE, 6, {5}, "*CCNc1ccccc1"),
        (METHYL_MESYLATE, 7, {0}, "*C"),
        (BENZYL_MESYLATE, 7, {5}, "*Cc1ccccc1"),
        (PHENYL_ISOCYANATE, 8, {2}, "*C(=O)Nc1ccccc1"),
        (PHENYL_ISOTHIOCYANATE, 8, {2}, "*C(=S)Nc1ccccc1"),
    ],
)
def test_glutathionation_aziridine_sulfonate_isocyanate(
    smi, forest_rxn, expected_soms, product_smi
):
    can = MolToSmiles(MolFromSmiles(smi))
    hits = _sites_for_rxn(
        Glutathionation(star_label="GSH"), can, forest_rxn, family="gsh"
    )
    assert hits, f"no products for forest rxn {forest_rxn} on {can}"
    soms = {next(iter(site)) for site, _ in hits}
    assert expected_soms <= soms
    for site, products in hits:
        assert len(site) == 1
        assert _smiles(products[0]) == product_smi


@pytest.mark.parametrize(
    "smi, forest_rxn, reason",
    [
        (AZETIDINE_PHENYL, 6, "azetidine"),
        (PYRROLIDINE, 6, "pyrrolidine"),
        ("CCNc1ccccc1", 6, "aniline"),
        ("c1ccncc1", 6, "pyridine"),
        (STYRENE_OXIDE, 6, "epoxide"),
        (SULFONAMIDE, 7, "sulfonamide"),
        (SULFONE, 7, "sulfone"),
        (SULFONIC_ACID, 7, "sulfonic acid"),
        ("CS(C)=O", 7, "sulfoxide"),
        (BENZONITRILE, 8, "nitrile"),
        (ACETANILIDE, 8, "amide"),
        ("O=C=O", 8, "CO2"),
        ("CC#N", 8, "acetonitrile"),
        (CARBODIIMIDE, 8, "carbodiimide"),
    ],
)
def test_glutathionation_lookalikes_miss_dedicated_smarts(smi, forest_rxn, reason):
    can = MolToSmiles(MolFromSmiles(smi))
    hits = _sites_for_rxn(
        Glutathionation(star_label="GSH"), can, forest_rxn, family="gsh"
    )
    assert hits == [], f"{reason} ({can}) unexpectedly matched forest rxn {forest_rxn}: {hits}"


def test_glutathionation_michael_on_beta_substituted_enone():
    """Port of the trailing case in ``tests/test_rules.py``."""

    mol = MolFromSmiles("CC=CC(=O)CO")
    hits = list(Glutathionation().metabolites(mol))
    assert hits
    soms = set()
    for por in hits:
        site = por.info["site"]
        site_fs = frozenset((site,) if isinstance(site, int) else site)
        assert len(site_fs) == 1
        soms.add(next(iter(site_fs)))
        assert "*" in MolToSmiles(por.products[0])
    assert 1 in soms
