"""Python-2 Metabolic Forest vs forest 0.2.3 (RDKit 2026) enumeration drift.

Py2 forest + old RDKit kept unsanitizable products and wrote aromatic ``can_smi``
strings that do not round-trip. Forest 0.2.3 fills valence caches, drops
RDKit-invalid products (whole fragment set), and emits kekulé canonical
SMILES.

Passing tests lock the 0.2.3 behavior. ``xfail(strict=True)`` tests encode
py2 strings or chemical overmatches we are not fixing yet; XPASS means the
generator changed. ``test_rdkit_warnings_do_not_print_to_stderr`` locks that
enumeration does not flood stderr with RDKit / sanitization noise.
"""

from __future__ import annotations

import warnings

import pytest
from rdkit import Chem

from xenosite.forest import PhaseOneRS, load_ruleset
from xenosite.forest.base import can_smi
from xenosite.forest.rules import QuinoneFormation
from xenosite.forest.utils import clean

DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"
APAP = "CC(=O)Nc1ccc(O)cc1"
STYRENE = "C=Cc1ccccc1"
SUDOXICAM = "CN1C(C(=O)Nc2nccs2)=C(O)c2ccccc2S1(=O)=O"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# py2 can_smi aromatic dealkylation (do not round-trip in RDKit 2026).
PY2_AROMATIC_DEALK_OH = "cccc(ccO)C(OCCN(C)C)c1=c-c=c-c=c-1"
PY2_AROMATIC_DEALK_OXO = "CN(C)CCOC(cccccc=O)c1=c-c=c-c=c-1"
PY2_PENTAVALENT_N_OXIDE = "CN(C)(O)CCOC(C1=CC=CC=C1)C1=CC=CC=C1"

# Forest 0.2.3 kekulé form of the same ring-opened dealkylation.
KEKULE_DEALK_OH = "C=CC(=CC=CO)C(OCCN(C)C)c1ccccc1"
N_OXIDE_ZWITTERION = "C[N+](C)([O-])CCOC(c1ccccc1)c1ccccc1"
N_DEMETHYL = "CNCCOC(c1ccccc1)c1ccccc1"
NAPQI = "CC(=O)N=C1C=CC(=O)C=C1"
STYRENE_VINYL_EPOXIDE = "c1ccc(C2CO2)cc1"
SUDOXICAM_ALDEHYDE = "CN1C(C=O)=C(O)c2ccccc2S1(=O)=O"
# Intact 2,5-thiazole quinone-imine (amide still attached).
SUDOXICAM_THIAZOLE_QUINONE_IMINE = "CN1C(C(=O)N=C2N=CC(=O)S2)=C(O)c2ccccc2S1(=O)=O"
# Thiazole fragment after amide C-N dealk.
SUDOXICAM_DEALK_QUINONE = "N=C1N=CC(=O)S1"


def _smiles_set(reactant: str, ruleset: str = "PhaseOneRS") -> set[str]:
    mol = Chem.MolFromSmiles(reactant)
    assert mol is not None, reactant
    out: set[str] = set()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rs = PhaseOneRS if ruleset == "PhaseOneRS" else load_ruleset(ruleset)
        for (_rule, _site), mols in rs.metabolites(mol, unique=True):
            for product in mols or []:
                smis = can_smi(rdmol=product)
                if smis:
                    out.add(smis[0])
    return out


def _quinone_formation_products(
    smiles: str,
) -> tuple[set[frozenset[int]], set[str]]:
    """Sites with at least one product, and canonical product SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    smis: set[str] = set()
    product_sites: set[frozenset[int]] = set()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for (_rule, site), mols in QuinoneFormation().metabolites(mol, unique=True):
            products = [p for p in (mols or []) if p is not None]
            if not products:
                continue
            product_sites.add(frozenset(site))
            for product in products:
                got = can_smi(rdmol=product)
                if got:
                    smis.add(got[0])
    return product_sites, smis


# --- pass: 0.2.3 is the intended generator ---------------------------------


def test_diphenhydramine_enumerates():
    """RDKit 2026 valence caches: Phase I no longer crashes on this suite molecule."""
    assert _smiles_set(DIPHENHYDRAMINE)


def test_diphenhydramine_products_round_trip():
    """Forest products parse; py2 aromatic ``can_smi`` strings often do not."""
    smis = _smiles_set(DIPHENHYDRAMINE)
    assert smis
    assert all(Chem.MolFromSmiles(s) is not None for s in smis)


def test_n_oxide_zwitterion_kept():
    """Valid N-oxide (charge-separated) is still emitted."""
    assert N_OXIDE_ZWITTERION in _smiles_set(DIPHENHYDRAMINE)


def test_n_demethyl_and_kekule_dealkylation_emitted():
    """N-dealkylation and ring-opened dealkylation use kekulé SMILES."""
    smis = _smiles_set(DIPHENHYDRAMINE)
    assert N_DEMETHYL in smis
    assert KEKULE_DEALK_OH in smis


def test_apap_napqi_quinone():
    """Acetaminophen still forms NAPQI (quinone-imine)."""
    smis = _smiles_set(APAP, "QF")
    assert NAPQI in smis


def test_styrene_vinyl_epoxide():
    mol = Chem.MolFromSmiles(STYRENE)
    smis = set()
    for (_rule, _site), mols in load_ruleset("SO.Epoxidation").metabolites(
        mol, unique=True
    ):
        for product in mols or []:
            got = can_smi(rdmol=product)
            if got:
                smis.add(got[0])
    assert STYRENE_VINYL_EPOXIDE in smis


def test_sudoxicam_valid_quinone_dealk_kept():
    """{6,8} (thiazole C2/C4) has no sanitizable products; {6,9} (C2/C5) does.

    Parent numbering: 6 = thiazole C2 (amide attachment), 8 = C4, 9 = C5.
    ``QuinoneFormation`` still *yields* the {6,8} SMARTS match, but ``clean()``
    empties that list — collect sites from non-empty lists, not iterator keys.

    {6,9} emits an intact 2,5-thiazole quinone-imine and a dealk pair: leftover
    benzothiazine-3-carbaldehyde (``Quinone=False``, amide carbon tagged
    ``dealk-noncarbon``) plus ``N=C1N=CC(=O)S1``. Both fragments parse, so
    ``clean()`` keeps the set. That aldehyde is a SMARTS overmatch (see xfail).
    """
    product_sites, smis = _quinone_formation_products(SUDOXICAM)
    assert frozenset({6, 8}) not in product_sites
    assert frozenset({6, 9}) in product_sites
    assert SUDOXICAM_THIAZOLE_QUINONE_IMINE in smis
    assert SUDOXICAM_DEALK_QUINONE in smis
    assert SUDOXICAM_ALDEHYDE in smis
    assert all(Chem.MolFromSmiles(s) is not None for s in smis)


def test_aspirin_products_round_trip():
    smis = _smiles_set(ASPIRIN)
    assert smis
    assert all(Chem.MolFromSmiles(s) is not None for s in smis)


def test_clean_does_not_keep_aldehyde_from_failed_quinone():
    """Py2 could keep the dealk leftover when the quinone fragment was invalid."""
    leftover = Chem.MolFromSmiles("CC=O")
    invalid = Chem.MolFromSmiles("[NH2]=c1nc(=O)cs1", sanitize=False)
    assert clean([leftover, invalid]) == []


# --- xfail: py2 strings / overmatches (strict: XPASS if the generator changes) ---


@pytest.mark.xfail(
    strict=True,
    reason="py2 kept pentavalent N-oxide CN(C)(O)…; forest 0.2.3 drops it as RDKit-invalid",
)
def test_py2_pentavalent_n_oxide_emitted():
    assert PY2_PENTAVALENT_N_OXIDE in _smiles_set(DIPHENHYDRAMINE)


@pytest.mark.xfail(
    strict=True,
    reason="py2 can_smi wrote aromatic ring-opened dealkylation; forest emits kekulé",
)
def test_py2_aromatic_dealkylation_oh_emitted():
    assert PY2_AROMATIC_DEALK_OH in _smiles_set(DIPHENHYDRAMINE)


@pytest.mark.xfail(
    strict=True,
    reason="py2 can_smi wrote unkekulized c1=c-c=c-c=c-1 dealkylation strings",
)
def test_py2_aromatic_dealkylation_oxo_emitted():
    assert PY2_AROMATIC_DEALK_OXO in _smiles_set(DIPHENHYDRAMINE)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "QuinoneFormation dealk ([#6R:1][#7,#8:2][#6:3]) matches the anilide "
        "C-N; leftover CHO is not a CYP N-dealkylation. Not tightening the "
        "SMARTS yet."
    ),
)
def test_sudoxicam_amide_dealk_aldehyde_not_emitted():
    """Amide carbonyl should not be treated as an N-dealkylation alkyl carbon.

    The dealk query is for aryl-O/N-alkyl (anisole/aniline → quinone + alkyl
    aldehyde). On sudoxicam it matches thiazole C2 – amide N – amide carbonyl,
    breaks C-N, and sanitizes the acyl fragment to benzothiazine-3-carbaldehyde.
    That is not a known sudoxicam metabolite and is not CYP N-dealkylation
    (which oxidizes α-CH2, not an amide carbonyl). The intact 2,5-thiazole
    quinone-imine at {6,9} is the plausible product and should remain.
    """
    _sites, smis = _quinone_formation_products(SUDOXICAM)
    assert SUDOXICAM_THIAZOLE_QUINONE_IMINE in smis
    assert SUDOXICAM_ALDEHYDE not in smis


# --- stderr must stay quiet during enumeration -------------------------------


def test_rdkit_warnings_do_not_print_to_stderr():
    """Enumeration must not print RDKit valence / kekulize noise to stderr.

    Pytest swallows warnings into its own log, so this runs a subprocess.
    ``rdBase.DisableLog`` covers RDKit's logger; dropped products are only
    logged at DEBUG and stay off stderr by default.
    """
    import subprocess
    import sys

    script = (
        "from rdkit import Chem\n"
        "from xenosite.forest import PhaseOneRS\n"
        f"mol = Chem.MolFromSmiles({DIPHENHYDRAMINE!r})\n"
        "list(PhaseOneRS.metabolites(mol, unique=True))\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        text=True,
    )
    err = proc.stderr
    assert "Explicit valence" not in err
    assert "Can't kekulize" not in err
    assert "Dropping RDKit-invalid" not in err
