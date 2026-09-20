"""RDKit 2026 / PhaseOne behavior locks (ported from forest drift suites).

Forest ``test_py2_rdkit_drift.py`` and ``test_rdkit_valence.py`` hit
``xenosite.forest``. These ports keep the *behavioral* asserts (enumerate,
round-trip, quiet stderr, key metabolites) on forest ``PhaseOne`` /
``QuinoneFormation`` so CI covers them without the forest API.
"""

from __future__ import annotations

import subprocess
import sys
import warnings

from rdkit import Chem

from xenosite.forest.phaseone import PhaseOneQF
from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import Epoxidation, QuinoneFormation
from xenosite.forest.rulesets import PhaseOne

DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"
APAP = "CC(=O)Nc1ccc(O)cc1"
STYRENE = "C=Cc1ccccc1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
N_OXIDE_ZWITTERION = "C[N+](C)([O-])CCOC(c1ccccc1)c1ccccc1"
N_DEMETHYL = "CNCCOC(c1ccccc1)c1ccccc1"
NAPQI = "CC(=O)N=C1C=CC(=O)C=C1"
STYRENE_VINYL_EPOXIDE = "c1ccc(C2CO2)cc1"
KEKULE_DEALK_OH = "C=CC(=CC=CO)C(OCCN(C)C)c1ccccc1"


def _smiles_set(reactant: str, ruleset=PhaseOne) -> set[str]:
    mol = MolFromSmiles(reactant)
    assert mol is not None, reactant
    out: set[str] = set()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for products, _info in ruleset.metabolize(mol):
            for product in products:
                out.add(product.xf.csmi)
    return out


def test_diphenhydramine_enumerates():
    assert _smiles_set(DIPHENHYDRAMINE)


def test_diphenhydramine_products_round_trip():
    smis = _smiles_set(DIPHENHYDRAMINE)
    assert smis
    assert all(Chem.MolFromSmiles(s) is not None for s in smis)


def test_n_oxide_zwitterion_kept():
    assert N_OXIDE_ZWITTERION in _smiles_set(DIPHENHYDRAMINE)


def test_n_demethyl_emitted():
    smis = _smiles_set(DIPHENHYDRAMINE)
    assert N_DEMETHYL in smis


def test_apap_napqi_quinone():
    smis = _smiles_set(APAP, PhaseOneQF)
    assert NAPQI in smis


def test_styrene_vinyl_epoxide():
    mol = MolFromSmiles(STYRENE)
    smis = {p.xf.csmi for _pl, _ in Epoxidation().metabolize(mol) for p in _pl}
    assert STYRENE_VINYL_EPOXIDE in smis


def test_aspirin_products_round_trip():
    smis = _smiles_set(ASPIRIN)
    assert smis
    assert all(Chem.MolFromSmiles(s) is not None for s in smis)


def test_quinone_formation_benzene_ortho_and_para():
    mol = MolFromSmiles("c1ccccc1")
    smis = {p.xf.csmi for _pl, _ in QuinoneFormation().metabolize(mol) for p in _pl}
    assert "O=C1C=CC=CC1=O" in smis
    assert "O=C1C=CC(=O)C=C1" in smis


def test_rdkit_warnings_do_not_print_to_stderr():
    script = (
        "from xenosite.forest.rdkit_api import MolFromSmiles\n"
        "from xenosite.forest.rulesets import PhaseOne\n"
        f"mol = MolFromSmiles({DIPHENHYDRAMINE!r})\n"
        "list(PhaseOne.metabolize(mol))\n"
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
