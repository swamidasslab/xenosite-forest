"""Hypothesis fuzz: PhaseOneRS BFS must not raise RDKit valence preconditions.

Catches regressions like GitHub issue #3, where resonance copies left empty
implicit-H caches and ``RunReactants`` aborted mid-enumeration.
"""

from __future__ import annotations

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from rdkit import Chem

from xenosite.forest import bfs

# Prior RDKit 2026 crashers + diverse Phase I substrates (docs / suite).
_CORPUS = (
    "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21",  # issue #3
    "CN(C)CCOC(c1ccccc1)c1ccccc1",  # diphenhydramine
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",  # ibuprofen
    "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",  # atenolol
    "C(=Cc1ccccc1)CN1CCN(C(c2ccccc2)c2ccccc2)CC1",  # cinnarizine
    "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",  # warfarin
    "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",  # trimethoprim
    "COc1cc2nc(SCc3ccccc3C)[nH]c2cc1OC",  # omeprazole-like
    "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F",  # fluconazole
    "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",  # sulfisoxazole
    "CCC1(c2ccccc2)C(=O)NC(=O)NC1=O",  # phenobarbital
    "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",  # chloramphenicol
    "CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2",  # atropine
    "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N21",  # penicillin G core
    "CN1C(C(=O)Nc2nccs2)=C(O)c2ccccc2S1(=O)=O",  # sudoxicam-like
    "CC(=O)Nc1ccc(O)cc1",  # APAP
    "C=Cc1ccccc1",  # styrene
    "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
    "c1ccccc1",
    "CCO",
    "CCN",
    "CCCl",
    "CCS",
    "CC=O",
    "C1OC1",
    "O=C1C=CC(=O)C=C1",
)

_ALKYL = ("C", "CC", "CCC", "C(C)C", "CCCC", "CC(C)C")
_RING = (
    "c1ccccc1",
    "c1ccncc1",
    "c1ccoc1",
    "c1ccsc1",
    "C1CCCCC1",
    "c1ccc2c(c1)CCO2",
    "c1ccc2[nH]ccc2c1",
)
_FUNC = ("", "O", "N", "Cl", "F", "C(=O)O", "C(=O)N", "OC", "NC", "C=O", "C#N")

_MAX_PRODUCTS = 250
_MAX_HEAVY = 36


@st.composite
def _generated_smiles(draw):
    """Build small, RDKit-valid organic SMILES from alkyl / ring / functional pieces."""
    kind = draw(st.sampled_from(("alkyl-func", "ring-func", "alkyl-ring", "fused-ish")))
    if kind == "alkyl-func":
        smi = draw(st.sampled_from(_ALKYL)) + draw(st.sampled_from(_FUNC))
    elif kind == "ring-func":
        ring = draw(st.sampled_from(_RING))
        func = draw(st.sampled_from(_FUNC))
        smi = func + ring if func and not func.startswith("C") else ring + func
    elif kind == "alkyl-ring":
        smi = draw(st.sampled_from(_ALKYL)) + draw(st.sampled_from(_RING))
    else:
        smi = (
            draw(st.sampled_from(_ALKYL))
            + "C(=O)N"
            + draw(st.sampled_from(_ALKYL))
            + draw(st.sampled_from(_RING))
        )
    mol = Chem.MolFromSmiles(smi)
    assume(mol is not None)
    return Chem.MolToSmiles(mol)


@st.composite
def phaseone_smiles(draw):
    """Corpus crashers plus generated substrates for Phase I BFS fuzzing."""
    if draw(st.booleans()):
        return draw(st.sampled_from(_CORPUS))
    return draw(_generated_smiles())


def _drain_bfs(smiles: str, *, depth: int = 2, limit: int = _MAX_PRODUCTS) -> int:
    n = 0
    for _ in bfs(smiles, ruleset="PhaseOneRS", depth=depth):
        n += 1
        if n >= limit:
            break
    return n


@given(smiles=phaseone_smiles())
@settings(
    max_examples=40,
    deadline=15_000,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_bfs_phaseone_depth2_no_rdkit_precondition(smiles: str):
    """Depth-2 PhaseOneRS BFS must finish without RDKit valence RuntimeErrors."""
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)
    # Any uncaught RuntimeError (calcImplicitValence, etc.) fails the example.
    _drain_bfs(smiles, depth=2)
