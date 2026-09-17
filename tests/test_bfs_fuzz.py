"""Hypothesis fuzz: randomly sample Full DFS two-step pathways.

DFS + shuffled branch order explores depth=2 without draining BFS. Example
database lives in ``.hypothesis/`` (gitignored) and is restored on CI.
"""

from __future__ import annotations

import random
from pathlib import Path

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from hypothesis.database import DirectoryBasedExampleDatabase
from rdkit import Chem

from xenosite.forest import dfs
from xenosite.forest.rules import Acetylation, Dehydrogenation

# Repo-root example DB so local and CI share the same on-disk cache location.
_HYPOTHESIS_DIR = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
_HYPOTHESIS_DB = DirectoryBasedExampleDatabase(str(_HYPOTHESIS_DIR))

# Prior RDKit 2026 crashers + diverse substrates (docs / suite).
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

# Fair sample of two-step paths; bound total yields so CI stays under timeout.
_SAMPLE_DEPTH2 = 50
_MAX_YIELDS = 400
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
def forest_smiles(draw):
    """Corpus crashers plus generated substrates for pathway fuzzing."""
    if draw(st.booleans()):
        return draw(st.sampled_from(_CORPUS))
    return draw(_generated_smiles())


def _sample_dfs_two_step(
    smiles: str,
    *,
    expand_star_conjugates: bool = False,
    seed: int = 0,
    sample: int | None = _SAMPLE_DEPTH2,
    max_paths: int | None = _MAX_YIELDS,
) -> tuple[int, int]:
    """Return (n_yields, n_depth2) after randomly sampling DFS pathways.

    ``sample`` / ``max_paths`` are optional caps (``None`` = unlimited).
    """
    rng = random.Random(seed)
    n = 0
    depth2 = 0
    for _, steps, _ in dfs(
        smiles,
        ruleset="Full",
        depth=2,
        expand_star_conjugates=expand_star_conjugates,
        shuffle_rng=rng,
        max_paths=max_paths,
    ):
        n += 1
        if steps and len(steps) >= 2:
            depth2 += 1
        if sample is not None and depth2 >= sample:
            break
    return n, depth2


@given(
    smiles=forest_smiles(),
    expand_star_conjugates=st.booleans(),
    seed=st.integers(0, 2**32 - 1),
)
@settings(
    max_examples=50,
    deadline=20_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_full_dfs_two_step_sample_no_rdkit_runtime_error(
    smiles: str, expand_star_conjugates: bool, seed: int
):
    """Randomly sample Full DFS depth-2 pathways; no RDKit RuntimeErrors."""
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)
    _sample_dfs_two_step(
        smiles,
        expand_star_conjugates=expand_star_conjugates,
        seed=seed,
    )


def test_full_dfs_two_step_issue3_with_star_expand():
    """Issue #3 parent: shuffled DFS reaches depth-2 with star expansion."""
    parent = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"
    n, depth2 = _sample_dfs_two_step(
        parent,
        expand_star_conjugates=True,
        seed=1,
        sample=20,
        max_paths=200,
    )
    assert n > 0
    assert depth2 >= 1


def test_dehydrogenation_on_star_acetyl_keeps_star():
    """Resonance reassembly must not strip conjugation ``*`` adducts."""
    parent = Chem.MolFromSmiles("CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21")
    acetyl = next(m for _, mols in Acetylation().metabolize(parent) for m in mols)
    assert any(a.GetSymbol() == "*" for a in acetyl.GetAtoms())
    n = sum(1 for _ in Dehydrogenation().metabolize(Chem.Mol(acetyl)))
    assert n > 0
    assert any(a.GetSymbol() == "*" for a in acetyl.GetAtoms())
