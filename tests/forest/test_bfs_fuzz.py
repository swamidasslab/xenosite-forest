"""Bounded DFS (and a small BFS smoke) over drug crashers + small organics.

No RDKit crash, no dotted product, traces present. Port of the hard
``tests/test_bfs_fuzz.py`` corpus (PhaseOne; live forest has no Full / star-expand).
DFS samples depth-2 without draining a BFS frontier on drug-sized mols.
A star acetyl is terminal here, so dehydrogenation does not expand it.

``canonical_emitted_sites`` is drawn with Hypothesis ``st.booleans()`` and
passed as kwargs on ``bfs`` / ``dfs`` (no env / custom flip helpers).
"""

from __future__ import annotations

import os as _os

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from rdkit import Chem

from xenosite.forest.find_path import bfs, dfs
from xenosite.forest.rules import Acetylation, Dehydrogenation, Hydroxylation
from xenosite.forest.rulesets import PhaseOne, RuleSet

# Prior RDKit crashers + suite anchors from forest bfs fuzz.

def _fuzz_examples(default: int) -> int:
    raw = _os.environ.get("XENOSITE_FUZZ_EXAMPLES")
    if raw:
        return int(raw)
    return default

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
    "c1ccccc1O",
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
_CAP = 25
_MAX_HEAVY = 36


@st.composite
def _generated_smiles(draw):
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
    if draw(st.booleans()):
        return draw(st.sampled_from(_CORPUS))
    return draw(_generated_smiles())


def _sample_dfs_depth2(smiles: str, *, canonical_emitted_sites: bool):
    """DFS reaches depth-2 without draining a BFS frontier (forest bfs fuzz)."""

    out = []
    depth2 = 0
    for product, _info in dfs(
        smiles,
        PhaseOne,
        depth=2,
        canonical_emitted_sites=canonical_emitted_sites,
    ):
        out.append(product)
        if product.xf.tracing.depth is not None and product.xf.tracing.depth >= 2:
            depth2 += 1
        if len(out) >= _CAP or depth2 >= 8:
            break
    return out


@given(smiles=forest_smiles(), canonical_emitted_sites=st.booleans())
@settings(
    max_examples=_fuzz_examples(16),
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_fuzz_dfs_depth2_is_connected_and_traced(
    smiles: str, canonical_emitted_sites: bool
):
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)
    products = _sample_dfs_depth2(
        smiles, canonical_emitted_sites=canonical_emitted_sites
    )
    # Tiny substrates may have no PhaseOne child; a miss is empty, not a crash.
    for product in products:
        assert "." not in product.xf.csmi
        assert product.xf.tracing.active
        assert product.xf.tracing.depth >= 1


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=20_000, derandomize=True)
def test_issue3_dfs_reaches_depth2(canonical_emitted_sites: bool):
    """Issue #3 parent: DFS under PhaseOne reaches a depth-2 metabolite."""

    parent = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"
    depth2 = 0
    for product, info in dfs(
        parent,
        PhaseOne,
        depth=2,
        canonical_emitted_sites=canonical_emitted_sites,
    ):
        if info is not None and product.xf.tracing.depth >= 2:
            depth2 += 1
            assert "." not in product.xf.csmi
        if depth2 >= 1:
            break
    assert depth2 >= 1


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=20_000, derandomize=True)
def test_bfs_smoke_on_small_substrates(canonical_emitted_sites: bool):
    """BFS still enumerates connected, traced products on tiny mols."""

    for smiles in ("CCO", "c1ccccc1", "CCN"):
        products = []
        for product, _info in bfs(
            smiles,
            PhaseOne,
            depth=1,
            canonical_emitted_sites=canonical_emitted_sites,
        ):
            products.append(product)
            if len(products) >= 5:
                break
        assert products
        for product in products:
            assert "." not in product.xf.csmi
            assert product.xf.tracing.active


def test_dehydrogenation_on_star_acetyl_is_terminal():
    parent = Chem.MolFromSmiles("CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21")
    acetyl = next(product for product, _info in Acetylation().metabolize(parent))
    assert any(atom.GetSymbol() == "*" for atom in acetyl.GetAtoms())
    children = list(Dehydrogenation().metabolize(acetyl))
    # Terminal conjugates are not expanded, including by dehydrogenation.
    assert children == []
    assert any(atom.GetSymbol() == "*" for atom in acetyl.GetAtoms())
    assert acetyl.xf.is_terminal


def test_bfs_does_not_expand_a_terminal():
    phenol = Chem.MolFromSmiles("Oc1ccccc1")
    acetyl, _info = next(Acetylation().metabolize(phenol))
    rules = RuleSet([Acetylation, Hydroxylation, Dehydrogenation], name="Mix")
    assert list(bfs(acetyl, rules, depth=2)) == []
    assert list(dfs(acetyl, rules, depth=2)) == []
