"""Bounded ``bfs`` / ``dfs`` over a small corpus.

No RDKit crash, no dotted product, traces present. A star acetyl survives
dehydrogenation, and a terminal conjugate is not expanded.
"""

from __future__ import annotations

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from rdkit import Chem

from xenosite.refactor_poc.find_path import bfs, dfs
from xenosite.refactor_poc.rules import Acetylation, Dehydrogenation, Hydroxylation
from xenosite.refactor_poc.rulesets import PhaseOne, RuleSet

_CORPUS = (
    "c1ccccc1",
    "CCO",
    "CCN",
    "CCCl",
    "C=Cc1ccccc1",
    "CC(=O)Nc1ccc(O)cc1",
    "c1ccccc1O",
)
_CAP = 25


def _sample(enumerator, smiles: str):
    out = []
    for product, info in enumerator(smiles, PhaseOne, depth=2):
        out.append((product, info))
        if len(out) >= _CAP:
            break
    return out


@given(smiles=st.sampled_from(_CORPUS), which=st.sampled_from(("bfs", "dfs")))
@settings(
    max_examples=6,
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_fuzz_enumerate_depth2_is_connected_and_traced(smiles: str, which: str):
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= 24)
    enumerator = bfs if which == "bfs" else dfs
    products = _sample(enumerator, smiles)
    assert products
    for product, _info in products:
        assert "." not in product.xf.csmi
        assert product.xf.tracing.active
        assert product.xf.tracing.depth >= 1


def test_dehydrogenation_on_star_acetyl_keeps_star():
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
