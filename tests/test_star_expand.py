"""Star-conjugate expansion option and DFS two-step sampling."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest import bfs, dfs
from xenosite.forest.rules import Acetylation
from xenosite.forest.utils import has_star_conjugate

ISSUE3_PARENT = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"


def _star_acetyl():
    parent = Chem.MolFromSmiles(ISSUE3_PARENT)
    return next(m for _, mols in Acetylation().metabolize(parent) for m in mols)


def test_has_star_conjugate():
    assert has_star_conjugate(Chem.MolFromSmiles("*NCC"))
    assert not has_star_conjugate(Chem.MolFromSmiles("CCN"))


def test_bfs_emits_star_acetyl_from_parent():
    """Star adducts are still returned as products (depth 1)."""
    saw = False
    for _, steps, mols in bfs(ISSUE3_PARENT, ruleset="Full", depth=1):
        if steps and steps[0][0] == "Acetylation":
            saw = True
            assert any(has_star_conjugate(m) for m in mols)
            break
    assert saw


def test_bfs_default_does_not_expand_star_conjugates():
    """Starting from a star adduct, default BFS yields nothing."""
    star = _star_acetyl()
    assert has_star_conjugate(star)
    rows = list(bfs(star, ruleset="Full", depth=1))
    assert rows == []


def test_bfs_expand_star_conjugates_metabolizes_star():
    star = _star_acetyl()
    n = 0
    for _, steps, _ in bfs(
        star, ruleset="Full", depth=1, expand_star_conjugates=True
    ):
        n += 1
        if n >= 1:
            assert steps
            break
    assert n >= 1


def test_dfs_samples_two_step_pathway_quickly():
    """DFS reaches a length-2 path without enumerating the full BFS frontier."""
    depth2 = None
    for smi, steps, _ in dfs(ISSUE3_PARENT, ruleset="Full", depth=2):
        if steps and len(steps) >= 2:
            depth2 = (smi, steps)
            break
    assert depth2 is not None


def test_dfs_default_skips_star_then_second_step():
    """Without expand_star_conjugates, DFS never continues past a star adduct."""
    n = 0
    for _, steps, _ in dfs(ISSUE3_PARENT, ruleset="Full", depth=2):
        n += 1
        if steps and len(steps) >= 2:
            assert steps[0][0] != "Acetylation"
        if steps and steps[0][0] == "Acetylation":
            assert len(steps) == 1
        if n >= 40:
            break
    assert n > 0
