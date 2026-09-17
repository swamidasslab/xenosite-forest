"""BFS options for star conjugates and Full-ruleset depth≥2 behavior."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest import bfs
from xenosite.forest.utils import has_star_conjugate

ISSUE3_PARENT = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"


def test_has_star_conjugate():
    assert has_star_conjugate(Chem.MolFromSmiles("*NCC"))
    assert not has_star_conjugate(Chem.MolFromSmiles("CCN"))


def test_bfs_default_does_not_expand_star_conjugates():
    """Star adducts are emitted, but not metabolized further at depth>1."""
    rows = list(bfs(ISSUE3_PARENT, ruleset="Full", depth=2))
    assert any(
        steps and steps[0][0] == "Acetylation" and len(steps) == 1 for _, steps, _ in rows
    )
    # No depth-2 path whose first step is Acetylation (star not expanded).
    assert not any(
        steps and len(steps) >= 2 and steps[0][0] == "Acetylation"
        for _, steps, _ in rows
    )


def test_bfs_expand_star_conjugates_allows_depth2_after_acetylation():
    rows = list(
        bfs(
            ISSUE3_PARENT,
            ruleset="Full",
            depth=2,
            expand_star_conjugates=True,
        )
    )
    assert any(
        steps and len(steps) >= 2 and steps[0][0] == "Acetylation"
        for _, steps, _ in rows
    )
