"""Nauty six-family orbits: ordered vs unordered; atom_bond identity."""

from __future__ import annotations

import importlib.util

import pytest

from xenosite.forest.graph_isomorphism import all_site_pair_orbits_nauty
from xenosite.forest.rdkitutil import MolFromSmiles

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pynauty") is None,
    reason="pynauty required for nauty six-family orbits",
)


def _mol(smi: str):
    mol = MolFromSmiles(smi)
    assert mol is not None
    return mol


def _partition(groups: list[list[tuple[int, int]]]) -> set[frozenset[tuple[int, int]]]:
    return {frozenset(g) for g in groups}


def test_atom_bond_ordered_equals_unordered_partition():
    """Types already distinguish atom vs bond — ordering does not split groups."""

    mol = _mol("c1ccccc1")
    families = all_site_pair_orbits_nauty(mol)
    assert _partition(families["atom_bond_unordered"]) == _partition(
        families["atom_bond_ordered"]
    )
    # Same list object / shared computation is allowed; at least equal content.
    assert families["atom_bond_unordered"] == families["atom_bond_ordered"]


def test_atom_atom_ordered_refines_or_differs_from_unordered():
    """Ordered pairs distinguish (a,b) from (b,a) when the automorphism does not."""

    mol = _mol("CCO")  # asymmetric ends
    families = all_site_pair_orbits_nauty(mol)
    unordered = _partition(families["atom_atom_unordered"])
    ordered = _partition(families["atom_atom_ordered"])
    # Every unordered pair {a,b} lifts to two ordered pairs in some ordered groups.
    assert all(len(g) >= 1 for g in families["atom_atom_ordered"])
    # On ethanol, (0,2) and (2,0) are different ordered pairs.
    flat_ord = {p for g in families["atom_atom_ordered"] for p in g}
    assert (0, 2) in flat_ord and (2, 0) in flat_ord
    flat_un = {p for g in families["atom_atom_unordered"] for p in g}
    assert (0, 2) in flat_un
    assert (2, 0) not in flat_un  # unordered stores a1 < a2
    assert unordered != ordered or len(flat_ord) > len(flat_un)


def test_benzene_atom_atom_unordered_three_orbits():
    mol = _mol("c1ccccc1")
    families = all_site_pair_orbits_nauty(mol)
    # ortho / meta / para
    assert len(families["atom_atom_unordered"]) == 3


def test_six_family_keys_present():
    mol = _mol("c1ccccc1")
    families = all_site_pair_orbits_nauty(mol)
    assert set(families) == {
        "atom_atom_unordered",
        "atom_atom_ordered",
        "atom_bond_unordered",
        "atom_bond_ordered",
        "bond_bond_unordered",
        "bond_bond_ordered",
    }
