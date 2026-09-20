"""Nauty six-family orbits: ordered vs unordered; atom_bond identity."""

from __future__ import annotations

import importlib.util

import pytest

from xenosite.refactor_poc.graph_isomorphism import (
    all_site_pair_orbits_nauty,
    atom_bond_generators_nauty,
    bond_atom_orbits_from_nauty_generators,
    bond_atom_pair_orbits_from_nauty_generators,
    endpoint_bond_atom_sites,
)
from xenosite.refactor_poc.rdkitutil import MolFromSmiles

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


def test_benzene_pair_of_bond_atom_sites_needs_joint_orbit():
    """Second-order: adjacent ≠ opposite pair of composite sites; rotation matches.

    Individual bond–atom orbits are transitive on benzene endpoints, so the old
    coarse key (two first-order orbit ids) collapses adjacent, opposite, and
    rotated-adjacent placements. The joint pair orbit must separate adjacent
    from opposite and identify the rotation.
    """

    mol = _mol("c1ccccc1")
    generators = atom_bond_generators_nauty(mol, include_stereo=True)

    def bond(i: int, j: int) -> int:
        b = mol.GetBondBetweenAtoms(i, j)
        assert b is not None
        return b.GetIdx()

    b01 = bond(0, 1)
    b12 = bond(1, 2)
    b23 = bond(2, 3)
    b34 = bond(3, 4)

    primitive = endpoint_bond_atom_sites(mol)
    _, primitive_orbit = bond_atom_orbits_from_nauty_generators(
        primitive, generators
    )

    assert primitive_orbit[(b01, 0)] == primitive_orbit[(b12, 1)]
    assert primitive_orbit[(b01, 0)] == primitive_orbit[(b34, 3)]

    def old_coarse_key(left, right):
        return tuple(
            sorted((primitive_orbit[left], primitive_orbit[right]))
        )

    A = tuple(sorted(((b01, 0), (b12, 1))))  # adjacent bonds
    B = tuple(sorted(((b01, 0), (b34, 3))))  # opposite bonds
    C = tuple(sorted(((b23, 2), (b34, 3))))  # rotation of A

    assert old_coarse_key(*A) == old_coarse_key(*B)
    assert old_coarse_key(*A) == old_coarse_key(*C)

    _, pair_orbit = bond_atom_pair_orbits_from_nauty_generators(
        primitive, generators, ordered=False
    )

    assert pair_orbit[A] != pair_orbit[B]
    assert pair_orbit[A] == pair_orbit[C]
