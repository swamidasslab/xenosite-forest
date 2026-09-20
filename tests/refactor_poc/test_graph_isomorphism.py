"""Pair-orbit validation: invariants + hand cases + optional pynauty oracle.

Validation matrix (what each test proves)
=========================================

Invariant-backed (no pynauty required)
--------------------------------------
- ``test_benzene_meta_and_para_share_ranks_but_not_atom_pair_orbits``:
  negative case — same atom ranks, different isotope pair_group (meta ≠ para).
- ``test_atom_signature_shape_and_unordered``:
  signature is ``AtomPairOrbitSignature``; key(a,b)==key(b,a); stable across calls.
- ``test_singleton_topeqiv_uses_trivial_pair_group``:
  either end size-1 → ``TRIVIAL_PAIR_GROUP``; no smiles/nauty table materialization.
- ``test_multi_multi_materializes_forest_cache_once``:
  both ends multi → forest nested cache; second call reuses same table object;
  ``pair_group`` is sequential ``PairGroupId``.
- ``test_benzene_bond_pairs_split_adjacent_skip_opposite``:
  hand bond-pair orbit sizes on benzene (3, 6, 6).
- ``test_benzene_bond_atom_isotope_separates_endpoint_classes``:
  bond–atom hand classes (endpoint / adjacent / opposite).
- ``test_naphthalene_hand_atom_pairs_distinct``:
  chemically clear 1–2 / 1–4 / 1–5 / 1–8 style pairs get distinct isotope ids
  when they are different orbits.
- ``test_partition_property_fixed_group_pair``:
  for fixed (ga,gb), every unordered pair gets a pair_group_id; equal ids iff
  isotope keys agree.
- ``test_bond_pair_and_bond_atom_signature_shapes``:
  NamedTuple types / unordered bond pairs / bond–atom ordered by kind.
- ``test_cip_sort_key_shapes_atom_before_bond``:
  atom / bond / pair / group membership are sorted tuples; atom before bond.
- ``test_quinone_unique_edit_keeps_ortho_and_para``:
  integration — QuinoneFormation emits both ortho and para products on benzene.

Oracle-backed (pynauty; skip only if import fails after optional install)
------------------------------------------------------------------------
- ``test_oracle_isotope_vs_pynauty_partitions``:
  atom–atom / bond–bond / bond–atom orbit partitions agree on several mols.
- ``test_oracle_canonical_pair_group_ids_agree``:
  when partitions agree, sequential CIP-numbered ``PairGroupId`` values agree
  across isotope and nauty tables.
- ``test_nauty_cache_fills_all_modes_up_front``:
  first multi-multi query fills nested nauty tables for all three modes.

Dispatcher-branch (depends on whether pynauty is installed)
-----------------------------------------------------------
- ``test_dispatcher_uses_isotope_pair_group_when_pynauty_absent``:
  skipped when pynauty is present (dispatcher prefers nauty).
"""

from __future__ import annotations

import importlib.util
from collections import defaultdict
from itertools import combinations

import pytest

from xenosite.refactor_poc.graph_isomorphism import (
    TRIVIAL_PAIR_GROUP,
    _nested_tables_from_groups,
    atom_pair_orbit_isotope,
    atom_pair_orbit_key,
    atom_site_cip_key,
    bond_atom_orbit_isotope,
    bond_atom_orbit_key,
    bond_pair_orbit_isotope,
    bond_pair_orbit_key,
    bond_site_cip_key,
    orbit_group_cip_key,
    orbit_membership,
    site_pair_cip_key,
    site_pair_orbits_nauty,
    site_pair_orbits_smiles,
)
from xenosite.refactor_poc.rdkitutil import Mol, MolFromSmiles, cip_ids
from xenosite.refactor_poc.records import (
    AtomPairOrbitSignature,
    BondAtomOrbitSignature,
    BondPairOrbitSignature,
    PairGroupId,
    TopoGroupId,
)
from xenosite.refactor_poc.rules import QuinoneFormation, _incident_orders

# Symmetric aromatics + a few larger poc/H2H substrates for oracle coverage.
_ORACLE_MOLSMILES = (
    "c1ccccc1",  # benzene
    "c1ccc2ccccc2c1",  # naphthalene
    "Cc1ccccc1",  # toluene (mixed singleton + multi)
    "Oc1ccccc1",  # phenol
    "c1ccc2c(c1)OCO2",  # benzodioxole
    "COc1ccc(O)cc1",  # MeOPhOH (H2H)
)


def _mol(smiles: str) -> Mol:
    mol = MolFromSmiles(smiles)
    assert mol is not None
    return mol


def _rank_key(mol: Mol, left: int, right: int):
    ranks = mol.xf.topol_equiv
    mapped = {1: left, 2: right}
    return (
        tuple((mapno, ranks[idx]) for mapno, idx in sorted(mapped.items())),
        _incident_orders(mol, ranks, mapped),
    )


def _has_pynauty() -> bool:
    return importlib.util.find_spec("pynauty") is not None


def _group_sizes(groups: dict) -> list[int]:
    return sorted(len(v) for v in groups.values())


def test_benzene_meta_and_para_share_ranks_but_not_atom_pair_orbits():
    """Negative: meta and para share ranks but not pair_group_id."""

    mol = _mol("c1ccccc1")
    meta, other_meta, para = (0, 2), (1, 3), (0, 3)
    assert _rank_key(mol, *meta) == _rank_key(mol, *para)
    assert atom_pair_orbit_isotope(mol, *meta) != atom_pair_orbit_isotope(mol, *para)
    assert atom_pair_orbit_isotope(mol, *meta) == atom_pair_orbit_isotope(
        mol, *other_meta
    )
    assert atom_pair_orbit_isotope(mol, *para) == atom_pair_orbit_isotope(
        mol, para[1], para[0]
    )
    groups: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, j in combinations(range(6), 2):
        groups[atom_pair_orbit_isotope(mol, i, j)].append((i, j))
    assert _group_sizes(groups) == [3, 6, 6]


def test_atom_signature_shape_and_unordered():
    mol = _mol("c1ccccc1")
    meta = atom_pair_orbit_key(mol, frozenset({0, 2}))
    para = atom_pair_orbit_key(mol, frozenset({0, 3}))
    assert isinstance(meta, AtomPairOrbitSignature)
    assert isinstance(para, AtomPairOrbitSignature)
    assert meta is not None and para is not None
    g = TopoGroupId(mol.xf.topol_equiv[0])
    assert meta.groups == (g, g) == para.groups
    assert isinstance(meta.pair_group, int)
    assert isinstance(para.pair_group, int)
    assert meta.pair_group != para.pair_group
    assert meta.pair_group >= 0 and para.pair_group >= 0
    assert atom_pair_orbit_key(mol, frozenset({2, 0})) == meta
    assert atom_pair_orbit_key(mol, frozenset({1, 3})) == meta
    assert atom_pair_orbit_key(mol, frozenset({0})) is None
    assert atom_pair_orbit_key(mol, frozenset({0, 2})) == meta


def test_singleton_topeqiv_uses_trivial_pair_group():
    mol = _mol("Cc1ccccc1")
    methyl = next(
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 3
    )
    ring = next(
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetTotalNumHs() == 1
    )
    sig = atom_pair_orbit_key(mol, frozenset({methyl, ring}))
    assert isinstance(sig, AtomPairOrbitSignature)
    ranks = mol.xf.topol_equiv
    assert sig.groups == tuple(
        sorted((TopoGroupId(ranks[methyl]), TopoGroupId(ranks[ring])))
    )
    assert sig.pair_group == TRIVIAL_PAIR_GROUP
    assert sig.pair_group == PairGroupId(-1)
    structure = mol._forest["cache"]
    assert "site_pair_orbits_smiles" not in structure
    assert "site_pair_orbits_nauty" not in structure


def test_multi_multi_materializes_forest_cache_once():
    mol = _mol("c1ccccc1")
    sig = atom_pair_orbit_key(mol, frozenset({0, 3}))
    assert isinstance(sig, AtomPairOrbitSignature)
    structure = mol._forest["cache"]
    if _has_pynauty():
        cache_key = "site_pair_orbits_nauty"
    else:
        cache_key = "site_pair_orbits_smiles"
    assert isinstance(sig.pair_group, int)
    assert sig.pair_group >= 0
    cached = structure[cache_key]
    assert "atom_atom" in cached
    assert sig.groups in cached["atom_atom"]
    assert (
        cached["atom_atom"][sig.groups][tuple(sorted((0, 3)))] == sig.pair_group
    )
    again = atom_pair_orbit_key(mol, frozenset({1, 4}))
    assert again == sig
    assert structure[cache_key] is cached


def test_cip_sort_key_shapes_atom_before_bond():
    """Concrete CIP tuple shapes; bonds are sorted endpoint tuples."""

    mol = _mol("c1ccccc1")
    cip = cip_ids(mol, include_stereo=False)
    atom_key = atom_site_cip_key(cip, 0)
    assert atom_key == ("atom", cip[0])
    bond = mol.GetBondWithIdx(0)
    ends = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    bond_key = bond_site_cip_key(mol, cip, 0)
    assert bond_key == ("bond", cip[ends[0]], cip[ends[1]])
    assert atom_key < bond_key  # atom before bond by tag convention

    pair_key = site_pair_cip_key(mol, cip, "atom_atom", (0, 3))
    assert pair_key == tuple(sorted((atom_site_cip_key(cip, 0), atom_site_cip_key(cip, 3))))

    ba_key = site_pair_cip_key(mol, cip, "bond_atom", (0, 3))
    assert ba_key == (bond_site_cip_key(mol, cip, 0), atom_site_cip_key(cip, 3))

    membership = orbit_membership(mol, cip, "atom_atom", [(3, 0), (1, 2)])
    assert isinstance(membership, tuple)
    assert all(isinstance(p, tuple) and len(p) == 2 for p in membership)
    assert membership == tuple(sorted(membership))
    group_key = orbit_group_cip_key(mol, cip, "atom_atom", membership)
    assert isinstance(group_key, tuple)
    assert all(isinstance(p, tuple) and len(p) == 2 for p in group_key)


def test_benzene_bond_pairs_split_adjacent_skip_opposite():
    mol = _mol("c1ccccc1")
    groups: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, j in combinations(range(6), 2):
        key = bond_pair_orbit_isotope(mol, i, j)
        groups[key].append((i, j))
        assert bond_pair_orbit_isotope(mol, j, i) == key
    assert _group_sizes(groups) == [3, 6, 6]


def test_benzene_bond_atom_isotope_separates_endpoint_classes():
    mol = _mol("c1ccccc1")
    endpoint = bond_atom_orbit_isotope(mol, 0, 0)
    adjacent = bond_atom_orbit_isotope(mol, 0, 2)
    opposite = bond_atom_orbit_isotope(mol, 0, 3)
    assert len({endpoint, adjacent, opposite}) == 3


def test_naphthalene_hand_atom_pairs_distinct():
    """Naphthalene: 1–2, 1–4, 1–5, 1–8 style pairs (indices on first ring + bridge)."""

    mol = _mol("c1ccc2ccccc2c1")
    keys = {
        j: atom_pair_orbit_isotope(mol, 0, j)
        for j in range(1, mol.GetNumAtoms())
    }
    assert len(set(keys.values())) >= 3


def test_partition_property_fixed_group_pair():
    mol = _mol("c1ccccc1")
    ranks = mol.xf.topol_equiv
    g = ranks[0]
    by_isotope: dict[str, list[tuple[int, int]]] = defaultdict(list)
    by_sig: dict[object, list[tuple[int, int]]] = defaultdict(list)
    for i, j in combinations(range(6), 2):
        assert ranks[i] == g == ranks[j]
        iso = atom_pair_orbit_isotope(mol, i, j)
        sig = atom_pair_orbit_key(mol, frozenset({i, j}))
        assert isinstance(sig, AtomPairOrbitSignature)
        assert sig.groups == (TopoGroupId(g), TopoGroupId(g))
        by_isotope[iso].append((i, j))
        by_sig[sig.pair_group].append((i, j))
    assert sum(len(v) for v in by_isotope.values()) == 15
    iso_sets = {frozenset(v) for v in by_isotope.values()}
    sig_sets = {frozenset(v) for v in by_sig.values()}
    assert iso_sets == sig_sets
    # Sequential ids within the mode (benzene has three atom–atom orbits).
    assert set(by_sig) == {PairGroupId(0), PairGroupId(1), PairGroupId(2)}


def test_bond_pair_and_bond_atom_signature_shapes():
    mol = _mol("c1ccccc1")
    bp = bond_pair_orbit_key(mol, frozenset({0, 3}))
    assert isinstance(bp, BondPairOrbitSignature)
    assert isinstance(bp.pair_group, int)
    assert bond_pair_orbit_key(mol, frozenset({3, 0})) == bp
    assert bond_pair_orbit_key(mol, frozenset({0})) is None
    ba = bond_atom_orbit_key(mol, 0, 3)
    assert isinstance(ba, BondAtomOrbitSignature)
    ba2 = bond_atom_orbit_key(mol, 0, 2)
    assert ba.groups[0] == ba2.groups[0]  # same bond group
    assert ba.pair_group != ba2.pair_group


def test_quinone_unique_edit_keeps_ortho_and_para():
    """Integration: benzene quinone formation keeps ortho and para products."""

    mol = _mol("c1ccccc1")
    qf = QuinoneFormation()
    by_product: dict[str, set[object]] = defaultdict(set)
    for info, products in qf.metabolites(mol):
        site = info["site"]
        if len(site) != 2:
            continue
        sig = atom_pair_orbit_key(mol, site)
        assert isinstance(sig, AtomPairOrbitSignature)
        for product in products:
            by_product[product.xf.csmi].add(sig.pair_group)
    assert "O=C1C=CC=CC1=O" in by_product  # ortho
    assert "O=C1C=CC(=O)C=C1" in by_product  # para
    assert by_product["O=C1C=CC=CC1=O"].isdisjoint(by_product["O=C1C=CC(=O)C=C1"])


def test_dispatcher_uses_isotope_pair_group_when_pynauty_absent():
    if _has_pynauty():
        pytest.skip("pynauty installed; dispatcher uses nauty pair_group")
    mol = _mol("c1ccccc1")
    sig = atom_pair_orbit_key(mol, frozenset({0, 3}))
    assert sig is not None
    assert isinstance(sig.pair_group, int)
    assert sig.pair_group >= 0
    cached = mol._forest["cache"]["site_pair_orbits_smiles"]
    assert (
        cached["atom_atom"][sig.groups][tuple(sorted((0, 3)))] == sig.pair_group
    )
    bsig = bond_pair_orbit_key(mol, frozenset({0, 3}))
    assert bsig is not None
    assert isinstance(bsig.pair_group, int)
    assert isinstance(bond_atom_orbit_key(mol, 0, 3).pair_group, int)


# ---------------------------------------------------------------------------
# Oracle (pynauty)
# ---------------------------------------------------------------------------


def _orbit_partition_isotope(mol: Mol, mode: str) -> set[frozenset[tuple[int, int]]]:
    groups: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if mode == "atom_atom":
        n = mol.GetNumAtoms()
        for i, j in combinations(range(n), 2):
            groups[atom_pair_orbit_isotope(mol, i, j)].append((i, j))
    elif mode == "bond_bond":
        n = mol.GetNumBonds()
        for i, j in combinations(range(n), 2):
            groups[bond_pair_orbit_isotope(mol, i, j)].append(
                tuple(sorted((i, j)))
            )
    else:
        for b in range(mol.GetNumBonds()):
            for a in range(mol.GetNumAtoms()):
                groups[bond_atom_orbit_isotope(mol, b, a)].append((b, a))
    return {frozenset(v) for v in groups.values()}


def test_oracle_isotope_vs_pynauty_partitions():
    pytest.importorskip("pynauty")
    for smiles in _ORACLE_MOLSMILES:
        mol = _mol(smiles)
        for mode in ("atom_atom", "bond_bond", "bond_atom"):
            if mode == "bond_bond" and mol.GetNumBonds() < 2:
                continue
            smiles_sets = _orbit_partition_isotope(mol, mode)
            nauty = site_pair_orbits_nauty(mol, modes=(mode,))[mode]  # type: ignore[arg-type]
            nauty_sets = {frozenset(g) for g in nauty}
            assert smiles_sets == nauty_sets, (smiles, mode)


def test_oracle_canonical_pair_group_ids_agree():
    """Partitions and sequential CIP PairGroupIds agree across backends."""

    pytest.importorskip("pynauty")
    for smiles in _ORACLE_MOLSMILES:
        mol = _mol(smiles)
        smiles_orbits = site_pair_orbits_smiles(mol)
        nauty_orbits = site_pair_orbits_nauty(mol)
        smiles_tables = _nested_tables_from_groups(mol, smiles_orbits)
        nauty_tables = _nested_tables_from_groups(mol, nauty_orbits)
        for mode in ("atom_atom", "bond_bond", "bond_atom"):
            if mode == "bond_bond" and mol.GetNumBonds() < 2:
                continue
            assert {frozenset(g) for g in smiles_orbits[mode]} == {
                frozenset(g) for g in nauty_orbits[mode]
            }, (smiles, mode)
            assert smiles_tables[mode] == nauty_tables[mode], (smiles, mode)
            # Memberships are sorted tuples; ids are 0..n-1 dense in the mode.
            ids = {
                gid
                for slice_map in smiles_tables[mode].values()
                for gid in slice_map.values()
            }
            if ids:
                assert ids == {PairGroupId(i) for i in range(len(ids))}


def test_nauty_cache_fills_all_modes_up_front():
    pytest.importorskip("pynauty")
    mol = _mol("c1ccccc1")
    sig = atom_pair_orbit_key(mol, frozenset({0, 3}))
    assert isinstance(sig, AtomPairOrbitSignature)
    cached = mol._forest["cache"]["site_pair_orbits_nauty"]
    assert set(cached) == {"atom_atom", "bond_bond", "bond_atom"}
    bp = bond_pair_orbit_key(mol, frozenset({0, 3}))
    assert isinstance(bp, BondPairOrbitSignature)
    assert bp.pair_group == cached["bond_bond"][bp.groups][tuple(sorted((0, 3)))]
    ba = bond_atom_orbit_key(mol, 0, 3)
    assert ba.pair_group == cached["bond_atom"][ba.groups][(0, 3)]


def test_unified_batch_helpers_still_cover_three_modes():
    mol = _mol("c1ccccc1")
    orbits = site_pair_orbits_smiles(mol)
    assert set(orbits) == {"atom_atom", "bond_bond", "bond_atom"}
    assert all(isinstance(g, tuple) for g in orbits["atom_atom"])
    assert _group_sizes({i: list(g) for i, g in enumerate(orbits["atom_atom"])}) == [
        3,
        6,
        6,
    ]
