"""Site-pair orbits under the molecular automorphism group.

Unified recipes compute ``atom_atom``, ``bond_bond``, and ``bond_atom``
together (share: marked-SMILES and optional ``pynauty``). Recipe 2
(self-substructure match) is not implemented.

Unique-edit pair signature (atom–atom; same shape for bond modes)::

    ((ga, gb), pair_group_id)

``ga, gb`` are topological group ids of the two ends (sorted for same-kind
pairs; bond then atom for ``bond_atom``). ``(ga, gb)`` selects a slice of
``(end_id, end_id) -> pair_group_id``. The signature is that outer group
tuple plus the looked-up ``pair_group_id``.

``pair_group_id`` is a sequential ``PairGroupId`` (int NewType) numbered
``0..n-1`` by sorting CIP-based orbit membership tuples. Singleton topeqiv
ends use ``TRIVIAL_PAIR_GROUP`` without materializing a table.

CIP sort-key shapes (atom before bond; all unordered collections are sorted
tuples)::

    atom site:  ("atom", cip)
    bond site:  ("bond", cip_lo, cip_hi)   # sorted endpoint CIPs
    site pair:  sorted (site_a, site_b) for same-kind; (bond, atom) for bond_atom
    group:      sorted tuple of site-pair keys → PairGroupId

Forest cache (on ``mol._forest["structure"]``, not on ephemeral ``Xf``):

- ``site_pair_orbits_nauty`` — nested tables filled up front on first need
  when backend is ``nauty`` (default when ``pynauty`` imports)::

      mode -> (ga, gb) -> {(a, b): pair_group_id, ...}

- ``site_pair_orbits_smiles`` — same shape; RDKit isotope path. Implemented and
  callable (``backend="smiles"`` / ``XENOSITE_PAIR_ORBIT_BACKEND=smiles``), but
  **not** the auto fallback when pynauty is absent — that default is ``none``
  (topeqiv + ``TRIVIAL_PAIR_GROUP``; product csmi is the safety net).

Public access is ``mol.xf.atom_pair_orbit_key`` / ``bond_pair_orbit_key`` /
``bond_atom_orbit_key`` / ``site_pair_orbits`` / ``pair_orbit_backend``. Free
functions here are the profileable implementations xf calls.
"""

from __future__ import annotations

import importlib
import os
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from itertools import combinations, product
from typing import Any, Final, Literal, cast

from xenosite.refactor_poc.rdkit_api import AssignStereochemistry, Mol, MolToSmiles
from xenosite.refactor_poc.rdkitutil import cip_ids
from xenosite.refactor_poc.records import (
    AtomPairOrbitSignature,
    AtomSiteCipKey,
    BondAtomOrbitSignature,
    BondPairOrbitSignature,
    BondSiteCipKey,
    NautyPairGroup,
    OrbitGroupCipKey,
    OrbitMembership,
    PairGroupId,
    PairOrbitByGroups,
    PairOrbitSignature,
    SiteCipKey,
    SitePairCipKey,
    SitePairOrbitTables,
    SmilesPairGroup,
    TopoGroupId,
)

PairMode = Literal["atom_atom", "bond_bond", "bond_atom"]
SiteKind = Literal["atom", "bond"]
MarkedSite = tuple[SiteKind, int]
# Dispatcher backends. ``smiles`` is RDKit isotope marking — implemented and
# callable, but not the default when pynauty is absent (does not improve speed).
PairOrbitBackend = Literal["nauty", "smiles", "none"]

PAIR_MODES: dict[PairMode, tuple[SiteKind, SiteKind]] = {
    "atom_atom": ("atom", "atom"),
    "bond_bond": ("bond", "bond"),
    "bond_atom": ("bond", "atom"),
}

_DEFAULT_MODES: tuple[PairMode, ...] = ("atom_atom", "bond_bond", "bond_atom")
_CACHE_NAUTY = "site_pair_orbits_nauty"
_CACHE_SMILES = "site_pair_orbits_smiles"
_ENV_BACKEND = "XENOSITE_PAIR_ORBIT_BACKEND"
_MEMBERSHIP_BASE = 4
# Documented trivial pair_group when either end is a singleton topeqiv group,
# or when backend is ``none`` for multi–multi (cheap path; no isotope tables).
# Reserved negative id so sequential orbit ids stay 0..n-1.
TRIVIAL_PAIR_GROUP: Final[PairGroupId] = PairGroupId(-1)

# Process-wide override. ``None`` → env var → auto (nauty if importable else none).
_backend_override: PairOrbitBackend | None = None


def get_pair_orbit_backend() -> PairOrbitBackend:
    """Resolved backend: override, else ``XENOSITE_PAIR_ORBIT_BACKEND``, else auto."""

    if _backend_override is not None:
        return _backend_override
    env = os.environ.get(_ENV_BACKEND, "").strip().lower()
    if env in ("nauty", "smiles", "none"):
        return cast(PairOrbitBackend, env)
    if _pynauty_available():
        return "nauty"
    return "none"


def set_pair_orbit_backend(backend: PairOrbitBackend | None) -> None:
    """Set process-wide backend, or ``None`` to clear (env / auto resume)."""

    global _backend_override
    if backend is not None and backend not in ("nauty", "smiles", "none"):
        raise ValueError(f"unknown pair_orbit_backend: {backend!r}")
    _backend_override = backend


# --- CIP sort keys (canonical ordering of orbit groups) ---


def atom_site_cip_key(cip: Sequence[int], atom_idx: int) -> AtomSiteCipKey:
    """Sort key for an atom site: ``("atom", cip)``."""

    return ("atom", int(cip[atom_idx]))


def bond_site_cip_key(mol: Mol, cip: Sequence[int], bond_idx: int) -> BondSiteCipKey:
    """Sort key for a bond site: ``("bond", cip_lo, cip_hi)``.

    Endpoints are a sorted index tuple first, then replaced by CIP ids.
    """

    bond = mol.GetBondWithIdx(bond_idx)
    ends = _sorted_pair(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return ("bond", int(cip[ends[0]]), int(cip[ends[1]]))


def site_cip_key(mol: Mol, cip: Sequence[int], site: MarkedSite) -> SiteCipKey:
    """CIP key for one marked site. Atom before bond by tag convention."""

    kind, idx = site
    if kind == "atom":
        return atom_site_cip_key(cip, idx)
    if kind == "bond":
        return bond_site_cip_key(mol, cip, idx)
    raise ValueError(f"Unknown site type: {kind}")


def site_pair_cip_key(
    mol: Mol, cip: Sequence[int], mode: PairMode, pair: tuple[int, int]
) -> SitePairCipKey:
    """CIP key for one index pair under ``mode``.

    Same-kind pairs are unordered → sorted site keys. ``bond_atom`` stays
    ``(bond_key, atom_key)``.
    """

    kind1, kind2 = PAIR_MODES[mode]
    left = site_cip_key(mol, cip, (kind1, pair[0]))
    right = site_cip_key(mol, cip, (kind2, pair[1]))
    if mode == "bond_atom":
        return (left, right)
    return cast(SitePairCipKey, tuple(sorted((left, right))))


def orbit_membership(
    mol: Mol,
    cip: Sequence[int],
    mode: PairMode,
    pairs: Iterable[tuple[int, int]],
) -> OrbitMembership:
    """Sorted tuple of index pairs, ordered by CIP site-pair keys."""

    unordered = mode != "bond_atom"
    normalized = [
        _sorted_pair(pair[0], pair[1]) if unordered else pair for pair in pairs
    ]
    # CIP primary; index pair breaks ties when ranks collide (e.g. benzene).
    return tuple(
        sorted(
            normalized,
            key=lambda p: (site_pair_cip_key(mol, cip, mode, p), p),
        )
    )


def orbit_group_cip_key(
    mol: Mol, cip: Sequence[int], mode: PairMode, membership: OrbitMembership
) -> OrbitGroupCipKey:
    """Group sort key: sorted tuple of CIP site-pair keys (membership order)."""

    return tuple(site_pair_cip_key(mol, cip, mode, pair) for pair in membership)


# --- Marked SMILES (recipe 1), named for profiling ---


def mark_site(mol: Mol, membership: list[int], site: MarkedSite, bit: int) -> None:
    """Mark an atom or all endpoints of a bond into ``membership``."""

    kind, idx = site
    if kind == "atom":
        membership[idx] |= bit
    elif kind == "bond":
        bond = mol.GetBondWithIdx(idx)
        membership[bond.GetBeginAtomIdx()] |= bit
        membership[bond.GetEndAtomIdx()] |= bit
    else:
        raise ValueError(f"Unknown site type: {kind}")


def marked_site_pair_smiles(
    mol: Mol, site1: MarkedSite, site2: MarkedSite, *, include_stereo: bool = True
) -> SmilesPairGroup:
    """Canonical key for one ordered site pair (site1=bit1, site2=bit2)."""

    marked = Mol(mol)
    membership = [0] * marked.GetNumAtoms()
    mark_site(marked, membership, site1, bit=1)
    mark_site(marked, membership, site2, bit=2)
    for atom in marked.GetAtoms():
        atom.SetIsotope(
            _MEMBERSHIP_BASE * atom.GetIsotope() + membership[atom.GetIdx()]
        )
    return MolToSmiles(marked, canonical=True, isomericSmiles=include_stereo)


def atom_pair_orbit_isotope(mol: Mol, left: int, right: int) -> SmilesPairGroup:
    """Recipe 1 for one unordered atom pair. Membership marks; min of both orientations."""

    return min(
        marked_site_pair_smiles(mol, ("atom", left), ("atom", right)),
        marked_site_pair_smiles(mol, ("atom", right), ("atom", left)),
    )


def bond_pair_orbit_isotope(mol: Mol, bond_a: int, bond_b: int) -> SmilesPairGroup:
    """Recipe 1 for one unordered bond pair."""

    return min(
        marked_site_pair_smiles(mol, ("bond", bond_a), ("bond", bond_b)),
        marked_site_pair_smiles(mol, ("bond", bond_b), ("bond", bond_a)),
    )


def bond_atom_orbit_isotope(mol: Mol, bond_idx: int, atom_idx: int) -> SmilesPairGroup:
    """Recipe 1 for one (bond, atom) pair. Types distinguish ends; no swap."""

    return marked_site_pair_smiles(mol, ("bond", bond_idx), ("atom", atom_idx))


def site_pair_orbits_smiles(
    mol: Mol,
    modes: Sequence[PairMode] = _DEFAULT_MODES,
    bond_indices: Iterable[int] | None = None,
    atom_indices: Iterable[int] | None = None,
    *,
    include_stereo: bool = True,
) -> dict[PairMode, list[OrbitMembership]]:
    """Unified RDKit marked-SMILES orbits. Each group is a CIP-sorted tuple."""

    bonds = _index_tuple(bond_indices, mol.GetNumBonds())
    atoms = _index_tuple(atom_indices, mol.GetNumAtoms())
    cip = cip_ids(mol, include_stereo=include_stereo)
    result: dict[PairMode, list[OrbitMembership]] = {}

    for mode in modes:
        kind1, kind2 = PAIR_MODES[mode]
        sites1: list[MarkedSite] = [
            (kind1, i) for i in (bonds if kind1 == "bond" else atoms)
        ]
        sites2: list[MarkedSite] = [
            (kind2, i) for i in (bonds if kind2 == "bond" else atoms)
        ]
        groups: dict[str, list[tuple[int, int]]] = defaultdict(list)

        if kind1 == kind2:
            for s1, s2 in combinations(sites1, 2):
                key = min(
                    marked_site_pair_smiles(mol, s1, s2, include_stereo=include_stereo),
                    marked_site_pair_smiles(mol, s2, s1, include_stereo=include_stereo),
                )
                groups[key].append(_sorted_pair(s1[1], s2[1]))
        else:
            for s1, s2 in product(sites1, sites2):
                key = marked_site_pair_smiles(
                    mol, s1, s2, include_stereo=include_stereo
                )
                groups[key].append((s1[1], s2[1]))

        result[mode] = [
            orbit_membership(mol, cip, mode, pairs) for pairs in groups.values()
        ]
    return result


# --- pynauty (recipe 3), named for profiling ---


def atom_bond_generators_nauty(
    mol: Mol, *, include_stereo: bool = True
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Generators as ``(atom_permutation, bond_permutation)``."""

    nauty = cast(Any, importlib.import_module("pynauty"))
    graph, n_atoms, n_bonds = _colored_graph(mol, nauty, include_stereo=include_stereo)
    raw_generators, _, _, _, _ = nauty.autgrp(graph)
    generators: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for perm in raw_generators:
        atom_map = tuple(int(perm[i]) for i in range(n_atoms))
        bond_map = tuple(int(perm[n_atoms + b]) - n_atoms for b in range(n_bonds))
        generators.append((atom_map, bond_map))
    return generators


def site_pair_orbits_nauty(
    mol: Mol,
    modes: Sequence[PairMode] = _DEFAULT_MODES,
    bond_indices: Iterable[int] | None = None,
    atom_indices: Iterable[int] | None = None,
    *,
    include_stereo: bool = True,
) -> dict[PairMode, list[OrbitMembership]]:
    """Unified nauty orbits. Each group is a CIP-sorted tuple of index pairs."""

    bonds = _index_tuple(bond_indices, mol.GetNumBonds())
    atoms = _index_tuple(atom_indices, mol.GetNumAtoms())
    generators = atom_bond_generators_nauty(mol, include_stereo=include_stereo)
    cip = cip_ids(mol, include_stereo=include_stereo)
    result: dict[PairMode, list[OrbitMembership]] = {}

    for mode in modes:
        if mode == "atom_atom":
            candidates = list(combinations(atoms, 2))

            def apply_generator(
                pair: tuple[int, int],
                atom_map: tuple[int, ...],
                bond_map: tuple[int, ...],
            ) -> tuple[int, int]:
                a1, a2 = pair
                return _sorted_pair(atom_map[a1], atom_map[a2])

        elif mode == "bond_bond":
            candidates = list(combinations(bonds, 2))

            def apply_generator(
                pair: tuple[int, int],
                atom_map: tuple[int, ...],
                bond_map: tuple[int, ...],
            ) -> tuple[int, int]:
                b1, b2 = pair
                return _sorted_pair(bond_map[b1], bond_map[b2])

        elif mode == "bond_atom":
            candidates = list(product(bonds, atoms))

            def apply_generator(
                pair: tuple[int, int],
                atom_map: tuple[int, ...],
                bond_map: tuple[int, ...],
            ) -> tuple[int, int]:
                b, a = pair
                return (bond_map[b], atom_map[a])

        else:
            raise ValueError(mode)

        candidate_set = set(candidates)
        parent = {x: x for x in candidates}

        def find(x: tuple[int, int]) -> tuple[int, int]:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: tuple[int, int], y: tuple[int, int]) -> None:
            x, y = find(x), find(y)
            if x != y:
                parent[y] = x

        for pair in candidates:
            for atom_map, bond_map in generators:
                image = apply_generator(pair, atom_map, bond_map)
                if image in candidate_set:
                    union(pair, image)

        buckets: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
        for pair in candidates:
            buckets[find(pair)].append(pair)
        result[mode] = [
            orbit_membership(mol, cip, mode, pairs) for pairs in buckets.values()
        ]

    return result


def atom_pair_orbit_pynauty(mol: Mol, left: int, right: int) -> NautyPairGroup:
    """Recipe 3 for one atom pair. Orbit as a CIP-sorted tuple of pairs."""

    groups = site_pair_orbits_nauty(mol, modes=("atom_atom",))["atom_atom"]
    needle = _sorted_pair(left, right)
    for group in groups:
        if needle in group:
            return group
    raise KeyError((left, right))


def bond_pair_orbit_pynauty(mol: Mol, bond_a: int, bond_b: int) -> NautyPairGroup:
    """Recipe 3 for one bond pair."""

    groups = site_pair_orbits_nauty(mol, modes=("bond_bond",))["bond_bond"]
    needle = _sorted_pair(bond_a, bond_b)
    for group in groups:
        if needle in group:
            return group
    raise KeyError((bond_a, bond_b))


def bond_atom_orbit_pynauty(
    mol: Mol, bond_idx: int, atom_idx: int
) -> NautyPairGroup:
    """Recipe 3 for one (bond, atom) pair."""

    groups = site_pair_orbits_nauty(mol, modes=("bond_atom",))["bond_atom"]
    needle = (bond_idx, atom_idx)
    for group in groups:
        if needle in group:
            return group
    raise KeyError(needle)


# --- Dispatch + forest cache ---
#
# pair_group_id is PairGroupId: sequential 0..n-1 from CIP-sorted membership
# tuples (same numbering for nauty and isotope when partitions agree).


def atom_pair_orbit_key(
    mol: Mol, site: frozenset[int]
) -> AtomPairOrbitSignature | None:
    """``AtomPairOrbitSignature(groups, pair_group)`` for a two-atom site."""

    if len(site) != 2:
        return None
    left, right = tuple(site)
    return cast(
        AtomPairOrbitSignature,
        _pair_orbit_signature(mol, "atom_atom", left, right),
    )


def bond_pair_orbit_key(
    mol: Mol, bonds: frozenset[int]
) -> BondPairOrbitSignature | None:
    """``BondPairOrbitSignature(groups, pair_group)`` for an unordered bond pair."""

    if len(bonds) != 2:
        return None
    left, right = tuple(bonds)
    return cast(
        BondPairOrbitSignature,
        _pair_orbit_signature(mol, "bond_bond", left, right),
    )


def bond_atom_orbit_key(
    mol: Mol, bond_idx: int, atom_idx: int
) -> BondAtomOrbitSignature:
    """``BondAtomOrbitSignature``. Hook for Dehydrogenation; not wired yet."""

    return cast(
        BondAtomOrbitSignature,
        _pair_orbit_signature(mol, "bond_atom", bond_idx, atom_idx),
    )


# Compat names used by the unique-edit signature before the rename.
pair_orbit_isotope = atom_pair_orbit_isotope
pair_orbit_pynauty = atom_pair_orbit_pynauty
pair_orbit_signature = atom_pair_orbit_key


def _pair_orbit_signature(
    mol: Mol, mode: PairMode, left: int, right: int
) -> PairOrbitSignature:
    unordered = mode != "bond_atom"
    pair = _sorted_pair(left, right) if unordered else (left, right)
    g_left = _end_group_id(mol, mode, left, end=0)
    g_right = _end_group_id(mol, mode, right, end=1)
    groups = _group_pair(mode, g_left, g_right)

    if not _both_ends_multi(mol, mode, left, right):
        return _make_signature(mode, groups, TRIVIAL_PAIR_GROUP)

    backend = get_pair_orbit_backend()
    if backend == "none":
        # Cheap path: ranks already on the signature; skip isotope tables.
        return _make_signature(mode, groups, TRIVIAL_PAIR_GROUP)
    if backend == "nauty":
        tables = _ensure_nauty_tables(mol)
    else:
        tables = _ensure_smiles_tables(mol)
    pair_group = tables[mode][groups][pair]
    return _make_signature(mode, groups, pair_group)


def ensure_site_pair_orbit_tables(
    mol: Mol, backend: PairOrbitBackend | None = None
) -> SitePairOrbitTables | None:
    """Materialize nested forest tables for ``backend`` (default: resolved).

    ``none`` returns ``None`` without building. Used by ``mol.xf.site_pair_orbits``.
    """

    chosen = get_pair_orbit_backend() if backend is None else backend
    if chosen == "none":
        return None
    if chosen == "nauty":
        return _ensure_nauty_tables(mol)
    return _ensure_smiles_tables(mol)


def _make_signature(
    mode: PairMode,
    groups: tuple[TopoGroupId, TopoGroupId],
    pair_group: PairGroupId,
) -> PairOrbitSignature:
    if mode == "atom_atom":
        return AtomPairOrbitSignature(groups, pair_group)
    if mode == "bond_bond":
        return BondPairOrbitSignature(groups, pair_group)
    return BondAtomOrbitSignature(groups, pair_group)


def _group_pair(
    mode: PairMode, g_left: TopoGroupId, g_right: TopoGroupId
) -> tuple[TopoGroupId, TopoGroupId]:
    if mode == "bond_atom":
        return (g_left, g_right)
    return cast(
        tuple[TopoGroupId, TopoGroupId],
        tuple(sorted((g_left, g_right))),
    )


def _end_group_id(mol: Mol, mode: PairMode, idx: int, *, end: int) -> TopoGroupId:
    if mode == "atom_atom":
        return TopoGroupId(mol.xf.topol_equiv[idx])
    if mode == "bond_bond":
        return TopoGroupId(_bond_group_ids(mol)[idx])
    if end == 0:
        return TopoGroupId(_bond_group_ids(mol)[idx])
    return TopoGroupId(mol.xf.topol_equiv[idx])


def _both_ends_multi(mol: Mol, mode: PairMode, left: int, right: int) -> bool:
    if mode == "atom_atom":
        return _atom_group_size(mol, left) > 1 and _atom_group_size(mol, right) > 1
    if mode == "bond_bond":
        return _bond_group_size(mol, left) > 1 and _bond_group_size(mol, right) > 1
    return _bond_group_size(mol, left) > 1 and _atom_group_size(mol, right) > 1


def _atom_group_size(mol: Mol, idx: int) -> int:
    ranks = mol.xf.topol_equiv
    rank = ranks[idx]
    return sum(1 for value in ranks.values() if value == rank)


def _bond_class_key(mol: Mol, bond_idx: int) -> tuple[int, int]:
    ranks = mol.xf.topol_equiv
    bond = mol.GetBondWithIdx(bond_idx)
    return _sorted_pair(ranks[bond.GetBeginAtomIdx()], ranks[bond.GetEndAtomIdx()])


def _bond_group_size(mol: Mol, bond_idx: int) -> int:
    key = _bond_class_key(mol, bond_idx)
    return sum(
        1 for b in range(mol.GetNumBonds()) if _bond_class_key(mol, b) == key
    )


def _bond_group_ids(mol: Mol) -> dict[int, int]:
    """Dense int ids for bond topological classes (from endpoint atom ranks)."""

    structure = _structure(mol)
    cached = structure.get("bond_topeqiv")
    if cached is not None:
        return cast(dict[int, int], cached)
    key_to_id: dict[tuple[int, int], int] = {}
    classes: dict[int, int] = {}
    for bond_idx in range(mol.GetNumBonds()):
        key = _bond_class_key(mol, bond_idx)
        if key not in key_to_id:
            key_to_id[key] = len(key_to_id)
        classes[bond_idx] = key_to_id[key]
    structure["bond_topeqiv"] = classes
    return classes


def _ensure_nauty_tables(mol: Mol) -> SitePairOrbitTables:
    structure = _structure(mol)
    cached = structure.get(_CACHE_NAUTY)
    if cached is not None:
        return cast(SitePairOrbitTables, cached)
    orbits = site_pair_orbits_nauty(mol)
    tables = _nested_tables_from_groups(mol, orbits, include_stereo=True)
    structure[_CACHE_NAUTY] = tables
    return tables


def _ensure_smiles_tables(mol: Mol) -> SitePairOrbitTables:
    structure = _structure(mol)
    cached = structure.get(_CACHE_SMILES)
    if cached is not None:
        return cast(SitePairOrbitTables, cached)
    # Only multi-group ends: singleton pairs never need pair_group materialization.
    atom_indices = tuple(
        i for i in range(mol.GetNumAtoms()) if _atom_group_size(mol, i) > 1
    )
    bond_indices = tuple(
        i for i in range(mol.GetNumBonds()) if _bond_group_size(mol, i) > 1
    )
    orbits = site_pair_orbits_smiles(
        mol, atom_indices=atom_indices, bond_indices=bond_indices
    )
    tables = _nested_tables_from_groups(mol, orbits, include_stereo=True)
    structure[_CACHE_SMILES] = tables
    return tables


def _nested_tables_from_groups(
    mol: Mol,
    orbits: Mapping[PairMode, Sequence[OrbitMembership]],
    *,
    include_stereo: bool = True,
) -> SitePairOrbitTables:
    """Build mode -> (ga, gb) -> {(a, b): PairGroupId}.

    Groups are numbered ``0..n-1`` by sorting CIP membership keys within each
    mode. Nauty and isotope agree on these ids when partitions agree.
    """

    cip = cip_ids(mol, include_stereo=include_stereo)
    tables: SitePairOrbitTables = {}
    for mode, groups in orbits.items():
        unordered = mode != "bond_atom"
        # Membership may already be sorted; re-normalize for a stable key.
        memberships = [
            orbit_membership(mol, cip, mode, group) for group in groups
        ]
        # CIP group key primary; membership tuple breaks ties when CIP collides.
        memberships.sort(
            key=lambda m: (orbit_group_cip_key(mol, cip, mode, m), m)
        )
        by_groups: PairOrbitByGroups = defaultdict(dict)
        for gid, membership in enumerate(memberships):
            pair_group = PairGroupId(gid)
            for pair in membership:
                idx_pair = (
                    _sorted_pair(pair[0], pair[1]) if unordered else pair
                )
                g_left = _end_group_id(mol, mode, idx_pair[0], end=0)
                g_right = _end_group_id(mol, mode, idx_pair[1], end=1)
                by_groups[_group_pair(mode, g_left, g_right)][idx_pair] = pair_group
        tables[mode] = dict(by_groups)
    return tables


def _structure(mol: Mol) -> dict[str, Any]:
    from xenosite.refactor_poc.rdkitutil import _structure as forest_structure

    return cast(dict[str, Any], forest_structure(mol))


def _index_tuple(indices: Iterable[int] | None, n: int) -> tuple[int, ...]:
    if indices is None:
        return tuple(range(n))
    return tuple(indices)


def _sorted_pair(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left <= right else (right, left)


def _colored_graph(
    mol: Mol, nauty: Any, *, include_stereo: bool = True
) -> tuple[Any, int, int]:
    marked = Mol(mol)
    if include_stereo:
        AssignStereochemistry(marked, cleanIt=True, force=True)

    n_atoms = marked.GetNumAtoms()
    n_bonds = marked.GetNumBonds()
    n_vertices = n_atoms + n_bonds
    adjacency: dict[int, set[int]] = {vertex: set() for vertex in range(n_vertices)}
    labels: list[tuple[object, ...]] = [()] * n_vertices

    for atom in marked.GetAtoms():
        index = atom.GetIdx()
        stereo = (
            atom.GetProp("_CIPCode")
            if include_stereo and atom.HasProp("_CIPCode")
            else None
        )
        labels[index] = (
            "atom",
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            atom.GetIsotope(),
            atom.GetNumRadicalElectrons(),
            atom.GetIsAromatic(),
            atom.GetTotalNumHs(),
            stereo,
        )

    for bond in marked.GetBonds():
        bond_vertex = n_atoms + bond.GetIdx()
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        labels[bond_vertex] = (
            "bond",
            str(bond.GetBondType()),
            bond.GetIsAromatic(),
            str(bond.GetStereo()) if include_stereo else None,
        )
        adjacency[begin].add(bond_vertex)
        adjacency[bond_vertex].add(begin)
        adjacency[end].add(bond_vertex)
        adjacency[bond_vertex].add(end)

    buckets: dict[tuple[object, ...], set[int]] = defaultdict(set)
    for vertex, label in enumerate(labels):
        buckets[label].add(vertex)
    graph = nauty.Graph(
        number_of_vertices=n_vertices,
        directed=False,
        adjacency_dict=adjacency,
        vertex_coloring=list(buckets.values()),
    )
    return graph, n_atoms, n_bonds


def _pynauty_available() -> bool:
    try:
        importlib.import_module("pynauty")
    except ImportError:
        return False
    return True
