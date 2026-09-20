"""Site-pair orbits under the molecular automorphism group.

**Source of truth:** :func:`all_site_pair_orbits_nauty` — six families via
colored pynauty graphs (atoms + bond-as-vertex). RDKit isotope recipes remain
callable for profiling but are not the unique-edit authority.

Families::

    atom_atom_unordered / atom_atom_ordered
    atom_bond_unordered / atom_bond_ordered   # identical partitions (types differ)
    bond_bond_unordered / bond_bond_ordered

Same-kind unordered uses combinations + sorted images; same-kind ordered uses
permutations. Asymmetric PatternInfo end roles → ordered family; symmetric →
unordered. ``atom_bond`` ordered/unordered share one partition.

Legacy :func:`site_pair_orbits_nauty` maps ``atom_atom`` / ``bond_bond`` to
the unordered families (unique-edit tables). Cross-kind ``atom_bond`` stays in
the six-family API only — not a Site pattern and not a UniqueOrbit mode.

Unique-edit pair signature (not a Site)::

    ((ga, gb), pair_group_id)

Same-kind pairs use ``swap_group`` ordered/unordered. Dehydrogenation and other
pair Sites use unordered atom–atom orbits on the two endpoint atoms.

Ordered vs unordered for ResonancePair unique-edit reads
:func:`resolved_swap_group` (When → PatternInfo → ``name``). Map-rank
embeddings and formula bags live here too; rule chemistry stays in
``rules`` / ``records``.

Public access is ``mol.xf.atom_pair_orbit_key`` / ``bond_pair_orbit_key`` /
``site_pair_orbits`` / ``pair_orbit_backend``, plus :func:`site_signature` /
:func:`pair_site_signature`.
"""


from __future__ import annotations

import importlib
import os
from collections import defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from itertools import (
    combinations,
    permutations,
    product,
)
from typing import Final, Literal, NamedTuple, Protocol, TypeAlias, cast

from xenosite.forest.rdkit_api import AssignStereochemistry, Mol, MolToSmiles
from xenosite.forest.rdkitutil import cip_ids
from xenosite.forest.records import (
    AtomPairOrbitSignature,
    AtomSiteCipKey,
    BondPairOrbitSignature,
    BondSiteCipKey,
    Effect,
    NautyPairGroup,
    OrbitGroupCipKey,
    OrbitMembership,
    PairGroupId,
    PairOrbitByGroups,
    PairOrbitSignature,
    PairRoleKey,
    PairSiteInfo,
    PairSiteSignature,
    PatternInfo,
    SiteCipKey,
    SitePairCipKey,
    SitePairOrbitTables,
    SiteSignature,
    SmilesPairGroup,
    Structure,
    TopoGroupId,
)

# Legacy same-kind mode names used by unique-edit / forest cache tables.
PairMode = Literal["atom_atom", "bond_bond"]
# Six-family nauty API (source of truth). atom_bond ordered ≡ unordered.
OrbitFamily = Literal[
    "atom_atom_unordered",
    "atom_atom_ordered",
    "atom_bond_unordered",
    "atom_bond_ordered",
    "bond_bond_unordered",
    "bond_bond_ordered",
]
SiteKind = Literal["atom", "bond"]
MarkedSite = tuple[SiteKind, int]
# Dispatcher backends. ``smiles`` is RDKit isotope marking — implemented and
# callable, but not the default when pynauty is absent (does not improve speed).
# Unique-edit prefers nauty; isotope is not the source of truth for groups.
PairOrbitBackend = Literal["nauty", "smiles", "none"]

# Re-export schema types owned by records (call sites historically imported here).
__all__ = (
    "SiteSignature",
    "PairSiteSignature",
    "PairOrbitBackend",
    "PairMode",
    "OrbitFamily",
)

PAIR_MODES: dict[PairMode, tuple[SiteKind, SiteKind]] = {
    "atom_atom": ("atom", "atom"),
    "bond_bond": ("bond", "bond"),
}

_DEFAULT_MODES: tuple[PairMode, ...] = ("atom_atom", "bond_bond")
_ORBIT_FAMILIES: tuple[OrbitFamily, ...] = (
    "atom_atom_unordered",
    "atom_atom_ordered",
    "atom_bond_unordered",
    "atom_bond_ordered",
    "bond_bond_unordered",
    "bond_bond_ordered",
)
_CACHE_NAUTY = "site_pair_orbits_nauty"
_CACHE_NAUTY_FAMILIES = "site_pair_orbits_nauty_families"
_CACHE_SMILES = "site_pair_orbits_smiles"
_CACHE_LEX_REPS = "lexical_orbit_representatives"
_ENV_BACKEND = "XENOSITE_PAIR_ORBIT_BACKEND"
_MEMBERSHIP_BASE = 4
# Documented trivial pair_group when either end is a singleton topeqiv group,
# or when backend is ``none`` for multi–multi (cheap path; no isotope tables).
# Reserved negative id so sequential orbit ids stay 0..n-1.
TRIVIAL_PAIR_GROUP: Final[PairGroupId] = PairGroupId(-1)

# Process-wide override. ``None`` → env var → auto (nauty if importable else none).
_backend_override: PairOrbitBackend | None = None

SitePairOrbitGroups = dict[OrbitFamily, list[list[tuple[int, int]]]]

# Concrete site shapes for lex-orbit emission (not unique-edit signatures).
OrbitKind = Literal["atom", "atom_atom", "bond_bond"]
AtomSite: TypeAlias = int
AtomPairSite: TypeAlias = tuple[int, int]
BondPairSite: TypeAlias = tuple[int, int]
OrbitCandidate: TypeAlias = AtomSite | AtomPairSite | BondPairSite
LexicalRepTable: TypeAlias = dict[OrbitCandidate, OrbitCandidate]


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


def canonical_emitted_sites_requested(kwargs: Mapping[str, object]) -> bool:
    """Opt-in lex-orbit site emission from explicit kwargs only.

    Default off. Pass ``canonical_emitted_sites=True`` on ``metabolize`` /
    ``metabolites`` / ``find_path`` / ``bfs`` / ``dfs`` (search splats down).
    No env or process-wide toggle.
    """

    return bool(kwargs.get("canonical_emitted_sites", False))


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

    Same-kind pairs are unordered → sorted site keys.
    """

    kind1, kind2 = PAIR_MODES[mode]
    left = site_cip_key(mol, cip, (kind1, pair[0]))
    right = site_cip_key(mol, cip, (kind2, pair[1]))
    return cast(SitePairCipKey, tuple(sorted((left, right))))


def orbit_membership(
    mol: Mol,
    cip: Sequence[int],
    mode: PairMode,
    pairs: Iterable[tuple[int, int]],
) -> OrbitMembership:
    """Sorted tuple of index pairs, ordered by CIP site-pair keys."""

    normalized = [_sorted_pair(pair[0], pair[1]) for pair in pairs]
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


class _NautyGraph(Protocol):
    """Minimal pynauty Graph surface used here."""


class _NautyModule(Protocol):
    def Graph(
        self,
        *,
        number_of_vertices: int,
        directed: bool,
        adjacency_dict: Mapping[int, set[int]],
        vertex_coloring: list[set[int]],
    ) -> _NautyGraph: ...

    def autgrp(
        self, graph: _NautyGraph
    ) -> tuple[list[Sequence[int]], object, object, object, object]: ...


def atom_bond_generators_nauty(
    mol: Mol, *, include_stereo: bool = True
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Generators as ``(atom_permutation, bond_permutation)``.

    Automorphism group via nauty/Traces (McKay & Piperno 2014). Induced
    orbits on ordered pairs are classically *orbitals* (Sharp 1999); see
    ``docs/forest/PAIR_ORBITS.md`` References.
    """

    nauty = cast(_NautyModule, importlib.import_module("pynauty"))
    graph, n_atoms, n_bonds = _colored_graph(mol, nauty, include_stereo=include_stereo)
    raw_generators, _, _, _, _ = nauty.autgrp(graph)
    generators: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for perm in raw_generators:
        atom_map = tuple(int(perm[i]) for i in range(n_atoms))
        bond_map = tuple(int(perm[n_atoms + b]) - n_atoms for b in range(n_bonds))
        generators.append((atom_map, bond_map))
    return generators


def all_site_pair_orbits_nauty(
    mol: Mol,
    *,
    include_stereo: bool = True,
    atom_indices: Iterable[int] | None = None,
    bond_indices: Iterable[int] | None = None,
) -> SitePairOrbitGroups:
    """Six-family nauty site-pair orbits (source of truth for unique-edit).

    Families::

      atom_atom_unordered : (a1, a2) with a1 < a2
      atom_atom_ordered   : (a1, a2) with a1 != a2
      atom_bond_unordered : (atom_idx, bond_idx)
      atom_bond_ordered   : same partition as unordered (types already differ)
      bond_bond_unordered : (b1, b2) with b1 < b2
      bond_bond_ordered   : (b1, b2) with b1 != b2

    ``atom_bond_ordered`` and ``atom_bond_unordered`` share one computation —
    atoms and bonds are not interchangeable same-kind sites, so ordering does
    not change the partition.

    Optional ``atom_indices`` / ``bond_indices`` restrict returned members;
    orbits are still computed on the full molecular graph.
    """

    generators = atom_bond_generators_nauty(mol, include_stereo=include_stereo)
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    all_atoms = tuple(range(n_atoms))
    all_bonds = tuple(range(n_bonds))
    selected_atoms = set(all_atoms) if atom_indices is None else set(atom_indices)
    selected_bonds = set(all_bonds) if bond_indices is None else set(bond_indices)

    same_kind_specs: dict[OrbitFamily, tuple[SiteKind, SiteKind, bool]] = {
        "atom_atom_unordered": ("atom", "atom", False),
        "atom_atom_ordered": ("atom", "atom", True),
        "bond_bond_unordered": ("bond", "bond", False),
        "bond_bond_ordered": ("bond", "bond", True),
    }

    output: SitePairOrbitGroups = {}

    for name, (kind1, kind2, ordered) in same_kind_specs.items():
        ids1 = all_atoms if kind1 == "atom" else all_bonds
        if ordered:
            candidates = list(permutations(ids1, 2))
        else:
            candidates = list(combinations(ids1, 2))
        output[name] = _orbit_groups_from_generators(
            candidates,
            generators,
            kind1=kind1,
            kind2=kind2,
            ordered=ordered,
            selected_atoms=selected_atoms,
            selected_bonds=selected_bonds,
        )

    # Cross-kind: one partition for both ordered and unordered keys.
    atom_bond_candidates = list(product(all_atoms, all_bonds))
    atom_bond_groups = _orbit_groups_from_generators(
        atom_bond_candidates,
        generators,
        kind1="atom",
        kind2="bond",
        ordered=True,  # product images stay (atom, bond); no sort
        selected_atoms=selected_atoms,
        selected_bonds=selected_bonds,
    )
    output["atom_bond_unordered"] = atom_bond_groups
    output["atom_bond_ordered"] = atom_bond_groups

    return output


def _orbit_groups_from_generators(
    candidates: list[tuple[int, int]],
    generators: list[tuple[tuple[int, ...], tuple[int, ...]]],
    *,
    kind1: SiteKind,
    kind2: SiteKind,
    ordered: bool,
    selected_atoms: set[int],
    selected_bonds: set[int],
) -> list[list[tuple[int, int]]]:
    """Union-find orbit partition under automorphism generators."""

    parent: dict[tuple[int, int], tuple[int, int]] = {p: p for p in candidates}

    def find(item: tuple[int, int]) -> tuple[int, int]:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(a: tuple[int, int], b: tuple[int, int]) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for pair in candidates:
        first, second = pair
        for atom_map, bond_map in generators:
            mapped_first = atom_map[first] if kind1 == "atom" else bond_map[first]
            mapped_second = (
                atom_map[second] if kind2 == "atom" else bond_map[second]
            )
            image = (mapped_first, mapped_second)
            if kind1 == kind2 and not ordered:
                image = _sorted_pair(image[0], image[1])
            union(pair, image)

    buckets: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    for pair in candidates:
        first, second = pair
        first_ok = (
            first in selected_atoms if kind1 == "atom" else first in selected_bonds
        )
        second_ok = (
            second in selected_atoms if kind2 == "atom" else second in selected_bonds
        )
        if first_ok and second_ok:
            buckets[find(pair)].append(pair)
    return list(buckets.values())


def site_pair_orbits_nauty(
    mol: Mol,
    modes: Sequence[PairMode] = _DEFAULT_MODES,
    bond_indices: Iterable[int] | None = None,
    atom_indices: Iterable[int] | None = None,
    *,
    include_stereo: bool = True,
) -> dict[PairMode, list[OrbitMembership]]:
    """Legacy same-kind view over :func:`all_site_pair_orbits_nauty`.

    Maps ``atom_atom`` / ``bond_bond`` → unordered families. Each group is a
    CIP-sorted tuple of index pairs. Cross-kind ``atom_bond`` lives only in
    the six-family API.
    """

    families = all_site_pair_orbits_nauty(
        mol,
        include_stereo=include_stereo,
        atom_indices=atom_indices,
        bond_indices=bond_indices,
    )
    cip = cip_ids(mol, include_stereo=include_stereo)
    result: dict[PairMode, list[OrbitMembership]] = {}
    for mode in modes:
        if mode == "atom_atom":
            raw = families["atom_atom_unordered"]
        elif mode == "bond_bond":
            raw = families["bond_bond_unordered"]
        else:
            raise ValueError(mode)
        result[mode] = [
            orbit_membership(mol, cip, mode, pairs) for pairs in raw
        ]
    return result


def atom_pair_orbit_pynauty(mol: Mol, left: int, right: int) -> NautyPairGroup:
    """Recipe 3 for one unordered atom pair."""

    groups = site_pair_orbits_nauty(mol, modes=("atom_atom",))["atom_atom"]
    needle = _sorted_pair(left, right)
    for group in groups:
        if needle in group:
            return group
    raise KeyError((left, right))


def atom_pair_orbit_pynauty_ordered(
    mol: Mol, left: int, right: int
) -> NautyPairGroup:
    """Ordered atom–atom orbit membership (asymmetric end roles)."""

    families = all_site_pair_orbits_nauty(mol)
    needle = (left, right)
    for group in families["atom_atom_ordered"]:
        if needle in group:
            return tuple(sorted(group))
    raise KeyError((left, right))


def bond_pair_orbit_pynauty(mol: Mol, bond_a: int, bond_b: int) -> NautyPairGroup:
    """Recipe 3 for one unordered bond pair."""

    groups = site_pair_orbits_nauty(mol, modes=("bond_bond",))["bond_bond"]
    needle = _sorted_pair(bond_a, bond_b)
    for group in groups:
        if needle in group:
            return group
    raise KeyError((bond_a, bond_b))


# --- Dispatch + forest cache ---
#
# pair_group_id is PairGroupId: sequential 0..n-1 from CIP-sorted membership
# tuples (same numbering for nauty and isotope when partitions agree).


def atom_pair_orbit_key(
    mol: Mol, site: frozenset[int]
) -> AtomPairOrbitSignature | None:
    """Unordered ``AtomPairOrbitSignature`` for a two-atom site."""

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
    """Unordered ``BondPairOrbitSignature`` for an unordered bond pair."""

    if len(bonds) != 2:
        return None
    left, right = tuple(bonds)
    return cast(
        BondPairOrbitSignature,
        _pair_orbit_signature(mol, "bond_bond", left, right),
    )


# --- Tagged constructors (ordered vs unordered visible at the call site) ---


def unordered_atom_pair_orbit(
    groups: tuple[TopoGroupId, TopoGroupId], pair_group: PairGroupId
) -> AtomPairOrbitSignature:
    """Atom–atom signature with ``ordered=False`` (swappable ends)."""

    return AtomPairOrbitSignature(groups, pair_group, False, ())


def ordered_atom_pair_orbit(
    groups: tuple[TopoGroupId, TopoGroupId],
    pair_group: PairGroupId,
    end_ranks: tuple[int, int],
) -> AtomPairOrbitSignature:
    """Atom–atom signature with ``ordered=True`` and name-order ``end_ranks``."""

    return AtomPairOrbitSignature(groups, pair_group, True, end_ranks)


def unordered_bond_pair_orbit(
    groups: tuple[TopoGroupId, TopoGroupId], pair_group: PairGroupId
) -> BondPairOrbitSignature:
    """Bond–bond signature with ``ordered=False``."""

    return BondPairOrbitSignature(groups, pair_group, False, ())


def ordered_bond_pair_orbit(
    groups: tuple[TopoGroupId, TopoGroupId],
    pair_group: PairGroupId,
    end_ranks: tuple[int, int],
) -> BondPairOrbitSignature:
    """Bond–bond signature with ``ordered=True`` and name-order ``end_ranks``."""

    return BondPairOrbitSignature(groups, pair_group, True, end_ranks)


# Compat names used by the unique-edit signature before the rename.
pair_orbit_isotope = atom_pair_orbit_isotope
pair_orbit_pynauty = atom_pair_orbit_pynauty
pair_orbit_signature = atom_pair_orbit_key


def _pair_orbit_signature(
    mol: Mol, mode: PairMode, left: int, right: int
) -> PairOrbitSignature:
    pair = _sorted_pair(left, right)
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
    """Low-level table lookup → unordered same-kind pair signature."""

    if mode == "atom_atom":
        return unordered_atom_pair_orbit(groups, pair_group)
    if mode == "bond_bond":
        return unordered_bond_pair_orbit(groups, pair_group)
    raise ValueError(mode)


def _group_pair(
    mode: PairMode, g_left: TopoGroupId, g_right: TopoGroupId
) -> tuple[TopoGroupId, TopoGroupId]:
    return cast(
        tuple[TopoGroupId, TopoGroupId],
        tuple(sorted((g_left, g_right))),
    )


def _end_group_id(mol: Mol, mode: PairMode, idx: int, *, end: int) -> TopoGroupId:
    if mode == "atom_atom":
        return TopoGroupId(mol.xf.topol_equiv[idx])
    if mode == "bond_bond":
        return TopoGroupId(_bond_group_ids(mol)[idx])
    raise ValueError(mode)


def _both_ends_multi(mol: Mol, mode: PairMode, left: int, right: int) -> bool:
    if mode == "atom_atom":
        return _atom_group_size(mol, left) > 1 and _atom_group_size(mol, right) > 1
    if mode == "bond_bond":
        return _bond_group_size(mol, left) > 1 and _bond_group_size(mol, right) > 1
    raise ValueError(mode)


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
    for mode_key, groups in orbits.items():
        mode: PairMode = mode_key
        unordered = True
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


def _structure(mol: Mol) -> Structure:
    from xenosite.forest.rdkitutil import _structure as forest_structure

    return forest_structure(mol)


def _index_tuple(indices: Iterable[int] | None, n: int) -> tuple[int, ...]:
    if indices is None:
        return tuple(range(n))
    return tuple(indices)


def _sorted_pair(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left <= right else (right, left)


def _colored_graph(
    mol: Mol, nauty: _NautyModule, *, include_stereo: bool = True
) -> tuple[_NautyGraph, int, int]:
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


# --- Unique-edit signature assembly (swap_group / map ranks / site keys) ---
#
# Theory-heavy helpers formerly in rules.py. PatternInfo / Effect are opaque
# Rule chemistry stays in rules; unique-edit orbits are atom–atom / bond–bond.


def resolved_swap_group(
    info: PatternInfo, effect: Effect | None = None
) -> str | None:
    """Resolved swap group: When override, else PatternInfo, else ``name``.

    Default is ``name`` so same-role pair ends are unordered without an
    explicit annotation. Set ``swap_group`` only when grouping differs from
    ``name`` (see HEURISTICS).
    """

    if effect is not None:
        when = effect.get("when")
        if when is not None and "swap_group" in when:
            group = when.get("swap_group")
            if group:
                return group
    group = info.get("swap_group")
    if group:
        return group
    name = info.get("name")
    return name if name else None


def ends_swappable(
    info1: PatternInfo,
    info2: PatternInfo,
    effect1: Effect | None = None,
    effect2: Effect | None = None,
) -> bool:
    """True → unordered pair orbit; False → ordered.

    Reads resolved ``swap_group`` (When → PatternInfo → ``name``). Unordered
    when both ends share the same non-empty group. Unequal → ordered.
    Canonical pattern order for ordered keys uses ``PatternInfo.name``.
    """

    g1 = resolved_swap_group(info1, effect1)
    g2 = resolved_swap_group(info2, effect2)
    return g1 is not None and g1 == g2


# Compat alias used in early unique-edit drafts.
end_roles_symmetric = ends_swappable


def map_rank_key(
    ranks: Mapping[int, int], mapped: Mapping[int, int]
) -> tuple[tuple[int, int], ...]:
    """Stable (mapno, topological-rank) embedding identity."""

    return tuple((mapno, ranks[idx]) for mapno, idx in sorted(mapped.items()))


def bond_rank_key(
    ranks: Mapping[int, int], site: Collection[int]
) -> tuple[int, ...]:
    """Unordered bond ends as sorted topological ranks (``site_kind="bond"``)."""

    return tuple(sorted(ranks[i] for i in site))


def formula_key(value: str | None) -> str:
    """Order-invariant formula bag for signature keys (``HCl`` ≡ ``ClH``)."""

    if not value:
        return ""
    return "".join(sorted(value))


def site_orbit(
    mol: Mol,
    mapped: Mapping[int, int],
    site: Collection[int],
) -> PairOrbitSignature | None:
    """Unique-edit orbit field for a SMARTS site (atom–atom when ``len==2``)."""

    fs = site if isinstance(site, frozenset) else frozenset(site)
    return mol.xf.atom_pair_orbit_key(fs)


def pair_orbit(
    mol: Mol,
    map1: Mapping[int, int],
    map2: Mapping[int, int],
    site_a: int,
    site_b: int,
    site: frozenset[int],
    info1: PatternInfo,
    info2: PatternInfo,
    *,
    effect1: Effect | None = None,
    effect2: Effect | None = None,
) -> PairOrbitSignature | None:
    """Orbit identity for a ResonancePairRule emission.

    :func:`ends_swappable` → unordered atom–atom; else ordered by PatternInfo
    ``name`` with name-order ``end_ranks``.
    """

    swappable = ends_swappable(info1, info2, effect1, effect2)
    key = mol.xf.atom_pair_orbit_key(site)
    if key is None:
        return None
    if swappable:
        return key  # already unordered from atom_pair_orbit_key
    # Asymmetric atom_atom: tag ordered + name-order end ranks (frozenset
    # orbit alone loses role order).
    ranks = mol.xf.topol_equiv
    name1 = info1.get("name") or ""
    name2 = info2.get("name") or ""
    end_ranks = (
        (ranks[site_a], ranks[site_b])
        if (name1, site_a) <= (name2, site_b)
        else (ranks[site_b], ranks[site_a])
    )
    return ordered_atom_pair_orbit(key.groups, key.pair_group, end_ranks)


def incident_orders(
    mol: Mol, ranks: Mapping[int, int], mapped: Mapping[int, int]
) -> tuple[tuple[int, int, float], ...]:
    """Bond orders touching the matched atoms, in rank space.

    Two Kekulé forms of the same site differ here. Two equivalent carbons
    do not, so they stay one edit.
    """

    idxs = set(mapped.values())
    bonds: list[tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in idxs and j not in idxs:
            continue
        a, b = sorted((ranks[i], ranks[j]))
        bonds.append((a, b, bond.GetBondTypeAsDouble()))
    return tuple(sorted(bonds))


def site_signature(
    context: Mol,
    work: Mol,
    mapped: Mapping[int, int],
    ranks: Mapping[int, int],
    site: Collection[int],
    rxn_num: int,
    effect: Effect,
    *,
    site_kind: str = "atom",
) -> SiteSignature:
    """Dedup key. Last field is a pair-orbit signature, or ``None`` for one atom.

    ``site_kind="bond"``: sorted site ranks (undirected ends — Epoxidation).
    ``"directed_bond"`` / ``"atom"``: directed MapRankKey (Dealkylation map 1
    is the oxygenated carbon). ``"atom_pair"`` is ResonancePair only and uses
    :func:`pair_site_signature` instead.
    """

    if site_kind == "bond":
        map_key = bond_rank_key(ranks, site)
    else:
        map_key = map_rank_key(ranks, mapped)
    return (
        map_key,
        incident_orders(work, ranks, mapped),
        rxn_num,
        effect.get("adds"),
        effect.get("removes"),
        bool(effect.get("cleaves")),
        bool(effect.get("dearomatizes")),
        site_orbit(context, mapped, site),
    )


def pair_site_signature(
    mol: Mol,
    ranks: Mapping[int, int],
    map1: Mapping[int, int],
    map2: Mapping[int, int],
    site_a: int,
    site_b: int,
    info1: PatternInfo,
    info2: PatternInfo,
    preview: PairSiteInfo,
) -> PairSiteSignature:
    """Unique-edit key for a pair-path emission (before mol edit).

    Swappable ends (shared ``swap_group``): roles + map embeddings sorted so
    match order does not matter. Ordered ends: canonical ``name`` order so
    (pat_lo@site, pat_hi@other) is stable under argument swap, while swapping
    which pattern sits on which atom stays distinct. Map ranks distinguish
    dealkylate embeddings that share path ends but cleave different partners.
    """

    effect = preview["options"]
    # Pair unique-edit site is always the two end atoms (not a nested Site union).
    site = frozenset({site_a, site_b})
    ends = preview.get("ends")
    effect1: Effect | None = ends[0] if ends else None
    effect2: Effect | None = ends[1] if ends else None
    name1 = info1.get("name") or ""
    name2 = info2.get("name") or ""
    maps1 = map_rank_key(ranks, map1)
    maps2 = map_rank_key(ranks, map2)
    swappable = ends_swappable(info1, info2, effect1, effect2)
    if swappable:
        roles = cast(
            tuple[PairRoleKey, PairRoleKey],
            tuple(
                sorted(
                    (
                        (resolved_swap_group(info1, effect1) or name1, maps1),
                        (resolved_swap_group(info2, effect2) or name2, maps2),
                    )
                )
            ),
        )
    elif (name1, site_a, maps1) <= (name2, site_b, maps2):
        roles = ((name1, maps1), (name2, maps2))
    else:
        roles = ((name2, maps2), (name1, maps1))
    return PairSiteSignature(
        roles=roles,
        site_ranks=tuple(sorted(ranks[i] for i in site)),
        path_end_ranks=tuple(sorted(ranks[i] for i in preview["path_ends"])),
        adds=formula_key(effect.get("adds")),
        removes=formula_key(effect.get("removes")),
        cleaves=bool(effect.get("cleaves")),
        dearomatizes=bool(effect.get("dearomatizes")),
        methide=bool(effect.get("methide")),
        orbit=pair_orbit(
            mol,
            map1,
            map2,
            site_a,
            site_b,
            site,
            info1,
            info2,
            effect1=effect1,
            effect2=effect2,
        ),
    )


# Underscore aliases for gradual call-site migration / tests.
_resolved_swap_group = resolved_swap_group
_ends_swappable = ends_swappable
_end_roles_symmetric = ends_swappable
_map_rank_key = map_rank_key
_formula_key = formula_key
_site_orbit = site_orbit
_pair_orbit = pair_orbit
_incident_orders = incident_orders
_site_signature = site_signature
_pair_site_signature = pair_site_signature


# --- Lexical orbit representatives (opt-in canonical emitted sites) ---


class LexicalOrbitRepresentatives(NamedTuple):
    """Concrete candidate → lex-smallest orbit member, by kind / orderedness.

    Built from atom orbits and same-kind nauty six-family groups. Cached on
    the **parent** forest ``cache`` under :data:`_CACHE_LEX_REPS` (not on
    products whose ``clear_structure`` wiped theirs).
    """

    atom: LexicalRepTable
    atom_atom_unordered: LexicalRepTable
    atom_atom_ordered: LexicalRepTable
    bond_bond_unordered: LexicalRepTable
    bond_bond_ordered: LexicalRepTable


def atom_orbit_groups_from_nauty_generators(
    n_atoms: int,
    generators: Iterable[tuple[Sequence[int], Sequence[int]]],
) -> list[list[int]]:
    """Partition atom indexes into automorphism orbits (union-find)."""

    parent_of = list(range(n_atoms))

    def find(item: int) -> int:
        while parent_of[item] != item:
            parent_of[item] = parent_of[parent_of[item]]
            item = parent_of[item]
        return item

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent_of[right] = left

    for atom_map, _bond_map in generators:
        for idx in range(n_atoms):
            union(idx, int(atom_map[idx]))

    buckets: dict[int, list[int]] = defaultdict(list)
    for idx in range(n_atoms):
        buckets[find(idx)].append(idx)
    return list(buckets.values())


def normalize_orbit_candidate(
    candidate: OrbitCandidate,
    *,
    kind: OrbitKind,
    ordered: bool,
) -> OrbitCandidate:
    """Normalize candidate representation before lookup / emission.

    ``ordered`` matters for same-kind pairs. Ignored for one atom.
    """

    if kind == "atom":
        return int(cast(AtomSite, candidate))

    if kind in ("atom_atom", "bond_bond"):
        left, right = cast(tuple[int, int], candidate)
        return (left, right) if ordered else cast(tuple[int, int], tuple(sorted((left, right))))

    raise ValueError(f"Unknown kind: {kind}")


def lexical_orbit_representatives(
    orbit_groups: Iterable[Iterable[OrbitCandidate]],
    *,
    kind: OrbitKind,
    ordered: bool = False,
) -> LexicalRepTable:
    """Build concrete candidate → canonical concrete representative.

    Each representative is the lexicographically lowest normalized candidate
    in its nauty orbit. Emitted sites become independent of discovery order,
    while remaining relative to the current RDKit atom/bond indexing.
    """

    representative_of: LexicalRepTable = {}
    for raw_members in orbit_groups:
        members = {
            normalize_orbit_candidate(member, kind=kind, ordered=ordered)
            for member in raw_members
        }
        if not members:
            continue
        representative = min(members)
        for member in members:
            representative_of[member] = representative
    return representative_of


def canonical_emitted_site(
    candidate: OrbitCandidate,
    representative_of: LexicalRepTable,
    *,
    kind: OrbitKind,
    ordered: bool = False,
) -> OrbitCandidate:
    """Convert a discovered candidate to the concrete site that should be emitted."""

    normalized = normalize_orbit_candidate(
        candidate, kind=kind, ordered=ordered
    )
    try:
        return representative_of[normalized]
    except KeyError as exc:
        raise KeyError(
            f"orbit candidate {normalized!r} (kind={kind}, ordered={ordered}) "
            "missing from representative table"
        ) from exc


def ensure_lexical_orbit_representatives(
    mol: Mol,
    *,
    parent: Mol | None = None,
) -> LexicalOrbitRepresentatives | None:
    """Materialize lex-rep tables from nauty atom / same-kind pair families.

    Cache lives on ``parent`` when given, otherwise on ``mol``. Use
    ``parent=reactant`` when remapping a product whose ``of_products`` /
    ``clear_structure`` wiped ``_forest["cache"]`` (including
    ``lexical_orbit_representatives``). Orbit computation always runs on the
    cache host — site indexes are in the reactant frame, not the product's.

    Returns ``None`` when nauty is unavailable (opt-in path then skips
    remapping and keeps discovery-order emission).
    """

    if get_pair_orbit_backend() != "nauty" or not _pynauty_available():
        return None

    host = parent if parent is not None else mol
    structure = _structure(host)
    cached = structure.get(_CACHE_LEX_REPS)
    if cached is not None:
        return cast(LexicalOrbitRepresentatives, cached)

    families = all_site_pair_orbits_nauty(host)
    generators = atom_bond_generators_nauty(host, include_stereo=True)
    atom_groups = atom_orbit_groups_from_nauty_generators(
        host.GetNumAtoms(), generators
    )

    tables = LexicalOrbitRepresentatives(
        atom=lexical_orbit_representatives(
            atom_groups, kind="atom", ordered=False
        ),
        atom_atom_unordered=lexical_orbit_representatives(
            families["atom_atom_unordered"], kind="atom_atom", ordered=False
        ),
        atom_atom_ordered=lexical_orbit_representatives(
            families["atom_atom_ordered"], kind="atom_atom", ordered=True
        ),
        bond_bond_unordered=lexical_orbit_representatives(
            families["bond_bond_unordered"], kind="bond_bond", ordered=False
        ),
        bond_bond_ordered=lexical_orbit_representatives(
            families["bond_bond_ordered"], kind="bond_bond", ordered=True
        ),
    )
    structure[_CACHE_LEX_REPS] = tables
    return tables


def _image_orbit_candidate(
    candidate: OrbitCandidate,
    atom_map: Sequence[int],
    bond_map: Sequence[int],
    *,
    kind: OrbitKind,
    ordered: bool,
) -> OrbitCandidate:
    """Apply an atom/bond automorphism to a concrete orbit candidate."""

    if kind == "atom":
        return int(atom_map[int(cast(AtomSite, candidate))])

    if kind in ("atom_atom", "bond_bond"):
        left, right = cast(tuple[int, int], candidate)
        if kind == "atom_atom":
            imaged = (atom_map[left], atom_map[right])
        else:
            imaged = (bond_map[left], bond_map[right])
        return normalize_orbit_candidate(imaged, kind=kind, ordered=ordered)

    raise ValueError(f"Unknown kind: {kind}")


def automorphism_to_representative(
    mol: Mol,
    candidate: OrbitCandidate,
    representative: OrbitCandidate,
    *,
    kind: OrbitKind,
    ordered: bool,
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    """Atom/bond maps sending ``candidate`` → ``representative``, or None."""

    normalized = normalize_orbit_candidate(
        candidate, kind=kind, ordered=ordered
    )
    target = normalize_orbit_candidate(
        representative, kind=kind, ordered=ordered
    )
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    identity_a = tuple(range(n_atoms))
    identity_b = tuple(range(n_bonds))
    if normalized == target:
        return identity_a, identity_b

    generators = atom_bond_generators_nauty(mol, include_stereo=True)
    if not generators:
        return None

    from collections import deque

    seen: set[tuple[int, ...]] = {identity_a}
    queue: deque[tuple[tuple[int, ...], tuple[int, ...]]] = deque(
        [(identity_a, identity_b)]
    )
    while queue:
        atom_map, bond_map = queue.popleft()
        for gen_a, gen_b in generators:
            # gen ∘ current
            next_a = tuple(gen_a[i] for i in atom_map)
            next_b = tuple(gen_b[i] for i in bond_map)
            if next_a in seen:
                continue
            seen.add(next_a)
            imaged = _image_orbit_candidate(
                normalized,
                next_a,
                next_b,
                kind=kind,
                ordered=ordered,
            )
            if imaged == target:
                return next_a, next_b
            queue.append((next_a, next_b))
    return None


def remap_mapped_atoms(
    mapped: Mapping[int, int], atom_map: Sequence[int]
) -> dict[int, int]:
    """Apply an atom automorphism to a SMARTS mapno → atom-index match."""

    return {mapno: int(atom_map[idx]) for mapno, idx in mapped.items()}


def is_canonical_orbit_candidate(
    candidate: OrbitCandidate,
    representative_of: LexicalRepTable,
    *,
    kind: OrbitKind,
    ordered: bool = False,
) -> bool:
    """True when ``candidate`` is already the lex representative of its orbit."""

    normalized = normalize_orbit_candidate(
        candidate, kind=kind, ordered=ordered
    )
    return representative_of.get(normalized) == normalized


def select_rep_table(
    tables: LexicalOrbitRepresentatives,
    *,
    kind: OrbitKind,
    ordered: bool,
) -> LexicalRepTable:
    """Pick the lex-rep lookup for one orbit kind / orderedness."""

    if kind == "atom":
        return tables.atom
    if kind == "atom_atom":
        return (
            tables.atom_atom_ordered if ordered else tables.atom_atom_unordered
        )
    if kind == "bond_bond":
        return (
            tables.bond_bond_ordered if ordered else tables.bond_bond_unordered
        )
    raise ValueError(f"Unknown kind: {kind}")


def _remap_site(
    site: Collection[int], atom_map: Mapping[int, int]
) -> frozenset[int] | tuple[int, ...]:
    """Apply ``atom_map`` to ``site``, preserving tuple vs frozenset."""

    remapped = tuple(int(atom_map[idx]) for idx in site)
    if isinstance(site, tuple):
        return remapped
    return frozenset(remapped)


def _copy_site(site: Collection[int]) -> frozenset[int] | tuple[int, ...]:
    if isinstance(site, tuple):
        return site
    return frozenset(site)


def _remap_match_via_auto(
    mol: Mol,
    mapped: Mapping[int, int],
    site: Collection[int],
    candidate: OrbitCandidate,
    representative: OrbitCandidate,
    *,
    kind: OrbitKind,
    ordered: bool,
) -> tuple[dict[int, int], frozenset[int] | tuple[int, ...]] | None:
    """Apply an automorphism taking ``candidate`` → ``representative`` to a match."""

    if candidate == representative:
        return dict(mapped), _copy_site(site)
    auto = automorphism_to_representative(
        mol, candidate, representative, kind=kind, ordered=ordered
    )
    if auto is None:
        return None
    atom_map, _bond_map = auto
    new_mapped = remap_mapped_atoms(mapped, atom_map)
    return new_mapped, _remap_site(site, atom_map)


def canonicalize_smarts_match(
    mol: Mol,
    mapped: Mapping[int, int],
    site: Collection[int],
    *,
    parent: Mol | None = None,
) -> tuple[dict[int, int], frozenset[int] | tuple[int, ...]] | None:
    """Remap a SMARTS match onto its lex orbit representative for emission.

    Call **after** ``filter_sites`` accepted the discovery site. Chemistry and
    the emitted ``site`` key use the returned match; stash discovery indexes
    on ``discovered_site`` when they differ. ``parent`` (when set) is the
    cache host for lex-rep tables — required when ``mol`` is a work/product
    copy whose forest cache was cleared. Returns ``None`` only when a
    representative exists but no automorphism was found (skip). With no
    nauty tables, returns the input unchanged.

    Preserves ``tuple`` vs ``frozenset`` for the site container (directed_bond
    emits ordered tuples; undirected bond emits frozensets).
    """

    host = parent if parent is not None else mol
    tables = ensure_lexical_orbit_representatives(mol, parent=host)
    if tables is None:
        return dict(mapped), _copy_site(site)

    if len(site) == 1:
        # One-atom unique-edit: lex-smallest atom in the automorphism orbit.
        site_atom = next(iter(site))
        representative = cast(
            AtomSite,
            canonical_emitted_site(
                site_atom, tables.atom, kind="atom", ordered=False
            ),
        )
        return _remap_match_via_auto(
            host,
            mapped,
            site,
            site_atom,
            representative,
            kind="atom",
            ordered=False,
        )

    if len(site) == 2:
        # One-pattern two-atom sites use unordered atom–atom orbits.
        left, right = sorted(site)
        candidate = (left, right)
        representative = cast(
            AtomPairSite,
            canonical_emitted_site(
                candidate,
                tables.atom_atom_unordered,
                kind="atom_atom",
                ordered=False,
            ),
        )
        return _remap_match_via_auto(
            host,
            mapped,
            site,
            candidate,
            representative,
            kind="atom_atom",
            ordered=False,
        )

    return dict(mapped), _copy_site(site)


def canonicalize_pair_match(
    mol: Mol,
    map1: Mapping[int, int],
    map2: Mapping[int, int],
    site_a: int,
    site_b: int,
    info1: PatternInfo,
    info2: PatternInfo,
    *,
    effect1: Effect | None = None,
    effect2: Effect | None = None,
    parent: Mol | None = None,
) -> tuple[dict[int, int], dict[int, int], int, int] | None:
    """Remap a ResonancePair match onto its lex orbit representative.

    Same contract as :func:`canonicalize_smarts_match`: remap after
    ``filter_sites`` accept for chemistry and the emitted ``site`` key.
    ``parent`` is the lex-rep cache host when ``mol``'s forest cache was
    cleared. ``None`` → skip embedding.
    """

    host = parent if parent is not None else mol
    tables = ensure_lexical_orbit_representatives(mol, parent=host)
    if tables is None:
        return dict(map1), dict(map2), site_a, site_b

    ordered = not ends_swappable(info1, info2, effect1, effect2)

    if ordered:
        name1 = info1.get("name") or ""
        name2 = info2.get("name") or ""
        candidate_aa: AtomPairSite = (
            (site_a, site_b) if (name1, site_a) <= (name2, site_b) else (site_b, site_a)
        )
    else:
        candidate_aa = cast(
            AtomPairSite,
            normalize_orbit_candidate(
                (site_a, site_b), kind="atom_atom", ordered=False
            ),
        )
    rep_table = select_rep_table(tables, kind="atom_atom", ordered=ordered)
    representative_aa = cast(
        AtomPairSite,
        canonical_emitted_site(
            candidate_aa, rep_table, kind="atom_atom", ordered=ordered
        ),
    )
    if candidate_aa == representative_aa:
        return dict(map1), dict(map2), site_a, site_b
    auto = automorphism_to_representative(
        mol,
        candidate_aa,
        representative_aa,
        kind="atom_atom",
        ordered=ordered,
    )
    if auto is None:
        return None
    atom_map, _bond_map = auto
    return (
        remap_mapped_atoms(map1, atom_map),
        remap_mapped_atoms(map2, atom_map),
        int(atom_map[site_a]),
        int(atom_map[site_b]),
    )
