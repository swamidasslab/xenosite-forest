"""Records a reader uses to learn the problem. Still dicts and tuples."""

from __future__ import annotations

from collections.abc import Mapping
from typing import (
    TYPE_CHECKING,
    Literal,
    NamedTuple,
    NewType,
    NotRequired,
    Protocol,
    TypeAlias,
    TypedDict,
    TypeGuard,
)

if TYPE_CHECKING:
    # rules.py imports Effect, PatternInfo, SiteInfo, and When from this module.
    from xenosite.forest.rdkit_api import Mol
    from xenosite.forest.rules import ReactionRule

_ATOMIC_NUMBER = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16}


class AtomRef(NamedTuple):
    """An atom a later step needs. It may not exist yet.

    ``idx`` is an index in the ``depth`` frame. Depth 0 is the input
    molecule. ``element`` is the atom this names.

    ``AtomRef(0, "O")`` is the oxygen added to the first atom of the input.
    ``AtomRef(0, "H", 3)`` is index 0 in the depth-3 frame, resolved once
    that atom has an H added.
    """

    idx: int
    element: str
    depth: int = 0

    def resolve(self, mol: Mol) -> int:
        """Current index of this atom. Fails until the atom exists."""

        forest = getattr(mol, "_forest", None)
        if not forest or "atom_trace" not in forest:
            raise KeyError("atom_trace")
        trace = forest["atom_trace"]
        number = _ATOMIC_NUMBER.get(self.element)
        if number is None:
            raise KeyError(self.element)
        wanted = []
        for transform_id, detail in trace.get("additions", {}).items():
            if int(detail.get("depth", 0)) != self.depth:
                continue
            site = detail.get("site")
            if site is None:
                continue
            if self.idx not in _flat_ints(site):
                continue
            wanted.append(transform_id)
        for record in trace.get("records", {}).values():
            if record.get("added_by") not in wanted:
                continue
            current = record["idx"][-1]
            atom = mol.GetAtomWithIdx(current)
            if atom.GetAtomicNum() == number:
                return current
        raise KeyError((self.idx, self.element, self.depth))


def _flat_ints(site: Site) -> set[int]:
    """Flatten a known-index :class:`Site`. Resolve AtomRef leaves first."""

    if isinstance(site, int):
        return {site}
    if isinstance(site, tuple):
        return set(site)
    found: set[int] = set()
    for item in site:
        if isinstance(item, int):
            found.add(item)
        else:
            found |= _flat_ints(item)
    return found


# Chemistry string brands (NewType → still ``str`` at runtime).
Smarts = NewType("Smarts", str)
"""Substructure match pattern for RDKit ``MolFromSmarts``.

Match-only: no reaction arrow. Used for endpoint queries and the reactant
side of a :class:`Smirks` after ``split(">>")``.
"""

Smirks = NewType("Smirks", str)
"""Reaction transform ``reactant>>product`` (SMIRKS / reaction SMARTS).

Parsed by RDKit ``ReactionFromSmarts``. Rule ``smirks`` fields hold these;
pair ``endpoints`` hold :class:`Smarts` instead.
"""

# Known atom indices only. AtomRef is a leaf, not a site by itself.
# site_kind → arm: atom → AtomSite; bond → BondSite; atom_pair → AtomPairSite
# (same frozenset[int], different semantics); directed_bond discovery →
# DirectedBondSite; pair-of-pairs → BondPairSite.
AtomSite: TypeAlias = int
DirectedBondSite: TypeAlias = tuple[int, ...]
BondSite: TypeAlias = frozenset[int]
AtomPairSite: TypeAlias = frozenset[int]
BondPairSite: TypeAlias = frozenset[frozenset[int]]
Site: TypeAlias = (
    AtomSite | DirectedBondSite | BondSite | AtomPairSite | BondPairSite
)


def _is_atom_site(value: object) -> TypeGuard[AtomSite]:
    """``AtomSite``: one known atom index. ``bool`` counts, same as ``int``."""

    return isinstance(value, int)


def _is_directed_bond_site(value: object) -> TypeGuard[DirectedBondSite]:
    """``DirectedBondSite``: map-ordered ends."""

    return isinstance(value, tuple) and all(isinstance(item, int) for item in value)


def _is_index_set_site(value: object) -> TypeGuard[BondSite]:
    """``BondSite`` or ``AtomPairSite``. Both are ``frozenset[int]``.

    An empty frozenset matches this, not :func:`_is_bond_pair_site`.
    """

    return isinstance(value, frozenset) and all(isinstance(item, int) for item in value)


def _is_bond_pair_site(value: object) -> TypeGuard[BondPairSite]:
    """``BondPairSite``: a frozenset of int-frozensets. Empty is an index set."""

    return (
        isinstance(value, frozenset)
        and bool(value)
        and all(
            isinstance(item, frozenset) and all(isinstance(inner, int) for inner in item)
            for item in value
        )
    )

# Same nesting as Site. Top-level single deferred atom is AtomRef only
# (bare int belongs to Site). Nested leaves may be int or AtomRef.
FutureAtomSite: TypeAlias = AtomRef
FutureDirectedBondSite: TypeAlias = tuple[int | AtomRef, ...]
FutureBondSite: TypeAlias = frozenset[int | AtomRef]
FutureAtomPairSite: TypeAlias = frozenset[int | AtomRef]
FutureBondPairSite: TypeAlias = frozenset[frozenset[int | AtomRef]]
FutureSite: TypeAlias = (
    FutureAtomSite
    | FutureDirectedBondSite
    | FutureBondSite
    | FutureAtomPairSite
    | FutureBondPairSite
)

AnySite: TypeAlias = Site | FutureSite


class Formula(TypedDict):
    """Heavy-atom counts, hydrogens included, and formal charge.

    Stays a dict. The forest stores it, and other code writes these keys.
    """

    counts: dict[str, int]
    charge: int


class When(TypedDict, total=False):
    """Constraint that picks one branch of a SMARTS OR once atoms are known.

    ``swap_group`` optionally overrides the resolved PatternInfo group for this
    branch only (pair unique-edit). See ``PatternInfo.swap_group``.
    """

    map: int
    z: int
    h: int
    # Optional override of PatternInfo.swap_group for this when-branch.
    # Absent / None → fall through to PatternInfo (then name). See HEURISTICS.
    swap_group: str


class Effect(TypedDict, total=False):
    """One concrete outcome.

    Declared on a possibility, then filled in from the matched atoms.
    ``symbol``, ``h``, and ``site_aromatic`` are the site atom.
    ``partner`` and ``partner_h`` are the atom the SMARTS OR was ambiguous about.
    A filter or a later resolver reads ``leave_count`` and ``breaks_ring``.
    Small-ion cleavage and N-dealkylation are the same ``cleaves`` bit plus
    different ``leave_count``. Ring opening is ``breaks_ring``, not a separate
    rule class in the search.

    ``delta_formula`` is the declared net change as element → signed count
    (zeros omitted). Derived from ``adds`` / ``removes`` bags when not set;
    When branches that disagree (halogen removal) each carry their own map.
    """

    adds: str
    removes: str
    delta_formula: dict[str, int]
    cleaves: bool
    leave_count: int | None
    breaks_ring: bool
    dearomatizes: bool
    methide: bool
    needs: str
    partner: str
    partner_h: int
    symbol: str
    h: int
    site_aromatic: bool
    when: When


# Span values: certain → bare; disagreeing possibilities → tuple.
_StrSpan: TypeAlias = str | tuple[str, ...]
_BoolSpan: TypeAlias = bool | tuple[bool, ...]
_IntSpan: TypeAlias = int | tuple[int, ...]
_LeaveSpan: TypeAlias = int | None | tuple[int | None, ...]


class _SpanCore(TypedDict):
    """Effect fields :func:`~xenosite.forest.rules._span` always writes.

    Every possibility from :func:`~xenosite.forest.rules.describe` carries
    these keys, so the collapsed span does too. A certain value stays bare;
    disagreeing branches become a tuple of the distinct values.
    """

    adds: _StrSpan
    removes: _StrSpan
    delta_formula: dict[str, int] | tuple[dict[str, int], ...]
    cleaves: _BoolSpan
    leave_count: _LeaveSpan
    breaks_ring: _BoolSpan
    dearomatizes: _BoolSpan
    methide: _BoolSpan
    needs: _StrSpan


class Span(_SpanCore, total=False):
    """Collapsed effect fields across a pattern's possibilities.

    Closed: only named keys. Optional keys appear when a possibility set them.
    """

    partner: _StrSpan
    partner_h: _IntSpan
    symbol: _StrSpan
    h: _IntSpan
    site_aromatic: _BoolSpan


class PatternInfo(TypedDict, total=False):
    """What a SMARTS pattern can do, before and after a match.

    ``span`` is the collapsed effect fields. A value may be bare or a tuple of
    the branches that disagree. ``name`` distinguishes this pattern from the
    others on the same rule.

    ``swap_group`` (optional): pair unique-edit ordered vs unordered. Resolved
    as When override → this field → ``name`` (``None`` / absent means use
    ``name``). Same non-empty resolved group on both ends → unordered
    (swappable); unequal → ordered, ends sorted by ``name``. Set explicitly
    only when grouping should differ from ``name``. A resolved :class:`When`
    may override via ``When["swap_group"]``. See HEURISTICS.
    """

    name: str
    possibilities: tuple[Effect, ...]
    span: Span
    edit: str
    site_map: int | tuple[int, ...]
    # Map numbers that receive an isotope label. More than one map can be
    # pinned. The label applies to the whole atom query, and those atoms
    # are written first so the match cannot land on a different atom.
    pin: tuple[int, ...]
    skip_same_rings: bool
    # Pair unique-edit group. Optional override when grouping differs from
    # ``name``; None / absent → use ``name``. Same resolved group on both
    # ends → unordered orbit; else ordered by ``name``. When may override.
    swap_group: str


class _SiteInfoCore(TypedDict):
    """Fields every ``filter_sites`` bag carries before an edit.

    ``site`` is the site key used for emission / unique-edit bookkeeping after
    accept. With opt-in ``canonical_emitted_sites``, that key is the lex
    orbit representative — **not** necessarily what ``filter_sites`` saw.
    When remapping ran, ``discovered_site`` holds the pre-canonical discovery
    indexes (escape hatch for callers that need discovery).

    For ``site_kind="directed_bond"``, public ``site`` is the undirected
    frozenset (only conversion: raw tuple → frozenset). ``discovered_site``
    is the raw discovery site with no conversion (the ordered map-order
    tuple). With ``canonical_emitted_sites``, ``site`` is the frozenset of
    the lex representative and ``discovered_site`` remains the raw directed
    discovery tuple.
    """

    site: Site
    # Leaf first; each containing RuleSet appends itself (same order as
    # ``addition["rules"]``). Walk for ``.name`` / identity.
    rule: list[ReactionRule]
    options: Effect
    discovered_site: NotRequired[Site]


class SmirksSiteInfo(_SiteInfoCore):
    """Bag from :meth:`SmirksReactionRule.metabolites` and
    :meth:`ResonanceRule.metabolites` for ``filter_sites``.
    """

    rxn_num: int
    pattern: PatternInfo


class PairSiteInfo(_SiteInfoCore):
    """Bag from :meth:`ResonancePairRule.pair_metabolites` for ``filter_sites``."""

    ends: tuple[Effect, Effect]
    end_atoms: tuple[int, int]
    end_maps: tuple[dict[int, int], dict[int, int]]
    path_ends: frozenset[int]


SiteInfo: TypeAlias = SmirksSiteInfo | PairSiteInfo

# metabolize yields SiteInfo unchanged. Product SMILES: ``product.xf.csmi``
# (emission identity: ``frozenset(p.xf.csmi for p in products)``).
ProductInfo: TypeAlias = SiteInfo


class TraceInfo(TypedDict, total=False):
    """Pattern snapshot on ``TraceAddition.info``.

    Rule objects become names. ``options`` and ``pattern`` are dropped; the
    resolved effect and pattern live on sibling TraceAddition fields.
    Closed: only named keys. All optional so a stub Addition may omit them.
    """

    site: Site
    discovered_site: Site
    rule: str | None
    rule_chain: tuple[str | None, ...]
    rxn_num: int
    ends: tuple[Effect, Effect]
    end_atoms: tuple[int, int]
    end_maps: tuple[dict[int, int], dict[int, int]]
    path_ends: frozenset[int]


SitesOn: TypeAlias = Literal["atom_hydrogen", "bonds", "atoms", "atom_pairs"]
# Emitted Site shape (atom-index frozenset). Class data on ReactionRule —
# not Generic[SiteT] (shared emit path + heterogeneous RuleSets erase the param;
# pyright cannot enforce frozenset size). Distinct from
# graph_isomorphism.SiteKind (atom/bond marks for orbits).
# ``"atom"`` singleton; ``"bond"`` undirected bond (two adjacent atoms);
# ``"directed_bond"`` bond Site with directed map ranks (Dealkylation);
# ``"atom_pair"`` ResonancePair path/ends only — never on plain SMARTS rules.
RuleSiteKind: TypeAlias = Literal["atom", "bond", "directed_bond", "atom_pair"]

# Keyword values :func:`~xenosite.forest.rules.describe` / ``branches`` accept.
EffectField: TypeAlias = str | bool | int | When | None


class EditCounters(Protocol):
    """Duck-typed int fields :func:`_bump` may increment during edits."""

    mol_edits: int
    rule_expansions: int
    sites_considered: int
    sites_skipped: int
    sanitize_dropped: int


class Addition(NamedTuple):
    """One transform, when a function returns it. Callers use attributes."""

    site: Site
    rules: tuple[ReactionRule, ...]
    info: TraceInfo
    effect: Effect
    name: str | None
    phase1: None
    depth: int
    pattern: PatternInfo | None = None
    discovered_site: Site | None = None


class McsResult(NamedTuple):
    """Every full-size embedding of one query, not only the best score."""

    embeddings: tuple[tuple[int, ...], ...]


class FragmentSplit(NamedTuple):
    """Pieces of one split. This module does not name a molecule type."""

    pieces: tuple[Mol, ...]


class AtomRecord(TypedDict, total=False):
    """One heavy-atom tag.

    ``idx`` and ``depth`` grow together. The last index is the atom now.
    ``added_by`` and ``removed_by`` are transform ids such as ``R1``.
    """

    idx: list[int]
    depth: list[int]
    added_by: str
    removed_by: str


class TraceAddition(TypedDict):
    """The dict stored at ``atom_trace["additions"][id]``.

    Same fields as :class:`Addition`. Stored as a dict because the trace
    writes it that way. ``info`` is the pattern snapshot with rule objects
    replaced by names, and ``options`` dropped. ``pattern`` is the
    :class:`PatternInfo` that fired, the same object the rule holds. It is
    not a rule and is not on ``rules``. ``phase1`` is reserved; only ``None``
    until its schema is decided.

    ``site`` is the emission site key (lex representative when opt-in
    canonicalization ran). ``discovered_site`` is optional: pre-canonical
    discovery indexes when they differ from ``site``.
    """

    site: Site
    rules: tuple[ReactionRule, ...]
    info: TraceInfo
    effect: Effect
    name: str | None
    phase1: None
    depth: int
    pattern: PatternInfo | None
    discovered_site: NotRequired[Site]


class AtomTrace(TypedDict):
    """Schema of ``_forest["atom_trace"]``. Not a second object beside the forest.

    ``records`` is the heavy atoms still in the molecule, keyed by tag.
    ``deletes`` is the same record after ``removed_by`` is set.
    ``transforms`` is the id order. ``next_transform`` is the next ``R`` number.
    ``depth`` is how many transforms this molecule is from the root.
    ``last_tag`` is the last tag integer issued.
    """

    records: dict[str, AtomRecord]
    deletes: dict[str, AtomRecord]
    transforms: list[str]
    additions: dict[str, TraceAddition]
    formula: Formula
    delta_formula: dict[str, Formula]
    depth: int
    last_tag: int
    next_transform: int


class InitializedAtomTrace(TypedDict):
    """The ``atom_trace`` dict ``xf.tracing._stamp`` writes.

    Same fields as :class:`AtomTrace`. Every key is present. Not a second schema.
    """

    records: dict[str, AtomRecord]
    deletes: dict[str, AtomRecord]
    transforms: list[str]
    additions: dict[str, TraceAddition]
    formula: Formula
    delta_formula: dict[str, Formula]
    depth: int
    last_tag: int
    next_transform: int


class KekuleParents(TypedDict, total=False):
    """Partial kekulé parents. Helpers fill this dict. A rule stores it.

    One live mol per distinct assignment of one conjugated system. Other
    systems on that mol stay aromatic. ``orders`` is the bond orders written
    on that parent, same index as ``parents``. ``systems`` maps the system's
    atoms to those indexes. ``by_order`` maps ``(bond, order)`` to the parent
    index that has that bond order.
    """

    parents: list[Mol]
    orders: list[dict[tuple[int, int], float]]
    systems: dict[frozenset[int], tuple[int, ...]]
    by_order: dict[tuple[tuple[int, int], float], int]


class EndParents(NamedTuple):
    """Kekulé parents for two atoms.

    ``same_system`` is true when both atoms are in one conjugated system.
    ``parents`` is then that system's assignments. When they are not,
    ``parents`` is each system's assignments, not a product of every system.
    """

    parents: tuple[Mol, ...]
    same_system: bool


# --- Pair-orbit unique-edit signatures (graph_isomorphism) ---
#
# Unique-edit keys, **not** :data:`Site` definitions. Sites name concrete
# atom indexes (atom vs atom_pair via :data:`RuleSiteKind`); these name isomorphism
# classes of an *edit*.
#
# Core shape: (groups=(ga, gb), pair_group_id).
# ``pair_group_id`` is a sequential int (``PairGroupId``) from CIP-sorted
# orbit membership tuples, or ``TRIVIAL_PAIR_GROUP`` when either end is a
# singleton topeqiv class. Never a bare atom/bond index.
#
# Ordered vs unordered is **tagged on the type**, not only which nauty family
# was queried. ``ordered`` is part of NamedTuple equality:
#   - ``ordered=False``: ends are swappable (shared resolved ``swap_group``);
#     ``end_ranks`` is ``()``.
#   - ``ordered=True``: ends have distinct roles; ``end_ranks`` holds site
#     topeqiv ranks in canonical ``PatternInfo.name`` order.
#
# Atom–atom / bond–bond only. Constructors: ``unordered_*`` / ``ordered_*``
# in ``graph_isomorphism``. (Directed bond–atom UniqueOrbit: DROPPED.md.)

TopoGroupId = NewType("TopoGroupId", int)
"""Topological equivalence class id (atom topeqiv or bond class). Not an atom index."""

PairGroupId = NewType("PairGroupId", int)
"""Orbit id numbered 0..n-1 by CIP-sorted membership. Not an atom index."""

# Sorted for atom–atom and bond–bond.
AtomGroupPair: TypeAlias = tuple[TopoGroupId, TopoGroupId]
BondGroupPair: TypeAlias = tuple[TopoGroupId, TopoGroupId]

# Recipe-level keys (partitioning / profiling). Tables and signatures use
# ``PairGroupId``. Orbit membership is always a sorted tuple of index pairs.
SmilesPairGroup: TypeAlias = str
NautyPairGroup: TypeAlias = tuple[tuple[int, int], ...]
OrbitMembership: TypeAlias = tuple[tuple[int, int], ...]

# CIP sort-key shapes (atom before bond by convention):
#   atom site:  ("atom", cip)
#   bond site:  ("bond", cip_lo, cip_hi)   # sorted endpoint CIPs
#   site pair:  sorted (site_a, site_b) for same-kind pairs
#   group:      sorted tuple of site-pair keys → numbered to PairGroupId
AtomSiteCipKey: TypeAlias = tuple[Literal["atom"], int]
BondSiteCipKey: TypeAlias = tuple[Literal["bond"], int, int]
SiteCipKey: TypeAlias = AtomSiteCipKey | BondSiteCipKey
SitePairCipKey: TypeAlias = tuple[SiteCipKey, SiteCipKey]
OrbitGroupCipKey: TypeAlias = tuple[SitePairCipKey, ...]


class AtomPairOrbitSignature(NamedTuple):
    """Atom–atom unique-edit orbit. ``ordered`` is part of equality."""

    groups: AtomGroupPair
    pair_group: PairGroupId
    ordered: Literal[True, False]
    end_ranks: tuple[int, ...]  # () if unordered; (ra, rb) name-order if ordered


class BondPairOrbitSignature(NamedTuple):
    """Bond–bond unique-edit orbit. ``ordered`` is part of equality."""

    groups: BondGroupPair
    pair_group: PairGroupId
    ordered: Literal[True, False]
    end_ranks: tuple[int, ...]  # () if unordered; (ra, rb) name-order if ordered


# Unique-edit orbit identity for one emission. Not a Site; see class docs.
PairOrbitSignature: TypeAlias = AtomPairOrbitSignature | BondPairOrbitSignature

# Map embedding identity: sorted (mapno, topological-rank) pairs.
MapRankKey: TypeAlias = tuple[tuple[int, int], ...]
# One end's role in a pair unique-edit key: (swap_group|name, map ranks).
PairRoleKey: TypeAlias = tuple[str, MapRankKey]


class PairSiteSignature(NamedTuple):
    """Unique-edit dedup key for a ResonancePair emission (before mol edit).

    ``roles`` is two ``(label, map_rank_key)`` ends: sorted when swappable
    (shared ``swap_group``), else canonical ``PatternInfo.name`` order.
    """

    roles: tuple[PairRoleKey, PairRoleKey]
    site_ranks: tuple[int, ...]
    path_end_ranks: tuple[int, ...]
    adds: str
    removes: str
    cleaves: bool
    dearomatizes: bool
    methide: bool
    orbit: PairOrbitSignature | None


# Dedup key for one SMARTS site across Kekulé forms / equivalent carbons.
# First field: directed ``MapRankKey`` for ``atom`` / ``directed_bond``, or
# sorted site ranks (undirected bond ends) for ``bond``.
# Last field: pair-orbit signature, or None for a one-atom site.
# Built in ``graph_isomorphism.site_signature``; declared here so records own
# the schema the same way ``PairSiteSignature`` does.
BondRankKey: TypeAlias = tuple[int, ...]
SiteMapKey: TypeAlias = MapRankKey | BondRankKey
SiteSignature: TypeAlias = tuple[
    SiteMapKey,
    tuple[tuple[int, int, float], ...],
    int,
    str | None,
    str | None,
    bool,
    bool,
    PairOrbitSignature | None,
]

# Forest nested tables: (ga, gb) -> {(end_a, end_b): pair_group_id}
PairOrbitSlice: TypeAlias = dict[tuple[int, int], PairGroupId]
PairOrbitByGroups: TypeAlias = dict[tuple[TopoGroupId, TopoGroupId], PairOrbitSlice]
SitePairOrbitTables: TypeAlias = dict[
    Literal["atom_atom", "bond_bond"], PairOrbitByGroups
]


class Structure(TypedDict, total=False):
    """Schema of ``_forest["cache"]``. Answers about one molecular structure.

    ``total`` is false because each key is filled in later. Structure-dependent
    / ephemeral: :func:`~xenosite.forest.rdkitutil.forest_copy` drops this
    key when the destination is a different structure, and keeps it by identity
    when copying onto the same structure.

    ``mcs_matches`` and ``mcs_targets`` are keyed by the target canonical SMILES.
    ``kekule_parents`` is the dict :class:`KekuleParents` helpers fill. The
    helpers do not touch ``_forest``. The rule stores the dict on the forest.
    """

    sanitized: int
    topol_equiv: dict[int, int]
    csmi: str
    formula: Formula
    smarts_matches: dict[Smarts, tuple[dict[int, int], ...]]
    resonance_bonds: tuple[dict[tuple[int, int], float], ...]
    kekule_parents: KekuleParents
    conjugated_systems: tuple[frozenset[int], ...]
    aromatic_systems: tuple[frozenset[int], ...]
    rings: dict[int, tuple[tuple[int, ...], ...]]
    mcs_matches: dict[str, McsResult]
    mcs_targets: dict[str, McsResult]
    site_pair_orbits_nauty: SitePairOrbitTables
    site_pair_orbits_smiles: SitePairOrbitTables
    bond_topeqiv: dict[int, int]
    cip_ids: tuple[int, ...]
    cip_ids_stereo: tuple[int, ...]
    lexical_orbit_representatives: object


class ImmutableForest(TypedDict, total=False):
    """``_forest["immutable"]`` payload. Only immutable types downward.

    Runtime value is a ``MappingProxyType`` wrapping this shape. Nested maps
    (e.g. ``start_labels``) are also ``MappingProxyType``. Allowed leaves:
    ``str``, ``int``, ``float``, ``bool``, ``None``, ``tuple``, ``frozenset``,
    and immutable mappings. Pyright sees ``Mapping`` rather than ``dict``.

    ``start_labels`` maps depth-0 atom index → CX ``atomLabel`` captured when
    the forest is first established. Stamp/check restamps those props onto the
    surviving start atoms so reaction finishing cannot drop them.
    """

    start_labels: Mapping[int, str]


class Forest(TypedDict, total=False):
    """The molecule's ``_forest``. ``total`` is false because each key is filled in later.

    Three layers:

    - ``immutable`` — frozen after write; :func:`~xenosite.forest.rdkitutil.forest_copy`
      shallow-copies the map (shared nested values).
    - ``cache`` — structure-dependent ephemeral answers (renamed from ``structure``).
    - everything else — mutable tree deep-copied without duplicating rules.
    """

    immutable: ImmutableForest
    cache: Structure
    atom_trace: AtomTrace
    is_terminal_product: bool


class UntracedForest(TypedDict, total=False):
    """A forest with no ``atom_trace`` key. ``cache`` is filled in when needed."""

    immutable: ImmutableForest
    cache: Structure
    is_terminal_product: bool


class TracingForest(TypedDict):
    """A forest whose trace is initialized. Same keys :class:`InitializedAtomTrace` names."""

    immutable: ImmutableForest
    cache: Structure
    atom_trace: InitializedAtomTrace

