"""Records a reader uses to learn the problem. Still dicts and tuples."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, NamedTuple, TypedDict

if TYPE_CHECKING:
    # rules.py imports Effect, PatternInfo, and When from this module.
    from xenosite.refactor_poc.rdkit_api import Mol
    from xenosite.refactor_poc.rules import ReactionRule

# One atom, a pair of atoms, or a set of those pairs.
# The frozenset[int] form has exactly two atom indices.
Site = int | frozenset[int] | frozenset[frozenset[int]]


class Formula(TypedDict):
    """Heavy-atom counts, hydrogens included, and formal charge.

    Stays a dict. The forest stores it, and other code writes these keys.
    """

    counts: dict[str, int]
    charge: int


class When(TypedDict, total=False):
    """Constraint that picks one branch of a SMARTS OR once atoms are known."""

    map: int
    z: int
    h: int


class Effect(TypedDict, total=False):
    """One concrete outcome.

    Declared on a possibility, then filled in from the matched atoms.
    ``symbol``, ``h``, and ``site_aromatic`` are the site atom.
    ``partner`` and ``partner_h`` are the atom the SMARTS OR was ambiguous about.
    A filter or a later resolver reads ``leave_count`` and ``breaks_ring``.
    Small-ion cleavage and N-dealkylation are the same ``cleaves`` bit plus
    different ``leave_count``. Ring opening is ``breaks_ring``, not a separate
    rule class in the search.
    """

    adds: str
    removes: str
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


class PatternInfo(TypedDict, total=False):
    """What a SMARTS pattern can do, before and after a match.

    ``span`` stays an open dict. It is the collapsed effect fields, and a
    value may be a bare result or a tuple of the branches that disagree.
    Naming that as one TypedDict would hide the disagreement.
    """

    possibilities: tuple[Effect, ...]
    span: dict[str, Any]
    edit: str
    site_map: int | tuple[int, ...]
    skip_same_rings: bool


class Addition(NamedTuple):
    """One transform, when a function returns it. Callers use attributes."""

    site: Site
    rules: tuple[ReactionRule, ...]
    info: dict[str, Any]
    effect: Effect
    name: str | None
    phase1: Any
    depth: int


class McsResult(NamedTuple):
    """Every full-size embedding of one query, not only the best score."""

    embeddings: tuple[tuple[int, ...], ...]


class FragmentSplit(NamedTuple):
    """Pieces of one split. This module does not name a molecule type."""

    pieces: tuple[Any, ...]


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
    writes it that way. ``info`` stays open: it is the pattern snapshot
    with rule objects replaced by names, and ``options`` dropped.
    ``phase1`` is reserved. Its schema is not decided.
    """

    site: Site
    rules: tuple[ReactionRule, ...]
    info: dict[str, Any]
    effect: Effect
    name: str | None
    phase1: Any
    depth: int


class AtomTrace(TypedDict, total=False):
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


class Structure(TypedDict, total=False):
    """Cache filled in later. ``total`` is false because each key is filled in later.

    ``mcs_matches`` and ``mcs_targets`` are keyed by the target canonical SMILES.
    ``kekule_parents`` is the dict :class:`KekuleParents` helpers fill. The
    helpers do not touch ``_forest``. The rule stores the dict on the forest.
    """

    sanitized: int
    topol_equiv: dict[int, int]
    is_terminal_product: bool
    csmi: str
    formula: Formula
    smarts_matches: dict[str, tuple[dict[int, int], ...]]
    resonance_bonds: tuple[dict[tuple[int, int], float], ...]
    kekule_parents: KekuleParents
    conjugated_systems: tuple[frozenset[int], ...]
    aromatic_systems: tuple[frozenset[int], ...]
    rings: dict[int, tuple[tuple[int, ...], ...]]
    mcs_matches: dict[str, McsResult]
    mcs_targets: dict[str, McsResult]


class Forest(TypedDict, total=False):
    """The molecule's ``_forest``. ``total`` is false because each key is filled in later."""

    structure: Structure
    atom_trace: AtomTrace
    is_terminal_product: bool
