"""Records a reader uses to learn the problem. Still dicts and tuples."""

from typing import Any, NamedTuple, TypedDict


class Formula(TypedDict):
    """Heavy-atom counts, hydrogens included, and formal charge.

    Stays a dict. The forest stores it, and other code writes these keys.
    """

    counts: dict[str, int]
    charge: int


class Addition(NamedTuple):
    """One transform. Callers that receive this use attributes."""

    site: tuple[int, ...]
    rules: tuple[str, ...]
    info: dict[str, Any]
    effect: dict[str, Any]
    name: str | None
    phase1: Any
    depth: int


class McsResult(NamedTuple):
    """Every full-size embedding of one query, not only the best score."""

    embeddings: tuple[tuple[int, ...], ...]


class FragmentSplit(NamedTuple):
    """Pieces of one split. This module does not name a molecule type."""

    pieces: tuple[Any, ...]


class Structure(TypedDict, total=False):
    """Cache of plain data. ``total`` is false because each key is filled in later."""

    sanitized: int
    topol_equiv: dict[int, int]
    is_terminal_product: bool
    csmi: str
    formula: Formula
    smarts_matches: dict[str, tuple[dict[int, int], ...]]
    resonance_bonds: tuple[dict[tuple[int, int], float], ...]
    conjugated_systems: tuple[frozenset[int], ...]
    aromatic_systems: tuple[frozenset[int], ...]
    rings: dict[int, tuple[tuple[int, ...], ...]]


class Forest(TypedDict, total=False):
    """The molecule's ``_forest``. ``total`` is false because each key is filled in later.

    ``atom_trace`` stays a plain dict. It is not its own type.
    """

    structure: Structure
    atom_trace: dict[str, Any]
    parent_atom_trace: dict[str, Any]
    is_terminal_product: bool
