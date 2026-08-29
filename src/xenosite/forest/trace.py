"""Public 1-based atom tracing for tagged metabolites.

Atom numbers (``AtomTrace``, SMILES ``:N``, ``GetAtomMapNum()``, Phase I
``1.h``) are **1-based**. ``0`` on a map number means unmapped / new, not atom
zero.

Depths (``t.depths``, ``map(start_depth=0)``, reaction count) are **0-based**.
Depth 0 is the original reactant. These are not atom numbers.

RDKit ``GetIdx()`` is internal and 0-based. Convert with :func:`atom_no` and
:func:`rdkit_idx` only; do not call a 1-based value ``idx``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from rdkit.Chem.rdchem import Mol


def atom_no(rdkit_idx: int) -> int:
    """0-based RDKit GetIdx() -> 1-based atom number (SMILES :N)."""
    return rdkit_idx + 1


def rdkit_idx(atom_no: int) -> int:
    """1-based atom number -> 0-based RDKit GetIdx()."""
    return atom_no - 1


class AtomTrace:
    """Atom-mapping history for a tagged metabolite.

    Atom numbers are 1-based; depths are 0-based.

    Load history from a mol that was produced (or initialized) without
    ``do_not_tag_atoms=True``. Map numbers on the mol equal ``origin()`` for
    atoms that came from depth 0.
    """

    def __init__(self, mol: Mol):
        from .base import AtomTracker

        self._mol = mol
        try:
            self._tags = AtomTracker.tags(mol)
        except KeyError as err:
            raise ValueError(
                "mol has no atom-index tags; run a reaction without "
                "do_not_tag_atoms=True"
            ) from err
        self._depths = tuple(AtomTracker.depths(self._tags))

    @property
    def depths(self) -> tuple[int, ...]:
        """0-based reaction steps present on this mol. Depth 0 is the reactant."""
        return self._depths

    def _idx_at(self, rec: dict, depth: int) -> int | None:
        if depth not in rec["depth"]:
            return None
        return rec["idx"][rec["depth"].index(depth)]

    def _rec_at(self, atom: int, depth: int) -> dict | None:
        idx = rdkit_idx(atom)
        for rec in self._tags.values():
            if self._idx_at(rec, depth) == idx:
                return rec
        return None

    def map(self, start_depth: int = 0, end_depth: int | None = None) -> dict[int, int]:
        """``{atom_no at start: atom_no at end}`` for atoms present at both depths.

        Keys and values are 1-based atom numbers. Depths are 0-based.
        """
        if end_depth is None:
            if not self._depths:
                return {}
            end_depth = self._depths[-1]
        out: dict[int, int] = {}
        for rec in self._tags.values():
            start = self._idx_at(rec, start_depth)
            end = self._idx_at(rec, end_depth)
            if start is not None and end is not None:
                out[atom_no(start)] = atom_no(end)
        return out

    def follow(self, atom: int, depth: int = 0) -> tuple[int | None, ...]:
        """That atom's 1-based number at every depth (``None`` if gone).

        ``atom`` is the 1-based number at ``depth`` (0-based). Depth 0 is the
        original reactant.
        """
        rec = self._rec_at(atom, depth)
        if rec is None:
            return tuple(None for _ in self._depths)
        out: list[int | None] = []
        for d in self._depths:
            idx = self._idx_at(rec, d)
            out.append(None if idx is None else atom_no(idx))
        return tuple(out)

    def origin(self, atom: int, depth: int | None = None) -> int | None:
        """Level-0 1-based number, or ``None`` if this atom was added later.

        Equals ``GetAtomMapNum()`` when maps are stamped from depth 0.
        ``atom`` is the 1-based number at ``depth`` (last depth if omitted).
        """
        if depth is None:
            if not self._depths:
                return None
            depth = self._depths[-1]
        rec = self._rec_at(atom, depth)
        if rec is None:
            return None
        idx = self._idx_at(rec, 0)
        return None if idx is None else atom_no(idx)

    def added(self, depth: int | None = None) -> frozenset[int]:
        """1-based current numbers that are new (no depth-0 origin).

        If ``depth`` is set, only atoms whose history starts at that depth.
        Otherwise, every atom present at the last depth with no origin.
        """
        if not self._depths:
            return frozenset()
        query_depth = self._depths[-1] if depth is None else depth
        out: set[int] = set()
        for rec in self._tags.values():
            if self._idx_at(rec, 0) is not None:
                continue
            cur = self._idx_at(rec, query_depth)
            if cur is None:
                continue
            if depth is not None and rec["depth"][0] != query_depth:
                continue
            out.add(atom_no(cur))
        return frozenset(out)

    def removed(self, depth: int | None = None) -> frozenset[int]:
        """1-based level-0 numbers that have disappeared by ``depth`` (last if omitted)."""
        if depth is None:
            if not self._depths:
                return frozenset()
            depth = self._depths[-1]
        out: set[int] = set()
        for rec in self._tags.values():
            orig = self._idx_at(rec, 0)
            if orig is None:
                continue
            if self._idx_at(rec, depth) is None:
                out.add(atom_no(orig))
        return frozenset(out)
