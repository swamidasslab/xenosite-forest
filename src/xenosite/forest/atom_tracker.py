"""Deprecated ``AtomTracker`` facade over ``mol.xf.tracing`` / ``mol.xf``.

Kept working for call-site compatibility after the promote to ``xenosite.forest``
swap. Prefer ``mol.xf`` / ``mol.xf.tracing`` for new code. Old forest
internals live under ``xenosite._archive_forest`` for a while.

Not a port of archived ``AtomTracker`` internals. Methods that have an ``xf``
equivalent proxy to it. Tag-dict shape matches forest records when a forest
dict is present (label → ``{idx, depth, ...}``).
"""

from __future__ import annotations

import itertools
import warnings
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from copy import deepcopy
from typing import Any, TypeAlias, cast

from xenosite.forest.rdkit_api import Mol
from xenosite.forest.records import AtomRecord, Forest

TagRecordMap: TypeAlias = dict[str, AtomRecord]
CompactTagMap: TypeAlias = dict[str, dict[int, int]]
TagInput: TypeAlias = Mol | Mapping[str, AtomRecord] | list[Mol | Mapping[str, AtomRecord]]

_DEPRECATION = (
    "xenosite.forest.AtomTracker is deprecated; use mol.xf / mol.xf.tracing "
    "instead (see docs/forest/XF.md). The facade remains for compatibility "
    "and will be removed in a future release."
)


def _warn_deprecated(*, stacklevel: int = 3) -> None:
    warnings.warn(_DEPRECATION, DeprecationWarning, stacklevel=stacklevel)


class AtomTracker:
    """Deprecated compatibility surface for archived ``AtomTracker`` callers.

    Prefer ``mol.xf`` / ``mol.xf.tracing`` (migration tutorial in
    ``docs/forest/XF.md``). Instantiation and public helpers emit
    :class:`DeprecationWarning`.
    """

    tag_name = "ATOM_INDEX_PATHS"
    last_tag_name = "LAST_TAG"
    previous_index_prop_name = "current_idx"
    atom_tag_prop_name = "_forestLabel"

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        _warn_deprecated(stacklevel=2)
        super().__init__()

    @staticmethod
    def topol_equiv(mol: Mol) -> dict[int, int]:
        """Topological equivalence classes via ``mol.xf.topol_equiv``."""

        _warn_deprecated()
        return dict(mol.xf.topol_equiv)

    @staticmethod
    def site_to_topol_site(
        site: tuple[str, Iterable[int]], topol_equiv: Mapping[int, int]
    ) -> tuple[str, tuple[int, ...]]:
        """Map ``(rule_name, atom_idxs)`` to ``(bare_name, sorted ranks)``."""

        _warn_deprecated()
        name = site[0]
        if isinstance(name, str):
            name = name.split("_", 1)[0]
        ranks = tuple(sorted(topol_equiv[x] for x in site[1]))
        return name, ranks

    def initialize_tags(self, mol: Mol) -> Mol:
        """Install a depth-0 atom trace (``xf.tracing._stamp``)."""

        # ``__init__`` already warned when constructed; classmethod-style
        # use without an instance still needs a signal — callers typically
        # construct then call, so skip a second warn here.
        return mol.xf.tracing._stamp()

    @classmethod
    def tags(
        cls,
        record: TagInput,
        depth: int | None = None,
        idx: int | None = None,
        strict: bool = True,
        compact: bool = False,
        **kwargs: Any,
    ) -> TagRecordMap | CompactTagMap | Iterator[object]:
        """Forest-shaped tag records, or a chain over a list of mols.

        For a forest ``Mol``, reads ``_forest["atom_trace"]["records"]`` when
        present (same dict shape the archived tracker used). Does not invent
        labels.
        """

        _warn_deprecated()
        if isinstance(record, list):
            return itertools.chain(
                *[
                    cls.tags(x, depth=depth, idx=idx, strict=strict, compact=compact)
                    for x in record
                ]
            )
        tags: TagRecordMap
        if isinstance(record, Mol):
            forest_obj = getattr(record, "_forest", None) or {}
            forest = cast(Forest, forest_obj)
            trace = forest.get("atom_trace")
            if trace is None or "records" not in trace:
                if strict:
                    raise KeyError("atom_trace")
                return {}
            tags = dict(trace["records"])
        elif isinstance(record, Mapping):
            tags = dict(record)
        else:
            raise ValueError("Must submit RDKit Mol or dict")

        out: TagRecordMap = tags
        if depth is not None:
            out = {
                tag: data
                for tag, data in out.items()
                if depth in data.get("depth", ())
            }
        if idx is not None:
            out = {
                tag: data for tag, data in out.items() if idx in data.get("idx", ())
            }
        if idx is None and depth is None:
            out = deepcopy(out)
        if compact:
            return cls.compact_tags(out)
        return out

    @staticmethod
    def compact_tags(
        record: Mapping[str, AtomRecord], adjust_root_by: int = 1
    ) -> CompactTagMap:
        """1-based depth→idx map (SMILES / map convention)."""

        return {
            k: {
                int(d): int(i) + adjust_root_by
                for d, i in zip(list(v.get("depth", [])), list(v.get("idx", [])))
            }
            for k, v in record.items()
        }

    @classmethod
    def depths(
        cls, record: Mol | Mapping[str, AtomRecord], strict: bool = True
    ) -> list[int]:
        """Sorted unique depths from tag records or ``xf.tracing.depths``."""

        if isinstance(record, Mol):
            tracing = record.xf.tracing
            if tracing.active:
                return list(tracing.depths())
            if strict:
                raise KeyError("atom_trace")
            return []
        if isinstance(record, Mapping):
            tags = record
        else:
            raise ValueError("Must submit RDKit Mol or dict")

        return sorted(
            set(
                itertools.chain(
                    *[list(x.get("depth", [])) for x in list(tags.values())]
                )
            )
        )

    @classmethod
    def metabolite_index_to_reversed_index_record(
        cls, metabolite: Mol, exact_depth: int = 2, strict: bool = True
    ) -> dict[int, list[int]] | None:
        """Forest-compatible reverse-index helper over compact tags."""

        depth_list = cls.depths(metabolite, strict=strict)
        if len(depth_list) < exact_depth:
            return None
        reversed_depth = list(reversed(depth_list[-exact_depth:]))
        idx_record = cast(
            CompactTagMap, cls.tags(metabolite, compact=True, strict=strict)
        )
        out: dict[int, list[int]] = defaultdict(list)
        for depth_to_idx in list(idx_record.values()):
            if set(reversed_depth) != set(depth_to_idx):
                continue
            metabolite_idx = depth_to_idx[reversed_depth[-1]]
            for depth in reversed_depth:
                out[metabolite_idx].append(depth_to_idx[depth])
        return out
