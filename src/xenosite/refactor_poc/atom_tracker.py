"""Thin ``AtomTracker`` facade over ``mol.xf.tracing`` / ``mol.xf``.

Migration intent (see LOG.md / HEURISTICS.md): when the POC replaces the
live ``xenosite.forest`` tree, this class should keep the familiar names
(``tags``, ``depths``, ``initialize_tags``, ``topol_equiv``,
``atom_tag_prop_name``) so call sites swap with minimal churn. Old forest
code moves to an archive dir for a while; this facade is the seam.

Not a port of forest ``AtomTracker`` internals. Methods that have an ``xf``
equivalent proxy to it. Tag-dict shape matches forest records when a POC
forest is present (label → ``{idx, depth, ...}``).
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
from typing import Any

from xenosite.refactor_poc.rdkit_api import Mol


class AtomTracker:
    """Compatibility surface for forest ``AtomTracker`` callers.

    Prefer ``mol.xf.tracing.*`` for new POC code. This class exists so a
    later forest←POC swap can keep import and method names stable.
    """

    tag_name = "ATOM_INDEX_PATHS"
    last_tag_name = "LAST_TAG"
    previous_index_prop_name = "current_idx"
    atom_tag_prop_name = "_forestLabel"

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__()

    @staticmethod
    def topol_equiv(mol: Mol) -> dict[int, int]:
        """Topological equivalence classes via ``mol.xf.topol_equiv``."""

        return dict(mol.xf.topol_equiv)

    @staticmethod
    def site_to_topol_site(
        site: tuple[str, Any], topol_equiv: Mapping[int, int]
    ) -> tuple[str, tuple[int, ...]]:
        """Map ``(rule_name, atom_idxs)`` to ``(bare_name, sorted ranks)``."""

        name = site[0]
        if isinstance(name, str):
            name = name.split("_", 1)[0]
        ranks = tuple(sorted(topol_equiv[x] for x in site[1]))
        return name, ranks

    def initialize_tags(self, mol: Mol) -> Mol:
        """Install a depth-0 atom trace (``xf.tracing._stamp``)."""

        return mol.xf.tracing._stamp()

    @classmethod
    def tags(
        cls,
        record: Mol | Mapping | list,
        depth: int | None = None,
        idx: int | None = None,
        strict: bool = True,
        compact: bool = False,
        **kwargs: Any,
    ) -> Any:
        """Forest-shaped tag records, or a chain over a list of mols.

        For a POC ``Mol``, reads ``_forest["atom_trace"]["records"]`` when
        present (same dict shape forest used). Does not invent labels.
        """

        if isinstance(record, list):
            return itertools.chain(
                *[
                    cls.tags(x, depth=depth, idx=idx, strict=strict, compact=compact)
                    for x in record
                ]
            )
        if isinstance(record, Mol):
            forest = getattr(record, "_forest", None) or {}
            trace = forest.get("atom_trace")
            if trace is None or "records" not in trace:
                if strict:
                    raise KeyError("atom_trace")
                return {}
            record = trace["records"]
        elif not isinstance(record, dict):
            raise ValueError("Must submit RDKit Mol or dict")

        out: dict = dict(record)
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
    def compact_tags(record: Mapping, adjust_root_by: int = 1) -> dict:
        """1-based depth→idx map (SMILES / map convention)."""

        return {
            k: {d: i + adjust_root_by for d, i in zip(v["depth"], v["idx"])}
            for k, v in record.items()
        }

    @classmethod
    def depths(cls, record: Mol | Mapping, strict: bool = True) -> list[int]:
        """Sorted unique depths from tag records or ``xf.tracing.depth``."""

        if isinstance(record, Mol):
            tracing = record.xf.tracing
            if tracing.active and tracing.depth is not None:
                # Prefer full record set when present so multi-depth paths agree
                # with forest AtomTracker.depths.
                try:
                    tags = cls.tags(record, strict=True)
                except KeyError:
                    return [int(tracing.depth)]
                record = tags
            else:
                if strict:
                    raise KeyError("atom_trace")
                return []
        elif not isinstance(record, dict):
            raise ValueError("Must submit RDKit Mol or dict")

        return sorted(
            set(itertools.chain(*[x["depth"] for x in list(record.values())]))
        )

    @classmethod
    def metabolite_index_to_reversed_index_record(
        cls, metabolite: Mol, exact_depth: int = 2, strict: bool = True
    ) -> dict | None:
        """Forest-compatible reverse-index helper over compact tags."""

        depth_list = cls.depths(metabolite, strict=strict)
        if len(depth_list) < exact_depth:
            return None
        reversed_depth = list(reversed(depth_list[-exact_depth:]))
        idx_record = cls.tags(metabolite, compact=True, strict=strict)
        out: dict = defaultdict(list)
        for depth_to_idx in list(idx_record.values()):
            if set(reversed_depth) != set(depth_to_idx):
                continue
            metabolite_idx = depth_to_idx[reversed_depth[-1]]
            for depth in reversed_depth:
                out[metabolite_idx].append(depth_to_idx[depth])
        return out
