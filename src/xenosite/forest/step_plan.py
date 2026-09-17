"""Partial-order plans over reaction steps (domain-agnostic).

``AtomRef`` / ``Step`` / ``Linearization`` form a reusable apply+resolve layer.
Creation bookkeeping lives on ``mol._forest["atom_refs"]`` (private), not in
AtomTracker tag strings.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterator, Sequence

from rdkit import Chem

from .base import AtomTracker, _copy_forest, _forest_state


# ---------------------------------------------------------------------------
# Rule registry (name -> ReactionRule instance)
# ---------------------------------------------------------------------------

_RULE_CACHE: dict = {}


def get_rule(name: str):
    """Return a shared ``ReactionRule`` instance for ``name``."""
    rule = _RULE_CACHE.get(name)
    if rule is None:
        from . import rules as forest_rules

        cls = getattr(forest_rules, name, None)
        if cls is None:
            raise KeyError("unknown forest rule %r" % (name,))
        rule = cls()
        _RULE_CACHE[name] = rule
    return rule


# ---------------------------------------------------------------------------
# Private atom_refs index on mol._forest
# ---------------------------------------------------------------------------


class AtomRefsIndex:
    """Maps (added_by, reactant-frame at) -> current GetIdx on a mol."""

    __slots__ = ("_entries",)

    def __init__(self, entries=None):
        # (rule_name, frozenset[int]) -> int
        self._entries = dict(entries or ())

    def copy(self) -> AtomRefsIndex:
        return AtomRefsIndex(self._entries)

    def set(self, added_by: str, at: frozenset, idx: int) -> None:
        self._entries[(added_by, frozenset(at))] = int(idx)

    def lookup(self, added_by: str, at: frozenset):
        return self._entries.get((added_by, frozenset(at)))

    def items(self):
        return self._entries.items()

    def remap_from_parent(self, parent_mol, child_mol) -> AtomRefsIndex:
        """Rebuild entries with child GetIdx via react_atom_idx from parent."""
        out = AtomRefsIndex()
        for key, parent_idx in self._entries.items():
            for atom in child_mol.GetAtoms():
                if (
                    atom.HasProp("react_atom_idx")
                    and int(atom.GetProp("react_atom_idx")) == parent_idx
                ):
                    out._entries[key] = atom.GetIdx()
                    break
        return out


def _atom_refs_index(mol) -> AtomRefsIndex:
    forest = _forest_state(mol)
    index = forest.get("atom_refs")
    if index is None:
        index = AtomRefsIndex()
        forest["atom_refs"] = index
    return index


def _install_product_forest(parent, product) -> AtomRefsIndex:
    """New ``_forest`` dict: share resonance, remap/copy atom_refs."""
    parent_forest = getattr(parent, "_forest", None) or {}
    child = {}
    if "resonance" in parent_forest:
        child["resonance"] = parent_forest["resonance"]
    prev = parent_forest.get("atom_refs")
    if prev is not None:
        child["atom_refs"] = prev.remap_from_parent(parent, product)
    else:
        child["atom_refs"] = AtomRefsIndex()
    product._forest = child
    return child["atom_refs"]


def _record_creations(step: Step, before, after, index: AtomRefsIndex) -> None:
    """Record atoms created by ``step`` onto ``index`` (keyed by origin site)."""
    origin_at = frozenset(
        ref.origin for ref in step.site if ref.origin is not None
    )
    if not origin_at:
        return

    old_on_after = set()
    for atom in after.GetAtoms():
        if atom.GetAtomMapNum() > 0 or atom.HasProp("react_atom_idx"):
            old_on_after.add(atom.GetIdx())
    new_atoms = [
        atom for atom in after.GetAtoms() if atom.GetIdx() not in old_on_after
    ]
    if not new_atoms:
        return

    chosen = None
    if len(new_atoms) == 1:
        chosen = new_atoms[0]
    else:
        targets = set()
        for origin in origin_at:
            try:
                targets.add(_resolve_origin(after, origin))
            except KeyError:
                continue
        for atom in new_atoms:
            for t in targets:
                if atom.GetIdx() != t and after.GetBondBetweenAtoms(
                    atom.GetIdx(), t
                ):
                    chosen = atom
                    break
            if chosen is not None:
                break
        if chosen is None:
            chosen = new_atoms[0]

    index.set(step.rule, origin_at, chosen.GetIdx())


def _resolve_origin(mol, origin: int) -> int:
    """Map reactant-frame GetIdx to current GetIdx (maps / react_atom_idx)."""
    origin = int(origin)
    has_maps = False
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() > 0:
            has_maps = True
            if atom.GetAtomMapNum() == origin + 1:
                return atom.GetIdx()
        if atom.HasProp("react_atom_idx"):
            has_maps = True
            if int(atom.GetProp("react_atom_idx")) == origin:
                return atom.GetIdx()
    if not has_maps and 0 <= origin < mol.GetNumAtoms():
        return origin
    raise KeyError("cannot resolve origin atom %s on mol" % (origin,))


def _ensure_apply_ready(mol) -> None:
    """Initialize atom maps on a fresh mol so later resolves can track origins."""
    if any(a.GetAtomMapNum() > 0 for a in mol.GetAtoms()):
        return
    if any(a.HasProp("react_atom_idx") for a in mol.GetAtoms()):
        return
    # Use any registered rule for AtomTracker mixin methods.
    try:
        probe = get_rule("Hydroxylation")
    except KeyError:
        return
    probe.initialize_tags(mol)
    AtomTracker.add_current_idx_as_atom_prop(
        mol, propname=AtomTracker.previous_index_prop_name
    )


def _fragment_retains_refs(product, later_refs) -> bool:
    """True if ``product`` still carries atoms needed for later AtomRefs."""
    for ref in later_refs:
        if ref.origin is not None:
            try:
                ref.resolve(product)
            except KeyError:
                return False
            continue
        if ref.added_by is not None:
            try:
                ref.resolve(product)
                continue
            except KeyError:
                for origin in ref.at or ():
                    try:
                        _resolve_origin(product, origin)
                    except KeyError:
                        return False
    return True


# ---------------------------------------------------------------------------
# AtomRef / Step / Linearization / StepPlan
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtomRef:
    """Reactant-stable atom identity (origin or created by a prior step)."""

    origin: int | None = None
    added_by: str | None = None
    at: frozenset | None = None

    def __post_init__(self):
        if self.at is not None and not isinstance(self.at, frozenset):
            object.__setattr__(self, "at", frozenset(self.at))
        if self.origin is not None and self.added_by is not None:
            raise ValueError("AtomRef cannot set both origin and added_by")
        if self.origin is None and self.added_by is None:
            raise ValueError("AtomRef needs origin or added_by")
        if self.added_by is not None and self.at is None:
            raise ValueError("added_by AtomRef requires at= reactant site")

    @classmethod
    def coerce(cls, value) -> AtomRef:
        if isinstance(value, AtomRef):
            return value
        if isinstance(value, int):
            return cls(origin=int(value))
        if isinstance(value, dict):
            return cls.from_json(value)
        raise TypeError("cannot coerce %r to AtomRef" % (value,))

    def to_json(self):
        if self.origin is not None:
            return self.origin
        return {
            "added_by": self.added_by,
            "at": sorted(self.at),
        }

    @classmethod
    def from_json(cls, data) -> AtomRef:
        if isinstance(data, int):
            return cls(origin=data)
        if "origin" in data and data["origin"] is not None:
            return cls(origin=int(data["origin"]))
        return cls(
            added_by=data["added_by"],
            at=frozenset(data["at"]),
        )

    def resolve(self, mol) -> int:
        """Current GetIdx on ``mol`` (maps / ``_forest['atom_refs']``)."""
        if self.origin is not None:
            return _resolve_origin(mol, self.origin)
        index = _atom_refs_index(mol)
        idx = index.lookup(self.added_by, self.at)
        if idx is None:
            raise KeyError(
                "cannot resolve AtomRef added_by=%r at=%r" % (self.added_by, self.at)
            )
        return idx


def _coerce_site(site) -> frozenset:
    return frozenset(AtomRef.coerce(x) for x in site)


def _site_to_json(site: frozenset):
    items = [ref.to_json() for ref in site]

    def sort_key(item):
        if isinstance(item, int):
            return (0, item)
        return (1, json.dumps(item, sort_keys=True))

    return sorted(items, key=sort_key)


def _site_from_json(items) -> frozenset:
    return frozenset(AtomRef.from_json(x) for x in items)


@dataclass(frozen=True)
class Step:
    """One named reaction at a site of :class:`AtomRef` (ints coerce to origin)."""

    rule: str
    site: frozenset

    def __post_init__(self):
        object.__setattr__(self, "site", _coerce_site(self.site))

    def resolve_site(self, mol) -> frozenset:
        """Opt-in: map this step's AtomRefs to current GetIdx on ``mol``."""
        return frozenset(ref.resolve(mol) for ref in self.site)

    def apply(self, mol, **kwargs) -> list:
        """Run the named Forest rule at the resolved site; update ``_forest``."""
        parent = mol
        mol = Chem.Mol(mol)
        _copy_forest(parent, mol)
        _ensure_apply_ready(mol)
        rule = get_rule(self.rule)
        site = self.resolve_site(mol)
        products = []
        for _outsite, frags in rule.metabolites_from_sites(
            mol,
            site,
            tag_atoms=True,
            only_emit_topologically_distinct_sites=False,
            **kwargs
        ):
            for frag in frags:
                if not frag:
                    continue
                index = _install_product_forest(mol, frag)
                _record_creations(self, mol, frag, index)
                products.append(frag)
        return products


@dataclass(frozen=True)
class Linearization:
    """Ordered steps; ``apply`` chains resolve → rule → ``atom_refs`` record."""

    steps: tuple

    def __post_init__(self):
        object.__setattr__(self, "steps", tuple(self.steps))

    def apply(self, mol, toward=None, **kwargs) -> list:
        """Apply steps in order.

        ``toward``: optional extra :class:`AtomRef` sequence used only for
        multi-fragment retention (e.g. final DH sites when applying a prep
        prefix). Not applied as reactions.
        """
        if not self.steps:
            return [Chem.Mol(mol)]
        root = Chem.Mol(mol)
        _copy_forest(mol, root)
        _ensure_apply_ready(root)
        currents = [root]
        toward_refs = list(toward or ())
        for i, step in enumerate(self.steps):
            later_refs = list(toward_refs)
            for later in self.steps[i + 1 :]:
                later_refs.extend(later.site)
            nxt = []
            for cur in currents:
                try:
                    products = step.apply(cur, **kwargs)
                except KeyError:
                    return []
                if not products:
                    return []
                if later_refs:
                    products = [
                        product
                        for product in products
                        if _fragment_retains_refs(product, later_refs)
                    ]
                nxt.extend(products)
            if not nxt:
                return []
            currents = nxt
        return currents


class StepPlan:
    """Partial order over :class:`Step` nodes.

    Edges in ``precedes`` are pairs of indices ``(i, j)`` meaning
    ``steps[i]`` must occur before ``steps[j]``. Within a layer from
    :meth:`layers`, order is undefined.
    """

    __slots__ = ("_steps", "_precedes")

    def __init__(
        self,
        steps: Sequence[Step],
        precedes: Sequence[tuple[int, int]] = (),
    ):
        self._steps = tuple(steps)
        self._precedes = tuple((int(a), int(b)) for a, b in precedes)
        n = len(self._steps)
        for a, b in self._precedes:
            if not (0 <= a < n and 0 <= b < n):
                raise ValueError("precedes index out of range")
            if a == b:
                raise ValueError("precedes cannot be reflexive")

    @classmethod
    def singleton(cls, rule: str, site) -> StepPlan:
        """One step, no ordering constraints."""
        return cls((Step(rule, site),), ())

    @classmethod
    def layers(cls, layers: Sequence[Sequence[Step]]) -> StepPlan:
        """Stack layers: every node in layer k precedes every node in layer k+1.

        Order within a layer is undefined (no edges).
        """
        steps: list[Step] = []
        ranges: list[tuple[int, int]] = []
        for layer in layers:
            layer = list(layer)
            if not layer:
                continue
            start = len(steps)
            steps.extend(layer)
            ranges.append((start, len(steps)))
        precedes: list[tuple[int, int]] = []
        for (a0, a1), (b0, b1) in zip(ranges, ranges[1:]):
            for i in range(a0, a1):
                for j in range(b0, b1):
                    precedes.append((i, j))
        return cls(steps, precedes)

    def __len__(self) -> int:
        return len(self._steps)

    def __eq__(self, other) -> bool:
        if not isinstance(other, StepPlan):
            return NotImplemented
        return self._steps == other._steps and self._precedes == other._precedes

    def __hash__(self) -> int:
        return hash((self._steps, self._precedes))

    def __repr__(self) -> str:
        return "StepPlan(steps=%r, precedes=%r)" % (self._steps, self._precedes)

    @property
    def steps(self) -> tuple[Step, ...]:
        return self._steps

    @property
    def precedes(self) -> tuple[tuple[int, int], ...]:
        return self._precedes

    def to_json(self) -> dict:
        return {
            "steps": [
                {"rule": s.rule, "site": _site_to_json(s.site)} for s in self._steps
            ],
            "precedes": [list(p) for p in self._precedes],
        }

    @classmethod
    def from_json(cls, data: dict) -> StepPlan:
        steps = [
            Step(item["rule"], _site_from_json(item["site"])) for item in data["steps"]
        ]
        precedes = [tuple(p) for p in data.get("precedes", ())]
        return cls(steps, precedes)

    @classmethod
    def from_mol(cls, mol) -> StepPlan:
        """Read the ``phase1_steps`` mol property."""
        if not mol.HasProp("phase1_steps"):
            raise ValueError("mol has no phase1_steps property")
        return cls.from_json(json.loads(mol.GetProp("phase1_steps")))

    def attach_to_mol(self, mol) -> None:
        """Stamp ``phase1_steps`` JSON onto ``mol``."""
        mol.SetProp("phase1_steps", json.dumps(self.to_json(), separators=(",", ":")))

    def apply_linearization(
        self, mol, order: Sequence[Step], toward=None, **kwargs
    ) -> list:
        """Apply an ordered sequence of steps (typically from ``iter_linearizations``)."""
        return Linearization(tuple(order)).apply(mol, toward=toward, **kwargs)

    def iter_linearizations(self) -> Iterator[tuple[Step, ...]]:
        """Lazily yield every total order consistent with ``precedes``."""
        n = len(self._steps)
        if n == 0:
            yield ()
            return

        successors = [[] for _ in range(n)]
        indegree = [0] * n
        for a, b in self._precedes:
            successors[a].append(b)
            indegree[b] += 1

        def search(ready: list[int], remaining_indegree: list[int], path: list[int]):
            if len(path) == n:
                yield tuple(self._steps[i] for i in path)
                return
            # Stable choice order for determinism among equal-priority nodes.
            for idx, node in enumerate(list(ready)):
                path.append(node)
                next_ready = ready[:idx] + ready[idx + 1 :]
                next_indegree = list(remaining_indegree)
                for succ in successors[node]:
                    next_indegree[succ] -= 1
                    if next_indegree[succ] == 0:
                        next_ready.append(succ)
                yield from search(sorted(next_ready), next_indegree, path)
                path.pop()

        initial = sorted(i for i, d in enumerate(indegree) if d == 0)
        if not initial and n:
            raise ValueError("StepPlan precedes graph has a cycle")
        yield from search(initial, indegree, [])
