"""Partial-order / nested pathway plans over reaction steps.

``AtomRef`` / ``Step`` / ``Linearization`` form a reusable apply+resolve layer.
``StepPlan`` is an ordered sequence of children (Seq); ``And``, ``Or``, and
``Deps`` subclass it. ``Deps`` is flat steps plus precedes edges (linearizations
are topological sorts). Creation bookkeeping lives on
``mol._forest["atom_refs"]`` (private), not in AtomTracker tag strings.
"""

from __future__ import annotations

import json
import math
from collections import deque
from dataclasses import dataclass
from functools import cache
from typing import Iterator, Sequence, Union

from rdkit import Chem

from .base import AtomTracker, _copy_forest, _forest_state
from .unstable import unstable

# Private storage key for StepPlan.attach_to_mol / from_mol / try_from_mol.
# Callers must use those APIs; do not read or write this prop directly.
_PHASE1_STEPS_PROP = "phase1_steps"

# Pathway nodes: Step (leaf) | StepPlan (Seq) | And | Or | Deps
PathwayNode = Union["Step", "StepPlan"]


# ---------------------------------------------------------------------------
# Rule registry (name -> ReactionRule instance)
# ---------------------------------------------------------------------------

_RULE_CACHE: dict = {}


def get_rule(name: str, pathways=()):
    """Return a shared ``ReactionRule`` instance for ``name``.

    ``pathways`` selects optional uncommon chemotypes (e.g. ``("methide",)`` on
    Dehydrogenation). Default empty keeps broad Phase I enumeration lean.
    """
    pathways = frozenset(pathways or ())
    key = (name, pathways) if pathways else name
    rule = _RULE_CACHE.get(key)
    if rule is None:
        from . import rules as forest_rules

        cls = getattr(forest_rules, name, None)
        if cls is None:
            raise KeyError("unknown forest rule %r" % (name,))
        if pathways and name == "Dehydrogenation":
            rule = cls(pathways=pathways)
        else:
            rule = cls()
        _RULE_CACHE[key] = rule
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

    def __str__(self) -> str:
        if self.origin is not None:
            return str(self.origin)
        at = ",".join(str(i) for i in sorted(self.at))
        return "%s[%s]" % (self.added_by, at)

    def __repr__(self) -> str:
        if self.origin is not None:
            return "AtomRef(%s)" % (self.origin,)
        return "AtomRef(%s[%s])" % (
            self.added_by,
            ",".join(str(i) for i in sorted(self.at)),
        )


def _atomref_sort_key(ref: AtomRef):
    if ref.origin is not None:
        return (0, ref.origin, "", ())
    return (1, -1, ref.added_by or "", tuple(sorted(ref.at or ())))


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
    """One named reaction at a site of :class:`AtomRef` (ints coerce to origin).

    ``pathways`` names optional uncommon chemotypes on the rule (e.g.
    ``("methide",)`` for quinone-methide Dehydrogenation). Empty by default so
    Phase I enumeration stays lean unless a plan or caller opts in.
    """

    rule: str
    site: frozenset
    pathways: frozenset = frozenset()

    def __post_init__(self):
        object.__setattr__(self, "site", _coerce_site(self.site))
        object.__setattr__(self, "pathways", frozenset(self.pathways or ()))

    def __str__(self) -> str:
        sites = ", ".join(
            str(ref) for ref in sorted(self.site, key=_atomref_sort_key)
        )
        base = "%s[%s]" % (self.rule, sites)
        if self.pathways:
            return "%s{%s}" % (base, ",".join(sorted(self.pathways)))
        return base

    def __repr__(self) -> str:
        if self.pathways:
            return "Step(%r, %s, pathways=%r)" % (
                self.rule,
                "{" + ", ".join(repr(r) for r in sorted(self.site, key=_atomref_sort_key)) + "}",
                set(self.pathways),
            )
        return "Step(%r, %s)" % (
            self.rule,
            "{" + ", ".join(repr(r) for r in sorted(self.site, key=_atomref_sort_key)) + "}",
        )

    def resolve_site(self, mol) -> frozenset:
        """Opt-in: map this step's AtomRefs to current GetIdx on ``mol``."""
        return frozenset(ref.resolve(mol) for ref in self.site)

    def apply(self, mol, **kwargs) -> list:
        """Run the named Forest rule at the resolved site; update ``_forest``."""
        parent = mol
        mol = Chem.Mol(mol)
        _copy_forest(parent, mol)
        _ensure_apply_ready(mol)
        rule = get_rule(self.rule, pathways=self.pathways)
        site = self.resolve_site(mol)
        products = []
        apply_kw = dict(kwargs)
        if self.pathways:
            apply_kw.setdefault("pathways", self.pathways)
        for _outsite, frags in rule.metabolites_from_sites(
            mol,
            site,
            tag_atoms=True,
            only_emit_topologically_distinct_sites=False,
            **apply_kw
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

    def apply(self, mol, toward=None, drop_last: int = 0, **kwargs) -> list:
        """Apply steps in order.

        ``toward``: optional extra :class:`AtomRef` sequence used only for
        multi-fragment retention (e.g. sites of steps not yet applied).
        Not applied as reactions.

        ``drop_last``: omit this many trailing steps from the run, but pass
        their sites as ``toward`` so fragments needed for them are kept.
        Equivalent to a prep-only replay when the final step is Dehydrogenation.
        """
        drop_last = int(drop_last)
        if drop_last < 0:
            raise ValueError("drop_last must be >= 0")
        if drop_last:
            if drop_last >= len(self.steps):
                return [Chem.Mol(mol)]
            toward_refs = list(toward or ())
            for step in self.steps[-drop_last:]:
                toward_refs.extend(step.site)
            return Linearization(self.steps[:-drop_last]).apply(
                mol, toward=toward_refs, drop_last=0, **kwargs
            )

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


# ---------------------------------------------------------------------------
# Nested pathway algebra: StepPlan (= Seq), And, Or
# ---------------------------------------------------------------------------


def _as_node(child):
    """Normalize a child to Step or StepPlan."""
    if isinstance(child, Step):
        return child
    if isinstance(child, StepPlan):
        return child
    raise TypeError("pathway child must be Step or StepPlan, got %r" % (type(child),))


def _node_str(node) -> str:
    return str(node)


class StepPlan:
    """Ordered sequence of pathway children (``Step`` or nested plans).

    This is the Seq form: children run left-to-right. Subclasses ``And`` and
    ``Or`` change only how :meth:`linearizations` expands children.
    """

    __slots__ = ("_children",)


    def __init__(self, children=(), precedes=()):
        if precedes:
            raise TypeError(
                "use Deps(steps, precedes) or StepPlan.from_steps_precedes(...); "
                "StepPlan(..., precedes=) no longer converts via Kahn layers"
            )
        if isinstance(children, (Step, StepPlan)) and not isinstance(
            children, (list, tuple)
        ):
            children = (children,)
        # Legacy flat list of Steps with no precedes → Seq of those steps
        self._children = tuple(_as_node(c) for c in children)

    @classmethod
    def singleton(cls, rule: str, site) -> StepPlan:
        return cls((Step(rule, site),))

    @classmethod
    def empty(cls) -> StepPlan:
        return Or(())

    @classmethod
    def layers(cls, layers) -> StepPlan:
        """Stack layers: each layer is an ``And`` (or a single ``Step``)."""
        parts = []
        for layer in layers:
            layer = tuple(layer)
            if not layer:
                continue
            if len(layer) == 1:
                parts.append(layer[0])
            else:
                parts.append(And(layer))
        if not parts:
            return cls.empty()
        if len(parts) == 1:
            return parts[0] if isinstance(parts[0], StepPlan) else cls((parts[0],))
        return cls(parts)

    @classmethod
    def from_steps_precedes(cls, steps, precedes=()) -> StepPlan:
        """Build a :class:`Deps` plan from flat steps + precedes edges."""
        return Deps(steps, precedes)

    @property
    def children(self) -> tuple:
        return self._children

    @property
    def expr(self) -> StepPlan:
        """Self (plans are the expression tree)."""
        return self

    def branches(self) -> list:
        """Top-level ``Or`` alternatives; default is ``[self]``."""
        return [self]

    def __bool__(self) -> bool:
        return bool(self._children)

    def __len__(self) -> int:
        return sum(
            1 if isinstance(c, Step) else len(c) for c in self._children
        )

    def __iter__(self):
        return iter(self._children)

    def __eq__(self, other) -> bool:
        if not isinstance(other, StepPlan):
            return NotImplemented
        return type(self) is type(other) and self._children == other._children

    def __hash__(self) -> int:
        return hash((type(self), self._children))

    def __repr__(self) -> str:
        return "%s(%r)" % (type(self).__name__, str(self))

    def __str__(self) -> str:
        if not self._children:
            return "%s()" % type(self).__name__
        parts = [_node_str(c) for c in self._children]
        if len(parts) == 1:
            return parts[0]
        return " → ".join(parts)

    @property
    def steps(self) -> tuple:
        """Leaf steps in declaration order (empty for multi-branch ``Or``)."""
        out = []
        for c in self._children:
            if isinstance(c, Step):
                out.append(c)
            else:
                out.extend(c.steps)
        return tuple(out)

    @property
    def precedes(self) -> tuple:
        """Legacy flat precedes for a pure Seq-of-And spine; else empty."""
        return _legacy_precedes(self)

    def to_json(self) -> dict:
        return pathway_to_json(self)

    @classmethod
    def from_json(cls, data: dict) -> StepPlan:
        node = pathway_from_json(data)
        if isinstance(node, Step):
            return StepPlan((node,))
        return node

    @classmethod
    def try_from_mol(cls, mol) -> StepPlan | None:
        if not mol.HasProp(_PHASE1_STEPS_PROP):
            return None
        return cls.from_json(json.loads(mol.GetProp(_PHASE1_STEPS_PROP)))

    @classmethod
    def from_mol(cls, mol) -> StepPlan:
        plan = cls.try_from_mol(mol)
        if plan is None:
            raise ValueError("mol has no attached StepPlan")
        return plan

    def attach_to_mol(self, mol) -> None:
        mol.SetProp(
            _PHASE1_STEPS_PROP, json.dumps(self.to_json(), separators=(",", ":"))
        )

    def linearizations(self) -> Iterator[Linearization]:
        """Lazily yield each total order of leaf :class:`Step`\\ s."""
        for order in self.iter_linearizations():
            yield Linearization(order)

    @unstable(name="StepPlan.n_linearizations")
    def n_linearizations(self) -> int:
        """Number of total orders — no enumeration.

        ``Seq`` / ``And`` / ``Or`` use closed recursive formulas (And needs a
        length→count map so variable-length ``Or`` children stay exact).
        ``Deps`` counts topological sorts by weakly connected components
        (bitmask DP per component + multinomial interleaving).
        """
        return _n_linearizations(self)

    def iter_linearizations(self) -> Iterator[tuple]:
        """Seq semantics: concatenate child linearizations in order."""
        yield from _seq_orders(self._children)

    def contains(self, recipe, *, by: str = "step") -> bool:
        """True iff ``recipe`` is exactly one linearization of this plan.

        Tree walk (Or / Seq splits / And interleavings) — does **not** expand
        every total order. Preferred over ``recipe in list(iter_linearizations())``.

        Parameters
        ----------
        recipe
            Sequence of :class:`Step` (``by="step"``) or rule-name strings
            (``by="rule"``). A lone :class:`Step` / :class:`Linearization` is
            accepted.
        by
            ``"step"`` — full Step equality; ``"rule"`` — compare ``Step.rule``
            only (useful for walk recipes that lack AtomRef sites).
        """
        if isinstance(recipe, Linearization):
            recipe = recipe.steps
        elif isinstance(recipe, Step):
            recipe = (recipe,)
        recipe = tuple(recipe)
        if by not in ("step", "rule"):
            raise ValueError("by must be 'step' or 'rule'")
        return _recipe_in_node(self, recipe, by=by)

    def __contains__(self, item) -> bool:
        if isinstance(item, (Linearization, Step, list, tuple)):
            return self.contains(item, by="step")
        return False


class And(StepPlan):
    """All children required; order unconstrained (all interleavings)."""

    def __str__(self) -> str:
        if not self._children:
            return "And()"
        parts = [_node_str(c) for c in self._children]
        if len(parts) == 1:
            return parts[0]
        return "(" + " & ".join(parts) + ")"

    def iter_linearizations(self) -> Iterator[tuple]:
        if not self._children:
            yield ()
            return
        from itertools import product

        child_orders = [list(_node_orders(c)) for c in self._children]
        for choice in product(*child_orders):
            yield from _interleave_orders(choice)


class Or(StepPlan):
    """Choose exactly one child pathway."""

    def __bool__(self) -> bool:
        return bool(self._children)

    def __len__(self) -> int:
        if not self._children:
            return 0
        return len(self.branches()[0])

    def __str__(self) -> str:
        if not self._children:
            return "Or()"
        parts = []
        for c in self._children:
            s = _node_str(c)
            if isinstance(c, StepPlan) and not isinstance(c, Or) and (
                "→" in s or "&" in s
            ):
                parts.append("(%s)" % s)
            else:
                parts.append(s)
        if len(parts) == 1:
            return parts[0]
        return " | ".join(parts)

    def branches(self) -> list:
        out = []
        for c in self._children:
            if isinstance(c, Step):
                out.append(StepPlan((c,)))
            elif isinstance(c, Or):
                out.extend(c.branches())
            else:
                out.append(c)
        return out

    @property
    def steps(self) -> tuple:
        if len(self._children) == 1:
            c = self._children[0]
            return (c,) if isinstance(c, Step) else c.steps
        return ()

    def iter_linearizations(self) -> Iterator[tuple]:
        if not self._children:
            return
        for child in self._children:
            yield from _node_orders(child)


@unstable
class Deps(StepPlan):
    """Flat steps + precedes edges; linearizations are topological sorts.

    Same role as ``And`` / ``Or``: a ``StepPlan`` specialization. Prefer this
    over Kahn→And/Seq when the dependency graph is the source of truth (guided
    emission, apply-replay correction). And/Or/Seq are series-parallel and
    cannot represent every precedes graph losslessly (different
    ``n_linearizations`` algorithms are a symptom of that gap).
    """

    __slots__ = ("_precedes",)

    def __init__(self, steps=(), precedes=()):
        steps = tuple(steps or ())
        nodes = []
        for s in steps:
            node = _as_node(s)
            if not isinstance(node, Step):
                raise TypeError("Deps steps must be Step instances, got %r" % (type(node),))
            nodes.append(node)
        precedes = tuple((int(a), int(b)) for a, b in (precedes or ()))
        n = len(nodes)
        for a, b in precedes:
            if not (0 <= a < n and 0 <= b < n) or a == b:
                raise ValueError("invalid precedes")
        if n > 0 and next(_topo_orders(nodes, precedes), None) is None:
            raise ValueError("cycle in precedes")
        self._children = tuple(nodes)
        # Always store the transitive reduction for stable output / JSON / keys.
        # Graphs are tiny; reduction cost is fine vs O(M+N) closure equality.
        self._precedes = (
            canonical_dependency_edges(n, precedes) if n else ()
        )

    @property
    def precedes(self) -> tuple:
        return self._precedes

    def same_linearizations(self, other: "Deps") -> bool:
        """True iff ``other`` admits exactly the same total orders.

        Requires **node identity** (same multiset of leaf :class:`Step`\\ s),
        then compares canonical precedes (stored transitive reductions).
        Prefer this over ``==`` only when declaration order of steps may differ;
        after construction, ``==`` also sees reduced edges.
        """
        if not isinstance(other, Deps):
            return False
        aligned = _align_deps_indices(self._children, other._children)
        if aligned is None:
            return False
        # Both sides already reduced; remap other's edges into self's index.
        edges_b = tuple(
            sorted((aligned[a], aligned[b]) for a, b in other._precedes)
        )
        return self._precedes == edges_b

    def __eq__(self, other) -> bool:
        """Same steps (order-sensitive) and same canonical precedes.

        Construction always reduces edges, so redundant input graphs compare
        equal once built. For order-insensitive step alignment use
        :meth:`same_linearizations`.
        """
        if not isinstance(other, Deps):
            return NotImplemented
        return self._children == other._children and self._precedes == other._precedes

    def __hash__(self) -> int:
        return hash((Deps, self._children, self._precedes))

    def __str__(self) -> str:
        if not self._children:
            return "Deps()"
        parts = [_node_str(c) for c in self._children]
        if not self._precedes:
            if len(parts) == 1:
                return parts[0]
            return "(" + " & ".join(parts) + ")"
        edges = ", ".join("%s≺%s" % (parts[a], parts[b]) for a, b in self._precedes)
        return "[%s | %s]" % ("; ".join(parts), edges)

    def iter_linearizations(self) -> Iterator[tuple]:
        yield from _topo_orders(self._children, self._precedes)


def or_of_plans(plans: Sequence) -> StepPlan:
    """One plan: empty ``Or``, a single plan, or ``Or`` of alternatives."""
    plans = tuple(
        p if isinstance(p, StepPlan) else StepPlan((p,)) for p in plans
    )
    if not plans:
        return Or(())
    if len(plans) == 1:
        return plans[0]
    return Or(plans)


def pathway_alternatives(plan) -> list:
    """Top-level ``Or`` branches as plans (compat helper)."""
    if plan is None:
        return []
    if isinstance(plan, StepPlan):
        return plan.branches()
    if isinstance(plan, Step):
        return [StepPlan((plan,))]
    return []


def pathway_linearizations(node) -> Iterator[Linearization]:
    if isinstance(node, StepPlan):
        yield from node.linearizations()
        return
    for order in _node_orders(node):
        yield Linearization(order)


def _node_orders(node) -> Iterator[tuple]:
    if isinstance(node, Step):
        yield (node,)
        return
    if isinstance(node, StepPlan):
        yield from node.iter_linearizations()
        return
    raise TypeError(type(node))


def _step_key(step: Step, by: str):
    if by == "rule":
        return step.rule
    return step


def _recipe_item_key(item, by: str):
    if by == "rule":
        if isinstance(item, Step):
            return item.rule
        return item
    if isinstance(item, Step):
        return item
    raise TypeError("recipe items must be Step when by='step'")


def _n_linearizations(node) -> int:
    """Total linearization count without materializing orders."""
    return sum(_lin_counts_by_length(node).values())


def _lin_counts_by_length(node) -> dict:
    """Map leaf-length → number of linearizations of that length."""
    if isinstance(node, Step):
        return {1: 1}
    if isinstance(node, Or):
        if not node._children:
            return {}
        out = {}
        for c in node._children:
            for length, count in _lin_counts_by_length(c).items():
                out[length] = out.get(length, 0) + count
        return out
    if isinstance(node, Deps):
        n = len(node._children)
        if n == 0:
            return {0: 1}
        return {n: count_transformation_orders(n, node._precedes)}
    if isinstance(node, And):
        return _and_lin_counts_by_length(node._children)
    if isinstance(node, StepPlan):
        return _seq_lin_counts_by_length(node._children)
    raise TypeError(type(node))


def _seq_lin_counts_by_length(children) -> dict:
    if not children:
        return {0: 1}
    acc = {0: 1}
    for child in children:
        nxt = {}
        child_map = _lin_counts_by_length(child)
        for la, ca in acc.items():
            for lb, cb in child_map.items():
                length = la + lb
                nxt[length] = nxt.get(length, 0) + ca * cb
        acc = nxt
    return acc


def _and_lin_counts_by_length(children) -> dict:
    if not children:
        return {0: 1}
    # Fold children left-to-right: interleave current bag with next child.
    acc = _lin_counts_by_length(children[0])
    for child in children[1:]:
        nxt = {}
        child_map = _lin_counts_by_length(child)
        for la, ca in acc.items():
            for lb, cb in child_map.items():
                ways = ca * cb * (
                    math.factorial(la + lb)
                    // (math.factorial(la) * math.factorial(lb))
                )
                length = la + lb
                nxt[length] = nxt.get(length, 0) + ways
        acc = nxt
    return acc


class ComponentTooLarge(RuntimeError):
    """Raised before DP if a weakly connected dependency component is too big."""

    def __init__(self, nodes, max_size):
        self.nodes = tuple(nodes)
        self.max_size = max_size
        super().__init__(
            "Dependency component has %d nodes (limit is %d): %s"
            % (len(nodes), max_size, list(nodes))
        )


def _align_deps_indices(steps_a, steps_b) -> dict | None:
    """Map indices in ``steps_b`` → indices in ``steps_a`` by Step equality.

    Returns ``None`` if the leaf multisets differ (node identity failed).
    Duplicate equal Steps (same rule+SOM) are allowed — expected mainly from
    pathological cycles, not normal emission — and matched greedily.
    """
    from collections import Counter

    if len(steps_a) != len(steps_b):
        return None
    if Counter(steps_a) != Counter(steps_b):
        return None
    used = [False] * len(steps_b)
    remap = {}
    for i, step in enumerate(steps_a):
        found = None
        for j, other in enumerate(steps_b):
            if used[j]:
                continue
            if other == step:
                found = j
                break
        if found is None:
            return None
        used[found] = True
        remap[found] = i
    return remap


def transitive_closure_masks(n, edges) -> tuple:
    """One bitmask per node: bit ``b`` set iff ``a`` must precede ``b``.

    Nodes are opaque index labels ``0..n-1`` — this does **not** check that
    two plans' steps are the same reactions/sites. Align node identity first
    (see :meth:`Deps.same_linearizations`).

    Raises ``ValueError`` if the graph contains a cycle.
    """
    outgoing = [[] for _ in range(n)]
    indegree = [0] * n
    for a, b in set(edges):
        if not (0 <= a < n and 0 <= b < n):
            raise ValueError("invalid edge: %r" % ((a, b),))
        if a == b:
            raise ValueError("cycle in precedes")
        outgoing[a].append(b)
        indegree[b] += 1

    queue = deque(v for v in range(n) if indegree[v] == 0)
    topo = []
    while queue:
        a = queue.popleft()
        topo.append(a)
        for b in outgoing[a]:
            indegree[b] -= 1
            if indegree[b] == 0:
                queue.append(b)
    if len(topo) != n:
        raise ValueError("cycle in precedes")

    closure = [0] * n
    for a in reversed(topo):
        for b in outgoing[a]:
            closure[a] |= (1 << b) | closure[b]
    return tuple(closure)


def canonical_dependency_edges(n, edges) -> tuple:
    """Unique minimal dependency edge set (transitive reduction) for this DAG.

    Two DAGs over the **same labeled nodes** ``0..n-1`` have exactly the same
    valid orderings iff this returns the same edge list. Node identity is the
    caller's responsibility. :class:`Deps` stores this form on construction
    for stable output.
    """
    closure = transitive_closure_masks(n, edges)
    canonical_edges = []
    for a in range(n):
        descendants = closure[a]
        reachable_through_another_node = 0
        remaining = descendants
        while remaining:
            bit = remaining & -remaining
            b = bit.bit_length() - 1
            remaining -= bit
            reachable_through_another_node |= closure[b]
        direct_children = descendants & ~reachable_through_another_node
        while direct_children:
            bit = direct_children & -direct_children
            b = bit.bit_length() - 1
            direct_children -= bit
            canonical_edges.append((a, b))
    return tuple(sorted(canonical_edges))


def count_transformation_orders(n, edges, *, max_component_size=20) -> int:
    """Count valid orders of ``n`` labeled nodes under precedes ``edges``.

    Nodes are index labels only. Splits into weakly connected components,
    counts each with memoized bitmask DP, then multiplies by multinomial
    interleavings across components. Returns 0 on a cycle. Raises
    :class:`ComponentTooLarge` if any component exceeds ``max_component_size``
    (``None`` disables the limit).
    """
    edge_set = set()
    for a, b in edges:
        if not (0 <= a < n and 0 <= b < n):
            raise ValueError("invalid dependency edge %r for n=%d" % ((a, b), n))
        if a != b:
            edge_set.add((a, b))

    undirected = [[] for _ in range(n)]
    for a, b in edge_set:
        undirected[a].append(b)
        undirected[b].append(a)

    component_id = [-1] * n
    components = []
    for start in range(n):
        if component_id[start] != -1:
            continue
        cid = len(components)
        component = []
        queue = deque([start])
        component_id[start] = cid
        while queue:
            u = queue.popleft()
            component.append(u)
            for v in undirected[u]:
                if component_id[v] == -1:
                    component_id[v] = cid
                    queue.append(v)
        if (
            max_component_size is not None
            and len(component) > max_component_size
        ):
            raise ComponentTooLarge(component, max_component_size)
        components.append(component)

    component_edges = [[] for _ in components]
    for a, b in edge_set:
        component_edges[component_id[a]].append((a, b))

    def count_component(nodes, local_edges):
        size = len(nodes)
        local_index = {node: i for i, node in enumerate(nodes)}
        prerequisites = [0] * size
        for a, b in local_edges:
            prerequisites[local_index[b]] |= 1 << local_index[a]
        all_done = (1 << size) - 1

        @cache
        def dp(done):
            if done == all_done:
                return 1
            total = 0
            for step in range(size):
                bit = 1 << step
                if not (done & bit) and (prerequisites[step] & ~done) == 0:
                    total += dp(done | bit)
            return total

        return dp(0)

    total_count = 1
    total_steps_so_far = 0
    for nodes, local_edges in zip(components, component_edges):
        component_count = count_component(nodes, local_edges)
        if component_count == 0:
            return 0
        component_size = len(nodes)
        total_count *= component_count
        total_count *= math.comb(
            total_steps_so_far + component_size, component_size
        )
        total_steps_so_far += component_size
    return total_count


def _possible_lengths(node) -> frozenset:
    """Possible leaf counts for ``node`` (Or branches may differ)."""
    if isinstance(node, Step):
        return frozenset({1})
    if isinstance(node, Or):
        if not node._children:
            return frozenset({0})
        out = set()
        for c in node._children:
            out |= _possible_lengths(c)
        return frozenset(out)
    if isinstance(node, Deps):
        return frozenset({len(node._children)})
    # Seq and And: sums of child lengths
    sets = [_possible_lengths(c) for c in node._children]
    if not sets:
        return frozenset({0})
    acc = {0}
    for s in sets:
        acc = {a + b for a in acc for b in s}
    return frozenset(acc)


def _recipe_in_node(node, recipe: tuple, *, by: str) -> bool:
    """Exact linearization membership without full expansion."""
    n = len(recipe)
    if n not in _possible_lengths(node):
        return False
    if isinstance(node, Step):
        return n == 1 and _recipe_item_key(recipe[0], by) == _step_key(node, by)
    if type(node) is Or:
        return any(_recipe_in_node(c, recipe, by=by) for c in node._children)
    if type(node) is And:
        return _recipe_in_and(node._children, recipe, by=by)
    if type(node) is Deps:
        return _recipe_in_deps(node, recipe, by=by)
    # StepPlan Seq
    return _recipe_in_seq(node._children, recipe, by=by)


def _recipe_in_deps(node: "Deps", recipe: tuple, *, by: str) -> bool:
    """True if ``recipe`` is a topo sort of ``node`` under ``by`` equality."""
    steps = node._children
    if len(recipe) != len(steps):
        return False
    # Map each recipe item to a unique step index.
    used = [False] * len(steps)
    index_of = []
    for item in recipe:
        key = _recipe_item_key(item, by)
        found = None
        for i, step in enumerate(steps):
            if used[i]:
                continue
            if _step_key(step, by) == key:
                found = i
                break
        if found is None:
            return False
        used[found] = True
        index_of.append(found)
    pos = {idx: p for p, idx in enumerate(index_of)}
    return all(pos[a] < pos[b] for a, b in node._precedes)


def _recipe_in_seq(children, recipe: tuple, *, by: str) -> bool:
    """Partition ``recipe`` into contiguous segments matching each child."""
    m = len(children)
    if m == 0:
        return len(recipe) == 0
    # dp[i][j]: children[:i] match recipe[:j]
    n = len(recipe)
    dp = [False] * (n + 1)
    dp[0] = True
    for child in children:
        lens = _possible_lengths(child)
        nxt = [False] * (n + 1)
        for j in range(n + 1):
            if not dp[j]:
                continue
            for L in lens:
                end = j + L
                if end > n:
                    continue
                if _recipe_in_node(child, recipe[j:end], by=by):
                    nxt[end] = True
        dp = nxt
    return dp[n]


def _recipe_in_and(children, recipe: tuple, *, by: str) -> bool:
    """``recipe`` is an interleaving of one linearization per child."""
    if not children:
        return len(recipe) == 0
    child0, *rest = children
    for L in _possible_lengths(child0):
        if L > len(recipe):
            continue
        from itertools import combinations

        for idxs in combinations(range(len(recipe)), L):
            sub = tuple(recipe[i] for i in idxs)
            if not _recipe_in_node(child0, sub, by=by):
                continue
            idx_set = set(idxs)
            complement = tuple(
                recipe[i] for i in range(len(recipe)) if i not in idx_set
            )
            if rest:
                if _recipe_in_and(rest, complement, by=by):
                    return True
            elif not complement:
                return True
    return False


def _seq_orders(children) -> Iterator[tuple]:
    if not children:
        yield ()
        return

    def concat(i: int, prefix: tuple):
        if i == len(children):
            yield prefix
            return
        for mid in _node_orders(children[i]):
            yield from concat(i + 1, prefix + mid)

    yield from concat(0, ())


def _interleave_orders(orders) -> Iterator[tuple]:
    if not orders:
        yield ()
        return
    total = sum(len(o) for o in orders)
    if total == 0:
        yield ()
        return

    def search(pos, path):
        if len(path) == total:
            yield tuple(path)
            return
        for i, o in enumerate(orders):
            if pos[i] < len(o):
                path.append(o[pos[i]])
                pos[i] += 1
                yield from search(pos, path)
                pos[i] -= 1
                path.pop()

    yield from search([0] * len(orders), [])


def _topo_orders(steps, precedes) -> Iterator[tuple]:
    """Yield every topological sort of ``steps`` under ``precedes`` edges."""
    n = len(steps)
    if n == 0:
        yield ()
        return
    successors = [[] for _ in range(n)]
    indegree = [0] * n
    for a, b in precedes:
        successors[a].append(b)
        indegree[b] += 1

    def search(remaining, path):
        if len(path) == n:
            yield tuple(steps[i] for i in path)
            return
        ready = [i for i in range(n) if remaining[i] == 0 and i not in path]
        for i in ready:
            remaining[i] = -1
            path.append(i)
            for succ in successors[i]:
                remaining[succ] -= 1
            yield from search(remaining, path)
            for succ in successors[i]:
                remaining[succ] += 1
            path.pop()
            remaining[i] = 0

    yield from search(list(indegree), [])


def _legacy_precedes(plan: StepPlan) -> tuple:
    """Recover precedes for Deps or a pure Seq of And/Step layers."""
    if type(plan) is Deps:
        return plan._precedes
    if type(plan) is not StepPlan:
        return ()
    layers = []
    for c in plan.children:
        if isinstance(c, Step):
            layers.append([c])
        elif type(c) is And:
            layers.append(list(c.steps))
        else:
            return ()
    steps = []
    ranges = []
    for layer in layers:
        start = len(steps)
        steps.extend(layer)
        ranges.append((start, len(steps)))
    precedes = []
    for (a0, a1), (b0, b1) in zip(ranges, ranges[1:]):
        for i in range(a0, a1):
            for j in range(b0, b1):
                precedes.append((i, j))
    return tuple(precedes)


def pathway_to_json(node) -> dict:
    if isinstance(node, Step):
        data = {"op": "step", "rule": node.rule, "site": _site_to_json(node.site)}
        if node.pathways:
            data["pathways"] = sorted(node.pathways)
        return data
    if isinstance(node, Deps):
        return {
            "op": "deps",
            "steps": [
                {
                    "rule": s.rule,
                    "site": _site_to_json(s.site),
                    **({"pathways": sorted(s.pathways)} if s.pathways else {}),
                }
                for s in node._children
            ],
            "precedes": [list(p) for p in node._precedes],
        }
    if isinstance(node, And):
        return {"op": "and", "children": [pathway_to_json(c) for c in node.children]}
    if isinstance(node, Or):
        return {"op": "or", "children": [pathway_to_json(c) for c in node.children]}
    if isinstance(node, StepPlan):
        # Seq
        return {"op": "seq", "children": [pathway_to_json(c) for c in node.children]}
    raise TypeError(type(node))


def pathway_from_json(data: dict):
    if not isinstance(data, dict):
        raise TypeError("pathway JSON must be a dict")
    op = data.get("op")
    if op in (None, "plan", "deps") and "steps" in data and "children" not in data:
        steps = []
        for item in data["steps"]:
            pathways = frozenset(item.get("pathways") or ())
            steps.append(
                Step(item["rule"], _site_from_json(item["site"]), pathways=pathways)
            )
        precedes = [tuple(p) for p in data.get("precedes", ())]
        return Deps(steps, precedes)
    if op == "step":
        pathways = frozenset(data.get("pathways") or ())
        return Step(
            data["rule"], _site_from_json(data["site"]), pathways=pathways
        )
    if op == "and":
        return And(tuple(pathway_from_json(c) for c in data["children"]))
    if op == "or":
        return Or(tuple(pathway_from_json(c) for c in data["children"]))
    if op == "seq":
        return StepPlan(tuple(pathway_from_json(c) for c in data["children"]))
    raise ValueError("unknown pathway op %r" % (op,))


# Back-compat alias: Seq was renamed — StepPlan *is* Seq.
Seq = StepPlan
