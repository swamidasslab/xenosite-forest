"""Guided reactant→product path finder (beside classic ``RuleSet.find_path``).

Uses :class:`~xenosite.forest.path_context.PathContext`, rule-owned hooks
(``is_redundant``, ``enumerate_for_path``, ``is_terminal_product``, …),
and instrumented counters. Hits are :class:`PathOutcome` (Required plan +
declarative cleavage-side Maybe) — Maybe is not searched or enumerated.
"""

from __future__ import annotations

import json
from collections import deque
from dataclasses import dataclass
from typing import Iterator, Optional

from rdkit import Chem

from .path_context import PathContext
from .phaseone import PhaseOneQF, PhaseOneRS
from .rulesets import RULESETS, RuleSet, load_ruleset
from .base import copy_mol
from .step_plan import (
    And,
    Deps,
    Linearization,
    Step,
    StepPlan,
    _origin_refs,
    canonical_dependency_edges,
    get_rule,
    pathway_linearizations,
)
from .utils import canon_smi, has_star_conjugate, unmapped_smiles
from .unstable import unstable

RULESETS["PhaseOneRS"] = PhaseOneRS
RULESETS["PhaseOneQF"] = PhaseOneQF

# Primary search stop — depth defaults to unbounded.
DEFAULT_MAX_EXPANSIONS = 200


@dataclass(frozen=True)
class CleavageSide:
    """Open side-metabolism bag for a fragment discarded at a bifurcation.

    ``site`` is the bifurcating cleavage site; ``opens`` are prior ring-open
    sites on the walk (oldest first). Together they define a segment of the
    original molecule. ``side`` is the canonical SMILES of the discarded
    fragment. Meaning: any metabolism on that segment may have occurred —
    we do not enumerate it. Preceding N-dealkylation at the same heteroatom
    (different site) is covered by this bag.
    """

    site: frozenset
    side: str
    opens: tuple = ()

    def span_sites(self) -> tuple:
        """Ring-open sites + bifurcating site (segment endpoints)."""
        return tuple(self.opens) + (self.site,)

    def __str__(self) -> str:
        def _fmt(s):
            return ", ".join(str(x) for x in sorted(s, key=str))

        if self.opens:
            parts = [_fmt(o) for o in self.opens] + [_fmt(self.site)]
            return "CleavageSide[%s]|%s" % ("; ".join(parts), self.side)
        return "CleavageSide[%s]|%s" % (_fmt(self.site), self.side)


def _site_key_set(site) -> frozenset:
    """Normalize a formation site to comparable atom keys (ints / AtomRef)."""
    if site is None:
        return frozenset()
    if isinstance(site, tuple) and len(site) >= 2 and isinstance(site[0], str):
        site = site[1]
    keys = set()
    for x in site if not isinstance(site, frozenset) else site:
        origin = getattr(x, "origin", None)
        if origin is not None:
            keys.add(int(origin))
            continue
        try:
            keys.add(int(x))
        except Exception:
            keys.add(str(x))
    return frozenset(keys)


@dataclass(frozen=True)
class MaybeFilter:
    """Conceptual bags of side metabolism that could accompany formation.

    Does **not** participate in :class:`~xenosite.forest.step_plan.StepPlan`
    algebra and is **not** searched. Bifurcating Required cleaves record
    :class:`CleavageSide` entries (optionally spanning prior ring-opens).
    """

    entries: tuple = ()

    def __bool__(self) -> bool:
        return bool(self.entries)

    def __str__(self) -> str:
        if not self.entries:
            return "MaybeFilter()"
        # Order-insensitive: And-equivalent walks accumulate the same bags.
        parts = sorted(str(e) for e in self.entries)
        return "MaybeFilter(" + "; ".join(parts) + ")"

    @classmethod
    def from_sides(cls, sides) -> MaybeFilter:
        return cls(entries=tuple(sides or ()))

    def sides(self) -> tuple:
        """Canonical SMILES of every discarded CleavageSide fragment."""
        return tuple(e.side for e in self.entries if isinstance(e, CleavageSide))

    def allows(self, rule_name: str = None, site=None, *, side: str = None) -> bool:
        """True if optional side metabolism is covered by a cleavage-side bag.

        - ``side`` SMILES: True if that discarded fragment is one of the bags.
        - ``site``: True if the site **shares atoms** with any site in a bag's
          span (``opens`` + bifurcate ``site``) but is **not identical** to the
          bifurcating formation site alone (e.g. demethyl vs TBA-forming dealk
          on the same nitrogen).
        """
        if not self.entries:
            return False
        if side is not None:
            side_c = _canon(side)
            return any(
                isinstance(e, CleavageSide) and e.side == side_c for e in self.entries
            )
        if site is None:
            return True
        site_keys = _site_key_set(site)
        if not site_keys:
            return False
        for e in self.entries:
            if not isinstance(e, CleavageSide):
                continue
            bif_keys = _site_key_set(e.site)
            if site_keys == bif_keys:
                continue  # formation cleavage itself — not Maybe
            for span in e.span_sites():
                if site_keys & _site_key_set(span):
                    return True
        return False


@unstable
@dataclass(frozen=True)
class PathOutcome:
    """One Required formation route plus its cleavage-side Maybe bag.

    Composition (not a StepPlan subclass): ``plan`` is Required-only;
    ``maybe`` never enters pathway algebra. Unpack as
    ``(smiles, steps, mols, plan, maybe)`` for back-compat.
    """

    plan: StepPlan
    maybe: MaybeFilter
    steps: tuple
    smiles: tuple
    mols: tuple = ()

    def allows(self, rule_name: str = None, site=None, *, side: str = None) -> bool:
        return self.maybe.allows(rule_name, site, side=side)

    def linearizations(self):
        return self.plan.linearizations()

    def contains(self, *args, **kwargs):
        return self.plan.contains(*args, **kwargs)

    def __iter__(self):
        yield self.smiles
        yield self.steps
        yield self.mols
        yield self.plan
        yield self.maybe

    def __len__(self) -> int:
        return 5

    def __getitem__(self, index):
        return (
            self.smiles,
            self.steps,
            self.mols,
            self.plan,
            self.maybe,
        )[index]

    def __str__(self) -> str:
        return "PathOutcome(plan=%s, maybe=%s)" % (self.plan, self.maybe)


@unstable
@dataclass
class PathSearchCounters:
    """Shared instrumentation for guided and classic ``find_path``.

    ``max_expansions`` caps :meth:`billed` work (linearizations + site applies),
    not :attr:`rule_expansions`. One Dealkylation ``enumerate_for_path`` can hide
    dozens of plan/linearization applies behind a single rule expansion.
    """

    rule_expansions: int = 0
    mol_edits: int = 0
    sites_considered: int = 0
    sites_skipped: int = 0
    nodes_enqueued: int = 0
    nodes_pruned: int = 0
    plans_expanded: int = 0
    linearizations_applied: int = 0
    site_applies: int = 0
    sanitize_dropped: int = 0
    budget_exhausted: bool = False
    _mirror: Optional[dict] = None

    def as_dict(self) -> dict:
        return {
            "rule_expansions": self.rule_expansions,
            "mol_edits": self.mol_edits,
            "sites_considered": self.sites_considered,
            "sites_skipped": self.sites_skipped,
            "nodes_enqueued": self.nodes_enqueued,
            "nodes_pruned": self.nodes_pruned,
            "plans_expanded": self.plans_expanded,
            "linearizations_applied": self.linearizations_applied,
            "site_applies": self.site_applies,
            "sanitize_dropped": self.sanitize_dropped,
            "billed": self.billed(),
            "budget_exhausted": self.budget_exhausted,
        }

    def billed(self) -> int:
        """Work charged against ``max_expansions``."""
        return int(self.linearizations_applied) + int(self.site_applies)

    def under_budget(self, max_expansions: Optional[int]) -> bool:
        """False once :meth:`billed` has reached ``max_expansions``."""
        if max_expansions is None:
            return True
        if self.billed() >= max_expansions:
            self.budget_exhausted = True
            self._sync()
            return False
        return True

    def _sync(self) -> None:
        if self._mirror is not None:
            self._mirror.update(self.as_dict())


def coerce_path_counters(counters, max_expansions=None) -> Optional[PathSearchCounters]:
    """Normalize ``None`` / dict / :class:`PathSearchCounters` for path search."""
    if counters is None:
        if max_expansions is None:
            return None
        return PathSearchCounters()
    if isinstance(counters, PathSearchCounters):
        return counters
    if isinstance(counters, dict):
        c = PathSearchCounters(
            rule_expansions=int(counters.get("rule_expansions", 0) or 0),
            mol_edits=int(counters.get("mol_edits", 0) or 0),
            sites_considered=int(counters.get("sites_considered", 0) or 0),
            sites_skipped=int(counters.get("sites_skipped", 0) or 0),
            nodes_enqueued=int(counters.get("nodes_enqueued", 0) or 0),
            nodes_pruned=int(counters.get("nodes_pruned", 0) or 0),
            plans_expanded=int(counters.get("plans_expanded", 0) or 0),
            linearizations_applied=int(
                counters.get("linearizations_applied", 0) or 0
            ),
            site_applies=int(counters.get("site_applies", 0) or 0),
            sanitize_dropped=int(counters.get("sanitize_dropped", 0) or 0),
            budget_exhausted=bool(counters.get("budget_exhausted", False)),
            _mirror=counters,
        )
        c._sync()
        return c
    raise TypeError(
        "counters must be PathSearchCounters, dict, or None; got %r"
        % (type(counters).__name__,)
    )


def _canon(mol) -> str:
    """Structure identity for path search — stereo is not relied on."""
    smi = canon_smi(unmapped_smiles(mol) if not isinstance(mol, str) else mol)
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return smi
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def _flatten_rule_instances(ruleset) -> list:
    out = []
    for item in ruleset.rules:
        if isinstance(item, RuleSet):
            out.extend(_flatten_rule_instances(item))
        else:
            out.append(item)
    return out


def _active_rules(ruleset) -> list:
    peers = _flatten_rule_instances(ruleset)
    return [r for r in peers if not r.is_redundant(peers)]


def _flatten_rules(ruleset) -> RuleSet:
    if isinstance(ruleset, RuleSet):
        return ruleset
    return load_ruleset(ruleset)


def _step_from_path_item(item) -> Step:
    if isinstance(item, Step):
        return item
    name, site = item[0], item[1]
    pathways = item[2] if len(item) > 2 else ()
    if (
        isinstance(site, tuple)
        and len(site) == 2
        and isinstance(site[0], str)
    ):
        site = site[1]
    return Step(name, site, pathways=pathways)


def _rule_name_is_cleavage(name: str) -> bool:
    try:
        return bool(get_rule(name).is_cleavage())
    except Exception:
        return False


def _step_origins(step: Step) -> set:
    return {
        int(ref.origin)
        for ref in step.site
        if getattr(ref, "origin", None) is not None
    }


def _step_depends_on(later: Step, earlier: Step) -> bool:
    """True if ``later``'s site references an atom created by ``earlier``."""
    earlier_origins = _step_origins(earlier)
    for ref in later.site:
        added_by = getattr(ref, "added_by", None)
        if not added_by or added_by[0] != earlier.rule:
            continue
        if set(added_by[1]) & earlier_origins:
            return True
    return False


def _rule_smarts_covers_site(rule, mol, site) -> bool:
    """True if rule query / reactant SMARTS covers ``site`` on ``mol``."""
    site = frozenset(int(x) for x in site)
    if not site:
        return False

    # ResonancePairRule / EditMol: query_smarts via match_queries.
    if hasattr(rule, "match_queries"):
        try:
            mq = rule.match_queries(mol)
        except Exception:
            mq = None
        if mq:
            covered = set()
            for hits in mq.values():
                for mapids, *_rest in hits:
                    if isinstance(mapids, dict):
                        covered.update(int(v) for v in mapids.values())
            # Pair sites: each endpoint appears in some query hit.
            if site <= covered:
                return True

    # Cleavage / SMARTS reactions: reactant-side site enumeration.
    if hasattr(rule, "iter_reactant_site_matches"):
        try:
            for matched in rule.iter_reactant_site_matches(mol):
                atoms = frozenset(int(x) for x in matched)
                if site == atoms or site <= atoms:
                    return True
        except Exception:
            pass

    # Raw reaction SMARTS LHS (mapid_site when present).
    mapids = list(getattr(rule, "mapid_site", None) or [])
    for smarts in getattr(rule, "smarts", None) or ():
        lhs = str(smarts).split(">>", 1)[0].strip()
        if lhs.startswith("(") and lhs.endswith(")"):
            lhs = lhs[1:-1]
        query = Chem.MolFromSmarts(lhs)
        if query is None:
            continue
        try:
            matches = mol.GetSubstructMatches(query)
        except Exception:
            continue
        for match in matches:
            if mapids:
                by_map = {}
                for qi, atom in enumerate(query.GetAtoms()):
                    mid = int(atom.GetAtomMapNum() or 0)
                    if mid and qi < len(match):
                        by_map[mid] = int(match[qi])
                try:
                    covered = frozenset(by_map[m] for m in mapids)
                except KeyError:
                    covered = frozenset(int(x) for x in match)
            else:
                covered = frozenset(int(x) for x in match)
            if site == covered or site <= covered:
                return True
    return False


def _step_smarts_matches(step: Step, mol) -> bool:
    """True if ``step``'s rule SMARTS matches at its resolved site on ``mol``."""
    try:
        site = step.resolve_site(mol)
    except Exception:
        return False
    try:
        rule = get_rule(step.rule, pathways=step.pathways)
    except Exception:
        return False
    return _rule_smarts_covers_site(rule, mol, site)


def _walk_mol_chain(steps, mols) -> list:
    """Normalize a walk to ``len(steps)+1`` mols (state before each step + final)."""
    steps = list(steps or ())
    chain = [m for m in (mols or ()) if m is not None]
    if not steps:
        return chain[:1] if chain else []
    if len(chain) >= len(steps) + 1:
        return list(chain[: len(steps) + 1])
    # Incomplete chain: keep what we have (caller may only pass reactant).
    return chain


def _plan_from_chain_smarts(nodes: list, mols: list) -> StepPlan:
    """Deps from SMARTS match-before/after on the walk metabolite chain.

    ``mols[k]`` is the structure *before* ``nodes[k]`` (and ``mols[-1]`` is the
    final product). Step ``j`` requires step ``i`` when ``j``'s SMARTS does not
    match on ``mols[i]`` but does match on ``mols[i+1]``.
    """
    n = len(nodes)
    if n <= 1:
        return StepPlan(tuple(nodes))
    if len(mols) < n + 1:
        return _plan_from_dep_layers(nodes)

    # match[j][k] = step j matches on mols[k]
    match = [
        [_step_smarts_matches(nodes[j], mols[k]) for k in range(n + 1)]
        for j in range(n)
    ]
    precedes = []
    for j in range(n):
        for i in range(j):
            if _step_depends_on(nodes[j], nodes[i]):
                precedes.append((i, j))
                continue
            # SMARTS: became matchable exactly across step i.
            if (not match[j][i]) and match[j][i + 1]:
                # Two cleavages may lose/regain sites across a split without a
                # real dependency (e.g. acetate hydrolysis then O-dealkylation).
                if _rule_name_is_cleavage(nodes[i].rule) and _rule_name_is_cleavage(
                    nodes[j].rule
                ):
                    continue
                precedes.append((i, j))
                continue
    if not precedes:
        # No enabling edges: unconstrained peers → Deps with empty precedes.
        return Deps(tuple(nodes))
    return StepPlan.from_steps_precedes(tuple(nodes), tuple(precedes))


def _plan_from_dep_layers(nodes: list) -> StepPlan:
    """Heuristic Deps when no walk mol chain is available.

    Order edges: AtomRef ``added_by`` deps; different non-cleavage rules keep
    walk Seq. Cleavage peers / cleave+prep without AtomRef dep stay unordered.
    """
    if not nodes:
        return StepPlan(())
    if len(nodes) == 1:
        return StepPlan(tuple(nodes))
    precedes = []
    for j in range(len(nodes)):
        for i in range(j):
            if _step_depends_on(nodes[j], nodes[i]):
                precedes.append((i, j))
                continue
            if nodes[i].rule == nodes[j].rule:
                continue
            if _rule_name_is_cleavage(nodes[i].rule) or _rule_name_is_cleavage(
                nodes[j].rule
            ):
                continue
            precedes.append((i, j))
    if not precedes:
        return Deps(tuple(nodes))
    return StepPlan.from_steps_precedes(tuple(nodes), tuple(precedes))


def _lin_reaches(reactant, steps, target_smi: str) -> bool:
    """True if applying ``steps`` in order yields ``target_smi``."""
    try:
        products = Linearization(tuple(steps)).apply(reactant)
    except Exception:
        return False
    if not products:
        return False
    want = _canon(target_smi)
    for product in products:
        try:
            if _canon(product) == want:
                return True
        except Exception:
            continue
    return False


def _plan_from_good_index_orders(nodes: list, good_perms: list) -> StepPlan:
    """Build Deps from index permutations that reach the target."""
    if not nodes:
        return StepPlan(())
    if len(nodes) == 1:
        return StepPlan(tuple(nodes))
    if not good_perms:
        return Deps(tuple(nodes))
    n = len(nodes)
    precedes = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if all(perm.index(i) < perm.index(j) for perm in good_perms):
                precedes.append((i, j))
    return Deps(tuple(nodes), tuple(precedes))


def _plan_lin_set_exact(plan: StepPlan, nodes: list, reactant, target_smi: str) -> bool:
    """True iff plan linearizations are exactly the leaf orders that reach T."""
    from itertools import permutations

    lin_set = set(plan.iter_linearizations())
    if not lin_set:
        return False
    if any(not _lin_reaches(reactant, order, target_smi) for order in lin_set):
        return False
    n = len(nodes)
    if n > 4:
        return True
    for perm in permutations(range(n)):
        order = tuple(nodes[i] for i in perm)
        reaches = _lin_reaches(reactant, order, target_smi)
        if reaches != (order in lin_set):
            return False
    return True


def _correct_plan_by_apply_replay(nodes: list, reactant, target_smi: str) -> StepPlan | None:
    """Last resort: rebuild Deps precedes from which leaf orders reach T.

    Corrects the dependency graph for emission — does not filter or drop a walk.
    """
    from itertools import permutations

    n = len(nodes)
    if n < 2 or n > 4:
        return None
    good = []
    for perm in permutations(range(n)):
        order = [nodes[i] for i in perm]
        if _lin_reaches(reactant, order, target_smi):
            good.append(perm)
    if not good:
        return None
    return _plan_from_good_index_orders(nodes, good)


def and_cleave_plan(steps, mols=None, reactant=None, target_smi: str | None = None) -> StepPlan:
    """Emit a pathway plan from a Required walk.

    When the walk metabolite chain ``mols`` is provided (``len == n_steps+1``),
    required order is inferred from rule SMARTS match-before/after on those
    structures (plus AtomRef ``added_by``), then materialized as :class:`Deps`.
    If that graph is not exact vs ``Linearization.apply`` (n≤4), replay
    metabolism as a last resort and emit a corrected ``Deps`` — never as a
    filter that drops the walk.
    ``reactant`` alone is a one-mol chain prefix. Without mols, AtomRef /
    cleavage heuristics only.
    """
    nodes = [_step_from_path_item(s) for s in (steps or ())]
    if not nodes:
        return StepPlan(())
    if len(nodes) == 1:
        return StepPlan(nodes)

    chain = _walk_mol_chain(nodes, mols)
    if not chain and reactant is not None:
        chain = [reactant]

    if len(chain) >= len(nodes) + 1:
        plan = _plan_from_chain_smarts(nodes, chain)
    else:
        plan = _plan_from_dep_layers(nodes)

    # Last resort: correct deps so the emitted plan matches reaching orders.
    root = chain[0] if chain else reactant
    tgt = target_smi
    if tgt is None and len(chain) >= len(nodes) + 1:
        try:
            tgt = _canon(chain[-1])
        except Exception:
            tgt = None
    if root is not None and tgt is not None and 2 <= len(nodes) <= 4:
        if not _plan_lin_set_exact(plan, nodes, root, tgt):
            corrected = _correct_plan_by_apply_replay(nodes, root, tgt)
            if corrected is not None:
                plan = corrected
    return plan


def _added_by_json_fp(site_obj: dict) -> tuple:
    """Fingerprint an AtomRef JSON site (new tuple form or legacy at=)."""
    ab = site_obj.get("added_by")
    depth = int(site_obj.get("depth", 0) or 0)
    if isinstance(ab, str):
        return (ab, tuple(sorted(site_obj.get("at") or ())), depth)
    if isinstance(ab, (list, tuple)) and len(ab) == 2:
        return (ab[0], tuple(sorted(ab[1] or ())), depth)
    return (ab, (), depth)


def _origin_json_fp(site_obj: dict) -> tuple:
    """Fingerprint an origin AtomRef JSON object (optional depth)."""
    return (int(site_obj.get("origin", -1)), int(site_obj.get("depth", 0) or 0))


def _site_json_fp(s) -> tuple:
    if isinstance(s, int):
        return ("o", s, 0)
    if isinstance(s, dict):
        if "origin" in s and s["origin"] is not None:
            o, d = _origin_json_fp(s)
            return ("o", o, d)
        return ("a",) + _added_by_json_fp(s)
    return ("x", repr(s))


def _canonical_plan_key(plan: StepPlan) -> tuple:
    """Order-insensitive fingerprint so And-equivalent plans collide."""
    data = plan.to_json() if plan is not None else {}
    return _canonicalize_plan_json(data)


def _canonicalize_plan_json(data):
    if data is None:
        return ()
    if isinstance(data, (int, str, bool)):
        return data
    if isinstance(data, list):
        return tuple(_canonicalize_plan_json(x) for x in data)
    if not isinstance(data, dict):
        return str(data)
    op = data.get("op")
    # Deps embeds steps as {rule, site} without op=step.
    if op == "step" or (
        op is None and "rule" in data and "site" in data and "children" not in data
        and "steps" not in data
    ):
        site = data.get("site", [])
        site_fp = tuple(sorted((_site_json_fp(s) for s in site), key=repr))
        pathways = tuple(sorted(data.get("pathways") or ()))
        return ("step", data.get("rule"), site_fp, pathways)
    if op == "and":
        kids = [_canonicalize_plan_json(c) for c in data.get("children", ())]
        return ("and", tuple(sorted(kids, key=repr)))
    if op == "or":
        kids = [_canonicalize_plan_json(c) for c in data.get("children", ())]
        return ("or", tuple(sorted(kids, key=repr)))
    if op == "seq":
        kids = [_canonicalize_plan_json(c) for c in data.get("children", ())]
        return ("seq", tuple(kids))
    if op == "deps" or (op in (None, "plan") and "steps" in data and "children" not in data):
        steps = [_canonicalize_plan_json(s) for s in data.get("steps", ())]
        # Order-insensitive: node identity (canon steps) + canonical edges.
        # Graph reduction alone does not encode Step identity.
        order = sorted(range(len(steps)), key=lambda i: repr(steps[i]))
        inv = {old: new for new, old in enumerate(order)}
        canon_steps = tuple(steps[i] for i in order)
        precedes = []
        for pair in data.get("precedes") or ():
            a, b = int(pair[0]), int(pair[1])
            if a in inv and b in inv:
                precedes.append((inv[a], inv[b]))
        try:
            canon_edges = canonical_dependency_edges(len(canon_steps), precedes)
        except Exception:
            canon_edges = tuple(sorted(set(precedes)))
        return ("deps", canon_steps, canon_edges)
    # Unknown / legacy
    return ("raw", json.dumps(data, sort_keys=True, default=str))


def _site_atoms(site) -> frozenset:
    """Return comparable atom keys (ints) from a walk / formation site."""
    if isinstance(site, tuple) and len(site) >= 2 and isinstance(site[0], str):
        site = site[1]
    return _site_key_set(site)


def _steps_key(steps) -> tuple:
    """Fingerprint Required walk sites for seen / hit dedup."""
    out = []
    for item in steps or ():
        name, site = item[0], item[1]
        out.append((name, frozenset(_site_key_set(site))))
    return tuple(out)


def _opens_key(opens) -> tuple:
    return tuple(frozenset(_site_key_set(s)) for s in (opens or ()))


def _pack_hit(smi_list, steps, mols, *, maybe=None) -> PathOutcome:
    tgt = None
    if smi_list:
        try:
            tgt = _canon(smi_list[-1])
        except Exception:
            tgt = None
    plan = and_cleave_plan(steps, mols=mols, target_smi=tgt)
    if maybe is None:
        maybe = MaybeFilter()
    elif not isinstance(maybe, MaybeFilter):
        maybe = MaybeFilter.from_sides(maybe)
    return PathOutcome(
        plan=plan,
        maybe=maybe,
        steps=tuple(steps or ()),
        smiles=tuple(smi_list or ()),
        mols=tuple(mols or ()),
    )


def _require_mol(piece, where):
    """Cohort pieces are molecules. Anything else means the pipeline broke."""
    if isinstance(piece, Chem.Mol):
        return piece
    raise TypeError(
        "%s: expected a molecule, got %s" % (where, type(piece).__name__)
    )


def _num_heavy(mol, where) -> int:
    try:
        return int(mol.GetNumHeavyAtoms())
    except Exception as err:
        raise TypeError("%s is not a molecule" % (where,)) from err


def _product_pieces(mol_or_list):
    if mol_or_list is None:
        return []
    if isinstance(mol_or_list, (list, tuple)):
        pieces = []
        for m in mol_or_list:
            pieces.extend(_product_pieces(m))
        return pieces
    try:
        frags = Chem.GetMolFrags(mol_or_list, asMols=True, sanitizeFrags=False)
    except Exception:
        # Unsanitized mol: keep the whole structure. A non-mol is a bug.
        return [_require_mol(mol_or_list, "product piece")]
    if len(frags) <= 1:
        return [_require_mol(mol_or_list, "product piece")]
    return list(frags)


def _expansion_cohorts(parent, products):
    """Split a reaction product list into cleavage cohorts."""
    pieces = []
    for item in products or []:
        pieces.extend(_product_pieces(item))
    if not pieces:
        return []
    if len(pieces) == 1:
        return [pieces]

    parent_ha = _num_heavy(parent, "expansion parent")
    if all(p.GetNumHeavyAtoms() >= parent_ha - 1 for p in pieces):
        return [[p] for p in pieces]

    if len(pieces) % 2 == 0:
        return [pieces[i : i + 2] for i in range(0, len(pieces), 2)]

    return [pieces]


def _cleavage_sides_for_kept(site, kept, cohort, opens=()) -> tuple:
    """Bags for discarded siblings when ``cohort`` is a multi-piece cleavage."""
    pieces = []
    for item in cohort:
        pieces.extend(_product_pieces(item))
    if len(pieces) < 2:
        return ()
    kept_smi = _canon(kept)
    site_f = frozenset(_site_key_set(site))
    open_f = tuple(frozenset(_site_key_set(o)) for o in (opens or ()))
    bags = []
    seen = set()
    for piece in pieces:
        smi = _canon(piece)
        if smi == kept_smi or smi in seen:
            continue
        seen.add(smi)
        bags.append(CleavageSide(site=site_f, side=smi, opens=open_f))
    return tuple(bags)


def _is_cleavage_rule(rule) -> bool:
    return bool(
        getattr(rule, "is_cleavage", lambda: False)()
        or getattr(rule, "cleave_alone", lambda: False)()
    )


def _advance_opens_and_bags(rule, site, kept, cohort, bags, opens):
    """Update bags / open_sites after expanding ``kept`` from ``cohort``."""
    site_f = frozenset(_site_key_set(site))
    if len(cohort) >= 2:
        new_bags = bags + _cleavage_sides_for_kept(site, kept, cohort, opens=opens)
        return new_bags, ()
    if _is_cleavage_rule(rule) and len(cohort) == 1:
        return bags, opens + (site_f,)
    return bags, opens


def _rules_for_mol(active, mol, product):
    """When T is smaller, expand cleavage peers only."""
    smaller_target = _num_heavy(product, "product") < _num_heavy(mol, "reactant")
    if not smaller_target:
        return list(active)

    peers = [
        r
        for r in active
        if getattr(r, "is_cleavage", lambda: False)()
        or getattr(r, "cleave_alone", lambda: False)()
    ]
    return peers if peers else list(active)


def _depth_remaining(depth_left) -> bool:
    return depth_left is None or depth_left > 0


def _consume_depth(depth_left, hops: int = 1):
    if depth_left is None:
        return None
    return depth_left - hops


@unstable
def find_path(
    reactant,
    product,
    ruleset="PhaseOneQF",
    depth: Optional[int] = None,
    expand_phase1_plans: bool = True,
    max_paths: Optional[int] = None,
    max_expansions: Optional[int] = DEFAULT_MAX_EXPANSIONS,
    counters: Optional[PathSearchCounters] = None,
    search: str = "bfs",
    **kwargs,
) -> Iterator[PathOutcome]:
    """Yield :class:`PathOutcome` formation routes.

    ``plan`` is Required-only. ``maybe`` holds declarative :class:`CleavageSide`
    bags (ring-open then bifurcate spans both sites). Maybe is not searched.

    Defaults: ``depth=None`` (unbounded hops); ``max_expansions`` caps billed
    work (:meth:`PathSearchCounters.billed` = linearizations + site applies),
    default :data:`DEFAULT_MAX_EXPANSIONS`.
    """
    if counters is None:
        counters = PathSearchCounters()

    if isinstance(reactant, str):
        reactant = Chem.MolFromSmiles(reactant)
    if isinstance(product, str):
        product = Chem.MolFromSmiles(product)
    if reactant is None or product is None:
        raise ValueError("invalid reactant or product")

    kwargs.pop("maybe_prefixes", None)
    kwargs.pop("max_maybe", None)
    kwargs.pop("maybe", None)

    rules = _flatten_rules(ruleset)
    ctx = PathContext.from_mols(reactant, product)
    active = _active_rules(rules)
    target_smi = _canon(product)

    from .path_context import unreachable_new_elements
    from .utils import path_counter_scope

    blocked = unreachable_new_elements(active, reactant, product)
    if blocked:
        # Impossible by element palette — do not burn budget.
        return

    if _canon(reactant) == target_smi:
        yield _pack_hit([target_smi], [], [copy_mol(reactant)])
        return

    seen_path_keys = set()
    n_hits = 0
    with path_counter_scope(counters):
        for hit in _guided_mol_search(
            reactant,
            product,
            target_smi,
            active,
            ctx,
            depth=depth,
            expand_phase1_plans=expand_phase1_plans,
            counters=counters,
            search=search,
            max_expansions=max_expansions,
            **kwargs,
        ):
            smi_list, steps, mols, bags = hit
            outcome = _pack_hit(smi_list, steps, mols, maybe=bags)
            # Deduplicate identical StepPlans (And-equivalent walks collapse).
            key = (
                target_smi,
                _canonical_plan_key(outcome.plan),
                str(outcome.maybe),
            )
            if key in seen_path_keys:
                continue
            seen_path_keys.add(key)
            yield outcome
            n_hits += 1
            if max_paths is not None and n_hits >= max_paths:
                return


find_path_guided = find_path


def _can_dearomatize(rule) -> bool:
    fn = getattr(rule, "can_dearomatize", None)
    if not callable(fn):
        return False
    try:
        return bool(fn())
    except Exception:
        return False


def _ref_atom_idxs(ref) -> frozenset:
    origin = getattr(ref, "origin", None)
    if origin is not None:
        try:
            return frozenset([int(origin)])
        except (TypeError, ValueError):
            pass
    added = getattr(ref, "added_by", None)
    if added and len(added) >= 2:
        out = []
        for idx in added[1]:
            try:
                out.append(int(idx))
            except (TypeError, ValueError):
                continue
        return frozenset(out)
    return frozenset()


def _expansion_site_groups(site, expr=None):
    """Atom-index groups for a hop site and, when present, each plan step."""
    groups = []
    atoms = _site_atoms(site)
    if atoms:
        groups.append(frozenset(int(i) for i in atoms))
    if expr is None:
        return groups
    try:
        lins = list(pathway_linearizations(expr))
    except Exception:
        lins = []
    steps = []
    if lins:
        for lin in lins:
            steps.extend(lin.steps)
    else:
        steps.extend(getattr(expr, "steps", ()) or ())
    for step in steps:
        idxs = set()
        for ref in getattr(step, "site", ()) or ():
            idxs |= _ref_atom_idxs(ref)
        if idxs:
            groups.append(frozenset(idxs))
    return groups


def _system_with_ends(systems, mol) -> set:
    """Aromatic-system atoms plus the atoms bonded to them.

    Quinone and dehydrogenation ends are the OH / NH atoms on the ring, which
    are not themselves aromatic. Atoms two bonds out (a methoxy carbon) are
    not included.
    """
    union = set().union(*systems) if systems else set()
    extended = set(union)
    for idx in union:
        try:
            atom = mol.GetAtomWithIdx(int(idx))
        except Exception:
            continue
        for nbr in atom.GetNeighbors():
            extended.add(nbr.GetIdx())
    return extended


def _site_on_systems(groups, systems, mol) -> bool:
    """True when the expansion site touches a system that must be dearomatized.

    Only the site is checked. A later plan step, such as removing a
    substituent, does not drop a site that is already on the system.
    """
    if not systems or not groups:
        return True
    return bool(set(groups[0]) & _system_with_ends(systems, mol))


def _guided_mol_search(
    reactant,
    product,
    target_smi,
    active,
    ctx,
    depth,
    expand_phase1_plans,
    counters,
    search="bfs",
    max_expansions=None,
    **kwargs,
):
    """Yield ``(smiles_path, step_path, mol_path, cleavage_side_bags)``."""
    start = copy_mol(reactant)
    # (mol, smi_path, step_path, mol_path, depth_left, bags, open_sites, boosted)
    # boosted: this node was reached by a match-boundary dearomatizing edit, so
    # its own boundary children jump ahead of older siblings.
    item0 = (start, [_canon(start)], [], [start], depth, (), (), False)
    if search == "dfs":
        frontier = [item0]
        pop = frontier.pop
    elif search == "bfs":
        frontier = deque([item0])
        pop = frontier.popleft
    else:
        raise ValueError("search must be 'bfs' or 'dfs', got %r" % (search,))

    seen = {(_canon(start), (), ())}
    counters.nodes_enqueued += 1

    from .path_context import dearomatization_systems, site_match_boundary_rank

    while frontier:
        if not counters.under_budget(max_expansions):
            return
        mol, smi_path, step_path, mol_path, depth_left, bags, opens, boosted = pop()
        if not _depth_remaining(depth_left):
            continue

        ctx = PathContext.from_mols(mol, product)
        systems = dearomatization_systems(mol, ctx)
        needs_dear = bool(systems)
        rules = _rules_for_mol(active, mol, product)
        if needs_dear:
            rules = [r for r in rules if _can_dearomatize(r)] + [
                r for r in rules if not _can_dearomatize(r)
            ]

        jump_buf = []
        rest_buf = []

        def enqueue(item, jump, _jump_buf=jump_buf, _rest_buf=rest_buf):
            (_jump_buf if jump else _rest_buf).append(item)

        for rule in rules:
            if not counters.under_budget(max_expansions):
                return
            if not rule.could_help(mol, product, ctx):
                counters.nodes_pruned += 1
                counters._sync()
                continue
            counters.rule_expansions += 1
            counters._sync()
            # A raise here is a broken rule, not a chemical prune.
            items = list(
                rule.enumerate_for_path(
                    mol,
                    ctx,
                    expand_phase1_plans=expand_phase1_plans,
                    **kwargs,
                )
            )
            prepared = []
            for kind, *payload in items:
                if kind == "plan":
                    expr, site = payload
                elif kind == "hop":
                    site = payload[0]
                    expr = None
                else:
                    prepared.append(((1, 1, 1, 1), kind, payload, False))
                    continue
                groups = _expansion_site_groups(site, expr)
                if (
                    needs_dear
                    and _can_dearomatize(rule)
                    and not _site_on_systems(groups, systems, mol)
                ):
                    counters.nodes_pruned += 1
                    counters._sync()
                    continue
                atoms = set()
                for group in groups:
                    atoms |= set(group)
                boundary = site_match_boundary_rank(ctx, atoms)
                rule_rank = 0 if needs_dear and _can_dearomatize(rule) else 1
                on_boundary = boundary[0] == 0 or boundary[1] == 0
                jump = (boosted or rule_rank == 0) and on_boundary and boundary[2] == 0
                prepared.append(((rule_rank,) + boundary, kind, payload, jump))
            prepared.sort(key=lambda row: row[0])
            for _rank, kind, payload, jump in prepared:
                if kind == "plan":
                    expr, site = payload
                    counters.plans_expanded += 1
                    counters.sites_considered += 1
                    counters._sync()
                    yield from _apply_plan_branches(
                        mol,
                        expr,
                        site,
                        rule,
                        product,
                        target_smi,
                        smi_path,
                        step_path,
                        mol_path,
                        depth_left,
                        bags,
                        opens,
                        enqueue,
                        jump,
                        seen,
                        counters,
                        max_expansions=max_expansions,
                    )
                elif kind == "hop":
                    site, products = payload
                    if not counters.under_budget(max_expansions):
                        return
                    counters.sites_considered += 1
                    counters.site_applies += 1
                    counters._sync()
                    for cohort in _expansion_cohorts(mol, products):
                        for child in cohort:
                            _require_mol(child, "hop product")
                            if not rule.child_may_reach(mol, child, product, ctx):
                                counters.nodes_pruned += 1
                                continue
                            counters.mol_edits += 1
                            counters._sync()
                            child_smi = _canon(child)
                            # Site idxs are GetIdx on ``mol`` at its current
                            # frame. Prefer depth-0 origins when the atom
                            # existed on the reactant; created atoms keep
                            # mid-frame depth.
                            stamped = _origin_refs(_site_atoms(site), mol)
                            new_steps = step_path + [(rule.name, stamped)]
                            new_bags, new_opens = _advance_opens_and_bags(
                                rule, stamped, child, cohort, bags, opens
                            )
                            if child_smi == target_smi:
                                yield (
                                    smi_path + [child_smi],
                                    new_steps,
                                    mol_path + [child],
                                    new_bags,
                                )
                                continue
                            if rule.is_terminal_product(child) or has_star_conjugate(
                                child
                            ):
                                continue
                            rem = _consume_depth(depth_left, 1)
                            if not _depth_remaining(rem):
                                continue
                            sk = (
                                child_smi,
                                _steps_key(new_steps),
                                _opens_key(new_opens),
                            )
                            if sk in seen:
                                counters.nodes_pruned += 1
                                continue
                            seen.add(sk)
                            enqueue(
                                (
                                    child,
                                    smi_path + [child_smi],
                                    new_steps,
                                    mol_path + [child],
                                    rem,
                                    new_bags,
                                    new_opens,
                                    jump,
                                ),
                                jump,
                            )
                            counters.nodes_enqueued += 1
                else:
                    raise ValueError(
                        "enumerate_for_path yielded %r; expected 'plan' or 'hop'"
                        % (kind,)
                    )
        for item in reversed(jump_buf):
            if search == "dfs":
                frontier.append(item)
            else:
                frontier.appendleft(item)
        if search == "dfs":
            frontier[:0] = rest_buf
        else:
            frontier.extend(rest_buf)


def _apply_plan_branches(
    mol,
    expr,
    site,
    rule,
    product,
    target_smi,
    smi_path,
    step_path,
    mol_path,
    depth_left,
    bags,
    opens,
    enqueue,
    jump,
    seen,
    counters,
    max_expansions=None,
):
    for lin in pathway_linearizations(expr):
        if not counters.under_budget(max_expansions):
            return
        counters.linearizations_applied += 1
        counters._sync()
        n_hops = len(lin.steps)
        if depth_left is not None and n_hops > depth_left:
            continue
        try:
            # Miss: this order does not fire on this mol. Not a broken search.
            products = lin.apply(mol)
        except Exception:
            continue
        step_label = [(s.rule, s.site) for s in lin.steps]
        ctx = PathContext.from_mols(mol, product)
        stamped_site = _origin_refs(_site_atoms(site), mol)
        for cohort in _expansion_cohorts(mol, products):
            for child in cohort:
                _require_mol(child, "plan product")
                if not rule.child_may_reach(mol, child, product, ctx):
                    counters.nodes_pruned += 1
                    continue
                counters.mol_edits += 1
                counters._sync()
                child_smi = _canon(child)
                new_steps = step_path + step_label
                new_bags, new_opens = _advance_opens_and_bags(
                    rule, stamped_site, child, cohort, bags, opens
                )
                new_smi = smi_path + [child_smi]
                new_mols = mol_path + [child]
                if child_smi == target_smi:
                    yield (new_smi, new_steps, new_mols, new_bags)
                    continue
                if rule.is_terminal_product(child) or has_star_conjugate(child):
                    continue
                rem = _consume_depth(depth_left, n_hops)
                if not _depth_remaining(rem):
                    continue
                sk = (child_smi, _steps_key(new_steps), _opens_key(new_opens))
                if sk in seen:
                    counters.nodes_pruned += 1
                    continue
                seen.add(sk)
                enqueue(
                    (child, new_smi, new_steps, new_mols, rem, new_bags, new_opens, jump),
                    jump,
                )
                counters.nodes_enqueued += 1
