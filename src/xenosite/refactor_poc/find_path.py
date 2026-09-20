"""Reactant-to-target search over the proof-of-concept rules.

One atom diff drives ``filter_rules`` and ``filter_sites``. A rule is a
SMARTS pattern plus ``describe()``; it does not grow this search. The yield
is a phase-I :class:`~xenosite.forest.step_plan.Deps` plus the cleavage
fragments that were not expanded.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from typing import NamedTuple

from xenosite.forest.step_plan import AtomRef as AddedRef
from xenosite.forest.step_plan import Deps, Step, StepPlan
from xenosite.refactor_poc.rdkitutil import (
    Atom,
    Mol,
    as_mol,
    cannonicalize_order,
    canon_smiles,
    copy_mol,
    is_tracing,
    mcs_matches,
    mcs_target_matches,
    sanitize_catch,
    split_fragments,
)
from xenosite.refactor_poc.records import AtomRef, Site, SiteInfo, _flat_ints
from xenosite.refactor_poc.rules import (
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
    _as_site,
    forest_trace,
    install_forest,
)
from xenosite.refactor_poc.rulesets import RuleSet

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------


class PathCounters:
    """Billed work for one search.

    ``mol_edits`` counts a ``RunReactants`` or a kekulé overlay that a filter
    has already accepted. ``billed`` is edits plus nodes visited. It is a
    measurement, not a second search mode.
    """

    def __init__(self):
        self.rule_expansions = 0
        self.sites_considered = 0
        self.sites_skipped = 0
        self.mol_edits = 0
        self.sanitize_dropped = 0
        self.nodes = 0

    @property
    def billed(self):
        return self.mol_edits + self.nodes


# ---------------------------------------------------------------------------
# Outcome
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CleavageSide:
    """Fragment discarded at a cleavage. Not searched.

    ``site`` is the cleavage site on the parent. ``opens`` are earlier
    ring-open sites on the same walk, oldest first. ``side`` is the
    canonical SMILES of the discarded fragment.
    """

    site: frozenset[int]
    side: str
    opens: tuple[frozenset[int], ...] = ()

    def span_sites(self) -> tuple[frozenset[int], ...]:
        return self.opens + (self.site,)


def _site_keys(site: Site | Sequence[int | AddedRef | AtomRef] | None) -> frozenset[int | str]:
    if site is None:
        return frozenset()
    if isinstance(site, int):
        return frozenset({site})
    keys: set[int | str] = set()
    for item in site:
        if isinstance(item, frozenset):
            keys |= _site_keys(item)
            continue
        origin = getattr(item, "origin", None)
        if origin is not None:
            keys.add(int(origin))
            continue
        if isinstance(item, AtomRef):
            keys.add(item.idx)
            continue
        if isinstance(item, int):
            keys.add(item)
            continue
        keys.add(str(item))
    return frozenset(keys)


@dataclass(frozen=True)
class Maybe:
    """Uncleared cleavage fragments. Not a step, and not searched."""

    entries: tuple = ()

    def __bool__(self):
        return bool(self.entries)

    def sides(self):
        return tuple(entry.side for entry in self.entries)

    def allows(self, rule_name=None, site=None, *, side=None):
        """True when ``side`` is a discarded fragment, or ``site`` overlaps one.

        The bifurcating cleavage site itself does not pass: that step is
        already in the plan. A different reaction on those atoms does.
        """

        if not self.entries:
            return False
        if side is not None:
            want = canon_smiles(side)
            return any(entry.side == want for entry in self.entries)
        if site is None:
            return True
        keys = _site_keys(site)
        if not keys:
            return False
        for entry in self.entries:
            if keys == _site_keys(entry.site):
                continue
            for span in entry.span_sites():
                if keys & _site_keys(span):
                    return True
        return False


@dataclass(frozen=True)
class PathOutcome:
    """Required phase-I plan, plus cleavage fragments left on :class:`Maybe`."""

    # Deps is @unstable; pyright does not see it as a class. StepPlan is the API.
    plan: StepPlan
    maybe: Maybe
    smiles: str

    def allows(
        self,
        rule_name: str | None = None,
        site: Site | Sequence[int | AddedRef | AtomRef] | None = None,
        *,
        side: Mol | str | None = None,
    ) -> bool:
        return self.maybe.allows(rule_name, site, side=side)


# ---------------------------------------------------------------------------
# Atom diff
# ---------------------------------------------------------------------------


def _hydrogens(atom: Atom):
    try:
        return atom.GetTotalNumHs()
    except RuntimeError:
        return 0


class AtomDiff:
    """Local reactant-to-target differences. This is the whole heuristic.

    A mapped atom can need an element, lose aromaticity, or lose hydrogens.
    An unmapped reactant atom, or a mapped bond that is absent in the target,
    is cleavage. ``filter_rules`` sees the molecule-level summary. ``filter_sites``
    sees one site.
    """

    def __init__(
        self,
        reactant: Mol,
        target: Mol,
        mapping: dict[int, int],
        needs_oxygen,
        needs_carbonyl,
        needs_alcohol,
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta,
        n_extra,
        bond_order_mismatches,
        bond_raises,
    ):
        self.reactant = reactant
        self.target = target
        self.mapping = mapping
        self.needs_oxygen = frozenset(needs_oxygen)
        self.needs_carbonyl = frozenset(needs_carbonyl)
        self.needs_alcohol = frozenset(needs_alcohol)
        self.cleaved = frozenset(cleaved)
        self.cleavage_bonds = set(cleavage_bonds)
        self.loses_aromaticity = frozenset(loses_aromaticity)
        self.h_delta = dict(h_delta)
        self.n_extra = n_extra
        self.bond_order_mismatches = bond_order_mismatches
        self.bond_raises = set(bond_raises)
        self.reactant_heavy = reactant.GetNumHeavyAtoms()
        self.target_heavy = target.GetNumHeavyAtoms()
        self.mappings: tuple[dict[int, int], ...] = (dict(mapping),)
        self._view_costs: tuple[int, ...] | None = None

    @property
    def target_smaller(self):
        return self.target_heavy < self.reactant_heavy

    @property
    def has_cleavage(self):
        return bool(self.cleaved or self.cleavage_bonds)

    @property
    def h_loss(self):
        return any(delta < 0 for delta in self.h_delta.values())

    def site_is_cleavage(self, atoms):
        """True when ``atoms`` is the bond that separates kept from gone.

        Touching a leaving-group atom is not enough. The site has to be
        the bridge itself, or one atom of that bridge.
        """

        atoms = set(atoms)
        if not atoms:
            return False
        for bond in self.cleavage_bonds:
            if bond <= atoms or atoms <= bond:
                return True
        return False

    def cost(self):
        """How much of the reactant still disagrees with the target.

        Used only to refuse a child that moved away. Not a ranker.
        When several placements are open, a step that gets closer on any
        one of them is not refused.
        """

        if self._view_costs is not None:
            return min(self._view_costs)
        return self._field_cost()

    def _field_cost(self):
        h_off = sum(1 for delta in self.h_delta.values() if delta)
        return (
            3 * len(self.cleaved)
            + 3 * self.n_extra
            + 2 * len(self.needs_oxygen)
            + len(self.loses_aromaticity)
            + h_off
            + 3 * len(self.cleavage_bonds)
            + self.bond_order_mismatches
        )


def _mapping_score(reactant: Mol, target: Mol, r_match, t_match):
    """Prefer a pairing that keeps rings on rings and bond orders put."""

    score = 0
    paired = {}
    for r_idx, t_idx in zip(r_match, t_match):
        paired[r_idx] = t_idx
        ra = reactant.GetAtomWithIdx(r_idx)
        ta = target.GetAtomWithIdx(t_idx)
        if ra.GetAtomicNum() != ta.GetAtomicNum():
            score -= 50
        if ra.GetIsAromatic() == ta.GetIsAromatic():
            score += 4
        else:
            score -= 12
        if ra.IsInRing() == ta.IsInRing():
            score += 3
        else:
            score -= 8
        if _hydrogens(ra) == _hydrogens(ta):
            score += 1
    for bond in reactant.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in paired or j not in paired:
            continue
        other = target.GetBondBetweenAtoms(paired[i], paired[j])
        if other is None:
            score -= 4
            continue
        if abs(other.GetBondTypeAsDouble() - bond.GetBondTypeAsDouble()) < 0.2:
            score += 1
        else:
            score -= 1
    return score


def _mappings(reactant: Mol, target: Mol) -> tuple[dict[int, int], ...]:
    """One alignment per placement. The cache's embeddings are all read.

    Reorderings of the same reactant atoms are one placement. The winning
    alignment is the one :func:`_mapping_score` likes. A shorter remainder
    is paired only with target embeddings of that same length.
    """

    reactant_hits = mcs_matches(reactant, target).embeddings
    target_hits = mcs_target_matches(reactant, target).embeddings
    by_size: dict[int, list[tuple[int, ...]]] = {}
    for hit in target_hits:
        by_size.setdefault(len(hit), []).append(hit)
    if not reactant_hits or not by_size:
        return ()

    best: dict[frozenset[int], tuple[int, dict[int, int]]] = {}
    for r_match in reactant_hits:
        mates = by_size.get(len(r_match))
        if not mates:
            continue
        key = frozenset(r_match)
        for t_match in mates:
            score = _mapping_score(reactant, target, r_match, t_match)
            held = best.get(key)
            aligned = {r: t for r, t in zip(r_match, t_match)}
            if held is None or score > held[0]:
                best[key] = (score, aligned)
    ranked = sorted(best.values(), key=lambda item: item[0], reverse=True)
    return tuple(item[1] for item in ranked)


def _merge_views(views: tuple[AtomDiff, ...]) -> AtomDiff:
    """One diff whose filters can see every placement.

    ``mapping`` stays the best-scoring alignment. ``cost`` is the minimum
    across placements, so a step that helps a second ring is not refused
    for failing to help the first.
    """

    primary = views[0]
    primary._view_costs = tuple(view._field_cost() for view in views)
    primary.mappings = tuple(dict(view.mapping) for view in views)
    if len(views) == 1:
        return primary
    primary.needs_oxygen = frozenset().union(*(view.needs_oxygen for view in views))
    primary.needs_carbonyl = frozenset().union(
        *(view.needs_carbonyl for view in views)
    )
    primary.needs_alcohol = frozenset().union(*(view.needs_alcohol for view in views))
    primary.cleaved = frozenset().union(*(view.cleaved for view in views))
    bonds: set[frozenset[int]] = set()
    raises: set[frozenset[int]] = set()
    for view in views:
        bonds.update(view.cleavage_bonds)
        raises.update(view.bond_raises)
    primary.cleavage_bonds = bonds
    primary.bond_raises = raises
    primary.loses_aromaticity = frozenset().union(
        *(view.loses_aromaticity for view in views)
    )
    h_delta: dict[int, int] = {}
    for view in views:
        for atom, delta in view.h_delta.items():
            held = h_delta.get(atom)
            if held is None or delta < held:
                h_delta[atom] = delta
    primary.h_delta = h_delta
    return primary


def atom_diff(reactant: Mol | str, target: Mol | str) -> AtomDiff:
    """Pair reactant atoms with target atoms and record the local change.

    Every cached placement is kept on ``mappings``. ``mapping`` is the
    best-scoring one. String inputs are parsed with ``MolFromSmiles``.
    Indexes then refer to that parse, which is stable for a given SMILES.
    """

    reactant_mol = as_mol(reactant)
    target_mol = as_mol(target)
    mappings = _mappings(reactant_mol, target_mol) or ({},)
    views = tuple(
        _diff_for(reactant_mol, target_mol, mapping) for mapping in mappings
    )
    return _merge_views(views)


def _diff_for(reactant: Mol, target: Mol, mapping: dict[int, int]) -> AtomDiff:
    image = set(mapping.values())

    needs_oxygen = set()
    needs_carbonyl = set()
    needs_alcohol = set()
    for atom in target.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in image:
            continue
        for neighbor in atom.GetNeighbors():
            bond = target.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond is None:
                continue
            for r_idx, t_idx in mapping.items():
                if t_idx != neighbor.GetIdx():
                    continue
                needs_oxygen.add(r_idx)
                if bond.GetBondTypeAsDouble() >= 1.5:
                    needs_carbonyl.add(r_idx)
                else:
                    needs_alcohol.add(r_idx)

    cleaved = set()
    for atom in reactant.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.GetIdx() not in mapping:
            cleaved.add(atom.GetIdx())

    cleavage_bonds = set()
    loses_aromaticity = set()
    h_delta = {}
    bond_order_mismatches = 0
    bond_raises = set()
    for r_idx, t_idx in mapping.items():
        ra = reactant.GetAtomWithIdx(r_idx)
        ta = target.GetAtomWithIdx(t_idx)
        if ra.GetIsAromatic() and not ta.GetIsAromatic():
            loses_aromaticity.add(r_idx)
        h_delta[r_idx] = _hydrogens(ta) - _hydrogens(ra)

    for bond in reactant.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in mapping or j not in mapping:
            # A cut that can reach the target joins a kept atom to a leaving
            # atom. Bonds inside the leaving group are not that cut.
            if (i in mapping) != (j in mapping):
                cleavage_bonds.add(frozenset((i, j)))
            continue
        other = target.GetBondBetweenAtoms(mapping[i], mapping[j])
        if other is None:
            cleavage_bonds.add(frozenset((i, j)))
            continue
        delta = other.GetBondTypeAsDouble() - bond.GetBondTypeAsDouble()
        if abs(delta) >= 0.2:
            bond_order_mismatches += 1
        if delta >= 0.2:
            bond_raises.add(frozenset((i, j)))

    n_extra = sum(
        1
        for atom in target.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetIdx() not in image
    )
    return AtomDiff(
        reactant,
        target,
        mapping,
        needs_oxygen,
        needs_carbonyl,
        needs_alcohol,
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta,
        n_extra,
        bond_order_mismatches,
        bond_raises,
    )


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------


def _span_values(span, key, default):
    if key not in span:
        return (default,)
    value = span[key]
    if isinstance(value, tuple):
        return value
    return (value,)


def _any_span(span, key, pred, default):
    return any(pred(value) for value in _span_values(span, key, default))


def _all_span(span, key, pred, default):
    values = _span_values(span, key, default)
    return bool(values) and all(pred(value) for value in values)


def _effect_adds_oxygen(effect):
    return "O" in (effect.get("adds") or "") or "O" in (effect.get("needs") or "")


def _pattern_could_help(info, diff):
    """``filter_rules`` sees ``span`` before any match."""

    span = info.get("span") or {}
    can_cleave = _any_span(span, "cleaves", bool, False)
    # The kept piece is smaller. A non-cleaving pattern cannot get there,
    # and editing the intact parent is the cost this filter exists to avoid.
    if diff.target_smaller and not can_cleave:
        return False
    if (
        _all_span(span, "adds", lambda value: "O" in (value or ""), "")
        and not can_cleave
        and not diff.needs_oxygen
    ):
        return False
    if _all_span(span, "cleaves", bool, False) and not diff.has_cleavage:
        return False
    if (
        _all_span(span, "dearomatizes", bool, False)
        and not diff.loses_aromaticity
    ):
        return False
    drops_h_only = (
        not can_cleave
        and not _any_span(span, "adds", lambda value: bool(value), "")
        and _any_span(span, "removes", lambda value: bool(value) and "H" in value, "")
    )
    if drops_h_only and not diff.h_loss and not diff.loses_aromaticity:
        return False
    return True


def _alkyl_bond_raises(mol: Mol, atom_idx, diff):
    """True when an exocyclic C-C bond at ``atom_idx`` is higher in the target."""

    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 6 or neighbor.GetIsAromatic():
            continue
        if frozenset((atom_idx, neighbor.GetIdx())) in diff.bond_raises:
            return True
    return False


def _site_could_help(site: Site, info: SiteInfo, diff: AtomDiff, mol: Mol) -> bool:
    """``filter_sites`` sees one resolved effect and the local atom diff."""

    effect = info["options"]
    atoms = _flat_ints(site)
    if effect.get("cleaves"):
        return diff.site_is_cleavage(atoms)

    if "ends" in info:
        ends = info["ends"]
        end_atoms = info["end_atoms"]
        for atom, end in zip(end_atoms, ends):
            if _effect_adds_oxygen(end) and atom not in diff.needs_oxygen:
                return False
            # An alkyl partner turns the ring bond into an exocyclic double
            # bond (methide). Skip it unless that C-C bond is higher in the target.
            if (end.get("partner") or "") == "C" and not _alkyl_bond_raises(
                mol, atom, diff
            ):
                return False
    elif _effect_adds_oxygen(effect) and not effect.get("dearomatizes"):
        oxygen_sites = [atom for atom in atoms if atom in diff.needs_oxygen]
        if not oxygen_sites:
            return False
        # The local change is a carbonyl on a ring that stops being aromatic.
        # A bare hydroxylation does not do that; the dearomatizing edit does.
        if all(
            atom in diff.needs_carbonyl and atom in diff.loses_aromaticity
            for atom in oxygen_sites
        ):
            return False

    path_ends: frozenset[int] | tuple[()] = (
        info["path_ends"] if "path_ends" in info else ()
    )
    if effect.get("dearomatizes"):
        scope = atoms | set(path_ends)
        if not (scope & set(diff.loses_aromaticity)):
            return False

    removes = effect.get("removes") or ""
    if (
        isinstance(removes, str)
        and "H" in removes
        and not _effect_adds_oxygen(effect)
        and not effect.get("cleaves")
    ):
        scope = atoms | set(path_ends)
        loses_h = any(diff.h_delta.get(atom, 0) < 0 for atom in scope)
        if not loses_h and not (scope & set(diff.loses_aromaticity)):
            return False
    return True


def _filters(diff, enabled, mol: Mol):
    if not enabled:
        return (lambda rule, info: True), (lambda site, info: True)

    def filter_rules(rule, info):
        return _pattern_could_help(info, diff)

    def filter_sites(site, info):
        return _site_could_help(site, info, diff, mol)

    return filter_rules, filter_sites


def _rule_can_cleave(rule: object) -> bool:
    patterns: list[tuple[str, dict[str, object]]] = []
    for group in (getattr(rule, "smarts", None), getattr(rule, "endpoints", None)):
        if not group:
            continue
        patterns.extend(group)
    for _smarts, info in patterns:
        span = info.get("span") or {}
        if not isinstance(span, dict):
            continue
        if _any_span(span, "cleaves", bool, False):
            return True
    return False


# ---------------------------------------------------------------------------
# Steps and dependency
# ---------------------------------------------------------------------------


def _atom_ref(mol: Mol, idx):
    """Origin index, or ``added_by`` when this atom was created by a step.

    ``added_by`` on the record is a transform id (``R1``). The rule name and
    site are on ``atom_trace["additions"]``.
    """

    atom = mol.GetAtomWithIdx(idx)
    if not is_tracing(mol) or not atom.HasProp("forestLabel"):
        return idx
    trace = mol._forest["atom_trace"]
    record = trace["records"].get(atom.GetProp("forestLabel"))
    if record is None:
        return idx
    added = record.get("added_by")
    if not added:
        return idx
    if isinstance(added, str):
        detail = trace["additions"].get(added)
        if detail is None:
            return idx
    else:
        detail = added
    name = detail.get("name")
    if name is None:
        rule = detail.get("rule")
        name = rule if isinstance(rule, str) else getattr(rule, "name", None)
    site = detail.get("site") or ()
    if isinstance(site, int):
        site = (site,)
    if name is None:
        return idx
    return AddedRef(added_by=(name, frozenset(site)))


class _PlanStep(NamedTuple):
    """One phase-I step before it is handed to :class:`Deps`.

    ``site`` is a flat tuple of atom notes for that step: a known index,
    a forest :class:`AddedRef`, or a records :class:`AtomRef` for an atom
    that does not exist yet. That is not a reaction :class:`FutureSite`.
    """

    rule: str
    site: tuple[int | AddedRef | AtomRef, ...]


def _step(mol: Mol, rule_name: str, site: Site) -> _PlanStep:
    if isinstance(site, int):
        atoms: tuple[int, ...] = (site,)
    elif isinstance(site, tuple):
        atoms = site
    else:
        atoms = tuple(sorted(_flat_ints(site)))
    return _PlanStep(rule_name, tuple(_atom_ref(mol, idx) for idx in atoms))


def _steps_for(mol: Mol, info: SiteInfo):
    """Phase-I steps for one accepted edit.

    QuinoneFormation is not itself a step. The hop stands in for the
    hydroxylations that supply each oxygen and the dehydrogenation that
    follows them.
    """

    rule_name = info["rule"].name
    if rule_name is None:
        return ()
    if rule_name == "QuinoneFormation":
        return _quinone_phase1(mol, info)
    return (_step(mol, rule_name, info["site"]),)


def _bonded(mol: Mol, idx, atomic_num):
    """Neighbor of ``idx`` with this atomic number, if the mol already has one."""

    atom = mol.GetAtomWithIdx(idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == atomic_num:
            return neighbor.GetIdx()
    return None


def _quinone_phase1(mol: Mol, info: SiteInfo):
    """Hydroxylations that supply missing oxygens, then one dehydrogenation.

    An end that already carries oxygen keeps that atom. An end that
    ``needs`` oxygen becomes a Hydroxylation; the dehydrogenation site
    points at the oxygen that step would add. No ``end_maps`` field is
    required: the partner atom is the neighbor already on ``mol``.
    """

    if "ends" not in info:
        return (_step(mol, "Dehydrogenation", info["site"]),)
    ends = info["ends"]
    end_atoms = info["end_atoms"]
    hydroxylations = []
    dh_refs = []
    for end, atom in zip(ends, end_atoms):
        partner = end.get("partner") or ""
        adds_oxygen = "O" in (end.get("needs") or "") or (
            "O" in (end.get("adds") or "") and partner != "O"
        )
        if adds_oxygen:
            anchor = _atom_ref(mol, atom)
            hydroxylations.append(_PlanStep("Hydroxylation", (anchor,)))
            if isinstance(anchor, int):
                idx, depth = anchor, 0
            elif anchor.origin is not None:
                idx, depth = anchor.origin, anchor.depth
            else:
                idx, depth = atom, 0
            dh_refs.append(AtomRef(idx, "O", depth))
            continue
        atomic_num = {"O": 8, "N": 7, "C": 6, "S": 16}.get(partner)
        hetero = _bonded(mol, atom, atomic_num) if atomic_num else None
        if hetero is not None:
            dh_refs.append(_atom_ref(mol, hetero))
    if not dh_refs:
        return (_step(mol, "Dehydrogenation", info["site"]),)
    return tuple(hydroxylations) + (_PlanStep("Dehydrogenation", tuple(dh_refs)),)


def _anchor(item) -> int | None:
    if isinstance(item, AtomRef):
        return item.idx
    if isinstance(item, int):
        return item
    origin = getattr(item, "origin", None)
    if origin is None:
        return None
    return int(origin)


def _anchors(step) -> set[int]:
    return {anchor for anchor in (_anchor(item) for item in step.site) if anchor is not None}


def _forest_item(item, steps, later):
    """Hand :class:`Deps` a forest ref. The rule name is the earlier step's."""

    if not isinstance(item, AtomRef):
        return item
    for previous in steps[:later]:
        if item.idx in _anchors(previous):
            return AddedRef(added_by=(previous.rule, frozenset({item.idx})))
    return item.idx


def _deps(steps: Sequence[_PlanStep]) -> StepPlan:
    """A later step depends on an earlier one when its site names an atom that step added.

    An :class:`AtomRef` is that note: its ``idx`` is the atom the earlier
    step changed, and its ``element`` is what had to be added.
    """

    edges: list[tuple[int, int]] = []
    forest_steps: list[Step] = []
    for later, step in enumerate(steps):
        for item in step.site:
            if isinstance(item, AtomRef):
                for earlier, previous in enumerate(steps[:later]):
                    if item.idx in _anchors(previous):
                        edges.append((earlier, later))
                continue
            added = getattr(item, "added_by", None)
            if not added:
                continue
            rule_name, site = added
            wanted = frozenset(site)
            for earlier, previous in enumerate(steps):
                if previous.rule != rule_name:
                    continue
                if _anchors(previous) == wanted:
                    edges.append((earlier, later))
        forest_steps.append(
            Step(step.rule, frozenset(_forest_item(item, steps, later) for item in step.site))
        )
    # Deps is wrapped by @unstable; pyright does not treat the constructor as Deps/StepPlan.
    return Deps(tuple(forest_steps), edges)  # pyright: ignore[reportReturnType, reportCallIssue]


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------


def default_ruleset():
    """The poc catalog as one rule. Phase I is the ``Deps`` search yields."""

    return RuleSet(
        (Dealkylation, QuinoneFormation, Hydroxylation, Dehydrogenation),
        name="Poc",
    )


def _finish(parent: Mol, raw_products, info, counters):
    """Sanitize and trace each fragment. Failed sanitizes are dropped."""

    finished = []
    for raw in raw_products:
        for piece in split_fragments(raw).pieces:
            if sanitize_catch(piece):
                counters.sanitize_dropped += 1
                continue
            for atom in piece.GetAtoms():
                atom.SetAtomMapNum(0)
            forest_trace(parent, piece, info)
            ordered, smiles = cannonicalize_order(piece)
            finished.append((ordered, smiles))
    return finished


def _keep_fragment(finished, target: Mol):
    """The fragment closest to ``target``. The rest were cleaved off."""

    best = None
    best_cost = None
    for mol, smiles in finished:
        if smiles == canon_smiles(target):
            cost = -1
        else:
            cost = atom_diff(mol, target).cost()
        if best_cost is None or cost < best_cost:
            best = (mol, smiles)
            best_cost = cost
    if best is None:
        return None, []
    discarded = [item for item in finished if item[0] is not best[0]]
    return best, discarded


@dataclass
class _Walk:
    mol: Mol
    steps: tuple[_PlanStep, ...]
    sides: tuple[CleavageSide, ...]
    opens: tuple[frozenset[int], ...]


def _cleavage_site(value: object) -> frozenset[int]:
    """Known-index site as a frozenset for cleavage bookkeeping."""

    site = _as_site(value)
    if isinstance(site, int):
        return frozenset({site})
    if isinstance(site, tuple):
        return frozenset(site)
    sample = next(iter(site), None)
    if isinstance(sample, frozenset):
        return frozenset(_flat_ints(site))
    return frozenset(item for item in site if isinstance(item, int))


def find_path(
    reactant: Mol | str,
    target: Mol | str,
    ruleset=None,
    counters=None,
    *,
    use_filters=True,
    max_paths=1,
    max_nodes=800,
):
    """Yield phase-I plans that turn ``reactant`` into ``target``.

    ``ruleset`` is one :class:`~xenosite.refactor_poc.rules.RuleSet`. Filters
    come from :func:`atom_diff` and are applied to each child pattern; the
    set does not hide them. Pass ``use_filters=False`` to bill the same pair
    with every site edited. ``counters`` is optional; pass one in when the
    test needs ``billed``.
    """

    reactant = as_mol(reactant)
    target_mol = as_mol(target)
    reactant = copy_mol(reactant)
    reactant._forest = None
    reactant = install_forest(reactant)
    target_smiles = canon_smiles(target_mol)
    if ruleset is None:
        ruleset = default_ruleset()
    if counters is None:
        counters = PathCounters()

    queue = deque([_Walk(reactant, (), (), ())])
    seen = {canon_smiles(reactant)}
    found = 0

    while queue and found < max_paths and counters.nodes < max_nodes:
        walk = queue.popleft()
        counters.nodes += 1
        here = canon_smiles(walk.mol)
        if here == target_smiles:
            yield PathOutcome(
                plan=_deps(walk.steps),
                maybe=Maybe(walk.sides),
                smiles=here,
            )
            found += 1
            continue

        diff = atom_diff(walk.mol, target_mol)
        parent_cost = diff.cost()
        # See HEURISTICS.md before changing what these filters are allowed to see.
        filter_rules, filter_sites = _filters(diff, use_filters, walk.mol)
        # Cleavage children run first. The set still passes each pattern to
        # the filters; this only picks an order.
        order_key = None
        if diff.target_smaller or diff.has_cleavage:

            def _cleavage_first(rule):
                return (0 if _rule_can_cleave(rule) else 1, rule.name)

            order_key = _cleavage_first

        for por in ruleset.metabolites(
            walk.mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            counters=counters,
            order_key=order_key,
        ):
            finished = _finish(walk.mol, por.products, por.info, counters)
            if not finished:
                continue
            kept, discarded = _keep_fragment(finished, target_mol)
            if kept is None:
                continue
            child, child_smiles = kept
            closer = (
                child_smiles == target_smiles
                or atom_diff(child, target_mol).cost() < parent_cost
            )
            if not closer:
                continue
            if child_smiles in seen and child_smiles != target_smiles:
                continue
            seen.add(child_smiles)

            options = por.info["options"]
            cleaves = bool(options.get("cleaves"))
            if len(finished) == 1 and cleaves:
                opens = walk.opens + (_cleavage_site(por.info["site"]),)
                sides = walk.sides
            elif len(finished) > 1:
                opens = walk.opens
                sides = walk.sides + tuple(
                    CleavageSide(
                        site=_cleavage_site(por.info["site"]),
                        side=smiles,
                        opens=walk.opens,
                    )
                    for _mol, smiles in discarded
                )
            else:
                opens = walk.opens
                sides = walk.sides

            steps = walk.steps + _steps_for(walk.mol, por.info)
            queue.append(_Walk(child, steps, sides, opens))


@dataclass
class _Expand:
    mol: Mol
    depth: int
    info: dict[str, object] | None = None


def _keep_rule(rule, info) -> bool:
    return True


def _keep_site(site, info) -> bool:
    return True


def _enumerate(
    reactant: Mol | str,
    ruleset,
    *,
    filter_rules,
    filter_sites,
    depth: int,
    pop,
    lifo: bool,
):
    """Metabolites of one ruleset, up to ``depth``. No atom diff, no closer drop.

    ``pop`` is the frontier. A stack (``lifo``) expands the newest child
    before the rest of that generation, so a depth-2 sample does not have
    to finish the depth-1 frontier. A queue does not. Filters are the same
    callbacks the set forwards to each child, and they run on the mol
    being expanded.
    """

    if ruleset is None:
        ruleset = default_ruleset()
    start = as_mol(reactant)
    frontier: deque[_Expand] = deque([_Expand(start, 0)])
    seen = {canon_smiles(start)}
    while frontier:
        node = pop(frontier)
        if node.info is not None:
            yield node.mol, node.info
        if node.depth >= depth:
            continue
        children: list[_Expand] = []
        for product, info in ruleset.metabolize(
            node.mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
        ):
            smiles = info["csmi"]
            if not isinstance(smiles, str) or smiles in seen:
                continue
            seen.add(smiles)
            children.append(_Expand(product, node.depth + 1, info))
        if lifo:
            children.reverse()
        frontier.extend(children)


def bfs(
    reactant: Mol | str,
    ruleset=None,
    *,
    filter_rules=_keep_rule,
    filter_sites=_keep_site,
    depth: int = 1,
):
    """Breadth-first metabolites. A queue. See :func:`_enumerate`."""

    yield from _enumerate(
        reactant,
        ruleset,
        filter_rules=filter_rules,
        filter_sites=filter_sites,
        depth=depth,
        pop=deque.popleft,
        lifo=False,
    )


def dfs(
    reactant: Mol | str,
    ruleset=None,
    *,
    filter_rules=_keep_rule,
    filter_sites=_keep_site,
    depth: int = 1,
):
    """Depth-first metabolites. A stack. See :func:`_enumerate`."""

    yield from _enumerate(
        reactant,
        ruleset,
        filter_rules=filter_rules,
        filter_sites=filter_sites,
        depth=depth,
        pop=deque.pop,
        lifo=True,
    )
