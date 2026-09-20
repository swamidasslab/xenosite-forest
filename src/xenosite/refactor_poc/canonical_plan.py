"""Canonical elementary plans for composite reaction hops.

A rule's :meth:`~xenosite.refactor_poc.rules.ReactionRule.canonical_plan`
returns the elementary steps a search should record. Ordinary rules are
already elementary: the plan is that rule at the site. A composite hop
(quinone formation) expands to hydroxylation then dehydrogenation.
Epoxidation and N-dealkylation look ahead to their phase-I group at the
same site (stable / unstable oxygenation).

``Addition.phase1`` stays reserved; its schema is not decided. Search
reads this method instead of inventing that field.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple, Protocol

from xenosite.forest.step_plan import AtomRef as AddedRef
from xenosite.forest.step_plan import Deps, Step, StepPlan
from xenosite.refactor_poc.rdkitutil import Mol, is_tracing
from xenosite.refactor_poc.records import AtomRef, Effect, Site, SiteInfo, _flat_ints

PlanAtom = int | AddedRef | AtomRef


class _NamedGroup(Protocol):
    name: str | None
    longname: str | None


class CanonicalStep(NamedTuple):
    """One elementary step in a canonical plan.

    ``site`` is a flat tuple of atom notes: a known index, a forest
    :class:`AddedRef` for an atom a prior transform added, or a records
    :class:`AtomRef` for an atom that does not exist yet.
    """

    rule: str
    site: tuple[PlanAtom, ...]


def plan_atom_note(mol: Mol, idx: int) -> int | AddedRef:
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


def _site_atoms(site: Site) -> tuple[int, ...]:
    if isinstance(site, int):
        return (site,)
    if isinstance(site, tuple):
        return site
    return tuple(sorted(_flat_ints(site)))


def identity_canonical_plan(
    rule_name: str, mol: Mol, site: Site
) -> tuple[CanonicalStep, ...]:
    """The rule is already elementary: one step at ``site``."""

    return (
        CanonicalStep(
            rule_name,
            tuple(plan_atom_note(mol, idx) for idx in _site_atoms(site)),
        ),
    )


def group_look_ahead_plan(
    group: _NamedGroup, mol: Mol, site: Site
) -> tuple[CanonicalStep, ...]:
    """One step naming a phase-I group ruleset at ``site``.

    Epoxidation looks ahead to StableOxygenation; N-dealkylation to
    UnstableOxygenation. The plan is not a multi-step invention: same site,
    group longname (or short name).
    """

    label = group.longname or group.name
    if not label:
        return ()
    return identity_canonical_plan(label, mol, site)


def _bonded(mol: Mol, idx: int, atomic_num: int) -> int | None:
    """Neighbor of ``idx`` with this atomic number, if the mol already has one."""

    atom = mol.GetAtomWithIdx(idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == atomic_num:
            return neighbor.GetIdx()
    return None


_HALOGEN = frozenset({"F", "Cl", "Br", "I", "At"})
_HALOGEN_Z = {"F": 9, "Cl": 17, "Br": 35, "I": 53, "At": 85}


def _end_needs_oxygen(end: Effect) -> bool:
    partner = end.get("partner") or ""
    return "O" in (end.get("needs") or "") or (
        "O" in (end.get("adds") or "") and partner != "O"
    )


def _oxygen_ref(anchor: PlanAtom, fallback_idx: int) -> AtomRef:
    """Dehydrogenation note for the oxygen a prep step will add at ``anchor``."""

    if isinstance(anchor, int):
        return AtomRef(anchor, "O", 0)
    if isinstance(anchor, AtomRef):
        return AtomRef(anchor.idx, "O", anchor.depth)
    if isinstance(anchor, AddedRef):
        origin = anchor.origin
        if origin is not None:
            return AtomRef(origin, "O", int(anchor.depth))
    return AtomRef(fallback_idx, "O", 0)


def hydroxylation_then_dehydrogenation(
    mol: Mol,
    ends: Sequence[Effect],
    end_atoms: Sequence[int],
) -> tuple[CanonicalStep, ...]:
    """Preps that supply missing oxygens, then one dehydrogenation.

    An end that already carries oxygen (or N/C partner) keeps that atom. A
    bare carbon that needs oxygen becomes a Hydroxylation. A carbon-halogen
    end that gains oxygen becomes OxidativeDehalogenation (reads ``partner``).
    The dehydrogenation site points at the oxygen that prep would add.
    """

    preps: list[CanonicalStep] = []
    dh_refs: list[PlanAtom] = []
    for end, atom in zip(ends, end_atoms):
        partner = end.get("partner") or ""
        if _end_needs_oxygen(end):
            anchor = plan_atom_note(mol, atom)
            if partner in _HALOGEN:
                halo = _bonded(mol, atom, _HALOGEN_Z[partner])
                if halo is None:
                    continue
                preps.append(
                    CanonicalStep(
                        "OxidativeDehalogenation",
                        (anchor, plan_atom_note(mol, halo)),
                    )
                )
            else:
                preps.append(CanonicalStep("Hydroxylation", (anchor,)))
            dh_refs.append(_oxygen_ref(anchor, atom))
            continue
        atomic_num = {"O": 8, "N": 7, "C": 6, "S": 16}.get(partner)
        hetero = _bonded(mol, atom, atomic_num) if atomic_num else None
        if hetero is not None:
            dh_refs.append(plan_atom_note(mol, hetero))
    if not dh_refs:
        return ()
    return tuple(preps) + (CanonicalStep("Dehydrogenation", tuple(dh_refs)),)


def quinone_canonical_plan(mol: Mol, info: SiteInfo) -> tuple[CanonicalStep, ...]:
    """Canonical elementary plan for one quinone (or related) hop."""

    if "ends" not in info:
        return identity_canonical_plan("Dehydrogenation", mol, info["site"])
    plan = hydroxylation_then_dehydrogenation(mol, info["ends"], info["end_atoms"])
    if plan:
        return plan
    return identity_canonical_plan("Dehydrogenation", mol, info["site"])


def _anchor(item: PlanAtom) -> int | None:
    if isinstance(item, AtomRef):
        return item.idx
    if isinstance(item, int):
        return item
    origin = getattr(item, "origin", None)
    if origin is None:
        return None
    return int(origin)


def _anchors(step: CanonicalStep) -> set[int]:
    return {
        anchor
        for anchor in (_anchor(item) for item in step.site)
        if anchor is not None
    }


def _forest_item(item: PlanAtom, steps: Sequence[CanonicalStep], later: int) -> PlanAtom:
    """Hand :class:`Deps` a forest ref. The rule name is the earlier step's.

    ``added_by`` sites the whole earlier step (forest creation keys), not only
    the AtomRef origin atom.
    """

    if not isinstance(item, AtomRef):
        return item
    for previous in steps[:later]:
        anchors = _anchors(previous)
        if item.idx in anchors:
            return AddedRef(added_by=(previous.rule, frozenset(anchors)))
    return item.idx


def as_deps(steps: Sequence[CanonicalStep]) -> StepPlan:
    """Turn a canonical plan into a forest :class:`Deps` for apply / search.

    A later step depends on an earlier one when its site names an atom that
    step added. An :class:`AtomRef` is that note: its ``idx`` is the atom the
    earlier step changed, and its ``element`` is what had to be added.
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
            Step(
                step.rule,
                frozenset(_forest_item(item, steps, later) for item in step.site),
            )
        )
    # Deps is wrapped by @unstable; pyright does not treat the constructor as Deps/StepPlan.
    return Deps(tuple(forest_steps), edges)  # pyright: ignore[reportReturnType, reportCallIssue]
