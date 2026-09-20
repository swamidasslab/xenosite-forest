"""Custom rules that force guided-search and MCS-cleavage branches.

Production Phase I rules emit phase1 plans, so the hop loop, prune exits, and
several cleavage-geometry decisions never run. These rules exist only to drive
those paths.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite._archive_forest import PathSearchCounters, RuleSet, find_path
from xenosite._archive_forest.base import ReactionRule, SmartsReactionRule
from xenosite._archive_forest.guided_path import _canon
from xenosite._archive_forest.path_context import (
    CLEAVE,
    PathContext,
    cleavage_site_bond,
    cleavage_site_deep_exterior,
    cleavage_site_deep_interior,
)
from xenosite._archive_forest.step_plan import Step, StepPlan


def _mol(smi: str):
    return Chem.MolFromSmiles(smi)


class ScriptedRule(ReactionRule):
    """Yields whatever ``emit(mol, ctx)`` returns. Flags force prune exits."""

    def __init__(
        self,
        name,
        emit,
        *,
        could_help=True,
        terminal=False,
        reach=True,
        redundant=False,
    ):
        super().__init__(name=name)
        self._emit = emit
        self._could = could_help
        self._terminal = terminal
        self._reach = reach
        self._redundant = redundant

    def could_help(self, mol, target, ctx) -> bool:
        if callable(self._could):
            return bool(self._could(mol, target, ctx))
        return bool(self._could)

    def is_terminal_product(self, mol) -> bool:
        if callable(self._terminal):
            return bool(self._terminal(mol))
        return bool(self._terminal)

    def child_may_reach(self, parent, child, target, ctx) -> bool:
        if callable(self._reach):
            return bool(self._reach(parent, child, target, ctx))
        return bool(self._reach)

    def is_redundant(self, peer_rules) -> bool:
        return bool(self._redundant)

    def enumerate_for_path(self, mol, ctx, **kwargs):
        yield from self._emit(mol, ctx)


class _BondCleave(SmartsReactionRule):
    """C–C cut. ``phase1_steps`` skips RunReactants; search still calls sites_toward."""

    phase1_equivalent = True
    smarts = (("[#6:1]-[#6:2]>>([#6:1].[#6:2])", {"formula_hint": CLEAVE}),)
    mapid_site = [1, 2]

    def phase1_steps(self, mol, site, **kwargs):
        return StepPlan.singleton(self.name, frozenset(site))


class _HeteroCleave(SmartsReactionRule):
    """C–N / C–O cut, so ring-open N is kept and ring-open O can be pruned."""

    phase1_equivalent = True
    smarts = (
        ("[#6:1]-[#7,#8:2]>>([#6:1].[#7,#8:2])", {"formula_hint": CLEAVE}),
    )
    mapid_site = [1, 2]

    def phase1_steps(self, mol, site, **kwargs):
        return StepPlan.singleton(self.name, frozenset(site))


class _SkipBond(SmartsReactionRule):
    """Site is the 1,3-pair, which is not a bond. cleavage_site_bond returns None."""

    phase1_equivalent = True
    smarts = (
        ("[#6:1]-[#6:2]-[#6:3]>>([#6:1]-[#6:2].[#6:3])", {"formula_hint": CLEAVE}),
    )
    mapid_site = [1, 3]

    def phase1_steps(self, mol, site, **kwargs):
        return StepPlan.singleton(self.name, frozenset(site))


def _run(rule, reactant, product, **kwargs):
    counters = PathSearchCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset=RuleSet([rule] if not isinstance(rule, list) else rule, name="Script"),
            counters=counters,
            max_expansions=kwargs.pop("max_expansions", 20),
            **kwargs,
        )
    )
    return hits, counters


def _two_hop(mol, ctx):
    smi = _canon(mol)
    if smi == "CCO":
        yield ("hop", frozenset({0}), [_mol("CC=O")])
    elif smi == "CC=O":
        yield ("hop", frozenset({0}), [_mol("CC(=O)O")])


def test_scripted_hop_reaches_on_bfs_and_dfs():
    rule = ScriptedRule("HopToAcid", _two_hop)
    for search in ("bfs", "dfs"):
        hits, counters = _run(rule, "CCO", "CC(=O)O", depth=3, search=search)
        assert hits, search
        assert [s[0] for s in hits[0].steps] == ["HopToAcid", "HopToAcid"]
        assert counters.nodes_enqueued >= 2
        assert counters.site_applies >= 2


def test_scripted_hop_prune_exits():
    def emit(mol, ctx):
        ald = _mol("CC=O")
        yield ("hop", frozenset({0}), [ald, ald])
        yield ("hop", frozenset({1}), [_mol("*C")])

    rule = ScriptedRule(
        "PruneHop",
        emit,
        terminal=lambda mol: _canon(mol) == "*C",
    )
    # Duplicate aldehyde is pruned after the first enqueue; star stops early.
    hits, counters = _run(rule, "CCO", "CC(=O)O", depth=2)
    assert hits == []
    assert counters.nodes_enqueued >= 2
    assert counters.nodes_pruned >= 1

    # One hop left, but no depth remains to enqueue it.
    _, shallow = _run(rule, "CCO", "CC(=O)O", depth=1)
    assert shallow.nodes_enqueued == 1
    assert shallow.site_applies >= 1

    _, capped = _run(rule, "CCO", "CC(=O)O", depth=0)
    assert capped.rule_expansions == 0

    def one(mol, ctx):
        yield ("hop", frozenset({0}), [_mol("CC=O")])

    _, blocked = _run(
        ScriptedRule("BlockHop", one, reach=False),
        "CCO",
        "CC(=O)O",
        depth=2,
    )
    assert blocked.nodes_pruned >= 1
    assert blocked.mol_edits == 0


def test_scripted_budget_could_help_and_redundant():
    def two(mol, ctx):
        yield ("hop", frozenset({0}), [_mol("CC=O")])
        yield ("hop", frozenset({1}), [_mol("C=O")])

    hits, counters = _run(
        ScriptedRule("BudgetHop", two),
        "CCO",
        "CC(=O)O",
        depth=2,
        max_expansions=1,
    )
    assert hits == []
    assert counters.budget_exhausted
    assert counters.site_applies == 1

    def boom(mol, ctx):
        raise RuntimeError("enumerate failed")

    with pytest.raises(RuntimeError, match="enumerate failed"):
        _run(ScriptedRule("BoomHop", boom), "CCO", "CC=O", depth=2)

    called = []

    def note(mol, ctx):
        called.append(_canon(mol))
        return iter(())

    dead = ScriptedRule("DeadHop", note, could_help=False)
    extra = ScriptedRule("AlsoDead", note, redundant=True)
    _, counters = _run([dead, extra], "CCO", "CC=O", depth=2)
    assert called == []
    assert counters.nodes_pruned >= 1
    assert counters.rule_expansions == 0


def test_broken_pipeline_raises():
    """Non-molecules, unknown kinds, and unknown search modes are bugs."""

    def junk(mol, ctx):
        yield ("hop", frozenset({0}), [object()])

    with pytest.raises(TypeError, match="molecule"):
        _run(ScriptedRule("JunkHop", junk), "CCO", "CC=O", depth=1)

    def weird(mol, ctx):
        yield ("nope", None)

    with pytest.raises(ValueError, match="plan"):
        _run(ScriptedRule("WeirdKind", weird), "CCO", "CC=O", depth=1)

    with pytest.raises(ValueError, match="search"):
        _run(ScriptedRule("HopToAcid", _two_hop), "CCO", "CC=O", search="best")


def test_scripted_plan_depth_failure_and_enqueue():
    def too_deep(mol, ctx):
        plan = StepPlan(
            (Step("Hydroxylation", frozenset({0})), Step("Hydroxylation", frozenset({1})))
        )
        yield ("plan", plan, frozenset({0}))

    _, counters = _run(ScriptedRule("DeepPlan", too_deep), "CC", "CCC", depth=1)
    assert counters.linearizations_applied >= 1
    assert counters.mol_edits == 0

    def bad(mol, ctx):
        yield ("plan", StepPlan((Step("NotARule", frozenset({0})),)), frozenset({0}))

    hits, counters = _run(ScriptedRule("BadPlan", bad), "CC", "CCC", depth=2)
    assert hits == []
    assert counters.linearizations_applied >= 1

    def hydroxylate(mol, ctx):
        if _canon(mol) == "CC":
            yield (
                "plan",
                StepPlan((Step("Hydroxylation", frozenset({0})),)),
                frozenset({0}),
            )

    hits, counters = _run(
        ScriptedRule("OhPlan", hydroxylate, reach=True, terminal=False),
        "CC",
        "CCC",
        depth=3,
    )
    assert hits == []
    assert counters.mol_edits >= 1
    assert counters.nodes_enqueued >= 2


def test_cleavage_rules_drive_mcs_geometry():
    bond = _BondCleave(name="BondCleave")

    # Long amine: far C–C cuts are deep exterior and dropped; the MCS end is kept.
    amine = _mol("NCCCCCCCCCCCCCCCCCCCC")
    methylamine = _mol("CN")
    ctx = PathContext.from_mols(amine, methylamine)
    kept = set(map(frozenset, bond.sites_toward(amine, methylamine, ctx)))
    exterior = [
        frozenset({b.GetBeginAtomIdx(), b.GetEndAtomIdx()})
        for b in amine.GetBonds()
        if cleavage_site_deep_exterior(amine, frozenset({b.GetBeginAtomIdx(), b.GetEndAtomIdx()}), ctx)
    ]
    assert exterior
    assert kept
    assert not (set(exterior) & kept)

    # tert-Butylbenzene: ring-interior C–C bonds are dropped; the alkyl cuts stay.
    tbu = _mol("CC(C)(C)c1ccccc1")
    benzene = _mol("c1ccccc1")
    ctx = PathContext.from_mols(tbu, benzene)
    kept = set(map(frozenset, bond.sites_toward(tbu, benzene, ctx)))
    interior = [
        frozenset({b.GetBeginAtomIdx(), b.GetEndAtomIdx()})
        for b in tbu.GetBonds()
        if cleavage_site_deep_interior(tbu, frozenset({b.GetBeginAtomIdx(), b.GetEndAtomIdx()}), ctx)
    ]
    assert interior
    assert kept
    assert not (set(interior) & kept)

    # 1,3-site is not a bond, so the bridge test is skipped and the site is kept.
    hexane = _mol("CCCCCC")
    sites = _SkipBond(name="SkipBond").sites_toward(hexane, _mol("CC"), PathContext.from_mols(hexane, _mol("CC")))
    assert sites
    assert all(cleavage_site_bond(hexane, s) is None for s in sites)

    # Piperidine N–C ring-open is never pruned.
    hetero = _HeteroCleave(name="HeteroCleave")
    pip = _mol("C1CCNCC1")
    pip_sites = hetero.sites_toward(pip, _mol("CCN"), PathContext.from_mols(pip, _mol("CCN")))
    assert pip_sites
    assert any(
        any(pip.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in s) for s in pip_sites
    )

    # Search actually consults the rule (plans), even though apply cannot resolve it.
    hits, counters = _run(bond, amine, methylamine, depth=2, max_expansions=4)
    assert hits == []
    assert counters.plans_expanded >= 1
    assert counters.linearizations_applied >= 1
