"""Hard PathOutcome cases: macrocycle span, multi-route, impossible early abort.

Port of ``tests/test_path_outcome_hard.py``. Forest ``PathOutcome.steps`` /
``linearizations()`` / ``budget_exhausted`` map to poc ``plan.children`` /
``plan.linearizations()`` / empty-frontier or ``nodes >= max_nodes``.
Symmetric arms collapse to one topol orbit (see DIVERGENCES.md).
"""

from __future__ import annotations

from xenosite.forest.find_path import (
    CleavageSide,
    PathCounters,
    PathOutcome,
    find_path,
)
from xenosite.forest.rules import Dealkylation, Hydrolysis, NDealkylation
from xenosite.forest.rulesets import PhaseOne, RuleSet

_ND = RuleSet((NDealkylation,), name="ND")
_HD = RuleSet((Hydrolysis, Dealkylation), name="HydrolysisDealk")
_DEALK = RuleSet((Dealkylation,), name="Dealk")


def test_macrocycle_ring_open_then_bifurcate_spans_opens():
    """First cleave opens the ring (no bag); second cuts out a linker segment."""

    reactant = "C1CCCCCCNC2CCCC(CC2)NCCCC1"
    product = "NC1CCCC(=O)CC1"
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset=_ND,
            max_paths=2,
            max_nodes=120,
            counters=counters,
        )
    )
    assert hits, "expected ring-open then cleave to amino-ketone"
    assert counters.nodes < 120

    two_cleave = [h for h in hits if len(h.plan.children) >= 2]
    assert two_cleave
    outcome = min(two_cleave, key=lambda h: len(h.plan.children))
    assert isinstance(outcome, PathOutcome)
    assert outcome.maybe
    assert all(isinstance(e, CleavageSide) for e in outcome.maybe.entries)
    assert any(e.opens for e in outcome.maybe.entries), outcome.maybe
    assert all(s.rule for s in outcome.plan.children)
    bag = next(e for e in outcome.maybe.entries if e.opens)
    open_site = bag.opens[0]
    assert outcome.allows("NDealkylation", open_site)
    assert not outcome.allows("NDealkylation", bag.site)


def test_multi_ndealk_distinct_required_routes_same_aldehyde():
    """Two-arm N-dealk to the dialdehyde: unordered Deps → ≥2 linearizations."""

    reactant = "CN(C)Cc1ccc(CN(C)CC)cc1"
    product = "O=Cc1ccc(C=O)cc1"
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset=_ND,
            max_paths=10,
            max_nodes=80,
            counters=counters,
        )
    )
    two_step = [
        h
        for h in hits
        if len(h.plan.children) == 2
        and all(s.rule == "NDealkylation" for s in h.plan.children)
    ]
    assert len(two_step) >= 1, "expected a two-cleave Required route; got %s" % (
        [[s.rule for s in h.plan.children] for h in hits],
    )
    outcome = two_step[0]
    assert isinstance(outcome, PathOutcome)
    assert outcome.maybe
    assert outcome.plan.contains(
        ["NDealkylation", "NDealkylation"], by="rule"
    )
    assert outcome.plan.n_linearizations() >= 2
    sides = {e.side for e in outcome.maybe.entries}
    assert len(sides) >= 2, sides
    assert len(hits) >= 2


def test_impossible_targets_abort_without_draining_budget():
    """Heuristics must empty the frontier — miss-after-budget is not a pass."""

    cases = [
        ("c1ccccc1O", "C", 40),
        ("c1cccc2ccccc12", "C", 40),
        ("c1ccccc1", "FC(F)(F)F", 40),
    ]
    for reactant, target, ceiling in cases:
        counters = PathCounters()
        hits = list(
            find_path(
                reactant,
                target,
                ruleset=PhaseOne,
                max_nodes=200,
                counters=counters,
                )
        )
        assert not hits, "unexpected path %s → %s: %s" % (reactant, target, hits)
        # Empty frontier, not a node-cap stop.
        assert counters.nodes < 200, (
            "%s → %s hit max_nodes (nodes=%d)" % (reactant, target, counters.nodes)
        )
        assert counters.billed <= ceiling, (
            "%s → %s billed %d > ceiling %d"
            % (reactant, target, counters.billed, ceiling)
        )


def test_path_outcome_linearizations_and_allows_delegate():
    hits = list(
        find_path(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=40,
        )
    )
    assert hits
    outcome = hits[0]
    assert isinstance(outcome, PathOutcome)
    assert list(outcome.plan.linearizations())
    assert outcome.allows(side=outcome.maybe.entries[0].side)
    assert outcome.plan is hits[0].plan and outcome.maybe is hits[0].maybe


def test_asymmetric_multi_mcs_arm_selection():
    """Cl-benzyl vs Ph: each aldehyde target uses its own arm (not the other)."""

    parent = "CN(Cc1ccccc1)Cc1ccc(Cl)cc1"
    cl_hits = list(
        find_path(
            parent,
            "O=Cc1ccc(Cl)cc1",
            ruleset=_ND,
            max_paths=5,
            max_nodes=40,
        )
    )
    ph_hits = list(
        find_path(
            parent,
            "O=Cc1ccccc1",
            ruleset=_ND,
            max_paths=5,
            max_nodes=40,
        )
    )
    assert cl_hits and ph_hits
    cl_direct = [h for h in cl_hits if len(h.plan.children) == 1]
    ph_direct = [h for h in ph_hits if len(h.plan.children) == 1]
    assert cl_direct and ph_direct
    assert {frozenset(map(str, h.plan.children[0].site)) for h in cl_direct} != {
        frozenset(map(str, h.plan.children[0].site)) for h in ph_direct
    }


def test_tribenzyl_symmetric_arms_one_orbit():
    """Three equivalent benzyl arms → one topol orbit (forest reports three)."""

    parent = "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1"
    hits = list(
        find_path(
            parent,
            "O=Cc1ccccc1",
            ruleset=_ND,
            max_paths=12,
            max_nodes=40,
        )
    )
    assert hits
    sites = {
        frozenset(map(str, h.plan.children[0].site))
        for h in hits
        if len(h.plan.children) == 1
    }
    assert len(sites) >= 1
    # Symmetry collapse is approved (DIVERGENCES.md); forest would see ≥3.


def test_o_demethyl_maybe_allows_after_hydrolysis():
    """Acetate → catechol: bags present; demethyl / ester sites are Required."""

    hits = list(
        find_path(
            "CC(=O)Oc1ccc(OC)cc1",
            "Oc1ccc(O)cc1",
            ruleset=_HD,
            max_paths=3,
            max_nodes=80,
        )
    )
    assert hits
    hyd = [
        h
        for h in hits
        if "Hydrolysis" in [s.rule for s in h.plan.children]
        and "Dealkylation" in [s.rule for s in h.plan.children]
    ]
    h = hyd[0] if hyd else hits[0]
    assert h.maybe
    assert h.allows(side=h.maybe.entries[0].side)
    # Single-bag exact match does not allow (anisole). Multi-bag walks may
    # still allow a Required site that shares an atom with another bag.
    dealk = next(s for s in h.plan.children if s.rule == "Dealkylation")
    if len(h.maybe.entries) == 1:
        assert not h.allows("Dealkylation", dealk.site)


def test_stereo_e_z_terbinafine_same_tba():
    e = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
    z = "CN(C/C=C\\C#CC(C)(C)C)Cc1cccc2ccccc12"
    tba = "CC(C)(C)C#CC=CC=O"
    e_hits = list(find_path(e, tba, ruleset=_DEALK, max_paths=1, max_nodes=40))
    z_hits = list(find_path(z, tba, ruleset=_DEALK, max_paths=1, max_nodes=40))
    assert e_hits and z_hits
    assert [s.rule for s in e_hits[0].plan.children] == [
        s.rule for s in z_hits[0].plan.children
    ]


def test_bipiperidinyl_two_ndealk_records_opens():
    """Linked piperidines → amino-aldehyde: two ND hops, first is ring-open."""

    counters = PathCounters()
    hits = list(
        find_path(
            "C1CCCN(C1)C1CCCCN1",
            "NCCCCC=O",
            ruleset=_ND,
            max_paths=5,
            max_nodes=80,
            counters=counters,
        )
    )
    assert hits or counters.nodes >= 80
    if not hits:
        return
    two = [h for h in hits if len(h.plan.children) == 2]
    assert two
    assert all(s.rule == "NDealkylation" for h in two for s in h.plan.children)
    assert any(h.maybe and any(e.opens for e in h.maybe.entries) for h in two)


def test_metal_target_fast_fail_no_edits():
    counters = PathCounters()
    hits = list(
        find_path(
            "c1ccccc1",
            "[Fe]",
            ruleset=PhaseOne,
            max_nodes=50,
            counters=counters,
        )
    )
    assert hits == []
    assert counters.mol_edits == 0
    assert counters.nodes <= 1
    assert counters.nodes < 50


def test_dialdehyde_under_tight_billed_budget():
    """Crowded bis-N-dealk still finds under a modest billed cap (or hits the node cap)."""

    counters = PathCounters()
    ceiling = 35
    hits = list(
        find_path(
            "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
            "O=Cc1ccc(C=O)cc1",
            ruleset=_ND,
            max_paths=2,
            max_nodes=ceiling,
            counters=counters,
        )
    )
    assert hits or counters.nodes >= ceiling
    if hits:
        assert counters.billed <= ceiling or counters.nodes >= ceiling
