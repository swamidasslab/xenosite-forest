"""Hard PathOutcome cases: macrocycle span, multi-route, impossible early abort."""

from __future__ import annotations

from xenosite.forest import PathOutcome, PathSearchCounters, find_path
from xenosite.forest.guided_path import CleavageSide, _site_key_set, _steps_key


def test_macrocycle_ring_open_then_bifurcate_spans_opens():
    """First cleave opens the ring (no bag); second cuts out a linker segment.

    Product is Forest-reachable amino-ketone (N-dealk leaves carbonyl), not diol.
    """
    reactant = "C1CCCCCCNC2CCCC(CC2)NCCCC1"
    product = "NC1CCCC(=O)CC1"
    counters = PathSearchCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset="PhaseOneRS",
            depth=3,
            max_paths=2,
            max_expansions=120,
            counters=counters,
        )
    )
    assert hits, "expected ring-open then cleave to amino-ketone"
    assert not counters.budget_exhausted

    two_cleave = [h for h in hits if len(h.steps) >= 2]
    assert two_cleave
    outcome = min(two_cleave, key=lambda h: len(h.steps))
    assert isinstance(outcome, PathOutcome)
    assert outcome.maybe
    assert all(isinstance(e, CleavageSide) for e in outcome.maybe.entries)
    # Bifurcating bag spans prior ring-open site(s).
    assert any(e.opens for e in outcome.maybe.entries), outcome.maybe
    # Formation sites are Required only in steps/plan.
    assert all(s[0] for s in outcome.steps)
    # Shared N with an open or bifurcate site (other than exact bif site) may pass.
    bag = next(e for e in outcome.maybe.entries if e.opens)
    open_site = bag.opens[0]
    # Site on the ring-open span (not the bifurcate formation site) passes.
    assert outcome.allows("NDealkylation", open_site)
    assert not outcome.allows("NDealkylation", bag.site)


def test_multi_ndealk_distinct_required_routes_same_aldehyde():
    """Two topologically distinct Required step paths to the same dialdehyde."""
    reactant = "CN(C)Cc1ccc(CN(C)CC)cc1"
    product = "O=Cc1ccc(C=O)cc1"
    counters = PathSearchCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset="ND",
            depth=3,
            max_paths=10,
            max_expansions=80,
            counters=counters,
        )
    )
    two_step = [
        h
        for h in hits
        if len(h.steps) == 2 and all(s[0] == "NDealkylation" for s in h.steps)
    ]
    assert len(two_step) >= 2, "expected ≥2 Required two-cleave routes; got %s" % (
        [h.plan for h in hits],
    )
    keys = {_steps_key(h.steps) for h in two_step}
    assert len(keys) >= 2, "Required step sites must differ across routes"
    maybes = {str(h.maybe) for h in two_step}
    assert len(maybes) >= 2, "each route should carry its own Maybe"
    for h in two_step:
        assert isinstance(h, PathOutcome)
        assert [s[0] for s in h.steps] == ["NDealkylation", "NDealkylation"]
        assert h.maybe


def test_impossible_targets_abort_without_budget_exhaustion():
    """Heuristics must empty the frontier — miss-after-budget is not a pass."""
    cases = [
        ("c1ccccc1O", "C", 40),  # phenol → methane
        ("c1cccc2ccccc12", "C", 40),  # naphthalene → methane
        ("c1ccccc1", "FC(F)(F)F", 40),  # benzene → CF4
    ]
    for r, t, ceiling in cases:
        counters = PathSearchCounters()
        hits = list(
            find_path(
                r,
                t,
                ruleset="PhaseOneRS",
                depth=None,
                max_expansions=200,
                counters=counters,
            )
        )
        assert not hits, "unexpected path %s → %s: %s" % (r, t, hits)
        assert not counters.budget_exhausted, (
            "%s → %s exhausted budget (exp=%d) — heuristics failed to abort"
            % (r, t, counters.rule_expansions)
        )
        assert counters.rule_expansions <= ceiling, (
            "%s → %s expansions %d > ceiling %d"
            % (r, t, counters.rule_expansions, ceiling)
        )


def test_path_outcome_linearizations_and_allows_delegate():
    hits = list(
        find_path(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            ruleset="PhaseOneRS",
            depth=2,
            max_paths=1,
            max_expansions=40,
        )
    )
    assert hits
    o = hits[0]
    assert isinstance(o, PathOutcome)
    assert list(o.linearizations())
    assert o.allows(side=o.maybe.entries[0].side)
    smiles, steps, mols, plan, maybe = o
    assert plan is o.plan and maybe is o.maybe


def test_asymmetric_multi_mcs_arm_selection():
    """Cl-benzyl vs Ph: each aldehyde target uses its own arm (not the other)."""
    parent = "CN(Cc1ccccc1)Cc1ccc(Cl)cc1"
    cl_hits = list(
        find_path(
            parent,
            "O=Cc1ccc(Cl)cc1",
            ruleset="ND",
            depth=2,
            max_paths=5,
            max_expansions=40,
        )
    )
    ph_hits = list(
        find_path(
            parent,
            "O=Cc1ccccc1",
            ruleset="ND",
            depth=2,
            max_paths=5,
            max_expansions=40,
        )
    )
    assert cl_hits and ph_hits
    cl_direct = [h for h in cl_hits if len(h.steps) == 1]
    ph_direct = [h for h in ph_hits if len(h.steps) == 1]
    assert cl_direct and ph_direct
    assert {frozenset(map(str, h.steps[0][1])) for h in cl_direct} != {
        frozenset(map(str, h.steps[0][1])) for h in ph_direct
    }


def test_tribenzyl_three_distinct_required_sites():
    parent = "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1"
    hits = list(
        find_path(
            parent,
            "O=Cc1ccccc1",
            ruleset="ND",
            depth=2,
            max_paths=12,
            max_expansions=40,
        )
    )
    assert hits
    sites = {frozenset(map(str, h.steps[0][1])) for h in hits if len(h.steps) == 1}
    assert len(sites) >= 3


def test_o_demethyl_maybe_allows_after_hydrolysis():
    """Acetate → catechol: bags present; demethyl site is Required Dealkylation."""
    hits = list(
        find_path(
            "CC(=O)Oc1ccc(OC)cc1",
            "Oc1ccc(O)cc1",
            ruleset="PhaseOneRS",
            depth=3,
            max_paths=3,
            max_expansions=80,
        )
    )
    assert hits
    h = hits[0]
    assert "Hydrolysis" in [s[0] for s in h.steps]
    assert "Dealkylation" in [s[0] for s in h.steps]
    assert h.maybe
    assert h.allows(side=h.maybe.entries[0].side)
    # Formation Dealkylation site is Required — not a Maybe bag site.
    dealk_site = next(s[1] for s in h.steps if s[0] == "Dealkylation")
    assert not h.allows("Dealkylation", dealk_site)


def test_stereo_e_z_terbinafine_same_tba():
    e = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
    z = "CN(C/C=C\\C#CC(C)(C)C)Cc1cccc2ccccc12"
    tba = "CC(C)(C)C#CC=CC=O"
    e_hits = list(
        find_path(e, tba, depth=2, max_paths=1, max_expansions=40)
    )
    z_hits = list(
        find_path(z, tba, depth=2, max_paths=1, max_expansions=40)
    )
    assert e_hits and z_hits
    assert [s[0] for s in e_hits[0].steps] == [s[0] for s in z_hits[0].steps]


def test_bipiperidinyl_two_ndealk_records_opens():
    """Linked piperidines → amino-aldehyde: two ND hops, first is ring-open."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "C1CCCN(C1)C1CCCCN1",
            "NCCCCC=O",
            ruleset="ND",
            depth=4,
            max_paths=5,
            max_expansions=80,
            counters=counters,
        )
    )
    assert hits or counters.budget_exhausted
    if not hits:
        return
    two = [h for h in hits if len(h.steps) == 2]
    assert two
    assert all(s[0] == "NDealkylation" for h in two for s in h.steps)
    assert any(h.maybe and any(e.opens for e in h.maybe.entries) for h in two)


def test_metal_target_fast_fail_no_budget():
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "c1ccccc1",
            "[Fe]",
            ruleset="PhaseOneRS",
            depth=3,
            max_expansions=50,
            counters=counters,
        )
    )
    assert hits == []
    assert counters.rule_expansions == 0
    assert counters.billed() == 0
    assert not counters.budget_exhausted


def test_dialdehyde_under_tight_billed_budget():
    """Crowded bis-N-dealk still finds under a modest billed cap (or flags budget)."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
            "O=Cc1ccc(C=O)cc1",
            ruleset="ND",
            depth=3,
            max_paths=2,
            max_expansions=35,
            counters=counters,
        )
    )
    assert hits or counters.budget_exhausted
    if hits:
        assert not counters.budget_exhausted or counters.billed() <= 35
