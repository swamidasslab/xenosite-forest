"""PathOutcome.maybe / allows: cleavage bags, not searched prefix chemistry.

Port of ``tests/test_maybe_path.py``. Forest ``MaybeFilter`` / ``CleavageSide``
map to forest ``Maybe`` / ``CleavageSide``. Forest helper ``_cleavage_sides_for_kept``
is not ported — bags are asserted on ``find_path`` outcomes.
"""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.find_path import (
    CleavageSide,
    Maybe,
    PathCounters,
    PathOutcome,
    canon_smiles,
    find_path,
)
from xenosite.forest.rules import Dealkylation, EpoxideOpening, Hydrolysis
from xenosite.forest.rulesets import PhaseOne, RuleSet

_HYD = RuleSet((Hydrolysis,), name="Hyd")
_DEALK = RuleSet((Dealkylation,), name="Dealk")


def test_cleavage_side_and_maybe_allows():
    side = CleavageSide(site=frozenset({1, 2}), side=canon_smiles("CCO"))
    empty = Maybe()
    assert not empty
    assert empty.entries == ()

    filled = Maybe(entries=(side,))
    assert filled
    assert filled.entries[0].side == canon_smiles("CCO")
    # Same nitrogen, different site (demethyl) passes; formation site does not.
    assert filled.allows("NDealkylation", frozenset({0, 1}))
    assert not filled.allows("NDealkylation", frozenset({1, 2}))
    assert filled.allows(side="CCO")
    assert not filled.allows("Hydroxylation", frozenset({9, 10}))


def test_hydrolysis_emits_cleavage_side_bag():
    """Ester → acid: discarded alcohol is a CleavageSide bag.

    Atom-diff filters refuse this Hydrolysis hop (all sites skipped); bags are
    still required when the chemistry runs (``use_filters=False``).
    """

    hits = list(
        find_path(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            ruleset=_HYD,
            max_paths=3,
            max_nodes=40,
            use_filters=False,
        )
    )
    assert hits
    with_bag = [h for h in hits if h.maybe]
    assert with_bag, "expected CleavageSide on bifurcating hydrolysis"
    outcome = with_bag[0]
    assert isinstance(outcome, PathOutcome)
    assert [s.rule for s in outcome.plan.children] == ["Hydrolysis"]
    assert all(isinstance(e, CleavageSide) for e in outcome.maybe.entries)
    assert any(e.side == canon_smiles("CC(C)(C)O") for e in outcome.maybe.entries)
    assert outcome.plan.contains(["Hydrolysis"], by="rule")


def test_dealkylation_emits_cleavage_side_bag():
    hits = list(
        find_path(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            ruleset=_DEALK,
            max_paths=4,
            max_nodes=40,
        )
    )
    assert hits
    with_bag = [h for h in hits if h.maybe]
    assert with_bag
    maybe = with_bag[0].maybe
    assert all(isinstance(e, CleavageSide) for e in maybe.entries)
    assert len(maybe.entries) >= 1
    parent_heavy = Chem.MolFromSmiles("CN(C)Cc1ccccc1").GetNumHeavyAtoms()
    assert all(
        Chem.MolFromSmiles(e.side).GetNumHeavyAtoms() < parent_heavy
        for e in maybe.entries
    )


def test_epoxide_opening_no_cleavage_side_bag():
    """Ring-open alternatives are parent-sized — no CleavageSide bag."""

    hits = list(
        find_path(
            "c1ccccc1C1OC1",
            "OCCc1ccccc1",
            ruleset=PhaseOne,
            max_paths=5,
            max_nodes=40,
        )
    )
    assert hits, "expected EpoxideOpening path"
    for outcome in hits:
        if outcome.plan.children and outcome.plan.children[0].rule == "EpoxideOpening":
            assert not outcome.maybe
            assert outcome.maybe.entries == ()


def test_multiple_cleavages_accumulate_bags_each_step():
    """Two sequential bifurcations → ≥1 CleavageSide per cleavage hop."""

    reactant = "CC(=O)Oc1ccc(OC(C)=O)cc1"
    product = "Oc1ccc(O)cc1"
    counters = PathCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset=_HYD,
            max_paths=8,
            max_nodes=100,
            counters=counters,
        )
    )
    assert hits, "expected path via two hydrolyses to hydroquinone"
    two_step = [
        h
        for h in hits
        if sum(1 for s in h.plan.children if s.rule == "Hydrolysis") >= 2
    ]
    assert two_step, "expected ≥2 Hydrolysis steps on a hit"
    outcome = min(two_step, key=lambda h: len(h.plan.children))
    n_hydro = sum(1 for s in outcome.plan.children if s.rule == "Hydrolysis")
    assert len(outcome.maybe.entries) >= n_hydro, (
        "each bifurcating cleavage must add a CleavageSide; "
        "hydrolysis=%d bags=%s" % (n_hydro, outcome.maybe)
    )
    assert all(isinstance(e, CleavageSide) for e in outcome.maybe.entries)
    assert all(e.side for e in outcome.maybe.entries)


def test_required_identity_empty_maybe():
    hits = list(
        find_path(
            "c1ccccc1",
            "c1ccccc1",
            max_paths=1,
            max_nodes=2,
        )
    )
    assert hits
    assert not hits[0].maybe


def test_tba_cleavage_side_and_preceding_ndealk_passes_maybe():
    """Terbinafine → TBA: formation Dealkylation; demethyl on same N passes Maybe."""

    terb = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
    tba = "CC(C)(C)C#CC=CC=O"
    counters = PathCounters()
    hits = list(
        find_path(
            terb,
            tba,
            ruleset=PhaseOne,
            max_paths=3,
            max_nodes=80,
            counters=counters,
        )
    )
    assert hits
    assert counters.nodes < 80

    direct = [
        h for h in hits if [s.rule for s in h.plan.children] == ["Dealkylation"]
    ]
    assert direct, "expected Dealkylation → TBA"
    outcome = direct[0]
    assert outcome.plan.contains(["Dealkylation"], by="rule")
    assert outcome.maybe
    assert all(isinstance(e, CleavageSide) for e in outcome.maybe.entries)
    assert any("cccc" in e.side for e in outcome.maybe.entries)

    demethyl_site = frozenset({0, 1})
    formation_site = frozenset({1, 2})
    assert outcome.allows("NDealkylation", demethyl_site), (
        "preceding demethyl on same N must pass Maybe; got %s" % (outcome.maybe,)
    )
    assert not outcome.allows("NDealkylation", formation_site)
    assert outcome.allows(side=outcome.maybe.entries[0].side)


def test_hydrolysis_metabolize_bifurcates_epoxide_does_not():
    """Live Hydrolysis metabolize is multi-piece; EpoxideOpening is not."""

    mol = Chem.MolFromSmiles("c1ccccc1C(=O)OC(C)(C)C")
    hyd = list(Hydrolysis().metabolize(mol))
    assert hyd
    # Bifurcation: acid and alcohol are separate yields (unique_csmi).
    csmi = {p.xf.csmi for p, _ in hyd}
    assert canon_smiles("O=C(O)c1ccccc1") in csmi
    assert canon_smiles("CC(C)(C)O") in csmi

    mol_e = Chem.MolFromSmiles("c1ccccc1C1OC1")
    eo = list(EpoxideOpening().metabolize(mol_e))
    assert eo
    # Single parent-sized product per site — not a discarded sibling bag.
    for product, _info in eo:
        assert product.GetNumHeavyAtoms() >= mol_e.GetNumHeavyAtoms() - 1
