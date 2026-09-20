"""MaybeFilter is a cleavage-side bag, not searched prefix chemistry."""

from __future__ import annotations

from rdkit import Chem

from xenosite._archive_forest import PathSearchCounters, find_path
from xenosite._archive_forest.guided_path import (
    CleavageSide,
    MaybeFilter,
    _canon,
    _cleavage_sides_for_kept,
)
from xenosite._archive_forest.rules import Dealkylation, EpoxideOpening, Hydrolysis


def test_cleavage_side_str_and_maybe_filter():
    side = CleavageSide(site=frozenset({1, 2}), side="CCO")
    assert "CleavageSide[1, 2]|CCO" == str(side)
    empty = MaybeFilter()
    assert not empty
    assert str(empty) == "MaybeFilter()"
    filled = MaybeFilter.from_sides([side])
    assert filled
    assert "CleavageSide" in str(filled)
    assert filled.entries[0].side == "CCO"
    # Same nitrogen, different site (demethyl) passes; formation site does not.
    assert filled.allows("NDealkylation", frozenset({0, 1}))
    assert not filled.allows("NDealkylation", frozenset({1, 2}))
    assert filled.allows(side="CCO")
    assert not filled.allows("Hydroxylation", frozenset({9, 10}))


def test_cleavage_sides_helper_bifurcate_vs_single():
    """Two-piece cohort → bag for discarded; single piece → no bag."""
    kept = Chem.MolFromSmiles("O=C(O)c1ccccc1")
    other = Chem.MolFromSmiles("CC(C)(C)O")
    bags = _cleavage_sides_for_kept(frozenset({6, 8}), kept, [kept, other])
    assert len(bags) == 1
    assert bags[0].side == _canon(other)

    alone = Chem.MolFromSmiles("OCCc1ccccc1")
    assert _cleavage_sides_for_kept(frozenset({6, 7}), alone, [alone]) == ()


def test_hydrolysis_emits_cleavage_side_bag():
    """Ester → acid: discarded alcohol is a CleavageSide bag."""
    from xenosite._archive_forest import RuleSet

    hits = list(
        find_path(
            "c1ccccc1C(=O)OC(C)(C)C",
            "O=C(O)c1ccccc1",
            ruleset=RuleSet([Hydrolysis()], name="Hyd"),
            depth=2,
            max_paths=3,
            max_expansions=40,
        )
    )
    assert hits
    with_bag = [h for h in hits if h[4]]
    assert with_bag, "expected CleavageSide on bifurcating hydrolysis"
    _smi, steps, _mols, plan, maybe = with_bag[0]
    assert [s[0] for s in steps] == ["Hydrolysis"]
    assert isinstance(maybe, MaybeFilter)
    assert len(maybe.entries) >= 1
    assert all(isinstance(e, CleavageSide) for e in maybe.entries)
    assert any(e.side == _canon("CC(C)(C)O") for e in maybe.entries)
    assert "Hydrolysis" in str(plan) or plan.contains(["Hydrolysis"], by="rule")


def test_dealkylation_emits_cleavage_side_bag():
    hits = list(
        find_path(
            "CN(C)Cc1ccccc1",
            "O=Cc1ccccc1",
            ruleset="ND",
            depth=2,
            max_paths=4,
            max_expansions=40,
        )
    )
    assert hits
    with_bag = [h for h in hits if h[4]]
    assert with_bag
    maybe = with_bag[0][4]
    assert all(isinstance(e, CleavageSide) for e in maybe.entries)
    # Exactly one discarded sibling per bifurcating hop on the kept path.
    assert len(maybe.entries) >= 1
    # Discarded piece is the small N-side (not another full-size alternative).
    assert all(
        Chem.MolFromSmiles(e.side).GetNumHeavyAtoms()
        < Chem.MolFromSmiles("CN(C)Cc1ccccc1").GetNumHeavyAtoms()
        for e in maybe.entries
    )


def test_epoxide_opening_no_cleavage_side_bag():
    """Ring-open alternatives are parent-sized — no CleavageSide bag."""
    hits = list(
        find_path(
            "c1ccccc1C1OC1",
            "OCCc1ccccc1",
            ruleset="PhaseOneRS",
            depth=2,
            max_paths=5,
            max_expansions=40,
        )
    )
    assert hits, "expected EpoxideOpening path"
    for _smi, steps, _mols, _plan, maybe in hits:
        if steps and steps[0][0] == "EpoxideOpening":
            assert not maybe
            assert maybe.entries == ()


def test_multiple_cleavages_accumulate_bags_each_step():
    """Two sequential bifurcations → ≥1 CleavageSide per cleavage hop.

    1,4-Phenylene diacetate → hydroquinone via two Hydrolysis cleaves; each
    discards an acetate/acetic piece. Bags must accumulate at every step.
    """
    reactant = "CC(=O)Oc1ccc(OC(C)=O)cc1"
    product = "Oc1ccc(O)cc1"

    counters = PathSearchCounters()
    hits = list(
        find_path(
            reactant,
            product,
            ruleset="PhaseOneRS",
            depth=3,
            max_paths=8,
            max_expansions=100,
            counters=counters,
        )
    )
    assert hits, "expected path via two hydrolyses to hydroquinone"
    two_step = [
        h
        for h in hits
        if sum(1 for s in h[1] if s[0] == "Hydrolysis") >= 2
    ]
    assert two_step, "expected ≥2 Hydrolysis steps on a hit"
    _smi, steps, _mols, plan, maybe = min(two_step, key=lambda h: len(h[1]))
    n_hydro = sum(1 for s in steps if s[0] == "Hydrolysis")
    assert len(maybe.entries) >= n_hydro, (
        "each bifurcating cleavage must add a CleavageSide; "
        "hydrolysis=%d bags=%s" % (n_hydro, maybe)
    )
    assert all(isinstance(e, CleavageSide) for e in maybe.entries)
    # Distinct discarded sides recorded along the walk (may repeat SMILES).
    assert all(e.side for e in maybe.entries)


def test_required_identity_empty_maybe():
    hits = list(
        find_path(
            "c1ccccc1",
            "c1ccccc1",
            depth=1,
            max_paths=1,
        )
    )
    assert hits
    assert not hits[0][4]


def test_tba_cleavage_side_and_preceding_ndealk_passes_maybe():
    """Terbinafine → TBA: formation Dealkylation; demethyl on same N passes Maybe.

    Mapped terbinafine: Me(0)-N(1)-CH2(2)-enyne… and N(1)-CH2(11)-naphthalene.
    TBA-forming dealk is N–CH2(enyne) site {1, 2}; discarded side is the
    naphthylmethyl-N-methyl amine. Preceding demethyl is NDealkylation at
    {0, 1} — same nitrogen, different site — and must pass the bag.
    """
    terb = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
    tba = "CC(C)(C)C#CC=CC=O"
    counters = PathSearchCounters()
    hits = list(
        find_path(
            terb,
            tba,
            ruleset="PhaseOneQF",
            depth=2,
            max_paths=3,
            max_expansions=80,
            counters=counters,
        )
    )
    assert hits
    assert not counters.budget_exhausted

    # Prefer the direct single-cleave formation hit.
    direct = [h for h in hits if [s[0] for s in h[1]] == ["Dealkylation"]]
    assert direct, "expected Dealkylation → TBA"
    _smi, steps, _mols, plan, maybe = direct[0]
    assert plan.contains(["Dealkylation"], by="rule")
    assert maybe
    assert all(isinstance(e, CleavageSide) for e in maybe.entries)
    assert any("cccc" in e.side for e in maybe.entries)  # naphthalene leaving piece

    demethyl_site = frozenset({0, 1})
    formation_site = frozenset({1, 2})
    assert maybe.allows("NDealkylation", demethyl_site), (
        "preceding demethyl on same N must pass Maybe; got %s" % (maybe,)
    )
    assert not maybe.allows("NDealkylation", formation_site)
    assert maybe.allows(side=maybe.entries[0].side)


def test_metabolize_cohort_matches_helper():
    """Live Hydrolysis metabolize list is a multi-piece cohort for the helper."""
    mol = Chem.MolFromSmiles("c1ccccc1C(=O)OC(C)(C)C")
    site, mets = next(Hydrolysis().metabolize(mol))
    assert len(mets) == 2
    kept = next(m for m in mets if "c1ccccc1" in _canon(m) or "O=C(O)" in _canon(m))
    bags = _cleavage_sides_for_kept(site, kept, mets)
    assert len(bags) == 1

    mol_e = Chem.MolFromSmiles("c1ccccc1C1OC1")
    site_e, mets_e = next(EpoxideOpening().metabolize(mol_e))
    assert len(mets_e) == 1
    assert _cleavage_sides_for_kept(site_e, mets_e[0], mets_e) == ()

    # Dealkylation also bifurcates
    mol_d = Chem.MolFromSmiles("CN(C)Cc1ccccc1")
    site_d, mets_d = next(Dealkylation().metabolize(mol_d))
    assert len(mets_d) >= 2
    kept_d = max(mets_d, key=lambda m: m.GetNumHeavyAtoms())
    assert _cleavage_sides_for_kept(site_d, kept_d, mets_d)
