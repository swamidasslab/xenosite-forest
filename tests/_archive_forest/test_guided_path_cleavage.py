"""Cleavage breaks formula heuristics — cleave-first search + And post-process."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from xenosite._archive_forest import (
    ADD_O,
    CLEAVE,
    Deps,
    PathSearchCounters,
    RuleSet,
    find_path,
)
from xenosite._archive_forest.guided_path import CleavageSide, and_cleave_plan, _canon
from xenosite._archive_forest.path_context import PathContext, formula_compatible
from xenosite._archive_forest.rules import Dealkylation, Hydrolysis, Hydroxylation
from xenosite._archive_forest.step_plan import And, StepPlan


def _smi(s):
    return Chem.MolFromSmiles(s)


def _n_cleavage_bags(maybe):
    return sum(1 for e in maybe.entries if isinstance(e, CleavageSide))


def test_formula_hints_before_after_half_molecule_cleavage():
    """On a large ester, CLEAVE fits and ADD_O does not; after cleave, ADD_O fits."""
    parent = _smi("CCCCCCCCCCCC(=O)Oc1ccccc1")
    target = _smi("Oc1ccccc1O")
    fragment = _smi("Oc1ccccc1")

    assert CalcMolFormula(parent).startswith("C18")
    assert formula_compatible(CLEAVE, parent, target)
    assert not formula_compatible(ADD_O, parent, target)

    assert formula_compatible(ADD_O, fragment, target)
    assert not formula_compatible(CLEAVE, fragment, target)

    ctx = PathContext.from_mols(parent, target)
    assert Hydrolysis().could_help(parent, target, ctx)
    assert not Hydroxylation().could_help(parent, target, ctx)

    ctx_frag = PathContext.from_mols(fragment, target)
    assert Hydroxylation().could_help(fragment, target, ctx_frag)


def test_long_ester_cleave_then_hydroxylation_and_plan():
    """Half-molecule: long-chain phenyl ester → catechol phenol.

    Search hydrolyzes first (formula-sound), then hydroxylates. Ring carbons
    project to depth-0 origins, so Deps leaves the two steps unordered.
    """
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "CCCCCCCCCCCC(=O)Oc1ccccc1",
            "Oc1ccccc1O",
            ruleset=RuleSet([Hydrolysis(), Hydroxylation()], name="cleave_oh"),
            depth=3,
            maybe_prefixes=False,
            max_paths=2,
            max_expansions=40,
            counters=counters,
        )
    )
    assert hits
    assert not counters.budget_exhausted
    smiles, steps, _mols, plan, maybe = hits[0]
    assert _n_cleavage_bags(maybe) >= 1  # discarded fatty acid
    assert _canon(smiles[-1]) == _canon("Oc1ccccc1O")
    assert [s[0] for s in steps] == ["Hydrolysis", "Hydroxylation"]
    assert isinstance(plan, Deps)
    assert plan.contains(["Hydrolysis", "Hydroxylation"], by="rule")
    assert plan.contains(["Hydroxylation", "Hydrolysis"], by="rule")


def test_anisole_acetate_two_cleaves_and_plan():
    """Acetate + methoxy → hydroquinone: two cleaves commute via Deps."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "CC(=O)Oc1ccc(OC)cc1",
            "Oc1ccc(O)cc1",
            ruleset=RuleSet([Hydrolysis(), Dealkylation()], name="two_cleave"),
            depth=3,
            maybe_prefixes=False,
            max_paths=2,
            max_expansions=40,
            counters=counters,
        )
    )
    assert hits
    smiles, steps, _mols, plan, maybe = hits[0]
    # Hydrolysis + Dealkylation each bifurcate → bags accumulate each hop.
    assert _n_cleavage_bags(maybe) >= 2
    assert _canon(smiles[-1]) == _canon("Oc1ccc(O)cc1")
    assert {s[0] for s in steps} == {"Hydrolysis", "Dealkylation"}
    assert isinstance(plan, Deps)
    # Sites project to depth-0 when atoms existed on the reactant, so the two
    # cleaves remain unordered (both linearizations are members).
    assert plan.contains(["Hydrolysis", "Dealkylation"], by="rule")
    assert plan.contains(["Dealkylation", "Hydrolysis"], by="rule")


def test_biphenyl_ester_cleaves_to_anisole_half():
    """Molecule cleaved in half: keep the methoxy-phenyl fragment + demethylate."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "c1ccccc1C(=O)Oc1ccc(OC)cc1",
            "Oc1ccc(O)cc1",
            ruleset=RuleSet([Hydrolysis(), Dealkylation()], name="half"),
            depth=3,
            maybe_prefixes=False,
            max_paths=2,
            max_expansions=80,
            counters=counters,
        )
    )
    assert hits
    smiles, steps, _mols, plan, maybe = hits[0]
    assert _n_cleavage_bags(maybe) >= 2
    assert _canon(smiles[-1]) == _canon("Oc1ccc(O)cc1")
    assert "Hydrolysis" in [s[0] for s in steps]
    assert "Dealkylation" in [s[0] for s in steps]
    assert isinstance(plan, Deps)
    # Cleavage first in the walk; Deps allows demethyl before hydrolysis.
    assert plan.contains([s[0] for s in steps], by="rule")
    assert plan.contains(["Dealkylation", "Hydrolysis"], by="rule") or plan.contains(
        ["Hydrolysis", "Dealkylation"], by="rule"
    )


def test_and_cleave_plan_empty_and_singleton():
    assert list(and_cleave_plan([]).iter_linearizations()) == [()]
    plan = and_cleave_plan([("Hydroxylation", frozenset({0}))])
    assert len(plan.steps) == 1
    assert plan.n_linearizations() == 1


def test_and_cleave_plan_pure_seq_when_no_cleavage():
    from xenosite._archive_forest.step_plan import Deps

    plan = and_cleave_plan(
        [("Hydroxylation", frozenset({0})), ("Dehydrogenation", frozenset({0, 1}))]
    )
    assert isinstance(plan, Deps)
    assert plan.contains(["Hydroxylation", "Dehydrogenation"], by="rule")
    assert not plan.contains(["Dehydrogenation", "Hydroxylation"], by="rule")


def test_and_cleave_plan_independent_preps_become_and_then_seq():
    """Two OHs with no mutual AtomRef dep → unordered; DH after both."""
    from xenosite._archive_forest.step_plan import AtomRef, Deps

    oh0 = ("Hydroxylation", frozenset({AtomRef(0)}))
    oh1 = ("Hydroxylation", frozenset({AtomRef(1)}))
    dh = (
        "Dehydrogenation",
        frozenset(
            {
                AtomRef(added_by=("Hydroxylation", frozenset({0}))),
                AtomRef(added_by=("Hydroxylation", frozenset({1}))),
            }
        ),
    )
    plan = and_cleave_plan([oh0, oh1, dh])
    assert isinstance(plan, Deps)
    # Interleaving of the two OHs is allowed before DH.
    assert plan.contains(
        ["Hydroxylation", "Hydroxylation", "Dehydrogenation"], by="rule"
    )
    from xenosite._archive_forest.guided_path import _canonical_plan_key

    assert _canonical_plan_key(plan) == _canonical_plan_key(
        and_cleave_plan([oh1, oh0, dh])
    )


def test_find_path_dedupes_identical_stepplans():
    """And-equivalent cleavage orderings emit one PathOutcome."""
    from xenosite._archive_forest.guided_path import _canonical_plan_key

    hits = list(
        find_path(
            "CN(C)Cc1ccc(CN(C)C)cc1",
            "O=Cc1ccc(C=O)cc1",
            ruleset="ND",
            depth=3,
            max_paths=20,
            max_expansions=80,
        )
    )
    assert hits
    keys = [_canonical_plan_key(h.plan) for h in hits]
    assert len(keys) == len(set(keys)), "duplicate StepPlans must not be emitted"


def test_hydrolysis_is_cleavage_hydroxylation_is_not():
    assert Hydrolysis().is_cleavage()
    assert Dealkylation().is_cleavage()
    assert not Hydroxylation().is_cleavage()
