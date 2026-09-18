"""Direct tests for error branches and helpers the integration suite skips."""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.guided_path import (
    CleavageSide,
    MaybeFilter,
    PathSearchCounters,
    _advance_opens_and_bags,
    _canon,
    _canonicalize_plan_json,
    _canonical_plan_key,
    _consume_depth,
    _depth_remaining,
    _expansion_cohorts,
    _flatten_rule_instances,
    _flatten_rules,
    _is_cleavage_rule,
    _pack_hit,
    _plan_from_dep_layers,
    _product_pieces,
    _rule_name_is_cleavage,
    _rule_smarts_covers_site,
    _rules_for_mol,
    _site_atoms,
    _site_key_set,
    _step_depends_on,
    _step_from_path_item,
    _steps_key,
    _walk_mol_chain,
    and_cleave_plan,
    coerce_path_counters,
)
from xenosite.forest.path_context import (
    ADD_O,
    NEUTRAL,
    FormulaHint,
    FormulaMatch,
    _site_atom_indices,
    any_hint_compatible,
    cleavage_site_bond,
    cleavage_site_may_reach,
    elements_may_add_from_hints,
    expand_formula_effects,
    formula_compatible,
    unreachable_new_elements,
)
from xenosite.forest.rules import Hydrolysis, Hydroxylation
from xenosite.forest.rulesets import RuleSet
from xenosite.forest.step_plan import (
    And,
    AtomRef,
    Deps,
    Linearization,
    Or,
    Step,
    StepPlan,
    count_transformation_orders,
    or_of_plans,
    pathway_alternatives,
    pathway_from_json,
    pathway_linearizations,
    pathway_to_json,
)


def _mol(smi):
    return Chem.MolFromSmiles(smi)


class _Adds:
    def __init__(self, adds):
        self._adds = adds

    def elements_may_add(self):
        return self._adds


class _AddsBoom:
    def elements_may_add(self):
        raise RuntimeError("no hints")


def test_formula_expand_error_branches():
    assert expand_formula_effects(None) == ()
    assert formula_compatible(None, _mol("C"), _mol("C")) is True

    spec = FormulaMatch(
        possible=(ADD_O,),
        resolve=lambda mol, match: NEUTRAL,
    )
    assert expand_formula_effects(spec, match=None) == (ADD_O,)
    assert expand_formula_effects(spec, _mol("C"), match=object()) == (NEUTRAL,)

    def one_arg(mol, match=None):
        raise TypeError("old")

    def old(mol):
        return ADD_O

    # TypeError on (mol, match), then the one-arg fallback also fails.
    assert expand_formula_effects(one_arg, _mol("C"), match=object()) == ()

    def raises(mol, match=None):
        raise RuntimeError("bad")

    assert expand_formula_effects(raises) == ()
    assert expand_formula_effects(lambda mol, match: None) == ()
    assert expand_formula_effects((ADD_O, NEUTRAL)) == (ADD_O, NEUTRAL)
    with pytest.raises(TypeError, match="unknown"):
        expand_formula_effects(3)

    # One-arg resolver still expands.
    assert expand_formula_effects(old, _mol("C"), match=object()) == (ADD_O,)


def test_formula_compatible_exception_and_negative_delta():
    class Boom:
        def GetAtoms(self):
            raise RuntimeError("no")

    assert FormulaHint(cleave=True).compatible(Boom(), Boom()) is True

    class NoHeavy:
        def GetAtoms(self):
            return []

        def GetNumHeavyAtoms(self):
            raise RuntimeError("ha")

    assert FormulaHint(cleave=True).compatible(NoHeavy(), NoHeavy()) is True
    # Target already has the oxygen a removal would take away.
    assert not FormulaHint(delta={"O": -1}).compatible(_mol("CO"), _mol("CO"))

    assert elements_may_add_from_hints(None) is None
    assert elements_may_add_from_hints([None]) is None
    assert elements_may_add_from_hints([NEUTRAL]) == frozenset()
    assert "O" in elements_may_add_from_hints([ADD_O])

    carbon = _mol("C")
    assert any_hint_compatible(None, carbon, carbon)
    assert any_hint_compatible([None], carbon, carbon)
    assert not any_hint_compatible([ADD_O], carbon, carbon)


def test_unreachable_elements_and_site_indices():
    carbon = _mol("C")
    iron = _mol("[Fe]")
    amine = _mol("N")
    assert "Fe" in unreachable_new_elements([_Adds(frozenset({"O"}))], carbon, iron)
    assert "Fe" in unreachable_new_elements([_Adds(None)], carbon, iron)
    assert "Fe" in unreachable_new_elements([_AddsBoom()], carbon, iron)
    assert "N" in unreachable_new_elements([_Adds(frozenset({"O"}))], carbon, amine)
    assert unreachable_new_elements([_Adds(frozenset())], carbon, carbon) == frozenset()
    assert unreachable_new_elements([], object(), object()) == frozenset()

    assert _site_atom_indices(("label", [1, 2])) == frozenset({1, 2})

    class Origin:
        origin = 4

    class Idx:
        idx = 2

    class Junk:
        origin = "no"
        idx = "no"

    assert _site_atom_indices([Origin(), Idx(), Junk(), "z"]) == frozenset({4, 2})
    assert cleavage_site_bond(_mol("CCO"), frozenset({0})) is None
    assert cleavage_site_may_reach(_mol("CCO"), frozenset({0, 1}), _mol("CO"), None)


def test_coerce_counters_and_small_path_helpers():
    assert coerce_path_counters(None) is None
    fresh = coerce_path_counters(None, max_expansions=3)
    assert isinstance(fresh, PathSearchCounters)
    again = coerce_path_counters(fresh)
    assert again is fresh
    mirrored = {}
    coerced = coerce_path_counters(
        {"rule_expansions": None, "budget_exhausted": True},
    )
    assert coerced.rule_expansions == 0
    assert coerced.budget_exhausted
    coerced._mirror = mirrored
    coerced._sync()
    assert "billed" in mirrored
    with pytest.raises(TypeError, match="counters"):
        coerce_path_counters(3)

    spent = PathSearchCounters(linearizations_applied=2)
    assert spent.under_budget(None)
    assert not spent.under_budget(2)
    assert spent.budget_exhausted

    assert _canon("not-a-smiles") == "not-a-smiles"
    stereo = _canon("C/C=C/C")
    assert "/" not in stereo and "\\" not in stereo

    inner = RuleSet([Hydrolysis()], name="inner")
    outer = RuleSet([Hydroxylation(), inner], name="outer")
    assert len(_flatten_rule_instances(outer)) == 2
    assert _flatten_rules(outer) is outer
    assert _flatten_rules("ND").name

    step = Step("Hydroxylation", {0})
    assert _step_from_path_item(step) is step
    wrapped = _step_from_path_item(("Hydroxylation", ("site", {1}), ("methide",)))
    assert wrapped.pathways == frozenset({"methide"})
    assert wrapped.site == frozenset({AtomRef(origin=1)})

    assert _rule_name_is_cleavage("Hydrolysis")
    assert not _rule_name_is_cleavage("NoSuchRule")

    later = Step(
        "Dehydrogenation",
        {AtomRef(added_by=("Hydroxylation", frozenset({0})))},
    )
    assert _step_depends_on(later, step)
    other = Step(
        "Dehydrogenation",
        {AtomRef(added_by=("Other", frozenset({0})))},
    )
    assert not _step_depends_on(other, step)

    assert _site_key_set(None) == frozenset()
    assert 0 in _site_key_set([AtomRef(origin=0), "x"])
    assert _site_atoms(("name", {3})) == frozenset({3})
    assert _steps_key([("Hydroxylation", {0})])
    walked = _walk_mol_chain([], [_mol("C")])
    assert len(walked) == 1
    assert _walk_mol_chain([step], [None]) == []
    assert _depth_remaining(None) and not _depth_remaining(0)
    assert _consume_depth(None) is None
    assert _consume_depth(4, 2) == 2


def test_plan_json_canonical_and_dep_layers():
    assert _canonical_plan_key(None)[0] == "raw"
    assert _canonicalize_plan_json(None) == ()
    assert _canonicalize_plan_json(3) == 3
    assert _canonicalize_plan_json(1.5) == "1.5"
    assert _canonicalize_plan_json([1, "a"]) == (1, "a")
    raw = _canonicalize_plan_json({"op": "mystery", "k": 1})
    assert raw[0] == "raw"

    a = {"op": "step", "rule": "A", "site": [0], "pathways": ["z", "a"]}
    b = {"op": "step", "rule": "B", "site": [{"added_by": "A", "at": [0]}]}
    and_key = _canonicalize_plan_json({"op": "and", "children": [b, a]})
    assert and_key[0] == "and"
    or_key = _canonicalize_plan_json({"op": "or", "children": [a, b]})
    assert or_key[0] == "or"
    seq_key = _canonicalize_plan_json({"op": "seq", "children": [a, b]})
    assert seq_key[0] == "seq"
    # Cycle falls back to the sorted edge set.
    cyc = _canonicalize_plan_json(
        {"op": "deps", "steps": [a, b], "precedes": [[0, 1], [1, 0]]}
    )
    assert cyc[0] == "deps"
    skipped = _canonicalize_plan_json(
        {"op": "plan", "steps": [a], "precedes": [[5, 6]]}
    )
    assert skipped[2] == ()

    oh = Step("Hydroxylation", {0})
    dh = Step("Dehydrogenation", {1})
    same = Step("Hydroxylation", {2})
    assert isinstance(_plan_from_dep_layers([]), StepPlan)
    assert len(_plan_from_dep_layers([oh]).steps) == 1
    ordered = _plan_from_dep_layers([oh, dh])
    assert ordered.precedes
    assert not _plan_from_dep_layers([oh, same]).precedes
    hyd = Step("Hydrolysis", {0, 1})
    assert not _plan_from_dep_layers([hyd, oh]).precedes

    packed = _pack_hit(["CCO"], [("Hydroxylation", frozenset({0}))], [_mol("CCO")])
    assert packed.smiles
    empty = _pack_hit([], [], [], maybe=[CleavageSide(frozenset({0}), "C")])
    assert empty.maybe
    assert and_cleave_plan([]).steps == ()


def test_product_pieces_cohorts_and_smarts_cover():
    assert _product_pieces(None) == []
    ethanol = _mol("CCO")
    bits = _mol("C.O")
    assert len(_product_pieces(bits)) == 2
    assert _product_pieces([ethanol, None]) == [ethanol]
    with pytest.raises(TypeError, match="molecule"):
        _product_pieces(object())

    assert _expansion_cohorts(ethanol, []) == []
    assert _expansion_cohorts(ethanol, [ethanol]) == [[ethanol]]
    # Same-size pieces are not a cleavage pair.
    same_size = _expansion_cohorts(ethanol, [ethanol, _mol("CCN")])
    assert len(same_size) == 2 and all(len(c) == 1 for c in same_size)
    small = [_mol("C"), _mol("O")]
    paired = _expansion_cohorts(ethanol, small)
    assert len(paired) == 1 and len(paired[0]) == 2
    odd = _expansion_cohorts(ethanol, [_mol("C"), _mol("O"), _mol("N")])
    assert len(odd) == 1
    with pytest.raises(TypeError, match="parent"):
        _expansion_cohorts(object(), small)

    class Cleaves:
        def is_cleavage(self):
            return True

        def cleave_alone(self):
            return False

    class Plain:
        def is_cleavage(self):
            return False

        def cleave_alone(self):
            return False

    assert _is_cleavage_rule(Cleaves())
    assert not _is_cleavage_rule(object())
    bags, opens = _advance_opens_and_bags(
        Cleaves(), frozenset({0}), ethanol, [ethanol, _mol("C")], (), ()
    )
    assert bags and opens == ()
    bags, opens = _advance_opens_and_bags(
        Cleaves(), frozenset({1}), ethanol, [ethanol], (), ()
    )
    assert opens
    assert _advance_opens_and_bags(Plain(), frozenset({1}), ethanol, [ethanol], (), ()) == (
        (),
        (),
    )

    assert _rules_for_mol([Plain(), Cleaves()], ethanol, ethanol)[0]
    only = _rules_for_mol([Plain(), Cleaves()], ethanol, _mol("C"))
    assert len(only) == 1
    assert _rules_for_mol([Plain()], ethanol, _mol("C"))
    with pytest.raises(TypeError, match="molecule"):
        _rules_for_mol([Plain()], object(), object())

    class Queries:
        def match_queries(self, mol):
            return {"q": [({1: 0},)]}

    assert _rule_smarts_covers_site(Queries(), ethanol, {0})
    assert not _rule_smarts_covers_site(Queries(), ethanol, set())

    class QueriesBoom:
        def match_queries(self, mol):
            raise RuntimeError("q")

        def iter_reactant_site_matches(self, mol):
            yield frozenset({0, 1})

    assert _rule_smarts_covers_site(QueriesBoom(), ethanol, {0})

    class SitesBoom:
        def iter_reactant_site_matches(self, mol):
            raise RuntimeError("sites")

        smarts = ("not smarts", "[C:1]>>[C:1]")
        mapid_site = [99]

    assert _rule_smarts_covers_site(SitesBoom(), ethanol, {0}) in (True, False)
    assert not _rule_smarts_covers_site(SitesBoom(), None, {0})


def test_step_plan_reprs_or_and_edges():
    with pytest.raises(ValueError):
        AtomRef(origin=1, added_by=("X", frozenset({0})))
    with pytest.raises(ValueError):
        AtomRef()
    with pytest.raises(ValueError):
        AtomRef(added_by="X")
    with pytest.raises(ValueError):
        AtomRef(added_by=("X",))
    ref = AtomRef(added_by=("Hydroxylation", [0, 1]))
    assert isinstance(ref.added_by[1], frozenset)
    assert AtomRef.from_json({"added_by": "Legacy", "at": [2]}).added_by == (
        "Legacy",
        frozenset({2}),
    )
    assert AtomRef.from_json({"added_by": ["New", [3]]}).added_by == (
        "New",
        frozenset({3}),
    )
    mid = AtomRef(added_by=("Hydroxylation", frozenset({1})), depth=2)
    assert mid.depth == 2
    assert mid.to_json() == {
        "added_by": ["Hydroxylation", [1]],
        "depth": 2,
    }
    assert AtomRef.from_json(mid.to_json()) == mid
    assert str(mid) == "Hydroxylation[1]@2"
    assert "Hydroxylation" in repr(ref)
    assert AtomRef.coerce(3).origin == 3
    assert AtomRef.coerce({"origin": 2}).origin == 2
    with pytest.raises(TypeError):
        AtomRef.coerce("no")
    with pytest.raises(ValueError):
        AtomRef(origin=0, depth=-1)

    tagged = Step("Dehydrogenation", {0}, pathways={"methide"})
    assert "methide" in repr(tagged) and "methide" in str(tagged)

    a, b, c = Step("A", {0}), Step("B", {1}), Step("C", {2})
    assert str(And(())) == "And()"
    assert str(And((a,))) == str(a)
    assert list(And(()).iter_linearizations()) == [()]
    assert str(Deps()) == "Deps()"
    assert str(Deps([a])) == str(a)

    empty_or = Or(())
    assert not empty_or and len(empty_or) == 0 and str(empty_or) == "Or()"
    assert list(empty_or.iter_linearizations()) == []
    seq = StepPlan((a, b))
    choice = Or((seq, c))
    assert "|" in str(choice) and "(" in str(choice)
    assert Or((a,)).steps == (a,)
    assert Or((a, b)).steps == ()
    assert Or((a, b)).n_linearizations() == 2

    assert isinstance(or_of_plans([]), Or)
    assert or_of_plans([a]).steps == (a,)
    assert isinstance(or_of_plans([a, b]), Or)
    assert pathway_alternatives(None) == []
    assert pathway_alternatives("nope") == []
    assert pathway_alternatives(a)
    assert list(pathway_linearizations(a))
    with pytest.raises(TypeError):
        StepPlan((3,))
    with pytest.raises(TypeError):
        list(pathway_linearizations(3))
    with pytest.raises(TypeError):
        pathway_to_json(3)
    with pytest.raises(TypeError):
        pathway_from_json([])

    assert StepPlan.layers([[]]).steps == ()
    layered = StepPlan.layers([[], [a, b], [c]])
    assert layered.precedes
    assert StepPlan((Or((a, b)),)).precedes == ()
    assert Deps([a, b], precedes=[(0, 1)]).precedes == ((0, 1),)

    plan = StepPlan((a,))
    assert plan.contains(a)
    assert plan.contains(Linearization((a,)))
    assert not plan.contains(["B"], by="rule")
    with pytest.raises(ValueError, match="by"):
        plan.contains(["A"], by="nope")
    with pytest.raises(TypeError):
        plan.contains(["A"])

    mol = _mol("C")
    assert Linearization(()).apply(mol)
    with pytest.raises(ValueError, match="drop_last"):
        Linearization((a,)).apply(mol, drop_last=-1)
    assert Linearization((a,)).apply(mol, drop_last=3)

    with pytest.raises(ValueError, match="invalid"):
        count_transformation_orders(2, [(0, 9)])
    assert not Deps([a, b]).same_linearizations(Deps([a]))
    assert "StepPlan" in repr(plan)
