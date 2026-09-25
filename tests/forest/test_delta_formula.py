"""PatternInfo delta_formula annotation + redundancy with adds/removes."""

from __future__ import annotations

from tests.forest.pattern_info_inventory import (
    PATTERNLESS_REACTION_RULE_BASES,
    discover_reaction_rule_classes,
    instantiate_rule,
    patterns_on,
)
from xenosite.forest.rules import (
    Hydroxylation,
    ReductiveDehalogenation,
    bag_counts,
    bag_delta_formula,
    describe,
)


def test_bag_delta_hydroxyl():
    assert bag_delta_formula("O", "H") == {"O": 1, "H": -1}
    assert bag_delta_formula("O", "O") == {}
    assert bag_counts("Cl") == {"Cl": 1}
    assert bag_counts("OH") == {"O": 1, "H": 1}


def test_describe_seals_delta_formula():
    info = describe(adds="O", removes="H", name="h")
    poss = info["possibilities"][0]
    assert poss["delta_formula"] == {"O": 1, "H": -1}
    assert info["span"]["delta_formula"] == {"O": 1, "H": -1}


def test_hydroxylation_patterns_carry_delta():
    for _smarts, info in Hydroxylation().smirks:
        for poss in info["possibilities"]:
            assert poss["delta_formula"] == {"O": 1, "H": -1}


def test_removes_partner_whens_disagree_on_delta():
    info = next(
        i for _s, i in ReductiveDehalogenation().smirks if i.get("name") == "cleave"
    )
    deltas = [tuple(sorted(p["delta_formula"].items())) for p in info["possibilities"]]
    flat = {item for d in deltas for item in d}
    assert ("Cl", -1) in flat
    assert ("F", -1) in flat
    # span collapses disagreeing maps to a tuple
    span_delta = info["span"]["delta_formula"]
    assert isinstance(span_delta, tuple)
    assert len(span_delta) >= 2


def test_delta_formula_matches_adds_removes_across_catalog():
    """adds/removes bags are redundant with sealed delta_formula."""

    mismatched = []
    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        rule = instantiate_rule(cls)
        for _group, _smarts, info in patterns_on(rule):
            for poss in info.get("possibilities") or ():
                expected = bag_delta_formula(
                    poss.get("adds") or "", poss.get("removes") or ""
                )
                got = dict(poss.get("delta_formula") or {})
                if got != expected:
                    mismatched.append(
                        (
                            cls.__name__,
                            info.get("name"),
                            poss.get("when"),
                            expected,
                            got,
                        )
                    )
    assert mismatched == []
