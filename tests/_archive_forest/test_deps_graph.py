"""Unit coverage for Deps graph algorithms (count / canon / identity)."""

from __future__ import annotations

import math

import pytest

from xenosite._archive_forest.step_plan import (
    ComponentTooLarge,
    Deps,
    Step,
    StepPlan,
    And,
    canonical_dependency_edges,
    count_transformation_orders,
    transitive_closure_masks,
)


def test_count_empty_and_singleton():
    assert count_transformation_orders(0, []) == 1
    assert count_transformation_orders(1, []) == 1
    assert Deps().n_linearizations() == 1
    assert Deps([Step("A", {0})]).n_linearizations() == 1


def test_count_free_nodes_is_factorial():
    assert count_transformation_orders(5, []) == math.factorial(5)


def test_count_two_components_interleave():
    # {0→1} and {2→3}: 2*2 * C(4,2) = 4*6 = 24? 
    # Component sizes 2,2: first contrib 1 * C(0+2,2)=1; then 1 * C(2+2,2)=6 → 6
    # Wait each chain has 1 lin: total = 1*1*C(4,2)=6. Yes (interleave two ordered pairs).
    assert count_transformation_orders(4, [(0, 1), (2, 3)]) == 6


def test_count_cycle_returns_zero():
    assert count_transformation_orders(2, [(0, 1), (1, 0)]) == 0


def test_count_component_too_large():
    n = 22
    # One connected chain → single component of size 22
    edges = [(i, i + 1) for i in range(n - 1)]
    with pytest.raises(ComponentTooLarge):
        count_transformation_orders(n, edges, max_component_size=20)
    # Unlimited still works (chain has 1 linearization)
    assert count_transformation_orders(n, edges, max_component_size=None) == 1


def test_transitive_closure_and_canonical():
    closure = transitive_closure_masks(3, [(0, 1), (1, 2), (0, 2)])
    assert closure[0] & (1 << 1) and closure[0] & (1 << 2)
    assert closure[1] & (1 << 2)
    assert canonical_dependency_edges(3, [(0, 1), (1, 2), (0, 2)]) == (
        (0, 1),
        (1, 2),
    )
    assert canonical_dependency_edges(3, [(0, 1), (1, 2)]) == canonical_dependency_edges(
        3, [(0, 1), (1, 2), (0, 2)]
    )


def test_transitive_closure_cycle_raises():
    with pytest.raises(ValueError, match="cycle"):
        transitive_closure_masks(2, [(0, 1), (1, 0)])
    with pytest.raises(ValueError, match="invalid|cycle"):
        transitive_closure_masks(2, [(0, 0)])


def test_deps_stores_canonical_precedes():
    a, b, c = Step("A", {0}), Step("B", {1}), Step("C", {2})
    d = Deps([a, b, c], precedes=[(0, 1), (1, 2), (0, 2)])
    assert d.precedes == ((0, 1), (1, 2))
    assert d.to_json()["precedes"] == [[0, 1], [1, 2]]


def test_deps_from_layers_precedes():
    a, b, c = Step("A", {0}), Step("B", {1}), Step("C", {2})
    layered = StepPlan.layers([[a, b], [c]])
    d = Deps(layered.steps, layered.precedes)
    assert isinstance(d, Deps)
    assert d.same_linearizations(Deps([a, b, c], precedes=[(0, 2), (1, 2)]))


def test_deps_json_roundtrip_op_deps():
    a, b = Step("A", {0}), Step("B", {1})
    d = Deps([a, b], precedes=[(0, 1)])
    restored = StepPlan.from_json(d.to_json())
    assert isinstance(restored, Deps)
    assert restored.same_linearizations(d)


def test_deps_invalid_edge():
    a, b = Step("A", {0}), Step("B", {1})
    with pytest.raises(ValueError):
        Deps([a, b], precedes=[(0, 5)])
    with pytest.raises(ValueError, match="cycle"):
        Deps([a, b], precedes=[(0, 1), (1, 0)])


def test_deps_contains_and_str():
    a, b, c = Step("A", {0}), Step("B", {1}), Step("C", {2})
    d = Deps([a, b, c], precedes=[(0, 2), (1, 2)])
    assert d.contains([a, b, c])
    assert d.contains([b, a, c])
    assert not d.contains([a, c, b])
    assert "≺" in str(d) or "A" in str(d)
    free = Deps([a, b])
    assert "&" in str(free)


def test_align_duplicate_steps_greedy():
    """Identical Steps (rare) still align via multiset match."""
    a1 = Step("A", {0})
    a2 = Step("A", {0})  # equal to a1
    b = Step("B", {1})
    d1 = Deps([a1, a2, b], precedes=[(0, 2), (1, 2)])
    d2 = Deps([a2, b, a1], precedes=[(0, 1), (2, 1)])
    # Same multiset + same constraints after align
    assert d1.same_linearizations(d2)


def test_and_n_linearizations_matches_deps_free():
    a, b, c = Step("A", {0}), Step("B", {1}), Step("C", {2})
    assert And((a, b, c)).n_linearizations() == Deps([a, b, c]).n_linearizations()
