"""Tests for :func:`xenosite.forest.unstable`."""

from __future__ import annotations

import warnings

import pytest

from xenosite.forest.unstable import UnstableWarning, _warned, unstable


@pytest.fixture(autouse=True)
def _clear_unstable_warned():
    _warned.clear()
    yield
    _warned.clear()


def test_unstable_function_warns_once():
    @unstable
    def sample():
        return 7

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UnstableWarning)
        assert sample() == 7
        assert sample() == 7
    unstable_only = [w for w in caught if issubclass(w.category, UnstableWarning)]
    assert len(unstable_only) == 1
    assert "unstable" in str(unstable_only[0].message).lower()


def test_unstable_class_warns_once():
    @unstable
    class Sample:
        def __init__(self, x):
            self.x = x

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UnstableWarning)
        assert Sample(1).x == 1
        assert Sample(2).x == 2
    unstable_only = [w for w in caught if issubclass(w.category, UnstableWarning)]
    assert len(unstable_only) == 1


def test_unstable_with_name_and_detail():
    @unstable(name="demo.api", detail="Expect churn.")
    def sample():
        return 1

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UnstableWarning)
        sample()
    assert len(caught) == 1
    msg = str(caught[0].message)
    assert "demo.api" in msg
    assert "Expect churn" in msg


def test_unstable_rejects_non_callable():
    with pytest.raises(TypeError, match="callable"):
        unstable(3)


def test_path_outcome_unpack_and_helpers():
    from xenosite.forest.guided_path import MaybeFilter, PathOutcome
    from xenosite.forest.step_plan import Step, StepPlan

    plan = StepPlan((Step("Hydroxylation", {0}),))
    outcome = PathOutcome(
        plan=plan,
        maybe=MaybeFilter(),
        steps=(("Hydroxylation", frozenset({0})),),
        smiles=("CCO", "CCO"),
        mols=(),
    )
    smiles, steps, mols, p, maybe = outcome
    assert smiles[-1] == "CCO"
    assert p is plan
    assert len(outcome) == 5
    assert outcome[3] is plan
    assert list(outcome.linearizations())
    assert outcome.contains(["Hydroxylation"], by="rule")
    assert "PathOutcome" in str(outcome)
    assert outcome.allows("Hydroxylation", frozenset({0})) in (True, False)
