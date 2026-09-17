"""Tests for Step / StepPlan partial-order helpers."""

from __future__ import annotations

import json

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest.step_plan import Step, StepPlan


def test_singleton_one_linearization():
    plan = StepPlan.singleton("Epoxidation", {0, 1})
    assert len(plan) == 1
    assert plan.steps[0] == Step("Epoxidation", frozenset({0, 1}))
    assert plan.precedes == ()
    assert list(plan.iter_linearizations()) == [
        (Step("Epoxidation", frozenset({0, 1})),)
    ]


def test_layers_two_prep_then_final_two_orders():
    h0 = Step("Hydroxylation", frozenset({0}))
    h3 = Step("Hydroxylation", frozenset({3}))
    dh = Step("Dehydrogenation", frozenset({0, 3}))
    plan = StepPlan.layers([[h0, h3], [dh]])
    assert len(plan) == 3
    orders = list(plan.iter_linearizations())
    assert len(orders) == 2
    assert set(orders) == {(h0, h3, dh), (h3, h0, dh)}


def test_json_round_trip():
    plan = StepPlan.layers(
        [
            [Step("Hydroxylation", {1}), Step("Hydroxylation", {2})],
            [Step("Dehydrogenation", {1, 2})],
        ]
    )
    restored = StepPlan.from_json(plan.to_json())
    assert restored == plan
    assert list(restored.iter_linearizations()) == list(plan.iter_linearizations())


def test_from_mol_and_attach():
    mol = MolFromSmiles("CCO")
    plan = StepPlan.singleton("Hydroxylation", {0})
    plan.attach_to_mol(mol)
    assert StepPlan.from_mol(mol) == plan
    raw = json.loads(mol.GetProp("phase1_steps"))
    assert raw["steps"][0]["rule"] == "Hydroxylation"


def test_from_mol_missing_prop():
    mol = MolFromSmiles("C")
    with pytest.raises(ValueError, match="phase1_steps"):
        StepPlan.from_mol(mol)


def test_empty_plan():
    plan = StepPlan((), ())
    assert list(plan.iter_linearizations()) == [()]


def test_cycle_raises():
    a = Step("A", {0})
    b = Step("B", {1})
    with pytest.raises(ValueError, match="cycle"):
        list(StepPlan((a, b), ((0, 1), (1, 0))).iter_linearizations())
