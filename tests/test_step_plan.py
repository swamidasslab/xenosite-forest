"""Tests for Step / StepPlan / AtomRef apply+resolve helpers."""

from __future__ import annotations

import json

import pytest
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest.step_plan import AtomRef, Step, StepPlan
from xenosite.forest.utils import unmapped_smiles


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
    dh = Step(
        "Dehydrogenation",
        frozenset(
            [
                AtomRef(added_by="Hydroxylation", at=frozenset({0})),
                AtomRef(added_by="Hydroxylation", at=frozenset({3})),
            ]
        ),
    )
    plan = StepPlan.layers([[h0, h3], [dh]])
    assert len(plan) == 3
    orders = list(plan.iter_linearizations())
    assert len(orders) == 2
    assert set(orders) == {(h0, h3, dh), (h3, h0, dh)}


def test_json_round_trip():
    plan = StepPlan.layers(
        [
            [Step("Hydroxylation", {1}), Step("Hydroxylation", {2})],
            [
                Step(
                    "Dehydrogenation",
                    {
                        AtomRef(added_by="Hydroxylation", at={1}),
                        AtomRef(added_by="Hydroxylation", at={2}),
                    },
                )
            ],
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


def test_compact_str():
    assert str(AtomRef(origin=3)) == "3"
    assert str(AtomRef(added_by="Hydroxylation", at={0})) == "Hydroxylation[0]"
    assert str(Step("Hydroxylation", {0})) == "Hydroxylation[0]"
    dh = Step(
        "Dehydrogenation",
        {
            AtomRef(added_by="Hydroxylation", at={0}),
            AtomRef(added_by="Hydroxylation", at={3}),
        },
    )
    assert str(dh) == "Dehydrogenation[Hydroxylation[0], Hydroxylation[3]]"
    plan = StepPlan.layers(
        [
            [Step("Hydroxylation", {0}), Step("Hydroxylation", {3})],
            [dh],
        ]
    )
    assert str(plan) == (
        "(Hydroxylation[0] & Hydroxylation[3]) → "
        "Dehydrogenation[Hydroxylation[0], Hydroxylation[3]]"
    )


def test_hydroxylation_apply_records_atom_ref():
    mol = MolFromSmiles("c1ccccc1")
    step = Step("Hydroxylation", {0})
    products = step.apply(mol)
    assert products
    product = products[0]
    ref = AtomRef(added_by="Hydroxylation", at=frozenset({0}))
    o_idx = ref.resolve(product)
    assert product.GetAtomWithIdx(o_idx).GetAtomicNum() == 8
    carbon_site = step.resolve_site(product)
    assert len(carbon_site) == 1


def test_benzene_prep_linearization_apply_agrees():
    mol = MolFromSmiles("c1ccccc1")
    h0 = Step("Hydroxylation", {0})
    h3 = Step("Hydroxylation", {3})
    dh = Step(
        "Dehydrogenation",
        frozenset(
            [
                AtomRef(added_by="Hydroxylation", at={0}),
                AtomRef(added_by="Hydroxylation", at={3}),
            ]
        ),
    )
    plan = StepPlan.layers([[h0, h3], [dh]])
    prep_smiles = [
        frozenset(unmapped_smiles(p) for p in products)
        for lin in plan.iter_as_linearizations()
        for products in [lin.apply(mol, drop_last=1)]
        if products
    ]
    assert prep_smiles and len(set(prep_smiles)) == 1
    got = Chem.MolToSmiles(Chem.MolFromSmiles(next(iter(next(iter(prep_smiles))))))
    assert got == Chem.MolToSmiles(Chem.MolFromSmiles("Oc1ccc(O)cc1"))
