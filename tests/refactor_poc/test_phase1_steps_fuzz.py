"""Canonical plans for quinone, epoxidation, and N-dealkylation.

``And`` / ``Or`` trees and plan JSON are not the poc surface. The check is
the elementary plan the rule reports: quinone ends in dehydrogenation after
the preps that supply oxygen; epoxidation and N-dealkylation are themselves.
"""

from __future__ import annotations

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from rdkit import Chem

from xenosite.refactor_poc.rules import (
    Epoxidation,
    NDealkylation,
    QuinoneFormation,
)

_CORPUS = (
    "c1ccccc1",
    "Oc1ccccc1",
    "Oc1ccc(O)cc1",
    "Nc1ccc(O)cc1",
    "C=C",
    "CCN",
)


def test_benzene_quinone_plan_ends_in_dehydrogenation():
    mol = Chem.MolFromSmiles("c1ccccc1")
    plans = []
    for _product, info in QuinoneFormation().metabolize(mol):
        steps = QuinoneFormation().canonical_plan(mol, info)
        if steps:
            plans.append(steps)
        if len(plans) >= 4:
            break
    assert plans
    for steps in plans:
        assert steps[-1].rule == "Dehydrogenation"
        assert any(step.rule == "Hydroxylation" for step in steps[:-1])


def test_epoxidation_plan_is_one_step():
    mol = Chem.MolFromSmiles("C=C")
    product, info = next(Epoxidation().metabolize(mol))
    steps = Epoxidation().canonical_plan(mol, info)
    assert [step.rule for step in steps] == ["Epoxidation"]
    assert "." not in product.xf.csmi


def test_ndealkylation_plan_is_one_step():
    mol = Chem.MolFromSmiles("CCN")
    product, info = next(NDealkylation().metabolize(mol))
    steps = NDealkylation().canonical_plan(mol, info)
    assert [step.rule for step in steps] == ["NDealkylation"]
    assert "." not in product.xf.csmi


@given(smiles=st.sampled_from(_CORPUS))
@settings(
    max_examples=6,
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_fuzz_quinone_plans_end_in_dehydrogenation(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    seen = 0
    for _product, info in QuinoneFormation().metabolize(mol):
        steps = QuinoneFormation().canonical_plan(mol, info)
        if not steps:
            continue
        seen += 1
        assert steps[-1].rule == "Dehydrogenation"
        assert "." not in _product.xf.csmi
        if seen >= 4:
            return
    assume(seen)
