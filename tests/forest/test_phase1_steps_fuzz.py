"""Canonical plans for quinone, epoxidation, and N-dealkylation.

``And`` / ``Or`` trees and plan JSON are not the forest surface. The check is
the elementary plan the rule reports: quinone ends in dehydrogenation after
the preps that supply oxygen; epoxidation and N-dealkylation are themselves.

``canonical_emitted_sites`` is drawn with Hypothesis ``st.booleans()``.
"""

from __future__ import annotations

import os as _os

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from rdkit import Chem

from xenosite.forest.rules import (
    Epoxidation,
    NDealkylation,
    QuinoneFormation,
)


def _fuzz_examples(default: int) -> int:
    raw = _os.environ.get("XENOSITE_FUZZ_EXAMPLES")
    if raw:
        return int(raw)
    return default

_CORPUS = (
    "c1ccccc1",
    "Oc1ccccc1",
    "Oc1ccc(O)cc1",
    "Nc1ccc(O)cc1",
    "C=C",
    "CCN",
)


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=20_000, derandomize=True)
def test_benzene_quinone_plan_ends_in_dehydrogenation(canonical_emitted_sites: bool):
    mol = Chem.MolFromSmiles("c1ccccc1")
    plans = []
    for _products, info in QuinoneFormation().metabolize(
        mol, canonical_emitted_sites=canonical_emitted_sites
    ):
        steps = QuinoneFormation().canonical_plan(mol, info)
        if steps:
            plans.append(steps)
        if len(plans) >= 4:
            break
    assert plans
    for steps in plans:
        assert steps[-1].rule == "Dehydrogenation"
        assert any(step.rule == "Hydroxylation" for step in steps[:-1])


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=10_000, derandomize=True)
def test_epoxidation_plan_is_one_step(canonical_emitted_sites: bool):
    mol = Chem.MolFromSmiles("C=C")
    products, info = next(
        Epoxidation().metabolize(
            mol, canonical_emitted_sites=canonical_emitted_sites
        )
    )
    product = products[0]
    steps = Epoxidation().canonical_plan(mol, info)
    assert [step.rule for step in steps] == ["Epoxidation"]
    assert "." not in product.xf.csmi


@given(canonical_emitted_sites=st.booleans())
@settings(max_examples=_fuzz_examples(4), deadline=10_000, derandomize=True)
def test_ndealkylation_plan_is_one_step(canonical_emitted_sites: bool):
    mol = Chem.MolFromSmiles("CCN")
    products, info = next(
        NDealkylation().metabolize(
            mol, canonical_emitted_sites=canonical_emitted_sites
        )
    )
    product = products[0]
    steps = NDealkylation().canonical_plan(mol, info)
    assert [step.rule for step in steps] == ["NDealkylation"]
    assert "." not in product.xf.csmi


@given(smiles=st.sampled_from(_CORPUS), canonical_emitted_sites=st.booleans())
@settings(
    max_examples=_fuzz_examples(6),
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_fuzz_quinone_plans_end_in_dehydrogenation(
    smiles: str, canonical_emitted_sites: bool
):
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    seen = 0
    for products, info in QuinoneFormation().metabolize(
        mol, canonical_emitted_sites=canonical_emitted_sites
    ):
        steps = QuinoneFormation().canonical_plan(mol, info)
        if not steps:
            continue
        seen += 1
        assert steps[-1].rule == "Dehydrogenation"
        assert all("." not in p.xf.csmi for p in products)
        if seen >= 4:
            return
    assume(seen)
