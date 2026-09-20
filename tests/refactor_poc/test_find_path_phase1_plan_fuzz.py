"""Create a phase-I or quinone product, then ``find_path`` it back.

A hit matches the target, is one component, and does not record an opaque
``QuinoneFormation`` step (the canonical plan is hydroxylation then
dehydrogenation). A miss is empty.
"""

from __future__ import annotations

from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from rdkit import Chem

from xenosite.refactor_poc.find_path import find_path
from xenosite.refactor_poc.rules import (
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
)
from xenosite.refactor_poc.rulesets import PhaseOne

_CORPUS = ("c1ccccc1", "Oc1ccccc1", "Oc1ccc(O)cc1", "Cc1ccccc1")


def _create_rules():
    return [QuinoneFormation(), Hydroxylation(), Dehydrogenation()]


def _candidates(mol, rules, seen: set[str]):
    out = []
    for rule in rules:
        for product, _info in rule.metabolize(mol):
            smi = product.xf.csmi
            if not smi or "." in smi or smi in seen:
                continue
            if product.GetNumHeavyAtoms() < max(4, mol.GetNumHeavyAtoms() // 3):
                continue
            out.append(smi)
            if len(out) >= 8:
                return out
    return out


def _assert_plan(outcome, target: str):
    assert "." not in outcome.smiles
    assert outcome.smiles == target
    names = [step.rule for step in outcome.plan.children]
    assert names
    assert "QuinoneFormation" not in names


@given(start=st.sampled_from(_CORPUS), data=st.data())
@settings(
    max_examples=4,
    deadline=25_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
def test_fuzz_created_target_hit_or_honest_miss(start: str, data):
    mol = Chem.MolFromSmiles(start)
    assume(mol is not None)
    seen = {mol.xf.csmi}
    pool = _candidates(mol, _create_rules(), seen)
    assume(pool)
    target = data.draw(st.sampled_from(pool))
    hits = list(
        find_path(start, target, ruleset=PhaseOne, max_nodes=150, max_paths=1)
    )
    if not hits:
        return
    _assert_plan(hits[0], target)


def test_benzene_quinone_plan_is_elementary():
    target = "O=C1C=CC(=O)C=C1"
    hits = list(
        find_path("c1ccccc1", target, ruleset=PhaseOne, max_nodes=400, max_paths=1)
    )
    assert hits
    _assert_plan(hits[0], Chem.MolFromSmiles(target).xf.csmi)
    names = [step.rule for step in hits[0].plan.children]
    assert names.count("Hydroxylation") == 2
    assert names.count("Dehydrogenation") == 1
