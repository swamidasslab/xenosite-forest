"""Cleavage walks: the kept piece matches the target; the other side is not searched.

``and_cleave_plan`` replay is not the poc surface. ``PathOutcome.maybe`` is
the discarded fragment.
"""

from __future__ import annotations

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from rdkit import Chem

from xenosite.refactor_poc.find_path import PathCounters, find_path
from xenosite.refactor_poc.rules import Dealkylation
from xenosite.refactor_poc.rulesets import PhaseOne

_ETHERS = ("COc1ccccc1", "CCOc1ccccc1", "COc1ccc(O)cc1")


def test_anisole_demethylation_keeps_phenol_and_records_the_side():
    counters = PathCounters()
    hits = list(
        find_path(
            "COc1ccccc1",
            "Oc1ccccc1",
            ruleset=PhaseOne,
            counters=counters,
            max_nodes=40,
            max_paths=1,
        )
    )
    assert hits
    outcome = hits[0]
    assert "." not in outcome.smiles
    assert outcome.smiles == Chem.MolFromSmiles("Oc1ccccc1").xf.csmi
    assert outcome.maybe
    side = outcome.maybe.sides()[0]
    assert side != outcome.smiles
    assert "." not in side
    # The discarded fragment is not another expanded node.
    assert counters.nodes == 2


@given(start=st.sampled_from(_ETHERS), data=st.data())
@settings(
    max_examples=4,
    deadline=20_000,
    derandomize=True,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)
def test_fuzz_dealkylation_hit_or_honest_miss(start: str, data):
    mol = Chem.MolFromSmiles(start)
    assume(mol is not None)
    products = []
    for product, info in Dealkylation().metabolize(mol):
        smi = product.xf.csmi
        if not smi or "." in smi or smi == mol.xf.csmi:
            continue
        if not info["options"].get("cleaves"):
            continue
        products.append(smi)
        if len(products) >= 6:
            break
    assume(products)
    target = data.draw(st.sampled_from(products))
    hits = list(
        find_path(start, target, ruleset=PhaseOne, max_nodes=60, max_paths=1)
    )
    if not hits:
        return
    assert hits[0].smiles == target
    assert "." not in hits[0].smiles
    for side in hits[0].maybe.sides():
        assert side != target
        assert "." not in side
