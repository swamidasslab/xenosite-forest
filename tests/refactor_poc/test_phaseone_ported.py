"""Port of ``tests/test_phaseone.py`` pathway cases onto poc ``find_path``.

Each case asks for the same reactant→product chemistry. Plan site sets stay
asserted when the forest suite named them; mark deferred bugs with
``xfail`` + ``regression``. Histidine validity uses poc ``PhaseOne`` /
``NitrogenOxidation``. Forest ``clean`` / ``bfs`` APIs and rule ``rxns=``
surgery are skipped (unfinished / different surface).
"""

from __future__ import annotations

import pytest

from xenosite.refactor_poc.find_path import find_path
from xenosite.refactor_poc.rdkit_api import MolFromSmiles, MolToSmiles
from xenosite.refactor_poc.rdkitutil import canon_smiles
from xenosite.refactor_poc.rules import NitrogenOxidation
from xenosite.refactor_poc.rulesets import PhaseOne

from .helpers import canon

HISTIDINE = "NC(Cc1cnc[nH]1)C(=O)O"


def _plan_rules(outcome) -> list[str]:
    return [step.rule for step in outcome.plan.children]


def _find(reactant: str, product: str, *, max_nodes: int = 200):
    return list(find_path(reactant, product, max_nodes=max_nodes, max_paths=3))


def test_dehydrogenation_path():
    hits = _find("CCO", "C=CO")
    assert hits
    assert hits[0].smiles == canon_smiles("C=CO")
    assert _plan_rules(hits[0]) == ["Dehydrogenation"]


def test_dealkylation1():
    hits = _find("CCN", "CCO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCO")
    assert _plan_rules(hits[0]) == ["Dealkylation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_dehydration1():
    hits = _find("CCO", "CC")
    assert hits
    assert hits[0].smiles == canon_smiles("CC")
    assert _plan_rules(hits[0]) == ["Dehydration"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_dephosphorylation1():
    hits = _find("COP(=O)(O)O", "CO")
    assert hits
    assert hits[0].smiles == canon_smiles("CO")
    assert _plan_rules(hits[0]) == ["Dephosphorylation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_epoxidation1():
    hits = _find("C=C", "C1OC1")
    assert hits
    assert hits[0].smiles == canon_smiles("C1CO1")
    assert _plan_rules(hits[0]) == ["Epoxidation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_epoxidation_opening1():
    hits = _find("C1OC1", "CCO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCO")
    assert _plan_rules(hits[0]) == ["EpoxideOpening"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_hydrogenation1():
    hits = _find("CC=O", "CCO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCO")
    assert _plan_rules(hits[0]) == ["Hydrogenation"]


def test_hydroxylation1():
    hits = _find("CC", "CCO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCO")
    assert _plan_rules(hits[0]) == ["Hydroxylation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_hydrolysis1():
    hits = _find("O=C(O)C", "CC=O")
    assert hits
    assert hits[0].smiles == canon_smiles("CC=O")
    assert _plan_rules(hits[0]) == ["Hydrolysis"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_nitrogen_oxidation1():
    hits = _find("CCN", "CCNO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCNO")
    assert _plan_rules(hits[0]) == ["NitrogenOxidation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_dehydration2():
    hits = _find("CCNO", "CCN")
    assert hits
    assert hits[0].smiles == canon_smiles("CCN")
    assert _plan_rules(hits[0]) == ["Dehydration"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_nitrogen_reduction1():
    reactant = "[O-]-[N+](C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O"
    product = "N(C1=CC=C(O1)C=NN2C(=O)NC(=O)C2)=O"
    hits = _find(reactant, product, max_nodes=400)
    assert hits
    assert hits[0].smiles == canon_smiles("O=Nc1ccc(C=NN2CC(=O)NC2=O)o1")
    assert _plan_rules(hits[0]) == ["NitrogenReduction"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_oxidative_dehalogenation():
    hits = _find("CCCl", "CCO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCO")
    assert _plan_rules(hits[0]) == ["OxidativeDehalogenation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_reductive_dehalogenation():
    hits = _find("CCCl", "CC")
    assert hits
    assert hits[0].smiles == canon_smiles("CC")
    assert _plan_rules(hits[0]) == ["ReductiveDehalogenation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_sulfur_oxidation():
    hits = _find("CCS", "CCSO")
    assert hits
    assert hits[0].smiles == canon_smiles("CCSO")
    assert _plan_rules(hits[0]) == ["SulfurOxidation"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_sulfur_reduction():
    hits = _find("CCSO", "CCS")
    assert hits
    assert hits[0].smiles == canon_smiles("CCS")
    assert _plan_rules(hits[0]) == ["SulfurReduction"]


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_histidine_phase1_does_not_emit_invalid_metabolites():
    mol = MolFromSmiles(HISTIDINE)
    assert mol is not None
    products = list(PhaseOne.metabolize(mol))
    assert products
    invalid = [
        (type(info["rule"]).__name__, info["site"], MolToSmiles(product))
        for product, info in products
        if MolFromSmiles(MolToSmiles(product)) is None
    ]
    assert invalid == []


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_histidine_nitrogen_oxidation_is_chemically_valid():
    mol = MolFromSmiles(HISTIDINE)
    assert mol is not None
    smiles = {canon(product) for product, _info in NitrogenOxidation().metabolize(mol)}

    assert canon("O=C(O)C(CC1=CN=CN1)NO") in smiles
    assert canon("O=NC(CC1=CN=CN1)C(=O)O") in smiles
    assert canon("NC(CC1=CN=CN1O)C(=O)O") in smiles
    assert canon("NC(CC1=C[N+]([O-])=CN1)C(=O)O") in smiles
    assert canon("NC(CC1=CN(O)=CN1)C(=O)O") not in smiles
    assert canon("NC(CC1=CN=CN1=O)C(=O)O") not in smiles
    assert all(MolFromSmiles(s) is not None for s in smiles)


# Skipped: test_rule_modification (rxns= surgery), forest clean() cases,
# duplicate dehydrogenation1/hydrogenation2 (same chemistry as above).
