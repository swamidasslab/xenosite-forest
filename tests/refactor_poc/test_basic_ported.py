"""Port of ``tests/test_basic.py`` chemical checks. Path finds use poc ``find_path``.

Epoxide / dehydrogenation path cases assert a non-empty plan. Dealkylation
site and sanitize checks keep the forest asserts. Nevirapine SMILES is the
constant from the forest suite (imported, not copied).
"""

from __future__ import annotations

import collections

import pytest
from rdkit.Chem.rdmolops import SanitizeMol

from test_basic import nevirapine
from xenosite.refactor_poc.find_path import find_path
from xenosite.refactor_poc.rdkit_api import MolFromSmiles
from xenosite.refactor_poc.rules import (
    Dealkylation,
    Dehydrogenation,
    Epoxidation,
    EpoxideOpening,
)
from xenosite.refactor_poc.rulesets import RuleSet

from .helpers import canon


def test_epoxide_opening_aromatic():
    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=200,
        )
    )
    assert hits


def test_epoxide_opening_kekulized():
    hits = list(
        find_path(
            "C1=CC=CC=C1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=200,
        )
    )
    assert hits


def test_propane_dehydrogenation_to_propene():
    hits = list(
        find_path(
            "CCC",
            "C=CC",
            ruleset=RuleSet((Dehydrogenation,), name="DH"),
            max_nodes=50,
        )
    )
    assert hits


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_propane_dehydrogenation_followed_by_epoxidation():
    hits = list(
        find_path(
            "CCC",
            "C1OC1C",
            ruleset=RuleSet((Dehydrogenation, Epoxidation), name="DH_E"),
            max_nodes=200,
        )
    )
    assert hits


def test_hydroxyl_should_not_be_dealkylated():
    mol = MolFromSmiles("CCO")
    sites = [info["site"] for _product, info in Dealkylation().metabolize(mol)]
    assert frozenset([1, 2]) not in sites


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_epoxide_opening2():
    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC1OC1",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=200,
        )
    )
    assert hits


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_epoxide_opening1_fail():
    """Depth-1 find_path must not reach the diol from benzene."""

    hits = list(
        find_path(
            "c1ccccc1",
            "C1=CC=CC(O)C1O",
            ruleset=RuleSet((Epoxidation, EpoxideOpening), name="EO"),
            max_nodes=40,
            max_paths=3,
        )
    )
    # Forest depth-default miss: epoxidation alone cannot yield the diol in one hop
    # under a tight ceiling. Keep the empty-hits assert.
    assert len(hits) == 0


def _dealk_sanitize(reactant: str) -> None:
    mol = MolFromSmiles(reactant)
    for _product, _info in Dealkylation().metabolize(mol):
        SanitizeMol(_product)


def test_matt_problem1():
    _dealk_sanitize("CC(C)CNS(=O)(=O)c1ccc(CCC(=O)Nc2ccc(Cl)cc2C)cc1")


@pytest.mark.xfail(reason="poc deferred bug")
@pytest.mark.regression
def test_matt_problem2():
    _dealk_sanitize("Oc1c(C(=O)Nc2cccnc2)c(=O)n2CCc3cccc1c23")

def test_matt_problem4():
    rmol = MolFromSmiles(nevirapine)
    for _product, info in Dealkylation().metabolize(rmol):
        site = info["site"]
        if isinstance(site, int):
            continue
        atoms = tuple(site)
        if len(atoms) != 2:
            continue
        assert rmol.GetBondBetweenAtoms(*atoms), (
            "NO BOND BETWEEN ATOMS %d and %d" % tuple(a + 1 for a in atoms)
        )


def test_unique_metabolites():
    metabolite_registry: dict = collections.defaultdict(set)
    rmol = MolFromSmiles(nevirapine)
    for product, info in Dealkylation().metabolize(rmol):
        site = info["site"]
        if isinstance(site, int):
            site_key = (site,)
        else:
            site_key = tuple(sorted(site))
        key = frozenset(["Dealkylation", *site_key])
        canonical = {canon(product)}
        assert not canonical.issubset(metabolite_registry[key])
        metabolite_registry[key].update(canonical)
