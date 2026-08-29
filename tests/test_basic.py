"""
Tests adapted from:
xenosite.metabolite.tests.Debug
xenosite.metabolite.tests.Debug2
"""

from rdkit.Chem.rdmolfiles import (
    MolFromSmiles,
    MolFromMolBlock,
    # MolToSmiles,
)
from rdkit.Chem.rdmolops import SanitizeMol
from xenosite.metabolite import rules
from xenosite.metabolite.rulesets import RuleSet
from xenosite.metabolite.base import can_smi
import collections
import pytest


def test_epoxide_opening_aromatic():
    reactant = MolFromSmiles("c1ccccc1")
    product = MolFromSmiles("C1=CC=CC(O)C1O")

    RS = RuleSet(rules=[rules.Epoxidation(), rules.EpoxideOpening()])

    path = next(RS.find_path(reactant, product, depth=2))
    assert len(path) > 0


def test_epoxide_opening_kekulized():
    reactant = MolFromSmiles("C1=CC=CC=C1")
    product = MolFromSmiles("C1=CC=CC(O)C1O")

    RS = RuleSet(rules=[rules.Epoxidation(), rules.EpoxideOpening()])

    path = next(RS.find_path(reactant, product, depth=2))
    assert len(path) > 0


def test_propane_dehydrogenation_to_propene():
    reactant = MolFromSmiles("CCC")
    product = MolFromSmiles("C=CC")

    RS = RuleSet(rules=[rules.Dehydrogenation()])

    path = next(RS.find_path(reactant, product))

    assert len(path) > 0


def test_propane_dehydrogenation_followed_by_epoxidation():
    reactant = MolFromSmiles("CCC")
    product = MolFromSmiles("C1OC1C")

    RS = RuleSet(rules=[rules.Dehydrogenation(), rules.Epoxidation()])

    path = next(RS.find_path(reactant, product, depth=2))
    assert len(path) > 0


def test_hydroxyl_should_not_be_dealkylated():
    mol = MolFromSmiles("CCO")
    R = rules.Dealkylation()
    sites = []

    for (rulename, site), metabolites in R.metabolites(mol):
        sites.append(site)

    assert frozenset([1, 2]) not in sites


def test_epoxide_opening2():
    RS = RuleSet(rules=[rules.Epoxidation(), rules.EpoxideOpening()])
    path = next(
        RS.find_path(MolFromSmiles("c1ccccc1"), end_mol=MolFromSmiles("C1=CC=CC1OC1"))
    )
    assert len(path) > 0


def test_epoxide_opening1_fail():
    reactant = MolFromSmiles("c1ccccc1")
    product = MolFromSmiles("C1=CC=CC(O)C1O")

    RS = RuleSet(rules=[rules.Epoxidation(), rules.EpoxideOpening()])
    path = list(RS.find_path(reactant, product))
    assert len(path) == 0


def dealk_tests(reactant):
    """Ensures all dealk metabolites are valid molecules"""
    mol = MolFromSmiles(reactant)
    R = rules.Dealkylation()
    for site, metabolites in R.metabolites(mol):
        for m in metabolites:
            SanitizeMol(m)


def test_matt_problem1():
    dealk_tests("CC(C)CNS(=O)(=O)c1ccc(CCC(=O)Nc2ccc(Cl)cc2C)cc1")


def test_matt_problem2():
    dealk_tests("Oc1c(C(=O)Nc2cccnc2)c(=O)n2CCc3cccc1c23")


@pytest.mark.xfail(reason="needs working ndealk model")
def test_matt_problem4(self):
    rmol = MolFromMolBlock(self.nevirapine)

    mg = RuleSet(rules=[rules.Dealkylation()])

    for (model, site), metabolites in mg.metabolites(rmol):
        for m in metabolites:
            assert rmol.GetBondBetweenAtoms(
                *site
            ), "NO BOND BETWEEN ATOMS %d and %d" % tuple([x + 1 for x in site])



def test_unique_metabolites():
    metabolite_registry = collections.defaultdict(set)

    rmol = MolFromSmiles(nevirapine)

    mg = RuleSet(rules=[rules.Dealkylation()])

    for (model, site), metabolites in mg.metabolites(rmol, unique=True):
        for m in metabolites:

            site = sorted([a  for a in site])

            key = frozenset([model] + site)
            canonical_metabolites = set(can_smi(rdmol=m))
            print(canonical_metabolites, metabolite_registry[key], key)

            assert not canonical_metabolites.issubset(metabolite_registry[key])

            metabolite_registry[key].update(canonical_metabolites)


nevirapine = "C1=CC2=C(C=NNC3=C2N=CC=C3)N=C1"