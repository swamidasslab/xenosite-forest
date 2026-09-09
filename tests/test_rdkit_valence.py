"""RDKit 2026 valence-cache regressions (resonance copies / RunReactants)."""

from __future__ import annotations

import warnings

import pytest
from rdkit import Chem

from xenosite.forest import PhaseOneRS
from xenosite.forest.base import SmartsReactionRule
from xenosite.forest.rules import Dehydrogenation
from xenosite.forest.utils import refresh_mol

# From the 327-molecule suite: crashed on RDKit 2026 until valence caches were filled.
DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"
IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
ETHANE = "CC"


def test_dehydrogenation_diphenhydramine_does_not_crash():
    mol = Chem.MolFromSmiles(DIPHENHYDRAMINE)
    n = sum(1 for _ in Dehydrogenation().metabolize(mol))
    assert n > 0


@pytest.mark.parametrize("smi", [DIPHENHYDRAMINE, IBUPROFEN])
def test_phaseone_unique_on_suite_crashers(smi):
    mol = Chem.MolFromSmiles(smi)
    rows = list(PhaseOneRS.metabolites(mol, unique=True))
    assert rows, smi


def test_phaseone_ethane_still_enumerates():
    mol = Chem.MolFromSmiles(ETHANE)
    rows = list(PhaseOneRS.metabolites(mol, unique=True))
    assert rows


def test_refresh_mol_allows_runreactants_on_resonance_copy():
    mol = Chem.MolFromSmiles(DIPHENHYDRAMINE)
    rule = Dehydrogenation()
    res = list(rule.resonance_structures(mol))
    assert len(res) > 1
    copy = res[1]
    refresh_mol(copy)
    prods = rule.rxns[1].RunReactants([copy])
    assert len(prods) >= 1


def test_metabolites_continues_after_runreactants_runtime_error():
    """A failed rxn must not abort later SMARTS on the same reactant."""

    class TwoRxns(SmartsReactionRule):
        name = "TwoRxns"
        smarts = [
            "[#6H3:1]>>[*:1]F",
            "[#6H3:1]>>[*:1]O",
        ]

    rule = TwoRxns()
    assert len(rule.rxns) == 2

    real0 = rule.rxns[0].RunReactants

    def boom(reactants):
        raise RuntimeError("simulated valence precondition")

    rule.rxns[0].RunReactants = boom

    mol = Chem.MolFromSmiles(ETHANE)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rows = list(rule.metabolites(mol))

    assert any("Skipping TwoRxns rxn 0" in str(w.message) for w in caught)
    assert rows, "second SMARTS should still yield after first RunReactants fails"
    # Restore in case the instance is reused (defensive).
    rule.rxns[0].RunReactants = real0
