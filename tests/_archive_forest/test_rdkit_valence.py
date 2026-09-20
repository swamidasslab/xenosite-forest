"""RDKit 2026 valence-cache regressions (resonance copies / RunReactants)."""

from __future__ import annotations

import logging

import pytest
from rdkit import Chem

from xenosite._archive_forest import PhaseOneRS, bfs
from xenosite._archive_forest.base import SmartsReactionRule
from xenosite._archive_forest.rules import Dehydrogenation
from xenosite._archive_forest.utils import refresh_mol

# From the 327-molecule suite: crashed on RDKit 2026 until valence caches were filled.
DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"
IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
ETHANE = "CC"
# GitHub issue #3: fused dihydrobenzofuran amide; forest 0.1.0 crashed after the
# first seven Dehydrogenation products on RDKit 2026.03.1.
ISSUE3_PARENT = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"


def test_dehydrogenation_diphenhydramine_does_not_crash():
    mol = Chem.MolFromSmiles(DIPHENHYDRAMINE)
    n = sum(1 for _ in Dehydrogenation().metabolize(mol))
    assert n > 0


def test_bfs_phaseone_issue3_parent_does_not_crash():
    """Issue #3: bfs(PhaseOneRS, depth=1) must finish, not stop after DH sites."""
    rows = list(bfs(ISSUE3_PARENT, ruleset="PhaseOneRS", depth=1))
    assert len(rows) > 7
    assert any(steps and steps[0][0] == "Dehydrogenation" for _, steps, _ in rows)
    assert any(steps and steps[0][0] != "Dehydrogenation" for _, steps, _ in rows)


@pytest.mark.parametrize("smi", [DIPHENHYDRAMINE, IBUPROFEN, ISSUE3_PARENT])
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


def test_metabolites_continues_after_runreactants_runtime_error(caplog):
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
    caplog.set_level(logging.DEBUG)
    rows = list(rule.metabolites(mol))

    assert any("Skipping TwoRxns rxn 0" in r.message for r in caplog.records)
    assert rows, "second SMARTS should still yield after first RunReactants fails"
    # Restore in case the instance is reused (defensive).
    rule.rxns[0].RunReactants = real0
