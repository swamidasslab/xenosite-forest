"""Focused RDKit valence crashers (not the whole forest valence suite).

Port of the crashers in ``tests/test_rdkit_valence.py``: diphenhydramine DH,
issue-3 PhaseOne bfs finish, and continue-after-failed RunReactants.
"""

from __future__ import annotations

from xenosite.forest.find_path import bfs
from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rdkitutil import run_reactants
from xenosite.forest.rules import Dehydrogenation, SmirksReactionRule, _describe
from xenosite.forest.rulesets import PhaseOne

DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c1ccccc1"
IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
ETHANE = "CC"
# GitHub issue #3: fused dihydrobenzofuran amide.
ISSUE3_PARENT = "CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21"


def test_dehydrogenation_diphenhydramine_does_not_crash():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    n = sum(1 for _ in Dehydrogenation().metabolize(mol))
    assert n > 0


def test_bfs_phaseone_issue3_parent_does_not_crash():
    """Issue #3: bfs(PhaseOne, depth=1) must finish past the first DH sites."""

    rows = list(bfs(ISSUE3_PARENT, ruleset=PhaseOne, depth=1))
    assert len(rows) > 7
    names = {info["rule"][0].name for _mol, info in rows}
    assert "Dehydrogenation" in names
    assert any(name != "Dehydrogenation" for name in names)


def test_phaseone_metabolize_on_suite_crashers():
    for smi in (DIPHENHYDRAMINE, IBUPROFEN, ISSUE3_PARENT, ETHANE):
        rows = list(PhaseOne.metabolize(MolFromSmiles(smi)))
        assert rows, smi


def test_metabolites_continues_after_runreactants_runtime_error(monkeypatch):
    """A failed rxn must not abort later SMIRKS on the same reactant."""

    class TwoRxns(SmirksReactionRule):
        name = "TwoRxns"
        smirks = (
            ("[#6H3:1]>>[*:1]F", _describe(adds="F", removes="H", name="f")),
            ("[#6H3:1]>>[*:1]O", _describe(adds="O", removes="H", name="o")),
        )

    rule = TwoRxns()
    real = run_reactants

    def boom(smarts: str, mol):
        # react_at uses isotope-pinned SMARTS; match the fluorine product arm.
        if "F" in smarts.split(">>", 1)[-1]:
            raise RuntimeError("simulated valence precondition")
        return real(smarts, mol)

    monkeypatch.setattr("xenosite.forest.rules.run_reactants", boom)
    rows = list(rule.metabolize(MolFromSmiles(ETHANE)))
    assert rows, "second SMIRKS should still yield after first RunReactants fails"
    assert {p.xf.csmi for _pl, _ in rows for p in _pl} == {"CCO"}
