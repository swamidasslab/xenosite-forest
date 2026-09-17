"""Degenerate Phase1-equivalent steps on Phase I / NDealkylation rules."""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest import StepPlan
from xenosite.forest.rules import (
    Acetylation,
    Epoxidation,
    Hydroxylation,
    NDealkylation,
    QuinoneFormation,
)
from xenosite.forest.utils import refresh_mol
from rdkit.Chem.rdmolfiles import MolToSmiles


def test_acetylation_phase1_steps_not_implemented():
    with pytest.raises(NotImplementedError):
        Acetylation().phase1_steps(MolFromSmiles("CCO"), frozenset({0}))


def test_epoxidation_degenerate_phase1_steps():
    mol = MolFromSmiles("C=C")
    plans = Epoxidation().phase1_steps(mol, frozenset({0, 1}))
    assert len(plans) == 1
    assert plans[0] == StepPlan.singleton("Epoxidation", frozenset({0, 1}))
    assert list(plans[0].iter_linearizations()) == [plans[0].steps]


def test_epoxidation_attach_phase1_steps():
    mol = MolFromSmiles("C=C")
    site, products = next(
        Epoxidation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    assert site[0] == "Epoxidation"
    plan = StepPlan.from_mol(products[0])
    assert plan == StepPlan.singleton("Epoxidation", frozenset(site[1]))


def test_ndealkylation_degenerate_phase1_steps():
    mol = MolFromSmiles("CCN")
    # Find an emitted N-dealkylation site, then ask phase1_steps for it.
    site, products = next(
        NDealkylation().metabolize(mol, tag_atoms=False)
    )
    plans = NDealkylation().phase1_steps(mol, site[1])
    assert plans == [StepPlan.singleton("NDealkylation", frozenset(site[1]))]


def test_ndealkylation_attach_matches_public():
    mol = MolFromSmiles("CCN")
    site, products = next(
        NDealkylation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    stamped = StepPlan.from_mol(products[0])
    public = NDealkylation().phase1_steps(mol, site[1])
    assert public and stamped == public[0]


def test_hydroxylation_degenerate():
    mol = MolFromSmiles("CC")
    site, _ = next(Hydroxylation().metabolize(mol, tag_atoms=False))
    plans = Hydroxylation().phase1_steps(mol, site[1])
    assert plans == [StepPlan.singleton("Hydroxylation", frozenset(site[1]))]


def test_quinone_not_yet_phase1_equivalent_flag():
    # QuinoneFormation overrides phase1_steps in a later commit; until then
    # phase1_equivalent stays False so base raises unless overridden.
    assert QuinoneFormation.phase1_equivalent is False
