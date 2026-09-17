"""Degenerate Phase1-equivalent steps on Phase I / NDealkylation rules."""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest import AtomRef, StepPlan
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


def test_quinone_apap_phase1_steps_single_dh():
    mol = MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    plans = QuinoneFormation().phase1_steps(mol, frozenset({4, 7}))
    assert plans
    # Prefer the single2double+single2double plan: only Dehydrogenation on N and O.
    dh_only = [
        p
        for p in plans
        if len(p) == 1 and p.steps[0].rule == "Dehydrogenation"
    ]
    assert dh_only
    assert dh_only[0].steps[0].site == frozenset(
        {AtomRef(origin=3), AtomRef(origin=8)}
    )


def test_quinone_benzene_addo_layers():
    mol = MolFromSmiles("c1ccccc1")
    # para carbons 0 and 3
    plans = QuinoneFormation().phase1_steps(mol, frozenset({0, 3}))
    assert plans
    layered = [p for p in plans if len(p) == 3]
    assert layered
    plan = layered[0]
    orders = list(plan.iter_linearizations())
    assert len(orders) == 2
    assert all(o[-1].rule == "Dehydrogenation" for o in orders)
    assert all(o[0].rule == "Hydroxylation" for o in orders)
    dh_site = orders[0][-1].site
    assert all(isinstance(r, AtomRef) and r.added_by == "Hydroxylation" for r in dh_site)


def test_quinone_attach_phase1_steps_matches_public():
    mol = MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    qf = QuinoneFormation()
    site, products = next(
        qf.metabolites_from_sites(
            mol, frozenset({4, 7}), attach_phase1_steps=True, tag_atoms=False
        )
    )
    stamped = [StepPlan.from_mol(p) for p in products if p.HasProp("phase1_steps")]
    assert stamped
    public = qf.phase1_steps(mol, frozenset({4, 7}))
    assert stamped[0] in public


def test_quinone_attach_phase1_steps_survives_tagging():
    """RenumberAtoms during tag/align must keep phase1_steps on products."""
    mol = MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    qf = QuinoneFormation()
    site, products = next(
        qf.metabolites_from_sites(
            mol, frozenset({4, 7}), attach_phase1_steps=True, tag_atoms=True
        )
    )
    stamped = [p for p in products if p.HasProp("phase1_steps")]
    assert stamped
    assert StepPlan.from_mol(stamped[0]) in qf.phase1_steps(mol, frozenset({4, 7}))


def test_quinone_not_yet_phase1_equivalent_flag():
    # QuinoneFormation overrides phase1_steps; flag stays False so metabolize
    # does not try the degenerate singleton attach path.
    assert QuinoneFormation.phase1_equivalent is False
