"""Degenerate Phase1-equivalent steps on Phase I / NDealkylation rules."""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from xenosite.forest import AtomRef, StepPlan
from xenosite.forest.base import ConjugatedSystems
from xenosite.forest.rules import (
    Acetylation,
    Dehydrogenation,
    Epoxidation,
    Hydroxylation,
    NDealkylation,
    QuinoneFormation,
)
from xenosite.forest.utils import canon_smi, refresh_mol, unmapped_smiles

from phase1_helpers import quinone_plan_reaches_product

# Same multi-system fixtures as resonance-cache tests (2 conjugated systems),
# plus an explicit 3-system alkyl-linked triphenyl.
NAPH_STYRYL = "C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C"
THREE_CONJ = "c1ccccc1CCCCc2ccccc2CCCCc3ccccc3"

# Small probes that collectively hit every QuinoneFormation query modification
# and every Dehydrogenation query / reaction SMARTS.
QUINONE_SMARTS_PROBES = (
    ("phenol", "Oc1ccccc1"),
    ("benzene", "c1ccccc1"),
    ("chlorobenzene", "Clc1ccccc1"),
    ("dimethylaniline", "CN(C)c1ccccc1"),
    ("anisole", "COc1ccccc1"),
    ("apap", "CC(=O)Nc1ccc(O)cc1"),
)

# Extra probes for Dehydrogenation reaction SMARTS not covered by quinone set.
# CS(O)C → [#16v4]-[Oh]; CCO → carbon–hetero DH.
DEHYDROGENATION_SMARTS_PROBES = (
    ("ethanol", "CCO"),
    ("dimethyl_sulfanol", "CS(O)C"),
)

SMARTS_COVERAGE_PROBES = QUINONE_SMARTS_PROBES + DEHYDROGENATION_SMARTS_PROBES


def _mod_names(modifications) -> frozenset[str]:
    if isinstance(modifications, str):
        return frozenset([modifications])
    return frozenset(modifications)


def _hit_query_modifications(rule, mol) -> frozenset[str]:
    """Modification names from ``rule.query_smarts`` that match ``mol``."""
    template = rule.standardize(mol)
    if template is None:
        return frozenset()
    hit: set[str] = set()
    for matches in rule.match_queries(template).values():
        for _maps, mods in matches:
            hit |= _mod_names(mods)
    return frozenset(hit)


def _hit_reaction_smarts(rule, mol) -> frozenset[str]:
    """Reaction SMARTS strings on ``rule`` for which ``RunReactants`` succeeds."""
    hit: set[str] = set()
    for smarts, rxn in zip(rule.smarts, rule.rxns):
        try:
            if rxn.RunReactants((mol,)):
                hit.add(smarts)
        except Exception:
            continue
    return frozenset(hit)


def test_phase1_smarts_collectively_hit():
    """Curated probes hit every QuinoneFormation and Dehydrogenation SMARTS."""
    qf = QuinoneFormation()
    dh = Dehydrogenation()
    qf_expected = frozenset(qf.valid_modifications)
    dh_query_expected = frozenset(dh.valid_modifications)
    dh_rxn_expected = frozenset(dh.smarts)
    assert qf_expected and dh_query_expected and dh_rxn_expected

    qf_hit: set[str] = set()
    dh_query_hit: set[str] = set()
    dh_rxn_hit: set[str] = set()
    for _id, smi in SMARTS_COVERAGE_PROBES:
        mol = MolFromSmiles(smi)
        qf_hit |= _hit_query_modifications(qf, mol)
        dh_query_hit |= _hit_query_modifications(dh, mol)
        dh_rxn_hit |= _hit_reaction_smarts(dh, mol)

    assert qf_expected <= qf_hit, "unhit quinone mods: %s" % sorted(
        qf_expected - qf_hit
    )
    assert dh_query_expected <= dh_query_hit, "unhit DH query mods: %s" % sorted(
        dh_query_expected - dh_query_hit
    )
    assert dh_rxn_expected <= dh_rxn_hit, "unhit DH reaction SMARTS: %s" % sorted(
        dh_rxn_expected - dh_rxn_hit
    )


@pytest.mark.parametrize(
    "probe_id,smiles",
    QUINONE_SMARTS_PROBES,
    ids=[p[0] for p in QUINONE_SMARTS_PROBES],
)
def test_quinone_smarts_probe_phase1_steps(probe_id, smiles):
    """Each quinone SMARTS probe: stamp matches public plan and reaches product."""
    mol = MolFromSmiles(smiles)
    assert _hit_query_modifications(QuinoneFormation(), mol)

    qf = QuinoneFormation()
    checked = 0
    for site, products in qf.metabolize(
        mol, attach_phase1_steps=True, tag_atoms=False
    ):
        for product in products:
            plan = StepPlan.try_from_mol(product)
            if plan is None:
                continue
            public = qf.phase1_steps(mol, site[1])
            assert plan in public
            quinone_plan_reaches_product(mol, plan, unmapped_smiles(product))
            checked += 1
            if checked >= 4:
                return
    assert checked, probe_id


@pytest.mark.parametrize(
    "probe_id,smiles",
    DEHYDROGENATION_SMARTS_PROBES,
    ids=[p[0] for p in DEHYDROGENATION_SMARTS_PROBES],
)
def test_dehydrogenation_smarts_probe_phase1_steps(probe_id, smiles):
    """DH coverage probes stamp degenerate phase1_steps and apply cleanly."""
    mol = MolFromSmiles(smiles)
    dh = Dehydrogenation()
    assert _hit_query_modifications(dh, mol) or _hit_reaction_smarts(dh, mol)

    site, products = next(
        dh.metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    plan = StepPlan.from_mol(products[0])
    public = dh.phase1_steps(mol, site[1])
    assert public and plan == public[0]
    assert next(plan.linearizations()).apply(mol)


def test_dehydrogenation_quinoid_smarts_hydroquinone_and_apap():
    """New #6H0–OH/NH query SMARTS enable Forest DH quinone products."""
    dh = Dehydrogenation()
    hq = MolFromSmiles("Oc1ccc(O)cc1")
    products = [
        unmapped_smiles(p)
        for _s, ps in dh.metabolize(hq, tag_atoms=False)
        for p in ps
    ]
    assert any(canon_smi(s) == canon_smi("O=C1C=CC(=O)C=C1") for s in products)

    apap = MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    products = [
        unmapped_smiles(p)
        for _s, ps in dh.metabolize(apap, tag_atoms=False)
        for p in ps
    ]
    assert any(canon_smi(s) == canon_smi("CC(=O)N=C1C=CC(=O)C=C1") for s in products)


def test_quinone_stepplan_full_forest_apply_apap():
    """With quinoid DH SMARTS, APAP DH-only StepPlan applies via Forest."""
    mol = MolFromSmiles("CC(=O)Nc1ccc(O)cc1")
    plan = next(
        p
        for p in QuinoneFormation().phase1_steps(mol, frozenset({4, 7}))
        if len(p) == 1
    )
    out = next(plan.linearizations()).apply(mol)
    assert out
    assert canon_smi(unmapped_smiles(out[0])) == canon_smi("CC(=O)N=C1C=CC(=O)C=C1")


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
    stamped = [
        plan for p in products if (plan := StepPlan.try_from_mol(p)) is not None
    ]
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
    stamped = [p for p in products if StepPlan.try_from_mol(p) is not None]
    assert stamped
    assert StepPlan.from_mol(stamped[0]) in qf.phase1_steps(mol, frozenset({4, 7}))


def test_quinone_not_yet_phase1_equivalent_flag():
    # QuinoneFormation overrides phase1_steps; flag stays False so metabolize
    # does not try the degenerate singleton attach path.
    assert QuinoneFormation.phase1_equivalent is False


@pytest.mark.parametrize(
    "smiles,n_conj",
    [
        (NAPH_STYRYL, 2),
        (THREE_CONJ, 3),
    ],
    ids=["two_conjugated_systems", "three_conjugated_systems"],
)
def test_phase1_steps_multi_conjugated_systems(smiles, n_conj):
    """phase1_steps / attach work on molecules with 2 and 3 conjugated systems."""
    mol = MolFromSmiles(smiles)
    assert len(ConjugatedSystems().systems(mol)) == n_conj

    qf = QuinoneFormation()
    site, products = next(
        qf.metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    stamped = [
        (p, plan)
        for p in products
        if (plan := StepPlan.try_from_mol(p)) is not None
    ]
    assert stamped
    product, plan = stamped[0]
    public = qf.phase1_steps(mol, site[1])
    assert public
    assert plan in public
    quinone_plan_reaches_product(mol, plan, unmapped_smiles(product))
