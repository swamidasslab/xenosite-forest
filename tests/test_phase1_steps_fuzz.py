"""Hypothesis fuzz: StepPlan linearizations via library apply API.

Example database lives in ``.hypothesis/`` (gitignored) and is restored on CI.
"""

from __future__ import annotations

from pathlib import Path

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from hypothesis.database import DirectoryBasedExampleDatabase
from rdkit import Chem

from xenosite.forest import StepPlan
from xenosite.forest.base import can_smi_set
from xenosite.forest.rules import Epoxidation, NDealkylation, QuinoneFormation
from xenosite.forest.utils import unmapped_smiles

_HYPOTHESIS_DIR = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
_HYPOTHESIS_DB = DirectoryBasedExampleDatabase(str(_HYPOTHESIS_DIR))

_CORPUS = (
    "CC(=O)Nc1ccc(O)cc1",
    "c1ccccc1",
    "Oc1ccc(O)cc1",
    "Oc1ccccc1",
    "Nc1ccc(O)cc1",
    "Clc1ccc(O)cc1",
    "COc1ccc(O)cc1",
    "Cc1ccc(O)cc1",
    "c1ccc2ccccc2c1",
    "C=C",
    "CCN",
    "C=Cc1ccccc1",
)

_MAX_HEAVY = 28
_MAX_STAMPED = 12


def _smi_sets(product_lists):
    return [
        frozenset(unmapped_smiles(m) for m in products)
        for products in product_lists
        if products
    ]


@st.composite
def corpus_smiles(draw):
    return draw(st.sampled_from(_CORPUS))


def test_benzene_prep_linearizations_agree():
    """Both OH orders for para-quinone prep yield the same hydroquinone."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    plan = next(
        p for p in QuinoneFormation().phase1_steps(mol, frozenset({0, 3})) if len(p) == 3
    )
    finals = _smi_sets(
        lin.apply(mol, drop_last=1) for lin in plan.iter_as_linearizations()
    )
    assert finals and len(set(finals)) == 1
    got = Chem.MolToSmiles(Chem.MolFromSmiles(next(iter(finals[0]))))
    assert got == Chem.MolToSmiles(Chem.MolFromSmiles("Oc1ccc(O)cc1"))


def test_epoxidation_linearization_replays_product():
    mol = Chem.MolFromSmiles("C=C")
    _, products = next(
        Epoxidation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    plan = StepPlan.from_mol(products[0])
    expected = can_smi_set(products)
    for lin in plan.iter_as_linearizations():
        result = lin.apply(mol)
        assert result
        assert can_smi_set(result) == expected


def test_ndealkylation_linearization_replays_product():
    mol = Chem.MolFromSmiles("CCN")
    _, products = next(
        NDealkylation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    plan = StepPlan.from_mol(products[0])
    expected = can_smi_set(products)
    assert any(
        expected == got or expected <= got or got <= expected
        for got in (
            can_smi_set(r)
            for lin in plan.iter_as_linearizations()
            for r in [lin.apply(mol)]
            if r
        )
    )


@given(smiles=corpus_smiles())
@settings(
    max_examples=40,
    deadline=30_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_quinone_linearizations_prep_agree_or_full_match(smiles: str):
    """Quinone plans: prep orders agree; full apply matches when all steps fire."""
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)
    assume(mol.GetNumHeavyAtoms() <= _MAX_HEAVY)

    qf = QuinoneFormation()
    stamped = 0
    for _site, products in qf.metabolites(mol, attach_phase1_steps=True):
        for product in products:
            if not product.HasProp("phase1_steps"):
                continue
            stamped += 1
            if stamped > _MAX_STAMPED:
                return
            plan = StepPlan.from_mol(product)
            target = unmapped_smiles(product)

            full_hits = _smi_sets(
                lin.apply(mol) for lin in plan.iter_as_linearizations()
            )
            for smis in full_hits:
                assert target in smis
            if full_hits:
                assert len(set(full_hits)) == 1

            if len(plan) > 1:
                prep_hits = _smi_sets(
                    lin.apply(mol, drop_last=1)
                    for lin in plan.iter_as_linearizations()
                )
                if prep_hits:
                    assert len(set(prep_hits)) == 1


@given(smiles=st.sampled_from(("C=C", "CC=C", "C=Cc1ccccc1", "CCN", "CCNC")))
@settings(
    max_examples=20,
    deadline=20_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_degenerate_bond_rules_replay(smiles: str):
    """Epoxidation / NDealkylation singleton plans apply to an emitted product set."""
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)

    for rule in (Epoxidation(), NDealkylation()):
        try:
            _site, products = next(
                rule.metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
            )
        except StopIteration:
            continue
        assume(products and products[0].HasProp("phase1_steps"))
        plan = StepPlan.from_mol(products[0])
        assert len(plan) == 1
        expected = can_smi_set(products)
        for lin in plan.iter_as_linearizations():
            result = lin.apply(mol)
            assert result
            got = can_smi_set(result)
            assert expected == got or expected <= got or got <= expected
