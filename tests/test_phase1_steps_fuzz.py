"""Hypothesis fuzz: StepPlan linearizations replay to the same metabolite.

Example database lives in ``.hypothesis/`` (gitignored) and is restored on CI.

QuinoneFormation finals are labeled ``Dehydrogenation``, but Forest's
``Dehydrogenation`` rule does not fire on typical aromatic phenols / hydroquinones.
For those plans we assert prep-step linearizations agree when they fire, and only
require full-path SMILES equality when every step (including the final DH) succeeds.
"""

from __future__ import annotations

from pathlib import Path

from hypothesis import HealthCheck, assume, given, settings, strategies as st
from hypothesis.database import DirectoryBasedExampleDatabase
from rdkit import Chem

from xenosite.forest import StepPlan
from xenosite.forest.base import AtomTracker, can_smi_set
from xenosite.forest.rules import (
    Dealkylation,
    Dehydrogenation,
    Epoxidation,
    Hydroxylation,
    NDealkylation,
    OxidativeDehalogenation,
    QuinoneFormation,
)
from xenosite.forest.utils import unmapped_smiles

_HYPOTHESIS_DIR = Path(__file__).resolve().parents[1] / ".hypothesis" / "examples"
_HYPOTHESIS_DIR.mkdir(parents=True, exist_ok=True)
_HYPOTHESIS_DB = DirectoryBasedExampleDatabase(str(_HYPOTHESIS_DIR))

_RULES = {
    "Hydroxylation": Hydroxylation(),
    "Dehydrogenation": Dehydrogenation(),
    "Dealkylation": Dealkylation(),
    "OxidativeDehalogenation": OxidativeDehalogenation(),
    "Epoxidation": Epoxidation(),
    "NDealkylation": NDealkylation(),
}

# Aromatic / quinone-prone substrates plus a few bond-rule parents.
_CORPUS = (
    "CC(=O)Nc1ccc(O)cc1",  # APAP
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


def _map_site(mol, origin_site):
    """Map reactant atom indices onto the current mol (maps / react_atom_idx)."""
    origin_site = frozenset(origin_site)
    mapping = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() > 0:
            mapping[atom.GetAtomMapNum() - 1] = atom.GetIdx()
        elif atom.HasProp("react_atom_idx"):
            mapping[int(atom.GetProp("react_atom_idx"))] = atom.GetIdx()
    if not mapping:
        return origin_site
    try:
        return frozenset(mapping[i] for i in origin_site)
    except KeyError:
        return None


def _apply_step(mol, step):
    """Return list of product-lists for ``step``, or None if the rule cannot fire."""
    rule = _RULES.get(step.rule)
    if rule is None:
        return None
    mapped = _map_site(mol, step.site)
    if mapped is None or len(mapped) != len(step.site):
        return None
    emissions = []
    for _site, products in rule.metabolites_from_sites(
        mol,
        mapped,
        tag_atoms=True,
        only_emit_topologically_distinct_sites=False,
    ):
        frags = [p for p in products if p]
        if frags:
            emissions.append(frags)
    return emissions or None


def _keep_fragments(products, remaining_origins):
    """Prefer fragments that still carry atoms needed for later steps."""
    if not remaining_origins:
        return products
    kept = []
    for product in products:
        mapped = _map_site(product, remaining_origins)
        if mapped is not None and len(mapped) == len(remaining_origins):
            kept.append(product)
    return kept or products


def replay_linearization(mol, order):
    """Apply Forest rules in ``order``. Return final mols, or None if a step fails."""
    root = Chem.Mol(mol)
    _RULES["Hydroxylation"].initialize_tags(root)
    AtomTracker.add_current_idx_as_atom_prop(
        root, propname=AtomTracker.previous_index_prop_name
    )
    currents = [root]
    for i, step in enumerate(order):
        remaining = frozenset().union(*(s.site for s in order[i + 1 :]))
        nxt = []
        for cur in currents:
            emissions = _apply_step(cur, step)
            if not emissions:
                return None
            for products in emissions:
                nxt.extend(_keep_fragments(products, remaining))
        if not nxt:
            return None
        currents = nxt
    return currents


def _prep_prefix(order):
    """Drop a trailing Dehydrogenation (quinone final); else return the full order."""
    if order and order[-1].rule == "Dehydrogenation" and len(order) > 1:
        return order[:-1]
    return order


def replay_prep_toward_quinone(mol, order):
    """Replay prep steps, keeping fragments that still carry the final DH sites."""
    prep = _prep_prefix(order)
    if prep == order:
        return replay_linearization(mol, order)
    dh_site = order[-1].site
    root = Chem.Mol(mol)
    _RULES["Hydroxylation"].initialize_tags(root)
    AtomTracker.add_current_idx_as_atom_prop(
        root, propname=AtomTracker.previous_index_prop_name
    )
    currents = [root]
    for i, step in enumerate(prep):
        # Always retain atoms needed for the final DH, plus any later prep sites.
        later = frozenset().union(*(s.site for s in prep[i + 1 :]), dh_site)
        nxt = []
        for cur in currents:
            emissions = _apply_step(cur, step)
            if not emissions:
                return None
            for products in emissions:
                nxt.extend(_keep_fragments(products, later))
        if not nxt:
            return None
        currents = nxt
    return currents


@st.composite
def corpus_smiles(draw):
    return draw(st.sampled_from(_CORPUS))


def test_benzene_prep_linearizations_agree():
    """Both OH orders for para-quinone prep yield the same hydroquinone."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    plans = QuinoneFormation().phase1_steps(mol, frozenset({0, 3}))
    layered = [p for p in plans if len(p) == 3]
    assert layered
    finals = []
    for order in layered[0].iter_linearizations():
        result = replay_prep_toward_quinone(mol, order)
        assert result is not None
        finals.append(frozenset(unmapped_smiles(m) for m in result))
    assert len(set(finals)) == 1
    # Kekule vs aromatic SMILES both OK; round-trip to one form.
    got = Chem.MolToSmiles(Chem.MolFromSmiles(next(iter(next(iter(finals))))))
    assert got == Chem.MolToSmiles(Chem.MolFromSmiles("Oc1ccc(O)cc1"))


def test_epoxidation_linearization_replays_product():
    mol = Chem.MolFromSmiles("C=C")
    site, products = next(
        Epoxidation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    plan = StepPlan.from_mol(products[0])
    expected = can_smi_set(products)
    for order in plan.iter_linearizations():
        result = replay_linearization(mol, order)
        assert result is not None
        assert can_smi_set(result) == expected


def test_ndealkylation_linearization_replays_product():
    mol = Chem.MolFromSmiles("CCN")
    site, products = next(
        NDealkylation().metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
    )
    plan = StepPlan.from_mol(products[0])
    expected = can_smi_set(products)
    matched = False
    for order in plan.iter_linearizations():
        result = replay_linearization(mol, order)
        if result is None:
            continue
        got = can_smi_set(result)
        if expected == got or expected <= got or got <= expected:
            matched = True
            break
    assert matched


@given(smiles=corpus_smiles())
@settings(
    max_examples=40,
    deadline=30_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
def test_quinone_linearizations_prep_agree_or_full_match(smiles: str):
    """Quinone plans: prep orders agree; full replay matches when DH can fire."""
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
            prep_endpoints = []
            full_hits = []
            for order in plan.iter_linearizations():
                full = replay_linearization(mol, order)
                if full is not None:
                    smis = frozenset(unmapped_smiles(m) for m in full)
                    full_hits.append(smis)
                    assert target in smis

                if _prep_prefix(order) != order:
                    mid = replay_prep_toward_quinone(mol, order)
                    if mid is not None:
                        prep_endpoints.append(
                            frozenset(unmapped_smiles(m) for m in mid)
                        )

            if prep_endpoints:
                assert len(set(prep_endpoints)) == 1
            if full_hits:
                assert len(set(full_hits)) == 1

@given(smiles=st.sampled_from(("C=C", "CC=C", "C=Cc1ccccc1", "CCN", "CCNC")))
@settings(
    max_examples=20,
    deadline=20_000,
    database=_HYPOTHESIS_DB,
    suppress_health_check=[HealthCheck.too_slow],
)
def test_degenerate_bond_rules_replay(smiles: str):
    """Epoxidation / NDealkylation singleton plans replay to an emitted product set."""
    mol = Chem.MolFromSmiles(smiles)
    assume(mol is not None)

    for rule in (Epoxidation(), NDealkylation()):
        try:
            site, products = next(
                rule.metabolize(mol, attach_phase1_steps=True, tag_atoms=False)
            )
        except StopIteration:
            continue
        assume(products and products[0].HasProp("phase1_steps"))
        plan = StepPlan.from_mol(products[0])
        assert len(plan) == 1
        expected = can_smi_set(products)
        for order in plan.iter_linearizations():
            result = replay_linearization(mol, order)
            assert result is not None
            got = can_smi_set(result)
            assert expected == got or expected <= got or got <= expected
