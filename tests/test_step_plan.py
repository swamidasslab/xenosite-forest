"""Tests for Step / StepPlan / AtomRef apply+resolve helpers."""

from __future__ import annotations

import math

import pytest
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest.step_plan import AtomRef, Step, StepPlan
from xenosite.forest.utils import unmapped_smiles


def test_singleton_one_linearization():
    plan = StepPlan.singleton("Epoxidation", {0, 1})
    assert len(plan) == 1
    assert plan.steps[0] == Step("Epoxidation", frozenset({0, 1}))
    assert plan.precedes == ()
    assert list(plan.iter_linearizations()) == [
        (Step("Epoxidation", frozenset({0, 1})),)
    ]


def test_deps_topo_sorts_not_kahn_layers():
    """Deps keeps Dealk free of OH→DH (unlike Kahn And(preps)→DH)."""
    from xenosite.forest.step_plan import Deps

    dealk = Step("Dealkylation", frozenset({0, 1}))
    oh = Step("Hydroxylation", frozenset({3}))
    dh = Step("Dehydrogenation", frozenset({1, 3}))
    # Only OH precedes DH.
    plan = Deps([dealk, oh, dh], precedes=[(1, 2)])
    orders = set(plan.iter_linearizations())
    assert len(orders) == 3
    assert plan.n_linearizations() == 3
    assert (oh, dh, dealk) in orders
    assert plan.contains([dealk, oh, dh], by="rule")
    assert plan.contains([oh, dh, dealk], by="rule")
    assert not plan.contains([dh, oh, dealk], by="rule")
    roundtrip = StepPlan.from_json(plan.to_json())
    assert isinstance(roundtrip, Deps)
    assert set(roundtrip.iter_linearizations()) == orders


def test_n_linearizations_matches_enumeration():
    from xenosite.forest.step_plan import And, Deps, Or

    a = Step("A", {0})
    b = Step("B", {1})
    c = Step("C", {2})
    d = Step("D", {3})

    seq = StepPlan((a, b, c))
    assert seq.n_linearizations() == 1
    assert seq.n_linearizations() == len(list(seq.iter_linearizations()))

    and_plan = And((a, b, c))
    assert and_plan.n_linearizations() == 6
    assert and_plan.n_linearizations() == len(list(and_plan.iter_linearizations()))

    # Seq of And(preps) then final — 2 orders.
    layered = StepPlan.layers([[a, b], [c]])
    assert layered.n_linearizations() == 2
    assert layered.n_linearizations() == len(list(layered.iter_linearizations()))

    # Or of different-length branches.
    or_plan = Or((StepPlan((a, b)), StepPlan((c,))))
    assert or_plan.n_linearizations() == 2
    assert or_plan.n_linearizations() == len(list(or_plan.iter_linearizations()))

    # And of Or (variable-length children) still exact.
    mixed = And((Or((a, StepPlan((b, c)))), d))
    assert mixed.n_linearizations() == len(list(mixed.iter_linearizations()))

    deps = Deps([a, b, c], precedes=[(0, 2), (1, 2)])
    assert deps.n_linearizations() == 2  # And(a,b)→c
    assert deps.n_linearizations() == len(list(deps.iter_linearizations()))

    free = Deps([a, b, c, d], precedes=[])
    assert free.n_linearizations() == 24

    # Independent components: many free nodes stay cheap (no 2^n on full n).
    many = Deps([Step("S", {i}) for i in range(12)], precedes=[])
    assert many.n_linearizations() == math.factorial(12)


def test_deps_same_linearizations_via_transitive_closure():
    """Lin-set identity is canonical edges — not ``==`` / raw precedes."""
    from xenosite.forest.step_plan import (
        Deps,
        canonical_dependency_edges,
    )

    a = Step("A", {0})
    b = Step("B", {1})
    c = Step("C", {2})
    chain = Deps([a, b, c], precedes=[(0, 1), (1, 2)])
    with_transitive = Deps([a, b, c], precedes=[(0, 1), (1, 2), (0, 2)])
    # Construction reduces edges → stable equal output.
    assert chain == with_transitive
    assert chain.precedes == ((0, 1), (1, 2))
    assert with_transitive.precedes == ((0, 1), (1, 2))
    assert chain.same_linearizations(with_transitive)
    assert canonical_dependency_edges(3, [(0, 1), (1, 2), (0, 2)]) == (
        (0, 1),
        (1, 2),
    )

    flipped = Deps([c, a, b], precedes=[(1, 2), (2, 0)])  # a≺b≺c
    assert chain.same_linearizations(flipped)

    layered = StepPlan.layers([[a, b], [c]])
    assert Deps(layered.steps, layered.precedes).same_linearizations(
        Deps([a, b, c], precedes=[(0, 2), (1, 2)])
    )

    free_dealk = Deps([a, b, c], precedes=[(1, 2)])  # only b≺c
    assert not free_dealk.same_linearizations(chain)

    # Graph algorithms ignore Step identity — same_linearizations must not.
    other_nodes = Deps(
        [Step("X", {0}), Step("Y", {1}), Step("Z", {2})],
        precedes=[(0, 1), (1, 2)],
    )
    assert canonical_dependency_edges(3, [(0, 1), (1, 2)]) == canonical_dependency_edges(
        3, [(0, 1), (1, 2)]
    )
    assert not chain.same_linearizations(other_nodes)


def test_layers_two_prep_then_final_two_orders():
    h0 = Step("Hydroxylation", frozenset({0}))
    h3 = Step("Hydroxylation", frozenset({3}))
    dh = Step(
        "Dehydrogenation",
        frozenset(
            [
                AtomRef(added_by=("Hydroxylation", frozenset({0}))),
                AtomRef(added_by=("Hydroxylation", frozenset({3}))),
            ]
        ),
    )
    plan = StepPlan.layers([[h0, h3], [dh]])
    assert len(plan) == 3
    orders = list(plan.iter_linearizations())
    assert len(orders) == 2
    assert set(orders) == {(h0, h3, dh), (h3, h0, dh)}


def test_json_round_trip():
    plan = StepPlan.layers(
        [
            [Step("Hydroxylation", {1}), Step("Hydroxylation", {2})],
            [
                Step(
                    "Dehydrogenation",
                    {
                        AtomRef(added_by=("Hydroxylation", {1})),
                        AtomRef(added_by=("Hydroxylation", {2})),
                    },
                )
            ],
        ]
    )
    restored = StepPlan.from_json(plan.to_json())
    assert restored == plan
    assert list(restored.iter_linearizations()) == list(plan.iter_linearizations())


def test_from_mol_and_attach():
    mol = MolFromSmiles("CCO")
    plan = StepPlan.singleton("Hydroxylation", {0})
    assert StepPlan.try_from_mol(mol) is None
    plan.attach_to_mol(mol)
    assert StepPlan.try_from_mol(mol) == plan
    assert StepPlan.from_mol(mol) == plan


def test_from_mol_missing_prop():
    mol = MolFromSmiles("C")
    with pytest.raises(ValueError, match="attached StepPlan"):
        StepPlan.from_mol(mol)


def test_empty_plan():
    plan = StepPlan((), ())
    assert list(plan.iter_linearizations()) == [()]


def test_cycle_raises():
    from xenosite.forest.step_plan import Deps

    a = Step("A", {0})
    b = Step("B", {1})
    with pytest.raises(ValueError, match="cycle"):
        Deps((a, b), ((0, 1), (1, 0)))


def test_compact_str():
    assert str(AtomRef(origin=3)) == "3"
    assert str(AtomRef(added_by=("Hydroxylation", {0}))) == "Hydroxylation[0]"
    assert str(Step("Hydroxylation", {0})) == "Hydroxylation[0]"
    dh = Step(
        "Dehydrogenation",
        {
            AtomRef(added_by=("Hydroxylation", {0})),
            AtomRef(added_by=("Hydroxylation", {3})),
        },
    )
    assert str(dh) == "Dehydrogenation[Hydroxylation[0], Hydroxylation[3]]"
    plan = StepPlan.layers(
        [
            [Step("Hydroxylation", {0}), Step("Hydroxylation", {3})],
            [dh],
        ]
    )
    assert str(plan) == (
        "(Hydroxylation[0] & Hydroxylation[3]) → "
        "Dehydrogenation[Hydroxylation[0], Hydroxylation[3]]"
    )


def test_hydroxylation_apply_records_atom_ref():
    mol = MolFromSmiles("c1ccccc1")
    step = Step("Hydroxylation", {0})
    products = step.apply(mol)
    assert products
    product = products[0]
    ref = AtomRef(added_by=("Hydroxylation", frozenset({0})))
    o_idx = ref.resolve(product)
    assert product.GetAtomWithIdx(o_idx).GetAtomicNum() == 8
    carbon_site = step.resolve_site(product)
    assert len(carbon_site) == 1


def test_added_by_and_origin_carry_frame_depth():
    """Created atoms lack depth-0; site idxs need the frame they were written in."""
    from xenosite.forest.step_plan import _frame_depth, _origin_refs

    benzene = MolFromSmiles("c1ccccc1")
    mono = Step("Hydroxylation", {0}).apply(benzene)[0]
    depth1 = _frame_depth(mono)
    assert depth1 == 1
    assert mono._forest["atom_trace"]["depth"] == 1

    # Second OH site is a carbon GetIdx on the mono-OH frame.
    carbons = [
        a.GetIdx()
        for a in mono.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetTotalNumHs() >= 1
    ]
    site_c = carbons[3] if len(carbons) > 3 else carbons[0]
    # Reactant-stable carbon → depth-0 origin; created O uses mid-frame.
    carbon_refs = _origin_refs([site_c], mono)
    c_ref = next(iter(carbon_refs))
    assert c_ref.depth == 0
    di = Step("Hydroxylation", carbon_refs).apply(mono)[0]
    depth2 = _frame_depth(di)
    assert depth2 == 2
    assert di._forest["atom_trace"]["depth"] == 2

    created = AtomRef(
        added_by=("Hydroxylation", frozenset([c_ref.origin])), depth=c_ref.depth
    )
    assert created.depth == 0
    o_idx = created.resolve(di)
    assert di.GetAtomWithIdx(o_idx).GetAtomicNum() == 8

    # New oxygens have no depth-0 identity — origin refs stay at di-OH frame.
    oxygens = [a.GetIdx() for a in di.GetAtoms() if a.GetAtomicNum() == 8]
    assert len(oxygens) == 2
    o_refs = _origin_refs(oxygens, di)
    assert all(r.depth == depth2 for r in o_refs)
    for r in o_refs:
        assert di.GetAtomWithIdx(r.resolve(di)).GetAtomicNum() == 8
        # That idx is not this oxygen's depth-0 identity. A hit is some
        # other atom; a miss is a KeyError.
        try:
            remapped = AtomRef(origin=r.origin, depth=0).resolve(di)
        except KeyError:
            continue
        assert remapped != r.resolve(di)
        assert di.GetAtomWithIdx(remapped).GetAtomicNum() != 8


def test_benzene_prep_linearization_apply_agrees():
    mol = MolFromSmiles("c1ccccc1")
    h0 = Step("Hydroxylation", {0})
    h3 = Step("Hydroxylation", {3})
    dh = Step(
        "Dehydrogenation",
        frozenset(
            [
                AtomRef(added_by=("Hydroxylation", {0})),
                AtomRef(added_by=("Hydroxylation", {3})),
            ]
        ),
    )
    plan = StepPlan.layers([[h0, h3], [dh]])
    prep_smiles = [
        frozenset(unmapped_smiles(p) for p in products)
        for lin in plan.linearizations()
        for products in [lin.apply(mol, drop_last=1)]
        if products
    ]
    assert prep_smiles and len(set(prep_smiles)) == 1


def test_created_atom_ref_follows_trace_label():
    """Dehydrogenation must not alias an older oxygen. A second hydroxylation
    at the same site resolves to the later oxygen on that mol, and the
    original frame does not resolve. A lying ``atom_refs`` index is ignored.
    """
    from xenosite.forest.base import AtomRefsIndex

    start = MolFromSmiles("COc1ccc(O)cc1")
    with pytest.raises(KeyError):
        AtomRef(added_by=("Hydroxylation", frozenset({0}))).resolve(start)

    catechol = Step("Hydroxylation", {4}).apply(start)[0]
    quinone = Step(
        "Dehydrogenation",
        {
            AtomRef(added_by=("Hydroxylation", frozenset({4}))),
            AtomRef(origin=6),
        },
    ).apply(catechol)[0]
    records = quinone._forest["atom_trace"]["records"]
    assert not any(
        rec.get("added_by") and rec["added_by"][0] == "Dehydrogenation"
        for rec in records.values()
    )
    frame = quinone._forest["atom_trace"]["depth"]
    oh_now = [
        rec["idx"][rec["depth"].index(frame)]
        for rec in records.values()
        if rec.get("added_by")
        and rec["added_by"][0] == "Hydroxylation"
        and 4 in rec["added_by"][1]
        and frame in rec["depth"]
    ]
    quinone._forest["atom_refs"] = AtomRefsIndex(
        {(("Hydroxylation", frozenset({4})), 0)}
    )
    followed = AtomRef(added_by=("Hydroxylation", frozenset({4}))).resolve(quinone)
    assert oh_now == [followed]
    assert followed != 0
    assert quinone.GetAtomWithIdx(followed).GetAtomicNum() == 8

    once = Step("Hydroxylation", {0}).apply(start)[0]
    twice = Step("Hydroxylation", {0}).apply(once)[0]
    frame = twice._forest["atom_trace"]["depth"]
    births = []
    for rec in twice._forest["atom_trace"]["records"].values():
        added = rec.get("added_by")
        if not added or added[0] != "Hydroxylation" or 0 not in added[1]:
            continue
        if frame not in rec["depth"]:
            continue
        births.append((rec["depth"][0], rec["idx"][rec["depth"].index(frame)]))
    assert len(births) == 2
    twice._forest["atom_refs"] = AtomRefsIndex(
        {(("Hydroxylation", frozenset({0})), min(births)[1])}
    )
    resolved = AtomRef(added_by=("Hydroxylation", frozenset({0}))).resolve(twice)
    assert resolved == max(births)[1]
    assert resolved != min(births)[1]
    once_frame = once._forest["atom_trace"]["depth"]
    once._forest["atom_refs"] = None
    once_idx = AtomRef(added_by=("Hydroxylation", frozenset({0}))).resolve(once)
    once_oh = [
        rec["idx"][rec["depth"].index(once_frame)]
        for rec in once._forest["atom_trace"]["records"].values()
        if rec.get("added_by")
        and rec["added_by"][0] == "Hydroxylation"
        and 0 in rec["added_by"][1]
        and once_frame in rec["depth"]
    ]
    assert once_oh == [once_idx]
