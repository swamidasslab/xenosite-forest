"""Reject cycles by ancestor dedup_smi; DH ends match target after application."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.find_path import (
    _dh_neighbors_match_any_view,
    _dh_product_ends_match_target,
    _is_dehydrogenation_effect,
    atom_diff,
)
from xenosite.forest.rules import (
    Dehydrogenation,
    Hydroxylation,
    _repeats_ancestor_dedup_smi,
    _describe,
)


def test_repeats_ancestor_dedup_smi_on_identity_hop():
    """Product whose dedup_smi equals a parent frame is a cycle."""

    parent = Chem.MolFromSmiles("CC")
    parent.xf.tracing._stamp()
    product = Chem.Mol(parent)
    product.xf.tracing._stamp()
    trace = product._forest["atom_trace"]
    trace["depth"] = 1
    trace["dedup_smi"] = ["CC", "CC"]
    assert _repeats_ancestor_dedup_smi(product)


def test_hydroxylation_product_does_not_repeat_parent():
    mol = Chem.MolFromSmiles("CC")
    finished, _ = next(Hydroxylation().metabolize(mol))
    product = finished[0]
    assert product.xf.tracing.dedup_smi == "CCO"
    assert not _repeats_ancestor_dedup_smi(product)


def test_is_dehydrogenation_effect_shape():
    assert _is_dehydrogenation_effect(
        _describe(removes="H", dearomatizes=True)["possibilities"][0]
    )
    assert not _is_dehydrogenation_effect(
        _describe(adds="O", removes="H")["possibilities"][0]
    )
    assert not _is_dehydrogenation_effect(
        _describe(removes="H", dearomatizes=True, cleaves=True)["possibilities"][0]
    )


def test_dh_neighbors_match_hydroquinone_to_quinone():
    """Para-hydroquinone → quinone: O ends keep the same heavy neighbors."""

    reactant = Chem.MolFromSmiles("Oc1ccc(O)cc1")
    target = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    diff = atom_diff(reactant, target)
    oxygens = [a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomicNum() == 8]
    assert len(oxygens) == 2
    for o in oxygens:
        assert _dh_neighbors_match_any_view(reactant, o, diff)


def test_dh_neighbors_mismatch_on_anisole_toward_quinone():
    """Anisole → quinone: methoxy O connectivity ≠ quinone carbonyl."""

    reactant = Chem.MolFromSmiles("COc1ccc(O)cc1")
    target = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    diff = atom_diff(reactant, target)
    oxygens = [a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomicNum() == 8]
    assert len(oxygens) == 2
    mismatches = [
        o for o in oxygens if not _dh_neighbors_match_any_view(reactant, o, diff)
    ]
    assert mismatches


def test_dh_product_ends_match_after_hydroquinone_dh():
    """Applied DH: product ends match the quinone target."""

    parent = Chem.MolFromSmiles("Oc1ccc(O)cc1")
    target = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    finished, info = next(Dehydrogenation().metabolize(parent))
    product = finished[0]
    assert product.xf.csmi == Chem.MolToSmiles(target)
    assert _dh_product_ends_match_target(
        parent, product, info["end_atoms"], target
    )


def test_dh_product_ends_refuse_when_connectivity_still_wrong():
    """Post-check uses product MCS — reactant mismatch alone is not the gate."""

    # Synthetic: check product (still anisole-like connectivity) vs quinone.
    parent = Chem.MolFromSmiles("COc1ccc(O)cc1")
    parent.xf.tracing._stamp()
    product = Chem.Mol(parent)
    product.xf.tracing._stamp()
    target = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    oxygens = [a.GetIdx() for a in parent.GetAtoms() if a.GetAtomicNum() == 8]
    assert not _dh_product_ends_match_target(parent, product, oxygens, target)
