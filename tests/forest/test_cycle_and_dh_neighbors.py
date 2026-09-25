"""Reject cycles by ancestor dedup_smi; DH sites need matching heavy neighbors."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.find_path import (
    _dh_neighbors_match_any_view,
    _is_dehydrogenation_effect,
    _site_could_help,
    atom_diff,
)
from xenosite.forest.rules import (
    Hydroxylation,
    _repeats_ancestor_dedup_smi,
    describe,
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
        describe(removes="H", dearomatizes=True)["possibilities"][0]
    )
    assert not _is_dehydrogenation_effect(
        describe(adds="O", removes="H")["possibilities"][0]
    )
    assert not _is_dehydrogenation_effect(
        describe(removes="H", dearomatizes=True, cleaves=True)["possibilities"][0]
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


def test_dh_filter_sites_refuses_mismatched_neighbors():
    reactant = Chem.MolFromSmiles("COc1ccc(O)cc1")
    target = Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")
    diff = atom_diff(reactant, target)
    oxygens = [a.GetIdx() for a in reactant.GetAtoms() if a.GetAtomicNum() == 8]
    bad = next(
        o for o in oxygens if not _dh_neighbors_match_any_view(reactant, o, diff)
    )
    other = next(o for o in oxygens if o != bad)
    end = describe(removes="H", dearomatizes=True, partner="O")["possibilities"][0]
    info = {
        "options": describe(removes="H", dearomatizes=True)["possibilities"][0],
        "ends": (end, end),
        "end_atoms": (bad, other),
        "path_ends": frozenset(),
        "site": frozenset({bad, other}),
    }
    assert not _site_could_help(info["site"], info, diff, reactant)
