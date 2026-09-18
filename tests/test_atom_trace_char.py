"""Characterization of AtomTracker tags / maps / atom_refs (pre-_forest refactor).

These snapshots are the regression oracle: after moving traces onto
``mol._forest``, this file must pass with no expected-value edits.
"""

from __future__ import annotations

from rdkit.Chem import AddHs
from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest import AtomTrace, rules
from xenosite.forest.base import AtomTracker
from xenosite.forest.bfs import bfs, dfs
from xenosite.forest.step_plan import AtomRef, Step, _atom_refs_index
from xenosite.forest.trace import atom_no


def _tag_snap(tags: dict):
    """Stable (tag, depths, idxs) tuples for frozen comparisons."""
    return tuple(
        sorted((int(k), tuple(v["depth"]), tuple(v["idx"])) for k, v in tags.items())
    )


def _map_snap(mol):
    return tuple(
        (a.GetIdx(), a.GetAtomicNum(), a.GetAtomMapNum())
        for a in mol.GetAtoms()
        if a.GetAtomicNum() != 1
    )


def test_cco_hydroxylation_tag_snapshot():
    mol = MolFromSmiles("CCO")
    _, products = next(rules.Hydroxylation().metabolites_from_sites(mol, [frozenset({0})]))
    p = products[0]
    tags = AtomTracker.tags(p)
    assert AtomTracker.depths(p) == [0, 1]
    assert _tag_snap(tags) == (
        (0, (0, 1), (0, 0)),
        (1, (0, 1), (1, 2)),
        (2, (0, 1), (2, 3)),
        (3, (1,), (1,)),
    )
    assert _map_snap(p) == ((0, 6, 1), (1, 8, 0), (2, 6, 2), (3, 8, 3))
    compact = AtomTracker.tags(p, compact=True)
    assert compact[0] == {0: 1, 1: 1}
    assert compact[3] == {1: 2}


def test_tags_filters_and_strict_false_snapshot():
    mol = MolFromSmiles("CCO")
    _, products = next(rules.Hydroxylation().metabolize(mol))
    p = products[0]
    by_depth = AtomTracker.tags(p, depth=1)
    assert all(1 in rec["depth"] for rec in by_depth.values())
    some_idx = next(iter(AtomTracker.tags(p).values()))["idx"][-1]
    by_idx = AtomTracker.tags(p, idx=some_idx)
    assert by_idx
    assert AtomTracker.tags(MolFromSmiles("CC"), strict=False) == {}


def test_initialize_tags_skips_explicit_hydrogens_snapshot():
    mol = AddHs(MolFromSmiles("C"))
    AtomTracker().initialize_tags(mol)
    tags = AtomTracker.tags(mol)
    assert _tag_snap(tags) == ((0, (0,), (0,)),)
    assert mol.GetAtomWithIdx(0).GetAtomicNum() == 6


def test_atom_trace_map_follow_origin_snapshot():
    mol = MolFromSmiles("CCO")
    _, products = next(
        rules.Hydroxylation().metabolites_from_sites(mol, [frozenset({0})])
    )
    p = products[0]
    t = AtomTrace(p)
    assert t.depths == (0, 1)
    assert t.map() == {1: 1, 2: 3, 3: 4}
    assert t.follow(1) == (1, 1)
    assert t.origin(2) is None
    assert t.added() == frozenset({2})
    assert t.removed() == frozenset()
    for atom in p.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        origin = t.origin(atom_no(atom.GetIdx()))
        if origin is None:
            assert atom.GetAtomMapNum() == 0
        else:
            assert origin == atom.GetAtomMapNum()


def test_do_not_tag_atoms_leaves_no_trace():
    mol = MolFromSmiles("CCO")
    _, products = next(rules.Hydroxylation().metabolize(mol, do_not_tag_atoms=True))
    p = products[0]
    assert AtomTracker.tags(p, strict=False) == {}
    assert all(a.GetAtomMapNum() == 0 for a in p.GetAtoms())
    import pytest

    with pytest.raises(ValueError):
        AtomTrace(p)


def test_bfs_and_dfs_path_products_carry_tags():
    kwargs = dict(
        ruleset="Full",
        depth=2,
        all_paths=False,
        outmols=True,
        phase1=True,
        ismi=True,
    )
    bfs_hits = list(bfs(["CC", "OCCO"], **kwargs))
    dfs_hits = list(dfs(["CC", "OCCO"], **kwargs))
    assert bfs_hits and dfs_hits
    for hits in (bfs_hits, dfs_hits):
        path_mols = hits[0][2]
        end = path_mols[-1]
        assert AtomTracker.tags(end)
        assert AtomTracker.depths(end) == list(range(len(path_mols)))
        tags = AtomTracker.tags(end)
        assert len(tags) >= end.GetNumHeavyAtoms()
        # Every heavy atom present at last depth is tagged.
        last = max(AtomTracker.depths(end))
        present = {
            rec["idx"][rec["depth"].index(last)]
            for rec in tags.values()
            if last in rec["depth"]
        }
        heavies = {
            a.GetIdx() for a in end.GetAtoms() if a.GetAtomicNum() not in (0, 1)
        }
        assert heavies <= present


def test_step_apply_hydroxylation_atom_refs_snapshot():
    """Frozen creation index from Step.apply (parity bar for shared recording)."""
    mol = MolFromSmiles("CCO")
    outs = Step("Hydroxylation", {AtomRef(origin=0)}).apply(mol)
    assert len(outs) == 1
    p = outs[0]
    index = _atom_refs_index(p)
    assert list(index.items()) == [(("Hydroxylation", frozenset({0})), 1)]
    # Tagged product still exposes the same public tag/map shape as metabolize.
    assert _tag_snap(AtomTracker.tags(p)) == (
        (0, (0, 1), (0, 0)),
        (1, (0, 1), (1, 2)),
        (2, (0, 1), (2, 3)),
        (3, (1,), (1,)),
    )


def test_metabolize_and_step_apply_share_atom_refs():
    mol = MolFromSmiles("CCO")
    _, metab = next(
        rules.Hydroxylation().metabolites_from_sites(mol, [frozenset({0})])
    )
    apply_outs = Step("Hydroxylation", {AtomRef(origin=0)}).apply(MolFromSmiles("CCO"))
    assert list(_atom_refs_index(metab[0]).items()) == list(
        _atom_refs_index(apply_outs[0]).items()
    )
    assert list(_atom_refs_index(metab[0]).items()) == [
        (("Hydroxylation", frozenset({0})), 1)
    ]


def test_atom_label_stable_across_step_string_equality():
    """Parent ``_forestLabel`` props persist on child atoms; new atoms get a new label once."""
    parent = MolFromSmiles("CC")
    AtomTracker().initialize_tags(parent)
    parent_labels = {
        a.GetProp(AtomTracker.atom_tag_prop_name)
        for a in parent.GetAtoms()
        if a.HasProp(AtomTracker.atom_tag_prop_name)
    }
    assert parent_labels == {"0", "1"}

    _, products = next(rules.Hydroxylation().metabolize(parent))
    child = products[0]
    child_by_label = {
        a.GetProp(AtomTracker.atom_tag_prop_name): a.GetIdx()
        for a in child.GetAtoms()
        if a.HasProp(AtomTracker.atom_tag_prop_name)
    }
    assert set(child_by_label) >= parent_labels
    # Surviving carbons keep the same label strings.
    assert "0" in child_by_label and "1" in child_by_label
    new_labels = set(child_by_label) - parent_labels
    assert len(new_labels) == 1
    # Labels are not rewritten on a second save of the same mol.
    before = {
        a.GetIdx(): a.GetProp(AtomTracker.atom_tag_prop_name)
        for a in child.GetAtoms()
        if a.HasProp(AtomTracker.atom_tag_prop_name)
    }
    AtomTracker()._save_tags(
        child, AtomTracker.tags(child), depth=child._forest["atom_trace"]["depth"]
    )
    after = {
        a.GetIdx(): a.GetProp(AtomTracker.atom_tag_prop_name)
        for a in child.GetAtoms()
        if a.HasProp(AtomTracker.atom_tag_prop_name)
    }
    assert before == after


def test_two_step_tag_history_public_shape():
    """Depth-2 products expose multi-depth history via tags() (public API)."""
    result = list(
        bfs(
            ["CC", "OCCO"],
            "Full",
            depth=2,
            all_paths=False,
            outmols=True,
            phase1=True,
            ismi=True,
        )
    )
    path = result[0][2]
    assert len(path) == 3
    for mol, depth in zip(path, range(3)):
        tags = AtomTracker.tags(mol, depth=depth)
        assert mol.GetNumHeavyAtoms() == len(tags)
    end_tags = AtomTracker.tags(path[-1])
    # At least one original carbon has a 3-deep history.
    assert any(len(rec["depth"]) == 3 for rec in end_tags.values())
    assert AtomTracker.depths(path[-1]) == [0, 1, 2]
