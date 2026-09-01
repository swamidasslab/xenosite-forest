"""Coverage for xenosite.forest.base helpers and reaction-rule plumbing."""

from __future__ import annotations

from collections import defaultdict
from unittest.mock import patch

import pytest
from rdkit.Chem import AddHs
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolFromSmarts, MolToSmiles

from xenosite.forest import rules
from xenosite.forest.base import (
    AromaticSystems,
    AtomTracker,
    ConjugatedSystems,
    EditMol,
    QueryMol,
    ReactionRule,
    Resonate,
    SmartsReactionRule,
    can_smi,
    can_smi_set,
)


# ---------------------------------------------------------------------------
# can_smi / can_smi_set
# ---------------------------------------------------------------------------


def test_can_smi_list_and_disconnected():
    mols = [MolFromSmiles("CCO"), MolFromSmiles("CC")]
    nested = can_smi(rdmol=mols)
    assert nested == [can_smi(rdmol=mols[0]), can_smi(rdmol=mols[1])]

    parts = can_smi(line="C.O")
    assert set(parts) == {"C", "O"}

    salt = MolFromSmiles("C.O")
    assert set(can_smi(rdmol=salt)) == {"C", "O"}


def test_can_smi_set_unifies_equivalent_mols():
    a = MolFromSmiles("CCO")
    b = MolFromSmiles("OCC")
    assert can_smi_set([a, b]) == frozenset({"CCO"})


# ---------------------------------------------------------------------------
# AtomTracker
# ---------------------------------------------------------------------------


def _tagged_hydroxylation_product(smiles="CCO"):
    mol = MolFromSmiles(smiles)
    _, products = next(rules.Hydroxylation().metabolize(mol))
    return mol, products[0]


def test_tags_list_strict_false_idx_compact_and_bad_type():
    mol, product = _tagged_hydroxylation_product()
    untagged = MolFromSmiles("CC")

    chained = list(AtomTracker.tags([product, untagged], strict=False))
    assert chained

    empty = AtomTracker.tags(untagged, strict=False)
    assert empty == {}

    with pytest.raises(ValueError, match="Must submit RDKit Mol or dict"):
        AtomTracker.tags("not a mol")

    recs = AtomTracker.tags(product)
    some_idx = next(iter(recs.values()))["idx"][-1]
    by_idx = AtomTracker.tags(product, idx=some_idx)
    assert by_idx
    assert all(some_idx in rec["idx"] for rec in by_idx.values())

    compact = AtomTracker.tags(product, compact=True)
    assert compact
    first = next(iter(compact.values()))
    assert min(first.values()) >= 1


def test_depths_rejects_non_mol():
    with pytest.raises(ValueError, match="Must submit RDKit Mol or dict"):
        AtomTracker.depths("not a mol")


def test_metabolite_index_to_reversed_index_record():
    _, product = _tagged_hydroxylation_product()
    rec = AtomTracker.metabolite_index_to_reversed_index_record(
        product, exact_depth=2
    )
    assert rec is not None
    assert rec
    for metabolite_idx, history in rec.items():
        assert metabolite_idx == history[-1]
        assert len(history) == 2

    assert (
        AtomTracker.metabolite_index_to_reversed_index_record(
            product, exact_depth=10
        )
        is None
    )


def test_initialize_tags_skips_explicit_hydrogens():
    mol = AddHs(MolFromSmiles("C"))
    AtomTracker().initialize_tags(mol)
    tags = AtomTracker.tags(mol)
    assert len(tags) == 1
    assert all(
        mol.GetAtomWithIdx(rec["idx"][0]).GetAtomicNum() != 1
        for rec in tags.values()
    )


def test_tag_copies_missing_previous_tags():
    mol = MolFromSmiles("CCO")
    _, dehydrated = next(rules.Dehydration().metabolize(mol))
    cc = dehydrated[0]
    _, hydroxylated = next(rules.Hydroxylation().metabolize(cc))
    product = hydroxylated[0]
    tags = AtomTracker.tags(product)
    depths = AtomTracker.depths(tags)
    assert 0 in depths and 1 in depths
    # Oxygen lost in dehydration is still in the tag record (copied forward).
    assert any(len(rec["depth"]) == 1 and rec["depth"] == [0] for rec in tags.values())


def test_next_tag_missing_and_non_integer():
    mol = MolFromSmiles("CC")
    assert AtomTracker._next_tag(mol) == 1

    mol.SetProp(AtomTracker.last_tag_name, "not-an-int")
    assert AtomTracker._next_tag(mol) == 1


def test_old_to_new_uses_react_atom_idx():
    mol = MolFromSmiles("CC")
    mol.GetAtomWithIdx(0).SetProp("react_atom_idx", "5")
    mapping = AtomTracker()._old_to_new_atom_indexes(mol)
    assert mapping[5] == 0


def test_stamp_origin_maps_untagged_empty_and_hydrogens():
    tracker = AtomTracker()
    untagged = MolFromSmiles("CC")
    tracker._stamp_origin_maps(untagged)

    empty = MolFromSmiles("CC")
    tracker._stamp_origin_maps(empty, tags={})

    mol = AddHs(MolFromSmiles("C"))
    tracker.initialize_tags(mol)
    tracker._stamp_origin_maps(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            assert atom.GetAtomMapNum() == 0
        else:
            assert atom.GetAtomMapNum() == 1


def test_reactant_aligned_order_emits_hydrogens_and_stragglers():
    tracker = AtomTracker()
    mol = AddHs(MolFromSmiles("CC"))
    origin_of = {0: 0, 1: 1}
    order = tracker._reactant_aligned_order(mol, origin_of)
    assert set(order) == set(range(mol.GetNumAtoms()))
    assert order[:2] == [0, 1]

    with patch("xenosite.forest.base.Chem.GetMolFrags", return_value=[]):
        leftover = tracker._reactant_aligned_order(mol, origin_of)
    assert leftover == list(range(mol.GetNumAtoms()))


def test_align_and_stamp_untagged_and_empty_tags():
    tracker = AtomTracker()
    untagged = MolFromSmiles("CC")
    assert tracker._align_and_stamp(untagged) is untagged

    tagged = MolFromSmiles("CC")
    tagged.SetProp(AtomTracker.tag_name, "{}")
    assert tracker._align_and_stamp(tagged) is tagged


def test_save_tags_converts_defaultdict():
    mol = MolFromSmiles("C")
    tags = defaultdict(dict)
    tags[0] = {"idx": [0], "depth": [0]}
    AtomTracker()._save_tags(mol, tags)
    loaded = AtomTracker.tags(mol)
    assert loaded[0]["idx"] == [0]


# ---------------------------------------------------------------------------
# Conjugated / aromatic systems
# ---------------------------------------------------------------------------


def test_add_orphan_double_bond_when_ring_already_assigned():
    mol = MolFromSmiles("C=CCCc1ccccc1")
    ring = next(iter(AromaticSystems().systems(mol)))
    sys_id = {0: set(ring)}
    ConjugatedSystems()._add_orphan_double_bond_systems(sys_id, mol)
    assert len(sys_id) >= 2
    vinyl = next(v for k, v in sys_id.items() if k != 0)
    assert vinyl & set(ring) == set()
    assert len(vinyl) >= 1

    empty = {}
    ConjugatedSystems()._add_orphan_double_bond_systems(empty, MolFromSmiles("C=C"))
    assert empty


def test_systems_kekulize_and_resonance_failures():
    mol = MolFromSmiles("c1ccccc1")
    with patch("xenosite.forest.base.Kekulize", side_effect=ValueError):
        systems = ConjugatedSystems().systems(mol)
    assert systems

    with patch("xenosite.forest.base.ResonanceMolSupplier", side_effect=ValueError):
        assert ConjugatedSystems().systems(mol) == [set()]


def test_fragment_on_bonds_kekulize_failure_and_empty_bonds():
    mol = MolFromSmiles("C=CC")
    cs = ConjugatedSystems()
    with patch("xenosite.forest.base.Kekulize", side_effect=ValueError):
        out, types = cs._fragment_on_bonds(MolFromSmiles("C=CC"), [(0, 1)])
    assert types
    assert out.GetNumAtoms() > mol.GetNumAtoms()

    intact, empty = cs._fragment_on_bonds(MolFromSmiles("CC"), [])
    assert empty == {}
    assert MolToSmiles(intact) == "CC"


def test_restore_bonds_updates_existing_bond():
    mol = MolFromSmiles("CC")
    out = ConjugatedSystems()._restore_bonds(
        mol, {frozenset([0, 1]): BondType.DOUBLE}
    )
    assert out.GetBondBetweenAtoms(0, 1).GetBondType() == BondType.DOUBLE


# ---------------------------------------------------------------------------
# QueryMol / EditMol
# ---------------------------------------------------------------------------


def test_bfs_atom_path_kekulize_failure():
    with patch("xenosite.forest.base.Kekulize", side_effect=ValueError):
        QueryMol()._bfs_atom_path(MolFromSmiles("c1ccccc1"), 0, 1)


def test_match_requires_query_or_smarts_and_builds_mapids():
    mol = MolFromSmiles("CCO")
    with pytest.raises(NotImplementedError, match="query mol or smarts"):
        next(EditMol().match(mol))

    query = MolFromSmarts("[#6:1][#8:2]")
    hits = list(EditMol().match(mol, query=query, queryidx2mapid=None))
    assert hits
    assert {1, 2}.issubset(hits[0])


def test_smarts_match_standardize_failure_and_bad_query():
    bad = MolFromSmiles("C#C(C)C", sanitize=False)
    assert QueryMol.smarts_match(bad, "[#6]") is None

    class _Mol:
        def HasProp(self, name):
            return True

        def GetSubstructMatches(self, query):
            raise RuntimeError("bad query")

    with pytest.raises(ValueError, match="Problem preparing query"):
        QueryMol.smarts_match(_Mol(), "[#6]")

    assert QueryMol.smarts_match("C=CC", ["[#6]=[#6]", "[#6]-[#6]"])
    mapped = QueryMol.smarts_match(
        "C=CC", {"S1": "[#6]=[#6]", "S2": "[#6]-[#6]"}
    )
    assert mapped["S1"] == frozenset({0, 1})
    assert mapped["S2"] == frozenset({1, 2})


def test_prepare_query_rejects_invalid_smarts():
    with pytest.raises(ValueError, match="Problem with this pattern"):
        EditMol()._prepare_query("[[[")


def test_match_queries_from_smiles_and_standardize_failure():
    editor = EditMol(query_smarts=[("alcohol", "[#6:1][#8:2]")])
    hits = editor.match_queries("CCO")
    assert hits

    bad = MolFromSmiles("C#C(C)C", sanitize=False)
    assert editor.match_queries(bad) == {}


def test_bonds_from_atom_path_runtime_error():
    class _Mol:
        def GetBondBetweenAtoms(self, *args):
            raise RuntimeError("bad idx")

    assert list(EditMol()._bonds_from_atom_path(_Mol(), [0, 1, 2])) == []


def test_swap_bonds_attach_path():
    mol = MolFromSmiles("C=CC=CCCl")
    EditMol().swap_bonds_along_path(mol, atoms=[0, 1, 2, 3, 4], attach_path=True)
    assert mol.GetProp("path") == "[0, 1, 2, 3, 4]"
    assert MolToSmiles(mol) == "CC=CC=CCl"


def test_editmol_modify_not_implemented():
    with pytest.raises(NotImplementedError):
        EditMol().modify(None, None, None)


# ---------------------------------------------------------------------------
# Resonate
# ---------------------------------------------------------------------------


def test_resfrags_resonance_supplier_fallbacks_and_none():
    mol = MolFromSmiles("c1ccccc1")
    with patch("xenosite.forest.base.ResonanceMolSupplier", side_effect=ValueError):
        frags = list(Resonate()._resfrags(mol))
    assert frags

    class _Supplier:
        def GetAtomConjGrpIdx(self, idx):
            return 0

        def GetNumConjGrps(self):
            return 1

        def __iter__(self):
            yield None
            yield mol

    with patch("xenosite.forest.base.ResonanceMolSupplier", return_value=_Supplier()):
        kept = list(Resonate()._resfrags(mol))
    assert kept
    assert all(item[0] is not None for item in kept)


# ---------------------------------------------------------------------------
# ReactionRule / SmartsReactionRule
# ---------------------------------------------------------------------------


def test_reaction_rule_longname_and_underscore_forbidden():
    rule = ReactionRule(name="Hydroxylation", longname="C-H hydroxylation")
    assert rule.longname == "C-H hydroxylation"
    with pytest.raises(AssertionError, match="Cannot have '_'"):
        ReactionRule(name="bad_name")


def test_reaction_rule_call_and_format_site():
    mol = MolFromSmiles("CCO")
    sites = list(rules.Hydroxylation()(mol))
    assert sites

    rule = rules.Hydroxylation()
    assert rule.format_site("Hydroxylation_SmartsReactionRuleRxn0") == "Hydroxylation"
    assert (
        rule.format_site(("Hydroxylation_x", frozenset({0})), just_rule_name=True)
        == "Hydroxylation_x"
    )
    formatted_list = rule.format_site(["Hydroxylation_x", "Dehydration_y"])
    assert formatted_list == ["Hydroxylation", "Dehydration"]
    assert rule.format_site(frozenset({0}), just_rule_name=True) == frozenset({0})


def test_metabolize_only_unique_largest_fragment_and_topol_retry():
    propane = MolFromSmiles("CCC")
    unique = list(
        rules.Hydroxylation().metabolize(propane, only_unique=True)
    )
    smiles = [can_smi_set(prods) for _, prods in unique]
    assert len(smiles) == len(set(smiles))

    dealk = list(
        rules.Dealkylation().metabolize(
            MolFromSmiles("CCN"), only_largest_fragment=True
        )
    )
    assert dealk
    for _, prods in dealk:
        assert len(prods) == 1

    real = AtomTracker.topol_equiv
    n = {"i": 0}

    def flaky(mol):
        n["i"] += 1
        if n["i"] == 2:
            raise RuntimeError("rank failed")
        return real(mol)

    with patch.object(AtomTracker, "topol_equiv", side_effect=flaky):
        out = list(rules.Hydroxylation().metabolize(MolFromSmiles("CCO")))
    assert out
    assert n["i"] >= 3


def test_reaction_rule_metabolites_not_implemented():
    with pytest.raises(NotImplementedError):
        ReactionRule(name="Dummy").metabolites(MolFromSmiles("CC"))


def test_metabolites_from_sites_deplete_and_cast_shapes():
    rule = rules.Hydroxylation()
    mol = MolFromSmiles("CCC")

    depleted = list(
        rule.metabolites_from_sites(
            mol, [frozenset({0}), frozenset({1})], deplete_sites=True
        )
    )
    assert depleted
    assert len(depleted) <= 2

    just = list(rule.metabolites_from_sites(MolFromSmiles("CCO"), 0, just_smiles=True))
    assert just
    assert all(isinstance(row, list) for row in just)

    unformatted = list(
        rules.Hydroxylation().metabolize(MolFromSmiles("CCO"), format_output_site=False)
    )
    assert unformatted
    assert unformatted[0][0][0].startswith("Hydroxylation")

    assert rule._cast_sites(0) == [("Hydroxylation", frozenset([0]))]
    assert rule._cast_sites((0, 1)) == [("Hydroxylation", frozenset((0, 1)))]
    assert rule._cast_sites(frozenset({0})) == [("Hydroxylation", frozenset({0}))]
    assert rule._cast_sites([0, 1]) == [("Hydroxylation", frozenset([0, 1]))]
    assert rule._cast_sites([["Hydroxylation", frozenset({0})]]) == [
        ("Hydroxylation", frozenset({0}))
    ]
    as_tuples = [("Hydroxylation", frozenset({1}))]
    assert rule._cast_sites(as_tuples) is as_tuples
    with pytest.raises(NotImplementedError, match="integer, a frozenset, or a list"):
        rule._cast_sites("nope")


def test_kekulize_copy_props_default_and_remove_props_list():
    mol = MolFromSmiles("c1ccccc1")
    rule = SmartsReactionRule(name="Dummy", rxns=["[#6:1]>>[#6:1]"])
    with patch("xenosite.forest.base.Kekulize", side_effect=ValueError):
        rule._kekulize(mol)

    copied = MolFromSmiles("CC")
    copied.GetAtomWithIdx(0).SetProp("react_atom_idx", "3")
    rule._copy_props(copied, props_to_copy=None)
    assert copied.GetAtomWithIdx(0).GetProp(AtomTracker.previous_index_prop_name) == "3"

    a = MolFromSmiles("CC")
    b = MolFromSmiles("CCO")
    a.GetAtomWithIdx(0).SetProp("old_mapno", "1")
    rule._remove_props([a, b])
    assert not a.GetAtomWithIdx(0).HasProp("old_mapno")


@pytest.mark.filterwarnings("ignore:.*found in sys.modules.*:RuntimeWarning")
def test_base_module_main_runs_doctest():
    import runpy

    with patch("doctest.testmod") as mocked:
        runpy.run_module("xenosite.forest.base", run_name="__main__")
    mocked.assert_called()
