"""Cache identity for the RDKit door. Callers read NamedTuple attributes."""

from typing import TYPE_CHECKING, cast

from rdkit import Chem

from xenosite.refactor_poc.rdkitutil import (
    copy_mol,
    ensure_forest,
    mcs_matches,
    resonance_bond_maps,
    rw_copy,
    sanitized_fragments,
    split_fragments,
)

if TYPE_CHECKING:
    from xenosite.refactor_poc.rdkit_api import (
        ForestMol,
        ForestNoTracingMol,
        NoForestMol,
        NoTracingMol,
        TracingMol,
    )

    def _needs_no_trace(mol: NoTracingMol) -> None:
        """A missing trace, whether or not a forest is present."""

    def _needs_forest(mol: ForestMol) -> None:
        """A forest, whether or not its trace is initialized."""

    def _subclass_relationships(
        bare: NoForestMol,
        untraced: ForestNoTracingMol,
        traced: TracingMol,
    ) -> None:
        """No-forest is no-trace. Forest-without-trace is both. Traced is a forest."""

        _needs_no_trace(bare)
        _needs_no_trace(untraced)
        _needs_forest(untraced)
        _needs_forest(traced)


def test_cached_answers_are_the_same_object():
    mol = ensure_forest(Chem.MolFromSmiles("CCC"))
    assert mol.xf.topol_equiv is mol.xf.topol_equiv
    assert mol.xf.csmi is mol.xf.csmi
    smarts = "[C:1][C:2]"
    assert mol.xf.smarts_matches(smarts) is mol.xf.smarts_matches(smarts)


def test_mcs_matches_cache_on_the_reactant():
    reactant = Chem.MolFromSmiles("c1ccccc1")
    target = ensure_forest(Chem.MolFromSmiles("c1ccccc1"))
    first = mcs_matches(reactant, target)
    second = mcs_matches(reactant, target)
    assert first is second
    assert len(first.embeddings) > 1
    held_reactant = ensure_forest(reactant)
    structure = held_reactant.xf.forest.get("cache")
    assert structure is not None
    held = (structure.get("mcs_matches") or {})[target.xf.csmi]
    assert held is first
    assert isinstance(held.embeddings, tuple)


def test_a_fresh_mol_does_not_reuse_parent_caches():
    parent = ensure_forest(Chem.MolFromSmiles("c1ccccc1O"))
    parent_csmi = parent.xf.csmi
    parent_maps = resonance_bond_maps(parent)

    product = Chem.Mol(parent)
    fresh = ensure_forest(product)
    structure = fresh.xf.forest.get("cache")
    assert structure is not None
    assert "csmi" not in structure
    assert "resonance_bonds" not in structure
    parent_structure = parent.xf.forest.get("cache")
    assert structure is not parent_structure
    assert fresh.xf.csmi == parent_csmi
    assert resonance_bond_maps(product) is not parent_maps

    child = copy_mol(parent)
    child_structure = ensure_forest(child).xf.forest.get("cache")
    assert child_structure is not None
    assert child_structure is parent_structure
    assert child_structure.get("csmi") == parent_csmi


def test_forest_stays_a_dict():
    """Forest is a TypedDict schema and still a plain dict at runtime."""

    from xenosite.refactor_poc.forest_copy import forest_copy

    mol = ensure_forest(Chem.MolFromSmiles("CC"))
    forest = mol.xf.forest
    # Intentional non-schema writes: the forest must remain a mutable dict.
    raw = cast(dict, forest)
    structure = raw.setdefault("cache", {})
    structure["csmi"] = "CC"
    raw["atom_trace"] = {"depth": 0}
    raw["not_a_schema_key"] = 1
    assert forest.get("cache", {}).get("csmi") == "CC"
    assert forest.get("atom_trace", {}).get("depth") == 0
    assert raw["not_a_schema_key"] == 1
    assert isinstance(forest, dict)

    formula = mol.xf.formula
    assert isinstance(formula, dict)
    assert formula["counts"]["C"] == 2
    assert formula["charge"] == 0
    assert formula is mol.xf.formula

    cloned = forest_copy(forest, same_structure=True)
    assert cloned.get("cache", {}).get("formula", {}).get("counts", {}).get("C") == 2
    assert cloned.get("atom_trace", {}).get("depth") == 0
    assert cloned["cache"] is forest["cache"]


def test_rw_copy_does_not_carry_the_structure_cache():
    parent = ensure_forest(Chem.MolFromSmiles("CCO"))
    _ = parent.xf.csmi
    child = rw_copy(parent)
    assert getattr(child, "_forest", None) is None
    copied = copy_mol(parent).xf.forest.get("cache")
    assert copied is not None
    assert copied is parent.xf.forest["cache"]
    assert copied.get("csmi") == parent.xf.csmi


def test_forestmol_bridge_and_wipe_honesty():
    """``xf.forestmol`` attaches; ``wipe_forest`` returns NoForestMol brand."""

    mol = Chem.MolFromSmiles("CC")
    assert mol.xf.has_forest is False
    held = mol.xf.forestmol
    assert held is mol and mol.xf.has_forest
    assert held.xf.forest is held._forest
    from xenosite.refactor_poc.rdkitutil import wipe_forest
    wiped = wipe_forest(held)
    assert wiped is mol and mol.xf.has_forest is False


def test_fragment_split_uses_pieces():
    mol = Chem.MolFromSmiles("CCO")
    split = split_fragments(mol)
    assert split.pieces == (mol,)
    fragments = sanitized_fragments(mol)
    assert len(fragments.pieces) == 1


def test_carry_forest_drops_sibling_records_and_clears_structure():
    from xenosite.refactor_poc.rdkitutil import GetMolFrags, carry_forest
    from xenosite.refactor_poc.rules import stamp_forest_labels

    parent = stamp_forest_labels(Chem.MolFromSmiles("C.O"))
    parent._forest["cache"]["csmi"] = "stale"
    frags = list(GetMolFrags(parent, asMols=True, sanitizeFrags=False))
    assert len(frags) == 2
    a = carry_forest(parent, frags[0])
    b = carry_forest(parent, frags[1])
    assert a._forest is not b._forest
    assert a._forest is not parent._forest
    assert a._forest["cache"] == {}
    assert b._forest["cache"] == {}
    for frag in (a, b):
        live = set(frag._forest["atom_trace"]["records"])
        on_mol = {
            atom.GetProp("forestLabel")
            for atom in frag.GetAtoms()
            if atom.HasProp("forestLabel")
        }
        assert live == on_mol
        # Siblings are absent from this fragment's live records (split-before-trace).
        assert not (set(parent._forest["atom_trace"]["records"]) - live) & live
