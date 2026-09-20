"""Cache identity for the RDKit door. Callers read NamedTuple attributes."""

import copy
from typing import TYPE_CHECKING, cast

from rdkit import Chem

from xenosite.refactor_poc.rdkitutil import (
    copy_mol,
    ensure_forest,
    get_csmi,
    get_forest,
    mcs_matches,
    molecule_formula,
    resonance_bond_maps,
    rw_copy,
    sanitized_fragments,
    smarts_matches,
    split_fragments,
    topol_equiv,
)

if TYPE_CHECKING:
    from xenosite.refactor_poc.rdkit_api import (
        ForestMol,
        ForestNoTracingMol,
        ForestTracingMol,
        NoForestMol,
        NoTracingMol,
    )

    def _needs_no_trace(mol: NoTracingMol) -> None:
        """A missing trace, whether or not a forest is present."""

    def _needs_forest(mol: ForestMol) -> None:
        """A forest, whether or not its trace is initialized."""

    def _subclass_relationships(
        bare: NoForestMol,
        untraced: ForestNoTracingMol,
        traced: ForestTracingMol,
    ) -> None:
        """No-forest is no-trace. Forest-without-trace is both. Traced is a forest."""

        _needs_no_trace(bare)
        _needs_no_trace(untraced)
        _needs_forest(untraced)
        _needs_forest(traced)


def test_cached_answers_are_the_same_object():
    mol = Chem.MolFromSmiles("CCC")
    assert topol_equiv(mol) is topol_equiv(mol)
    assert get_csmi(mol) is get_csmi(mol)
    smarts = "[C:1][C:2]"
    assert smarts_matches(mol, smarts) is smarts_matches(mol, smarts)


def test_mcs_matches_cache_on_the_reactant():
    reactant = Chem.MolFromSmiles("c1ccccc1")
    target = Chem.MolFromSmiles("c1ccccc1")
    first = mcs_matches(reactant, target)
    second = mcs_matches(reactant, target)
    assert first is second
    assert len(first.embeddings) > 1
    structure = get_forest(ensure_forest(reactant)).get("structure")
    assert structure is not None
    held = (structure.get("mcs_matches") or {})[get_csmi(target)]
    assert held is first
    assert isinstance(held.embeddings, tuple)


def test_a_fresh_mol_does_not_reuse_parent_caches():
    parent = Chem.MolFromSmiles("c1ccccc1O")
    parent_csmi = get_csmi(parent)
    parent_maps = resonance_bond_maps(parent)

    product = Chem.Mol(parent)
    fresh = ensure_forest(product)
    structure = get_forest(fresh).get("structure")
    assert structure is not None
    assert "csmi" not in structure
    assert "resonance_bonds" not in structure
    parent_structure = get_forest(ensure_forest(parent)).get("structure")
    assert structure is not parent_structure
    assert get_csmi(product) == parent_csmi
    assert resonance_bond_maps(product) is not parent_maps

    child = copy_mol(parent)
    child_structure = get_forest(ensure_forest(child)).get("structure")
    assert child_structure is not None
    assert child_structure.get("csmi") == parent_csmi


def test_forest_stays_a_dict():
    """Forest is a TypedDict schema and still a plain dict at runtime."""

    mol = Chem.MolFromSmiles("CC")
    forest = get_forest(ensure_forest(mol))
    # Intentional non-schema writes: the forest must remain a mutable dict.
    raw = cast(dict, forest)
    structure = raw.setdefault("structure", {})
    structure["csmi"] = "CC"
    raw["atom_trace"] = {"depth": 0}
    raw["not_a_schema_key"] = 1
    assert forest.get("structure", {}).get("csmi") == "CC"
    assert forest.get("atom_trace", {}).get("depth") == 0
    assert raw["not_a_schema_key"] == 1
    assert isinstance(forest, dict)

    formula = molecule_formula(mol)
    assert isinstance(formula, dict)
    assert formula["counts"]["C"] == 2
    assert formula["charge"] == 0
    assert formula is molecule_formula(mol)

    cloned = copy.deepcopy(forest)
    assert cloned.get("structure", {}).get("formula", {}).get("counts", {}).get("C") == 2
    assert cloned.get("atom_trace", {}).get("depth") == 0


def test_rw_copy_does_not_carry_the_structure_cache():
    parent = Chem.MolFromSmiles("CCO")
    get_csmi(parent)
    child = rw_copy(parent)
    assert getattr(child, "_forest", None) is None
    copied = get_forest(copy_mol(parent)).get("structure")
    assert copied is not None
    assert copied.get("csmi") == get_csmi(parent)


def test_fragment_split_uses_pieces():
    mol = Chem.MolFromSmiles("CCO")
    split = split_fragments(mol)
    assert split.pieces == (mol,)
    fragments = sanitized_fragments(mol)
    assert len(fragments.pieces) == 1
