"""ForestMol.xf facade: mint-on-read, forestmol bridge, tracing, terminal conjugates."""

import gc
import weakref

import pytest
from rdkit import Chem

from xenosite.forest.rdkitutil import (
    Xf,
    XfTracing,
    copy_mol,
    is_forest,
    wipe_forest,
)
from xenosite.forest.rules import Acetylation, Hydroxylation


def test_xf_is_minted_per_access_with_strong_ref():
    mol = Chem.MolFromSmiles("CCO")
    a = mol.xf
    b = mol.xf
    assert isinstance(a, Xf)
    assert a is not b
    assert a.mol is mol is b.mol
    assert Xf.__slots__ == ("_mol",)
    assert not hasattr(a, "__dict__")
    with pytest.raises(AttributeError, match="read-only"):
        mol.xf = a  # type: ignore[misc]


def test_xf_keeps_strong_ref_through_temporary_chain():
    csmi = Chem.MolFromSmiles("c1ccccc1O").xf.csmi
    assert csmi
    mol = Chem.MolFromSmiles("CC")
    xf = mol.xf
    ref = weakref.ref(mol)
    del mol
    gc.collect()
    assert ref() is not None  # kept alive by xf
    assert xf.csmi


def test_tracing_facade_is_minted_not_flat_on_xf():
    mol = Chem.MolFromSmiles("CCO")
    tracing = mol.xf.tracing
    assert isinstance(tracing, XfTracing)
    assert tracing.mol is mol
    assert not hasattr(mol.xf, "stamp")
    assert not hasattr(mol.xf, "trace")
    assert not hasattr(mol.xf, "depth")


def test_has_forest_and_forestmol_bridge():
    mol = Chem.MolFromSmiles("CCO")
    assert mol.xf.has_forest is False
    held = mol.xf.forestmol
    assert held is mol
    assert is_forest(held)
    assert mol.xf.has_forest is True
    wipe_forest(mol)
    assert mol.xf.has_forest is False
    assert mol.xf.forestmol is mol and mol.xf.has_forest


def test_copy_mol_does_not_share_xf_identity():
    parent = Chem.MolFromSmiles("CCO")
    _ = parent.xf.csmi
    child = copy_mol(parent)
    assert child.xf.mol is child
    assert child.xf.mol is not parent
    assert child.xf.has_forest


def test_csmi_on_demand_via_xf():
    mol = Chem.MolFromSmiles("c1ccccc1O")
    assert not mol.xf.has_forest
    first = mol.xf.csmi
    assert mol.xf.forest["cache"]["csmi"] == first
    assert mol.xf.csmi is first


def test_rings_and_conjugated_systems_on_xf():
    mol = Chem.MolFromSmiles("c1ccccc1")
    assert 0 in mol.xf.rings
    assert len(mol.xf.conjugated_systems) >= 1
    assert len(mol.xf.aromatic_systems) >= 1


def test_stamp_installs_tracing_without_exposing_labels():
    mol = Chem.MolFromSmiles("CCO")
    assert not mol.xf.tracing.active
    stamped = mol.xf.tracing._stamp()
    assert stamped is mol
    assert mol.xf.tracing.active
    assert mol.xf.tracing.depth == 0
    assert mol.xf.tracing.atom_origin(0) == 0
    assert "forestLabel" not in dir(mol.xf.tracing)


def test_tracing_migration_helpers_on_hydroxylation_product():
    """depths / added_indices / root_map / index_at replace AtomTracker patterns."""

    products, _info = next(Hydroxylation().metabolize(Chem.MolFromSmiles("CCO")))
    product = products[0]
    t = product.xf.tracing
    assert t.active
    assert t.depth == 1
    assert t.depths() == (0, 1)
    added = t.added_indices()
    assert added
    for i in added:
        assert t.atom_root(i) is None
        assert t.atom_added_by(i) is not None
        assert t.index_at(i, 0) is None
        assert t.index_at(i, 1) == i
    roots = t.root_map()
    assert roots
    assert not (set(roots) & added)
    for cur, root in roots.items():
        assert t.index_at(cur, 0) == root


def test_of_products_single_and_list():
    reactant = Chem.MolFromSmiles("CCO")
    reactant = reactant.xf.tracing._stamp()
    rule = Hydroxylation()
    work = copy_mol(reactant).xf.tracing._stamp()
    rule._clear_atom_maps(work)
    pieces = []
    site_info = None
    for por in rule.metabolites(work):
        site_info = por.info
        pieces.extend(por.products)
        if len(pieces) >= 2:
            break
    assert pieces and site_info is not None

    finished_one = work.xf._of_products(pieces[0], site_info)
    assert len(finished_one) == 1
    assert finished_one[0].xf.tracing.active
    assert finished_one[0].xf.tracing.depth == 1
    assert getattr(pieces[0], "_forest", None) is None

    finished_many = work.xf._of_products(pieces[:2], site_info)
    assert all(p.xf.tracing.active for p in finished_many)


def test_metabolize_uses_of_products_and_xf_csmi():
    reactant = Chem.MolFromSmiles("CCO").xf.tracing._stamp()
    products = list(Hydroxylation().metabolize(reactant))
    assert products
    product_list, info = products[0]
    product = product_list[0]
    assert product.xf.tracing.active
    assert product.xf.tracing.depth == 1
    assert "csmi" not in info
    assert product.xf.csmi


def test_conjugation_products_are_terminal_and_not_reexpanded():
    phenol = Chem.MolFromSmiles("c1ccccc1O").xf.tracing._stamp()
    products = list(Acetylation(as_star=True).metabolize(phenol))
    assert products
    for product_list, _info in products:
        for product in product_list:
            assert product.xf.is_terminal
            assert list(Acetylation(as_star=True).metabolize(product)) == []
            assert list(Hydroxylation().metabolize(product)) == []
