"""A rule may stamp a copy. It must not edit the chemistry it was given.

The parent gains a ``_forest`` when it had none. Every metabolite carries a
``_forest`` whose trace depth is one greater than that parent.
"""

from rdkit import Chem

from xenosite.refactor_poc.rules import (
    Dealkylation,
    Epoxidation,
    Hydroxylation,
    install_forest,
    molecule_formula,
)
from xenosite.refactor_poc.rulesets import RuleSet


def _chemistry(mol):
    return (
        Chem.MolToSmiles(mol),
        tuple(
            (
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
                atom.GetAtomMapNum(),
                atom.GetTotalNumHs(),
            )
            for atom in mol.GetAtoms()
        ),
        tuple(
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), str(bond.GetBondType()))
            for bond in mol.GetBonds()
        ),
    )


def _assert_parent_and_products(parent, products, parent_depth):
    trace = parent._forest["atom_trace"]
    assert trace["depth"] == parent_depth
    assert "formula" in trace
    for product, _info in products:
        child = product._forest["atom_trace"]
        assert child["depth"] == parent_depth + 1
        assert child["formula"] == molecule_formula(product)


def _assert_untouched(rule, smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert getattr(mol, "_forest", None) is None
    before = _chemistry(mol)
    products = list(rule.metabolize(mol))
    assert _chemistry(mol) == before
    assert products
    _assert_parent_and_products(mol, products, 0)
    return products


def test_hydroxylation_does_not_edit_its_input():
    _assert_untouched(Hydroxylation(), "CC")


def test_dealkylation_does_not_edit_its_input():
    _assert_untouched(Dealkylation(), "COc1ccccc1")


def test_epoxidation_does_not_edit_its_input():
    _assert_untouched(Epoxidation(), "C=CC")


def test_existing_parent_forest_is_kept():
    mol = Chem.MolFromSmiles("CC")
    mol = install_forest(mol)
    mol._forest["atom_trace"]["depth"] = 2
    before = _chemistry(mol)
    products = list(Hydroxylation().metabolize(mol))
    assert _chemistry(mol) == before
    assert mol._forest["atom_trace"]["depth"] == 2
    _assert_parent_and_products(mol, products, 2)


def test_ruleset_does_not_edit_its_input():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Poc")
    products = _assert_untouched(ruleset, "CC")
    product, info = products[0]
    trace = product._forest["atom_trace"]
    added = [
        by
        for record in trace["records"].values()
        if (by := record.get("added_by"))
    ]
    assert added
    transform_id = added[0]
    assert transform_id.startswith("R")
    addition = trace["additions"][transform_id]
    assert addition["name"] in {"Hydroxylation", "Dealkylation"}
    names = tuple(getattr(item, "name", item) for item in addition["rules"])
    assert "Poc" in names
    assert addition["name"] in names
    assert addition["depth"] == 0
    assert "site" in addition and "effect" in addition
    assert trace["formula"]["counts"]["C"] >= 1
    assert "H" in trace["formula"]["counts"]
    assert transform_id in trace["delta_formula"]
    assert trace["delta_formula"][transform_id]["counts"]
    assert trace["formula"] == molecule_formula(product)
