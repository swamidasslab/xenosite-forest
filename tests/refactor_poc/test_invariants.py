"""What :meth:`ReactionRule.metabolize` promises.

The same sentences are in the docstring of :meth:`ReactionRule.metabolize`
and :meth:`ReactionRule.metabolites`. A failure here means the docstring
and the library have drifted.
"""

from rdkit import Chem

from xenosite.refactor_poc.rules import (
    Dealkylation,
    Epoxidation,
    Hydroxylation,
    ReactionRule,
    install_forest,
    molecule_formula,
)
from xenosite.refactor_poc.rulesets import RuleSet

_ADDITION_FIELDS = ("site", "rules", "info", "effect", "name", "phase1", "depth")


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


def _pairs(rule, smiles="CC"):
    mol = Chem.MolFromSmiles(smiles)
    before = _chemistry(mol)
    products = list(rule.metabolize(mol))
    return mol, before, products


def test_input_chemistry_is_unchanged():
    rule = RuleSet((Hydroxylation, Dealkylation, Epoxidation), name="Poc")
    mol, before, products = _pairs(rule, "C=CC")
    assert products
    assert _chemistry(mol) == before


def test_a_parent_with_no_forest_keeps_one_afterward():
    mol, _before, _products = _pairs(Hydroxylation(), "CC")
    assert mol._forest["atom_trace"]["depth"] == 0
    assert "formula" in mol._forest["atom_trace"]


def test_an_existing_parent_depth_is_not_reset():
    mol = Chem.MolFromSmiles("CC")
    install_forest(mol)
    mol._forest["atom_trace"]["depth"] = 4
    products = list(Hydroxylation().metabolize(mol))
    assert mol._forest["atom_trace"]["depth"] == 4
    assert products[0][0]._forest["atom_trace"]["depth"] == 5


def test_each_product_is_one_depth_below_its_parent():
    mol, _before, products = _pairs(Hydroxylation(), "CC")
    parent_depth = mol._forest["atom_trace"]["depth"]
    for product, info in products:
        assert isinstance(product, Chem.Mol)
        assert isinstance(info, dict)
        child = product._forest["atom_trace"]
        assert child["depth"] == parent_depth + 1
        assert child["formula"] == molecule_formula(product)
        assert info["csmi"] == Chem.MolToSmiles(product, isomericSmiles=False)
        assert info["rule"].name == "Hydroxylation"


def test_canonical_smiles_are_not_repeated():
    _mol, _before, products = _pairs(
        RuleSet((Hydroxylation, Dealkylation), name="Poc"), "CC"
    )
    smiles = [info["csmi"] for _product, info in products]
    assert len(smiles) == len(set(smiles))


def test_a_new_atom_points_at_one_addition_record():
    mol, _before, products = _pairs(
        RuleSet((Hydroxylation,), name="Poc"), "CC"
    )
    parent_formula = mol._forest["atom_trace"]["formula"]
    product, _info = products[0]
    trace = product._forest["atom_trace"]
    ids = [
        record["added_by"]
        for record in trace["records"].values()
        if record.get("added_by")
    ]
    assert ids
    transform_id = ids[0]
    assert transform_id in trace["additions"]
    assert transform_id in trace["delta_formula"]
    addition = trace["additions"][transform_id]
    assert set(_ADDITION_FIELDS) <= set(addition)
    names = tuple(getattr(item, "name", item) for item in addition["rules"])
    assert names == ("Poc", "Hydroxylation")
    assert addition["name"] == "Hydroxylation"
    assert addition["depth"] == 0
    delta = trace["delta_formula"][transform_id]
    for element, change in delta["counts"].items():
        assert (
            trace["formula"]["counts"].get(element, 0)
            - parent_formula["counts"].get(element, 0)
            == change
        )
    assert (
        trace["formula"]["charge"] - parent_formula["charge"] == delta["charge"]
    )


def test_refusing_the_rule_edits_nothing_and_still_instruments_the_parent():
    mol = Chem.MolFromSmiles("CC")
    before = _chemistry(mol)
    products = list(
        Hydroxylation().metabolize(mol, filter_rules=lambda rule, info: False)
    )
    assert products == []
    assert _chemistry(mol) == before
    assert mol._forest["atom_trace"]["depth"] == 0


def test_refusing_every_site_edits_nothing():
    mol = Chem.MolFromSmiles("CC")
    before = _chemistry(mol)
    products = list(
        Hydroxylation().metabolize(mol, filter_sites=lambda site, info: False)
    )
    assert products == []
    assert _chemistry(mol) == before
    assert mol._forest["atom_trace"]["depth"] == 0


def test_a_rule_is_itself_when_iterated_and_its_name_has_no_underscore():
    rule = Hydroxylation()
    assert isinstance(rule, ReactionRule)
    assert list(rule) == [rule]
    assert "_" not in rule.name
    called = list(rule(Chem.MolFromSmiles("CC")))
    direct = list(rule.metabolize(Chem.MolFromSmiles("CC")))
    assert [info["csmi"] for _product, info in called] == [
        info["csmi"] for _product, info in direct
    ]


def test_a_ruleset_iterates_its_children():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Poc")
    assert [type(rule) for rule in ruleset] == [Hydroxylation, Dealkylation]
    assert "_" not in ruleset.name
