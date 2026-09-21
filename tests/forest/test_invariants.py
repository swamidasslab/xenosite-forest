"""What :meth:`ReactionRule.metabolize` promises.

The same sentences are in the docstring of :meth:`ReactionRule.metabolize`
and :meth:`ReactionRule.metabolites`. A failure here means the docstring
and the library have drifted.
"""

from rdkit import Chem

from xenosite.forest.rules import (
    Dealkylation,
    Epoxidation,
    Hydroxylation,
    ReactionRule,
)
from xenosite.forest.rulesets import RuleSet

_ADDITION_FIELDS = ("site", "rules", "info", "effect", "name", "phase1", "depth", "pattern")


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
    rule = RuleSet((Hydroxylation, Dealkylation, Epoxidation), name="Forest")
    mol, before, products = _pairs(rule, "C=CC")
    assert products
    assert _chemistry(mol) == before


def test_a_parent_with_no_forest_keeps_one_afterward():
    mol, _before, _products = _pairs(Hydroxylation(), "CC")
    assert mol._forest["atom_trace"]["depth"] == 0
    assert "formula" in mol._forest["atom_trace"]


def test_an_existing_parent_depth_is_not_reset():
    mol = Chem.MolFromSmiles("CC")
    mol = mol.xf.tracing._stamp()
    mol._forest["atom_trace"]["depth"] = 4
    products = list(Hydroxylation().metabolize(mol))
    assert mol._forest["atom_trace"]["depth"] == 4
    assert products[0][0][0]._forest["atom_trace"]["depth"] == 5


def test_each_product_is_one_depth_below_its_parent():
    mol, _before, products = _pairs(Hydroxylation(), "CC")
    parent_depth = mol._forest["atom_trace"]["depth"]
    for product_list, info in products:
        assert isinstance(product_list, list)
        assert len(product_list) == 1
        product = product_list[0]
        assert isinstance(product, Chem.Mol)
        assert isinstance(info, dict)
        child = product._forest["atom_trace"]
        assert child["depth"] == parent_depth + 1
        assert child["formula"] == product.xf.formula
        assert product.xf.csmi == Chem.MolToSmiles(product, isomericSmiles=False)
        assert "csmi" not in info
        rules = info["rule"]
        assert isinstance(rules, list)
        assert len(rules) == 1
        assert isinstance(rules[0], ReactionRule)
        assert rules[0].name == "Hydroxylation"


def test_canonical_smiles_are_not_repeated():
    """Emitted ``_unique_csmi_key`` values are unique under RuleSet metabolize.

    Per-rule keys still include pattern tokens. RuleSet drops only when a
    *different* child rule repeats an emission frozenset (INFO log, yield drop).
    """

    from xenosite.forest.rules import _unique_csmi_key

    _mol, _before, products = _pairs(
        RuleSet((Hydroxylation, Dealkylation), name="Forest"), "CC"
    )
    keys = [
        _unique_csmi_key(info, frozenset(p.xf.csmi for p in products))
        for products, info in products
    ]
    assert len(keys) == len(set(keys))


def test_a_new_atom_points_at_one_addition_record():
    mol, _before, products = _pairs(
        RuleSet((Hydroxylation,), name="Forest"), "CC"
    )
    parent_formula = mol._forest["atom_trace"]["formula"]
    product_list, _info = products[0]
    product = product_list[0]
    trace = product._forest["atom_trace"]
    ids = [
        by
        for record in trace["records"].values()
        if (by := record.get("added_by"))
    ]
    assert ids
    transform_id = ids[0]
    assert transform_id in trace["additions"]
    assert transform_id in trace["delta_formula"]
    addition = trace["additions"][transform_id]
    assert set(_ADDITION_FIELDS) <= set(addition)
    names = tuple(getattr(item, "name", item) for item in addition["rules"])
    assert names == ("Hydroxylation", "Forest")
    assert addition["name"] == "Hydroxylation"
    # Ethane carbons are h=3 → pattern ``h2`` (partitioned from ``h`` / h1).
    assert addition["pattern"] is Hydroxylation.smirks[1][1]
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
        Hydroxylation().metabolize(mol, filter_rules=lambda mol, rule, info: False)
    )
    assert products == []
    assert _chemistry(mol) == before
    assert mol._forest["atom_trace"]["depth"] == 0


def test_refusing_every_site_edits_nothing():
    mol = Chem.MolFromSmiles("CC")
    before = _chemistry(mol)
    products = list(
        Hydroxylation().metabolize(mol, filter_sites=lambda mol, site, info: False)
    )
    assert products == []
    assert _chemistry(mol) == before
    assert mol._forest["atom_trace"]["depth"] == 0


def test_a_rule_is_itself_when_iterated_and_its_name_has_no_underscore():
    rule = Hydroxylation()
    assert isinstance(rule, ReactionRule)
    assert list(rule) == [rule]
    assert rule.name is not None and "_" not in rule.name
    called = list(rule(Chem.MolFromSmiles("CC")))
    direct = list(rule.metabolize(Chem.MolFromSmiles("CC")))
    assert [
        frozenset(p.xf.csmi for p in products) for products, _info in called
    ] == [frozenset(p.xf.csmi for p in products) for products, _info in direct]


def test_a_ruleset_iterates_its_children():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Forest")
    assert [type(rule) for rule in ruleset] == [Hydroxylation, Dealkylation]
    assert ruleset.name is not None and "_" not in ruleset.name
