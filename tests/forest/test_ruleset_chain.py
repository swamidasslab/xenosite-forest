"""A ruleset stays on the chain after the rule that emitted the product."""

from rdkit import Chem

from xenosite.forest.phaseone import PhaseOneRS, StableOxygenation_PhaseOne
from xenosite.forest.rules import Hydroxylation, ReactionRule
from xenosite.forest.rulesets import RuleSet


def _addition(product):
    return product._forest["atom_trace"]["additions"]["R1"]


def _rules(product):
    return _addition(product)["rules"]


def test_leaf_info_rule_is_a_list_containing_the_emitting_rule():
    """``info["rule"]`` is ``list[ReactionRule]`` so callers can walk ``.name``."""

    leaf = Hydroxylation()
    _products, info = next(leaf.metabolize(Chem.MolFromSmiles("CC")))
    rules = info["rule"]
    assert isinstance(rules, list)
    assert rules == [leaf]
    assert all(isinstance(rule, ReactionRule) for rule in rules)
    assert [rule.name for rule in rules] == ["Hydroxylation"]


def test_ruleset_appends_itself_to_info_rule_after_the_leaf():
    """Live rulesets append ``self`` onto ``info["rule"]`` (same order as
    ``addition["rules"]``): leaf first, containing set last."""

    leaf = Hydroxylation()
    ruleset = RuleSet((leaf,), name="Outer")
    _products, info = next(ruleset.metabolize(Chem.MolFromSmiles("CC")))
    rules = info["rule"]
    assert isinstance(rules, list)
    assert type(rules[0]) is Hydroxylation
    assert rules[0] is leaf
    assert rules[-1] is ruleset
    assert len(rules) == 2
    assert [rule.name for rule in rules] == ["Hydroxylation", "Outer"]


def test_ruleset_appends_itself_after_the_leaf_on_the_same_molecule():
    leaf = Hydroxylation()
    ruleset = RuleSet((leaf,), name="Outer")
    emitted = []
    metabolize = leaf.metabolize

    def record(mol, **kwargs):
        for products, info in metabolize(mol, **kwargs):
            emitted.extend(products)
            yield products, info

    leaf.metabolize = record  # pyright: ignore[reportAttributeAccessIssue]
    products, _info = next(ruleset.metabolize(Chem.MolFromSmiles("CC")))
    product = products[0]
    assert product is emitted[0]
    rules = _rules(product)
    assert type(rules[0]) is Hydroxylation
    assert rules[-1] is ruleset
    assert len(rules) == 2


def test_nested_phaseone_ruleset_keeps_the_inner_set_on_the_chain():
    hydroxylation = StableOxygenation_PhaseOne.rules[0]
    assert type(hydroxylation) is Hydroxylation

    product = None
    info = None
    for products, emitted_info in PhaseOneRS.metabolize(
        Chem.MolFromSmiles("CC"),
        filter_rules=lambda mol, rule, info: type(rule) is Hydroxylation,
    ):
        product = products[0]
        info = emitted_info
        break

    assert product is not None
    assert info is not None
    addition = _addition(product)
    rules = addition["rules"]
    assert type(rules[0]) is Hydroxylation
    assert rules[1] is StableOxygenation_PhaseOne
    assert rules[-1] is PhaseOneRS
    assert len(rules) == 3
    # ``info["rule"]`` mirrors the addition chain: leaf, then each ruleset.
    assert isinstance(info["rule"], list)
    assert info["rule"] == list(rules)
    assert [rule.name for rule in info["rule"]] == ["Hydroxylation", "SO", None]
    assert StableOxygenation_PhaseOne.name == "SO"
    assert PhaseOneRS.name is None
    from xenosite.forest.phaseone import reaction_labels

    assert reaction_labels(addition) == ("Hydroxylation", "SO")
    assert addition["pattern"] is Hydroxylation.smirks[1][1]
    pattern = addition["pattern"]
    assert pattern is not None
    assert pattern.get("name") == "h2"
    assert addition["pattern"] not in rules


def test_unnamed_ruleset_stays_on_the_chain_and_emits_no_name():
    leaf = Hydroxylation()
    ruleset = RuleSet((leaf,), name="")
    assert ruleset.name is None

    products, _info = next(ruleset.metabolize(Chem.MolFromSmiles("CC")))

    product = products[0]
    addition = _addition(product)
    rules = addition["rules"]
    assert type(rules[0]) is Hydroxylation
    assert rules[-1] is ruleset
    from xenosite.forest.phaseone import reaction_labels

    assert reaction_labels(addition) == ("Hydroxylation",)
    assert addition["pattern"] is Hydroxylation.smirks[1][1]
