"""A ruleset stays on the chain after the rule that emitted the product."""

from rdkit import Chem

from xenosite.forest.phaseone import PhaseOneRS, StableOxygenation_PhaseOne
from xenosite.forest.rules import Hydroxylation
from xenosite.forest.rulesets import RuleSet


def _addition(product):
    return product._forest["atom_trace"]["additions"]["R1"]


def _rules(product):
    return _addition(product)["rules"]


def test_ruleset_appends_itself_after_the_leaf_on_the_same_molecule():
    leaf = Hydroxylation()
    ruleset = RuleSet((leaf,), name="Outer")
    emitted = []
    metabolize = leaf.metabolize

    def record(mol, **kwargs):
        for product, info in metabolize(mol, **kwargs):
            emitted.append(product)
            yield product, info

    leaf.metabolize = record  # pyright: ignore[reportAttributeAccessIssue]
    product, _info = next(ruleset.metabolize(Chem.MolFromSmiles("CC")))

    assert product is emitted[0]
    rules = _rules(product)
    assert type(rules[0]) is Hydroxylation
    assert rules[-1] is ruleset
    assert len(rules) == 2


def test_nested_phaseone_ruleset_keeps_the_inner_set_on_the_chain():
    hydroxylation = StableOxygenation_PhaseOne.rules[0]
    assert type(hydroxylation) is Hydroxylation

    product = None
    for product, _info in PhaseOneRS.metabolize(
        Chem.MolFromSmiles("CC"),
        filter_rules=lambda mol, rule, info: type(rule) is Hydroxylation,
    ):
        break

    assert product is not None
    addition = _addition(product)
    rules = addition["rules"]
    assert type(rules[0]) is Hydroxylation
    assert rules[1] is StableOxygenation_PhaseOne
    assert rules[-1] is PhaseOneRS
    assert len(rules) == 3
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

    product, _info = next(ruleset.metabolize(Chem.MolFromSmiles("CC")))
    addition = _addition(product)
    rules = addition["rules"]
    assert type(rules[0]) is Hydroxylation
    assert rules[-1] is ruleset
    from xenosite.forest.phaseone import reaction_labels

    assert reaction_labels(addition) == ("Hydroxylation",)
    assert addition["pattern"] is Hydroxylation.smirks[1][1]
