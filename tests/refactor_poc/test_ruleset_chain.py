"""A ruleset stays on the chain after the rule that emitted the product."""

from rdkit import Chem

from xenosite.refactor_poc.phaseone import PhaseOneRS, StableOxygenation_PhaseOne
from xenosite.refactor_poc.rules import Hydroxylation
from xenosite.refactor_poc.rulesets import RuleSet


def _rules(product):
    return product._forest["atom_trace"]["additions"]["R1"]["rules"]


def test_ruleset_appends_itself_after_the_leaf_on_the_same_molecule():
    leaf = Hydroxylation()
    ruleset = RuleSet((leaf,), name="Outer")
    emitted = []
    metabolize = leaf.metabolize

    def record(mol, **kwargs):
        for product, info in metabolize(mol, **kwargs):
            emitted.append(product)
            yield product, info

    leaf.metabolize = record
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
        filter_rules=lambda rule, info: type(rule) is Hydroxylation,
    ):
        break

    assert product is not None
    rules = _rules(product)
    assert type(rules[0]) is Hydroxylation
    assert rules[1] is StableOxygenation_PhaseOne
    assert rules[-1] is PhaseOneRS
    assert len(rules) == 3
