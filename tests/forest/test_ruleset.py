"""RuleSet is a rule that runs its children, and filters still see each child."""

from rdkit import Chem

from xenosite.forest.rules import Dealkylation, Hydroxylation, ReactionRule
from xenosite.forest.rulesets import PhaseOne, RuleSet


def test_ruleset_metabolize_runs_each_child():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Poc")
    assert isinstance(ruleset, ReactionRule)

    products = list(ruleset.metabolize(Chem.MolFromSmiles("CC")))
    names = {type(info["rule"]).__name__ for _product, info in products}

    assert names == {"Hydroxylation", "Dealkylation"}


def test_filter_rules_refuses_one_child_and_keeps_the_other():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Poc")
    seen = []

    def filter_rules(mol, rule, info):
        seen.append(type(rule).__name__)
        assert "span" in info
        return type(rule) is not Hydroxylation

    products = list(
        ruleset.metabolize(Chem.MolFromSmiles("CC"), filter_rules=filter_rules)
    )
    names = {type(info["rule"]).__name__ for _product, info in products}

    assert names == {"Dealkylation"}
    assert "Hydroxylation" in seen
    assert "Dealkylation" in seen
    assert "RuleSet" not in seen


def test_phaseone_lists_existing_rule_classes():
    assert isinstance(PhaseOne, RuleSet)
    assert [type(rule).__name__ for rule in PhaseOne] == [
        "Hydroxylation",
        "Epoxidation",
        "SulfurOxidation",
        "NitrogenOxidation",
        "Dehydrogenation",
        "QuinoneFormation",
        "Dephosphorylation",
        "EpoxideOpening",
        "Hydrolysis",
        "Dehydration",
        "Hydrogenation",
        "NitrogenReduction",
        "OxygenReduction",
        "ReductiveDehalogenation",
        "SulfurReduction",
        "Dealkylation",
        "OxidativeDehalogenation",
    ]
