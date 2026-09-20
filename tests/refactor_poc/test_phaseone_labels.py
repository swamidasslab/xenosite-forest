"""Legacy labels are the class names on the addition's rule chain."""

from rdkit import Chem

from xenosite.refactor_poc.phaseone import (
    Dehydrogenation_PhaseOne,
    PhaseOne,
    PhaseOneQF,
    PhaseOneRS,
    metabolize,
    reaction_labels,
)
from xenosite.refactor_poc.records import Addition
from xenosite.refactor_poc.rules import Epoxidation, Hydroxylation, QuinoneFormation
from xenosite.refactor_poc.rulesets import PhaseOne as CatalogPhaseOne
from xenosite.refactor_poc.rulesets import RuleSet


def _addition(ruleset, smiles, leaf):
    mol = Chem.MolFromSmiles(smiles)
    for product, _info in ruleset.metabolize(
        mol, filter_rules=lambda rule, info: type(rule) is leaf
    ):
        return product._forest["atom_trace"]["additions"]["R1"]
    raise AssertionError(leaf)


def test_phaseone_is_the_catalog_ruleset():
    assert PhaseOne is CatalogPhaseOne
    assert isinstance(PhaseOne, RuleSet)
    assert isinstance(PhaseOneRS, RuleSet)
    assert isinstance(PhaseOneQF, RuleSet)


def test_metabolize_calls_phaseone():
    mol = Chem.MolFromSmiles("CC")
    direct = {info["csmi"] for _product, info in PhaseOne.metabolize(mol)}
    wrapped = {info["csmi"] for _product, info in metabolize(mol)}
    assert wrapped == direct
    assert "CCO" in wrapped


def test_hydroxylation_label_is_the_chain_class_names():
    addition = _addition(PhaseOne, "CC", Hydroxylation)

    assert reaction_labels(addition) == ("RuleSet", "Hydroxylation")
    assert reaction_labels(addition) == tuple(
        type(rule).__name__ for rule in addition["rules"]
    )
    assert addition["phase1"] is None


def test_ruleset_class_name_is_on_the_chain_when_that_set_runs():
    addition = _addition(Dehydrogenation_PhaseOne, "CC", type(Dehydrogenation_PhaseOne.rules[0]))

    assert reaction_labels(addition) == ("RuleSet", "Dehydrogenation")
    assert addition["rules"][0].name == "DH"
    assert "DH" not in reaction_labels(addition)


def test_quinone_label_is_the_rule_on_the_chain():
    addition = _addition(PhaseOne, "c1ccccc1", QuinoneFormation)

    assert reaction_labels(addition) == ("RuleSet", "QuinoneFormation")
    assert "Hydroxylation" not in reaction_labels(addition)
    assert "Dehydrogenation" not in reaction_labels(addition)

    alone = _addition(QuinoneFormation(), "c1ccccc1", QuinoneFormation)
    assert reaction_labels(alone) == ("QuinoneFormation",)


def test_phaseoneqf_chain_includes_the_ruleset_class():
    addition = _addition(PhaseOneQF, "c1ccccc1", QuinoneFormation)

    assert reaction_labels(addition) == ("RuleSet", "QuinoneFormation")
    assert [type(rule).__name__ for rule in PhaseOneQF].count("QuinoneFormation") == 1


def test_epoxidation_label_is_not_a_look_ahead():
    addition = _addition(PhaseOne, "C=C", Epoxidation)

    assert reaction_labels(addition) == ("RuleSet", "Epoxidation")


def test_grouped_set_omits_deferred_rules():
    names = [type(rule).__name__ for rule in PhaseOneRS]
    assert "Epoxidation" in names
    assert "NDealkylation" not in names
    assert "Tautomerization" not in names
    assert "QuinoneFormation" not in names


def test_labels_read_a_named_addition_the_same_way():
    forest_addition = _addition(PhaseOne, "CC", Hydroxylation)
    addition = Addition(
        site=forest_addition["site"],
        rules=forest_addition["rules"],
        info={},
        effect={},
        name=forest_addition["name"],
        phase1=None,
        depth=forest_addition["depth"],
    )

    assert reaction_labels(addition) == reaction_labels(forest_addition)


def test_any_rules_on_the_chain_contribute_their_class_names():
    addition = Addition(
        site=(0,),
        rules=(PhaseOne, PhaseOneRS, Hydroxylation()),
        info={},
        effect={},
        name=None,
        phase1=None,
        depth=0,
    )

    assert reaction_labels(addition) == ("RuleSet", "RuleSet", "Hydroxylation")
