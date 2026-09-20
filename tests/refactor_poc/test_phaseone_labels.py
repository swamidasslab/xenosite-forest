"""Legacy labels are the initialization names on the addition's rule chain."""

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


def test_hydroxylation_label_is_the_chain_init_names():
    addition = _addition(PhaseOne, "CC", Hydroxylation)

    assert reaction_labels(addition) == ("Hydroxylation", "PhaseOne")
    assert type(addition["rules"][0]) is Hydroxylation
    assert addition["rules"][-1] is PhaseOne
    assert addition["phase1"] is None
    assert addition["pattern"] is Hydroxylation.smarts[0][1]
    assert addition["pattern"]["name"] == "h"
    assert all(not isinstance(rule, dict) for rule in addition["rules"])


def test_named_ruleset_reports_its_init_name():
    addition = _addition(Dehydrogenation_PhaseOne, "CC", type(Dehydrogenation_PhaseOne.rules[0]))

    assert reaction_labels(addition) == ("Dehydrogenation", "DH")
    assert type(addition["rules"][0]) is type(Dehydrogenation_PhaseOne.rules[0])
    assert addition["rules"][-1] is Dehydrogenation_PhaseOne
    assert addition["rules"][-1].name == "DH"


def test_quinone_label_is_the_rule_on_the_chain():
    addition = _addition(PhaseOne, "c1ccccc1", QuinoneFormation)

    assert reaction_labels(addition) == ("QuinoneFormation", "PhaseOne")
    assert "Hydroxylation" not in reaction_labels(addition)
    assert "Dehydrogenation" not in reaction_labels(addition)

    alone = _addition(QuinoneFormation(), "c1ccccc1", QuinoneFormation)
    assert reaction_labels(alone) == ("QuinoneFormation",)


def test_phaseoneqf_chain_includes_the_ruleset_name():
    addition = _addition(PhaseOneQF, "c1ccccc1", QuinoneFormation)

    assert reaction_labels(addition) == ("QuinoneFormation", "PhaseOneQF")
    assert type(addition["rules"][0]) is QuinoneFormation
    assert addition["rules"][-1] is PhaseOneQF
    assert [type(rule).__name__ for rule in PhaseOneQF].count("QuinoneFormation") == 1


def test_epoxidation_label_is_not_a_look_ahead():
    addition = _addition(PhaseOne, "C=C", Epoxidation)

    assert reaction_labels(addition) == ("Epoxidation", "PhaseOne")


def _leaf_names(rule):
    if isinstance(rule, RuleSet):
        names = []
        for child in rule:
            names.extend(_leaf_names(child))
        return names
    return [type(rule).__name__]


def test_grouped_set_omits_deferred_rules():
    assert [rule.longname for rule in PhaseOneRS] == [
        "Dehydrogenation",
        "Hydrolysis",
        "Reduction",
        "StableOxygenation",
        "UnstableOxygenation",
    ]
    names = _leaf_names(PhaseOneRS)
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
        pattern=forest_addition.get("pattern"),
    )

    assert reaction_labels(addition) == reaction_labels(forest_addition)


def test_unnamed_ruleset_stays_on_the_chain_but_emits_no_label():
    addition = Addition(
        site=(0,),
        rules=(PhaseOne, PhaseOneRS, Hydroxylation()),
        info={},
        effect={},
        name=None,
        phase1=None,
        depth=0,
    )

    assert PhaseOne.name == "PhaseOne"
    assert PhaseOneRS.name is None
    assert reaction_labels(addition) == ("PhaseOne", "Hydroxylation")
