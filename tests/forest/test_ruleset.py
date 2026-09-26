"""RuleSet is a rule that runs its children, and filters still see each child."""

import logging
import warnings

from rdkit import Chem

from xenosite.forest.records import PatternInfo
from xenosite.forest.rules import (
    Dealkylation,
    Hydroxylation,
    ReactionRule,
    RuleSiteKind,
    SiteDeduplicationWarning,
    SmirksReactionRule,
    _describe,
)
from xenosite.forest.rulesets import PhaseOne, RuleSet

_RULES_LOG = "xenosite.forest.rules"


class OverlapOhA(SmirksReactionRule):
    """Minimal hydroxylator for redundant-rules tests."""

    site_kind: RuleSiteKind = "atom"
    smirks: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h3:1]>>[*:1]O",
            _describe(adds="O", removes="H", name="oh", site_map=1),
        ),
    )


class OverlapOhB(SmirksReactionRule):
    """Same SMARTS as ``OverlapOhA``, different rule name."""

    site_kind: RuleSiteKind = "atom"
    smirks: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h3:1]>>[*:1]O",
            _describe(adds="O", removes="H", name="oh", site_map=1),
        ),
    )


def test_ruleset_metabolize_runs_each_child():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Forest")
    assert isinstance(ruleset, ReactionRule)

    products = list(ruleset.metabolize(Chem.MolFromSmiles("CC")))
    names = {type(info["rule"][0]).__name__ for _product, info in products}

    assert names == {"Hydroxylation", "Dealkylation"}


def test_ruleset_redundant_rules_info_names_both_rules(caplog):
    """Cross-rule same CSMI → INFO naming both rules; no Warning emission."""

    ruleset = RuleSet((OverlapOhA, OverlapOhB), name="OverlapSet")
    mol = Chem.MolFromSmiles("CC")

    with caplog.at_level(logging.INFO, logger=_RULES_LOG):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("error", SiteDeduplicationWarning)
            products = list(ruleset.metabolize(mol))

    assert [p.xf.csmi for pl, _ in products for p in pl] == ["CCO"]
    assert {type(info["rule"][0]).__name__ for _, info in products} == {"OverlapOhA"}
    assert not caught

    info_msgs = [
        r.getMessage()
        for r in caplog.records
        if r.name == _RULES_LOG and r.levelno == logging.INFO
    ]
    assert info_msgs, "expected INFO for overlapping rules"
    message = info_msgs[0]
    assert "OverlapOhA" in message
    assert "OverlapOhB" in message
    assert "CCO" in message
    assert "kept_rule" in message and "dropped_rule" in message


def test_ruleset_redundant_rules_logs_info(caplog):
    """Focused caplog entry for RuleSet CSMI = redundant rules (INFO only)."""

    ruleset = RuleSet((OverlapOhA, OverlapOhB), name="OverlapSet")
    with caplog.at_level(logging.INFO, logger=_RULES_LOG):
        list(ruleset.metabolize(Chem.MolFromSmiles("CC")))

    messages = " ".join(
        r.getMessage() for r in caplog.records if r.name == _RULES_LOG
    )
    assert "OverlapOhB" in messages
    assert "OverlapOhA" in messages
    assert "CCO" in messages


def test_filter_rules_refuses_one_child_and_keeps_the_other():
    ruleset = RuleSet((Hydroxylation, Dealkylation), name="Forest")
    seen = []

    def filter_rules(mol, rule, info):
        seen.append(type(rule).__name__)
        assert "span" in info
        return type(rule) is not Hydroxylation

    products = list(
        ruleset.metabolize(Chem.MolFromSmiles("CC"), filter_rules=filter_rules)
    )
    names = {type(info["rule"][0]).__name__ for _product, info in products}

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
