"""A rule set is a rule that runs the rules it contains.

Callers execute one set. ``filter_rules`` and ``filter_sites`` are passed to
each child, so the set does not hide that child's patterns. A pattern the
filter refuses is not run. A pattern it accepts is.
"""

from xenosite.refactor_poc.rdkitutil import Mol
from .rules import (
    Dealkylation,
    Dehydration,
    Dehydrogenation,
    Dephosphorylation,
    Epoxidation,
    EpoxideOpening,
    Hydrolysis,
    Hydroxylation,
    Hydrogenation,
    NitrogenOxidation,
    NitrogenReduction,
    OxidativeDehalogenation,
    OxygenReduction,
    QuinoneFormation,
    ReactionRule,
    ReductiveDehalogenation,
    SulfurOxidation,
    SulfurReduction,
)


class RuleSet(ReactionRule):
    """Container of rules, and itself a rule.

    Nested sets are flattened. Executing the outer set runs the leaves, and
    the filters still see each leaf pattern (``span`` before a match,
    resolved effect after). ``order_key`` only sorts children; it does not
    wrap them.
    """

    def __init__(self, rules=(), name=None, longname=None):
        contained = []
        for rule in rules:
            if isinstance(rule, type):
                rule = rule()
            if isinstance(rule, RuleSet):
                contained.extend(rule.rules)
            elif isinstance(rule, ReactionRule):
                contained.append(rule)
            else:
                raise TypeError(
                    "RuleSet holds rule classes or rules, not %r" % (rule,)
                )
        self.rules = tuple(contained)
        if name is None and len(self.rules) == 1:
            name = self.rules[0].name
        super().__init__(name=name or "RuleSet", longname=longname)

    def __iter__(self):
        yield from self.rules

    def metabolites(
        self,
        mol: Mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        order_key=None,
        **kwargs,
    ):
        rules = self.rules
        if order_key is not None:
            rules = tuple(sorted(rules, key=order_key))
        for rule in rules:
            yield from rule.metabolites(
                mol,
                filter_rules=filter_rules,
                filter_sites=filter_sites,
                **kwargs,
            )


# Phase I classes that already exist in rules.py. QuinoneFormation is included
# because that class exists. Forest rules with no class here are omitted.
PhaseOne = RuleSet(
    (
        Hydroxylation,
        Epoxidation,
        SulfurOxidation,
        NitrogenOxidation,
        Dehydrogenation,
        QuinoneFormation,
        Dephosphorylation,
        EpoxideOpening,
        Hydrolysis,
        Dehydration,
        Hydrogenation,
        NitrogenReduction,
        OxygenReduction,
        ReductiveDehalogenation,
        SulfurReduction,
        Dealkylation,
        OxidativeDehalogenation,
    ),
    name="PhaseOne",
    longname="Phase I",
)
