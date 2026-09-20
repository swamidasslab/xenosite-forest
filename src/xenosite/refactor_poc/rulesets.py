"""A rule set is a rule that runs the rules it contains.

Callers execute one set. ``filter_rules`` and ``filter_sites`` are passed to
each child, so the set does not hide that child's patterns. A pattern the
filter refuses is not run. A pattern it accepts is.
"""

from __future__ import annotations

from collections.abc import Generator, Iterator

from xenosite.refactor_poc.rdkit_api import ForestTracingMol
from xenosite.refactor_poc.rdkitutil import Mol
from xenosite.refactor_poc.records import ProductInfo

from .rules import (
    Dealkylation,
    Dehydration,
    Dehydrogenation,
    Dephosphorylation,
    Epoxidation,
    EpoxideOpening,
    FilterRules,
    FilterSites,
    Hydrogenation,
    Hydrolysis,
    Hydroxylation,
    NitrogenOxidation,
    NitrogenReduction,
    OxidativeDehalogenation,
    OxygenReduction,
    ProductsOfReaction,
    QuinoneFormation,
    ReactionRule,
    ReductiveDehalogenation,
    SulfurOxidation,
    SulfurReduction,
    ensure_tracing,
)


class RuleSet(ReactionRule):
    """Container of rules, and itself a rule.

    Nested sets stay nested. This set calls :meth:`metabolize` on each
    contained rule. That call's products already end with that rule. This
    set then appends itself on that product, after that rule. The molecule
    is not copied. The filters still see each leaf pattern (``span`` before
    a match, resolved effect after). ``order_key`` only sorts children; it
    does not wrap them.
    """

    def __init__(self, rules=(), name=None, longname=None):
        contained = []
        for rule in rules:
            if isinstance(rule, type):
                rule = rule()
            if isinstance(rule, ReactionRule):
                contained.append(rule)
            else:
                raise TypeError(
                    "RuleSet holds rule classes or rules, not %r" % (rule,)
                )
        self.rules = tuple(contained)
        if name is None and len(self.rules) == 1:
            name = self.rules[0].name
        if name == "":
            name = None
        # ReactionRule fills a class-name default when name is missing.
        # A ruleset keeps the given name, or no name. It still appends
        # itself on the chain either way.
        super().__init__(name=name, longname=longname)
        self.name = name
        if longname is None:
            self.longname = name
        else:
            self.longname = longname

    def __iter__(self) -> Iterator[ReactionRule]:
        yield from self.rules

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        order_key=None,
        **kwargs,
    ) -> Generator[ProductsOfReaction, None, None]:
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

    def metabolize(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        unique_csmi: bool = True,
        order_key=None,
        **kwargs,
    ) -> Generator[tuple[ForestTracingMol, ProductInfo], None, None]:
        """Run each contained rule, then append this set on that product.

        The contained rule, including a nested set, has already put itself
        last on the product. This set is the next rule. The product object
        is the one that rule yielded.
        """

        if mol is None:
            raise ValueError("mol is required")
        mol = ensure_tracing(mol)
        if self.is_terminal_product(mol):
            return

        rules = self.rules
        if order_key is not None:
            rules = tuple(sorted(rules, key=order_key))
        seen: set[str] = set()
        for rule in rules:
            for product, info in rule.metabolize(
                mol,
                filter_rules=filter_rules,
                filter_sites=filter_sites,
                unique_csmi=unique_csmi,
                **kwargs,
            ):
                trace = product._forest["atom_trace"]
                addition = trace["additions"][trace["transforms"][-1]]
                addition["rules"] = tuple(addition["rules"]) + (self,)
                if unique_csmi:
                    csmi = info["csmi"]
                    if csmi in seen:
                        continue
                    seen.add(csmi)
                yield product, info


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
