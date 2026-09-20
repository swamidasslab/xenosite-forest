"""A rule set is a rule that runs the rules it contains.

Callers execute one set. ``filter_rules`` and ``filter_sites`` are passed to
each child, so the set does not hide that child's patterns. A pattern the
filter refuses is not run. A pattern it accepts is.
"""

from __future__ import annotations

from collections.abc import Callable, Generator, Iterator
from typing import Any

from xenosite.forest.rdkit_api import TracingMol
from xenosite.forest.rdkitutil import Mol
from xenosite.forest.records import ProductInfo

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
    _accept_all_rules,
    _accept_all_sites,
    _report_redundant_rules_drop,
    _rule_dedup_name,
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

    def __init__(
        self,
        rules: (
            tuple[ReactionRule | type[ReactionRule], ...]
            | list[ReactionRule | type[ReactionRule]]
        ) = (),
        name: str | None = None,
        longname: str | None = None,
    ):
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
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        order_key: Callable[[ReactionRule], Any] | None = None,
        **kwargs: Any,
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
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        unique_csmi: bool = True,
        order_key: Callable[[ReactionRule], Any] | None = None,
        **kwargs: Any,
    ) -> Generator[tuple[TracingMol, ProductInfo], None, None]:
        """Run each contained rule, then append this set on that product.

        The contained rule, including a nested set, has already put itself
        last on the product. This set is the next rule. The product object
        is the one that rule yielded.
        """

        if mol is None:
            raise ValueError("mol is required")
        mol = mol.xf.tracing._stamp()
        if self.is_terminal_product(mol):
            return

        rules = self.rules
        if order_key is not None:
            rules = tuple(sorted(rules, key=order_key))
        # Children: yield off (``unique_csmi=False``) so alternate rules /
        # nested sets bubble up; leaf check still runs. This set's own yield
        # (caller ``unique_csmi``) is the cross-child CSMI layer — outermost
        # caller setting applies at each RuleSet that was asked to uniquify.
        # Cross-rule same CSMI → INFO (not SiteDeduplicationWarning).
        # Same rule, different PatternInfo tokens may still both emit.
        seen_csmi: dict[str, str] = {}
        for rule in rules:
            for product, info in rule.metabolize(
                mol,
                filter_rules=filter_rules,
                filter_sites=filter_sites,
                **kwargs,
                unique_csmi=False,
            ):
                trace = product._forest["atom_trace"]
                addition = trace["additions"][trace["transforms"][-1]]
                addition["rules"] = tuple(addition["rules"]) + (self,)
                if unique_csmi:
                    product_csmi = product.xf.csmi
                    rule_name = _rule_dedup_name(info["rule"])
                    kept = seen_csmi.get(product_csmi)
                    if kept is not None and kept != rule_name:
                        _report_redundant_rules_drop(
                            mol,
                            kept_rule=kept,
                            dropped_info=info,
                            product_csmi=product_csmi,
                        )
                        continue
                    if kept is None:
                        seen_csmi[product_csmi] = rule_name
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
