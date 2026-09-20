"""Shared inventory of PatternInfo possibilities for coverage + meta-tests.

Both ``test_pattern_info_coverage`` and the completeness meta-test call
:func:`iter_pattern_possibilities` so the parameterized rows cannot drift
from the rule-derived set.
"""

from __future__ import annotations

import inspect
from collections.abc import Iterator
from typing import Any, NamedTuple

from xenosite.refactor_poc import rules as rules_mod
from xenosite.refactor_poc import rulesets
from xenosite.refactor_poc.records import PatternInfo, When
from xenosite.refactor_poc.rules import (
    ReactionRule,
    ResonancePairRule,
    ResonanceRule,
    SmartsReactionRule,
)
from xenosite.refactor_poc.rulesets import RuleSet

# Bases / containers that inherit ReactionRule but are not themselves
# pattern-coverage targets. Keep this tight: every entry must still
# subclass ReactionRule, and every discovered subclass that is not listed
# here is required to appear in the pattern-info coverage inventory
# (vacuously if it has no patterns).
PATTERNLESS_REACTION_RULE_BASES: frozenset[type[ReactionRule]] = frozenset(
    {
        ReactionRule,
        SmartsReactionRule,
        ResonanceRule,
        ResonancePairRule,
        RuleSet,
    }
)


class PatternPossibility(NamedTuple):
    rule_cls: type[ReactionRule]
    group: str
    smarts: str
    pattern_name: str
    poss_i: int
    when: When | None


def when_key(when: When | None) -> tuple[int | None, int | None, int | None] | None:
    if when is None:
        return None
    return (when.get("map"), when.get("z"), when.get("h"))


def possibility_key(
    rule_cls: type[ReactionRule] | str,
    pattern_name: str,
    poss_i: int,
    when: When | None = None,
) -> tuple[str, str, int, tuple[int | None, int | None, int | None] | None]:
    name = rule_cls if isinstance(rule_cls, str) else rule_cls.__name__
    return (name, pattern_name, poss_i, when_key(when))


def discover_reaction_rule_classes() -> list[type[ReactionRule]]:
    """Unique ``ReactionRule`` subclasses declared under refactor_poc."""

    found: dict[type[ReactionRule], type[ReactionRule]] = {}
    for mod in (rules_mod, rulesets):
        for _name, obj in inspect.getmembers(mod, inspect.isclass):
            if not issubclass(obj, ReactionRule):
                continue
            found[obj] = obj
    return sorted(found.values(), key=lambda c: c.__name__)


def instantiate_rule(cls: type[ReactionRule]) -> ReactionRule:
    kwargs: dict[str, Any] = {}
    sig = inspect.signature(cls.__init__)
    if "as_star" in sig.parameters:
        kwargs["as_star"] = False
    return cls(**kwargs) if kwargs else cls()


def patterns_on(
    rule: ReactionRule,
) -> list[tuple[str, str, PatternInfo]]:
    """``(group, smarts, info)`` for ``smarts`` and ``endpoints``."""

    rows: list[tuple[str, str, PatternInfo]] = []
    for group in ("smarts", "endpoints"):
        for smarts, info in getattr(rule, group, ()) or ():
            rows.append((group, smarts, info))
    return rows


def iter_pattern_possibilities() -> Iterator[PatternPossibility]:
    """Every possibility on every non-whitelisted concrete ``ReactionRule``.

    Whitelisted bases are skipped entirely. Other classes must instantiate;
    a ``TypeError`` here is a real gap (add a whitelist entry with a comment,
    or fix construction), not a silent skip.
    """

    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        rule = instantiate_rule(cls)
        for group, smarts, info in patterns_on(rule):
            name = info.get("name") or ""
            for poss_i, poss in enumerate(info.get("possibilities") or ()):
                yield PatternPossibility(
                    rule_cls=cls,
                    group=group,
                    smarts=smarts,
                    pattern_name=name,
                    poss_i=poss_i,
                    when=poss.get("when"),
                )


def pattern_possibility_keys() -> set[
    tuple[str, str, int, tuple[int | None, int | None, int | None] | None]
]:
    return {
        possibility_key(row.rule_cls, row.pattern_name, row.poss_i, row.when)
        for row in iter_pattern_possibilities()
    }
