"""Shared inventory of PatternInfo possibilities for coverage + meta-tests.

``test_pattern_info_coverage``, the completeness meta-test, and within-rule
name uniqueness all call helpers here so parameterized rows and uniqueness
checks cannot drift from the rule-derived set.
"""

from __future__ import annotations

import inspect
from collections.abc import Iterator, Mapping
from typing import Any, NamedTuple

from xenosite.forest import rules as rules_mod
from xenosite.forest import rulesets
from xenosite.forest.records import PatternInfo, When
from xenosite.forest.rules import (
    ReactionRule,
    ResonancePairRule,
    ResonanceRule,
    SmartsReactionRule,
)
from xenosite.forest.rulesets import RuleSet

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
        rules_mod.TautomerRule,  # design stub; raises NotImplementedError
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
    """Unique ``ReactionRule`` subclasses declared under xenosite.forest."""

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


def _when_resolved_name(when: When | None) -> str | None:
    """Optional when→name label, if the schema ever carries one.

    ``when → name`` is only a HEURISTICS Schema proposal (Status: not
    decided). ``When`` has no ``name`` field today — do not invent it.
    Read via :class:`~collections.abc.Mapping` so a future optional key
    plugs into uniqueness without rewriting callers.
    """

    if when is None:
        return None
    raw = Mapping[str, Any](when).get("name")
    if isinstance(raw, str) and raw:
        return raw
    return None


def emitable_names_for_pattern(info: PatternInfo) -> list[str]:
    """Names one PatternInfo can emit (for ``unique_csmi`` / uniqueness).

    Always includes ``PatternInfo.name`` when set. Also folds any
    when-resolved names once that optional key appears on ``When``.
    """

    names: list[str] = []
    base = info.get("name") or ""
    if base:
        names.append(base)
    for poss in info.get("possibilities") or ():
        resolved = _when_resolved_name(poss.get("when"))
        if resolved:
            names.append(resolved)
    return names


def optional_when_name(when: When | None) -> str | None:
    """Optional emit name on a ``When`` branch.

    ``When`` has no ``name`` field yet (HEURISTICS Schema proposals,
    Status: not decided). Read the key if a future schema adds it so
    :func:`emitable_pattern_names` picks up resolved tokens without a
    new code path.
    """

    if when is None:
        return None
    raw = dict(when).get("name")
    if isinstance(raw, str) and raw:
        return raw
    return None


def emitable_pattern_names(info: PatternInfo) -> frozenset[str]:
    """Names this pattern may contribute to ``unique_csmi`` / traces.

    Always includes ``PatternInfo.name`` when set. Also includes every
    ``when``-branch name once when→name lands (optional ``When["name"]``).
    Empty when the pattern still lacks a name.
    """

    names: set[str] = set()
    base = info.get("name")
    if isinstance(base, str) and base:
        names.add(base)
    for poss in info.get("possibilities") or ():
        branch = optional_when_name(poss.get("when"))
        if branch is not None:
            names.add(branch)
    return frozenset(names)


def emitable_names_by_rule() -> dict[str, dict[str, list[str]]]:
    """``rule_cls_name -> emit_name -> [pattern sources]`` for concrete rules.

    Skips :data:`PATTERNLESS_REACTION_RULE_BASES`. A source string is
    ``"{group}:{PatternInfo.name|smarts}"`` so duplicate emit tokens point
    at the conflicting patterns. Cross-rule reuse of a name is fine.
    """

    by_rule: dict[str, dict[str, list[str]]] = {}
    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        rule = instantiate_rule(cls)
        rule_key = cls.__name__
        bucket = by_rule.setdefault(rule_key, {})
        for group, smarts, info in patterns_on(rule):
            label = info.get("name") or smarts
            source = f"{group}:{label}"
            for emit in emitable_pattern_names(info):
                bucket.setdefault(emit, []).append(source)
    return by_rule
