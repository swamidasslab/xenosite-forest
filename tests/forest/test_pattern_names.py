"""Every PatternInfo should carry a chemically meaningful name.

Auto-numbered leftovers from ``_assign_pattern_names`` ("1", "2", …) are
interim only. This test treats them as still empty until replaced.
Names must be unique within each rule; duplication across rules is fine.

Emit-name uniqueness also covers optional when→name resolution (HEURISTICS
Schema proposals): once a ``When`` carries ``name=``, that token joins the
per-rule uniqueness set via :func:`emitable_pattern_names`.
"""

from __future__ import annotations

from xenosite.forest.records import PatternInfo

from .pattern_info_inventory import (
    PATTERNLESS_REACTION_RULE_BASES,
    discover_reaction_rule_classes,
    emitable_names_by_rule,
    emitable_pattern_names,
    instantiate_rule,
    optional_when_name,
    patterns_on,
)


def _name_is_missing(info: PatternInfo) -> bool:
    name = info.get("name")
    if name is None or name == "":
        return True
    # Digit-only names are the interim auto-numbering, not chemical names.
    return str(name).isdigit()


def test_all_pattern_infos_have_names():
    missing: list[tuple[str, str | None, str | None]] = []
    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        rule = instantiate_rule(cls)
        for _group, _smarts, info in patterns_on(rule):
            if _name_is_missing(info):
                missing.append(
                    (cls.__name__, info.get("name"), info.get("edit"))
                )
    assert not missing, (
        f"{len(missing)} PatternInfo still lack chemically meaningful names: "
        f"{missing[:20]}{'…' if len(missing) > 20 else ''}"
    )


def test_emitable_names_unique_within_each_reaction_rule():
    """Every concrete ReactionRule: emit names unique within that rule.

    Collects :class:`PatternInfo` names plus any optional when-resolved
    names (see :func:`optional_when_name`). Cross-rule duplicates OK.
    """

    collisions: dict[str, dict[str, list[str]]] = {}
    for rule_name, by_emit in emitable_names_by_rule().items():
        dups = {
            emit: sources
            for emit, sources in by_emit.items()
            if len(sources) > 1
        }
        if dups:
            collisions[rule_name] = dups
    assert not collisions, (
        "duplicate emitable pattern names within a ReactionRule "
        f"(cross-rule reuse is OK): {collisions}"
    )


def test_emitable_pattern_names_includes_future_when_name():
    """Hook: when a When carries name=, it joins the emit set."""

    info: PatternInfo = {
        "name": "base",
        "possibilities": (
            {"when": {"map": 1, "z": 6, "h": 2, "name": "h2"}},  # type: ignore[typeddict-item]
            {"when": {"map": 1, "z": 6, "h": 3, "name": "h3"}},  # type: ignore[typeddict-item]
        ),
    }
    assert emitable_pattern_names(info) == frozenset({"base", "h2", "h3"})
    assert optional_when_name({"map": 1, "z": 6, "h": 2}) is None
    assert optional_when_name(
        {"map": 1, "z": 6, "h": 2, "name": "h2"}  # type: ignore[arg-type,typeddict-item]
    ) == "h2"
