"""Every PatternInfo should carry a chemically meaningful name.

Auto-numbered leftovers from ``_assign_pattern_names`` ("1", "2", …) are
interim only. This test treats them as still empty until replaced.
Names must be unique within each rule; duplication across rules is fine.
"""

from __future__ import annotations

import inspect
from collections import defaultdict

from xenosite.refactor_poc import rules as rules_mod
from xenosite.refactor_poc.records import PatternInfo
from xenosite.refactor_poc.rules import SmartsReactionRule


def _patterns_declared_on(cls: type) -> list[PatternInfo]:
    """PatternInfo dicts owned by ``cls``, not inherited from a base."""

    found: list[PatternInfo] = []
    for attr in ("smarts", "endpoints"):
        if attr not in cls.__dict__:
            continue
        for entry in cls.__dict__[attr] or ():
            if isinstance(entry, tuple) and len(entry) >= 2:
                info = entry[1]
                if isinstance(info, dict):
                    found.append(info)
    return found


def _iter_rule_patterns() -> list[tuple[str, PatternInfo]]:
    rows: list[tuple[str, PatternInfo]] = []
    for name, cls in inspect.getmembers(rules_mod, inspect.isclass):
        if not issubclass(cls, SmartsReactionRule):
            continue
        for info in _patterns_declared_on(cls):
            rows.append((name, info))
    return rows


def _name_is_missing(info: PatternInfo) -> bool:
    name = info.get("name")
    if name is None or name == "":
        return True
    # Digit-only names are the interim auto-numbering, not chemical names.
    return str(name).isdigit()


def test_all_pattern_infos_have_names():
    missing = [
        (rule_name, info.get("name"), info.get("edit"))
        for rule_name, info in _iter_rule_patterns()
        if _name_is_missing(info)
    ]
    assert not missing, (
        f"{len(missing)} PatternInfo still lack chemically meaningful names: "
        f"{missing[:20]}{'…' if len(missing) > 20 else ''}"
    )


def test_pattern_names_unique_within_rule():
    by_rule: dict[str, list[str]] = defaultdict(list)
    for rule_name, info in _iter_rule_patterns():
        name = info.get("name")
        if name is None or name == "":
            continue
        by_rule[rule_name].append(str(name))
    dups = {
        rule: names
        for rule, names in by_rule.items()
        if len(names) != len(set(names))
    }
    assert not dups, f"duplicate PatternInfo names within a rule: {dups}"
