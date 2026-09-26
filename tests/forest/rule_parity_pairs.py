"""Discover and pair Rust ↔ Python leaf rules for product parity.

Exceptions live on **rule attributes** (Python ``rust_parity_exception``,
Rust ``RuleSet.parity_exception``), not in test bodies. A missing twin or a
product mismatch is a failure unless that attribute is a non-empty reason.
"""

from __future__ import annotations

from xenosite.forest.rules import ReactionRule

from .pattern_info_inventory import (
    PATTERNLESS_REACTION_RULE_BASES,
    discover_reaction_rule_classes,
    instantiate_rule,
)


def python_leaf_classes() -> dict[str, type[ReactionRule]]:
    """Concrete pattern-bearing Python rules keyed by class / leaf name."""

    out: dict[str, type[ReactionRule]] = {}
    for cls in discover_reaction_rule_classes():
        if cls in PATTERNLESS_REACTION_RULE_BASES:
            continue
        # Prefer the live ``name`` after init (matches Rust leaf labels).
        rule = instantiate_rule(cls)
        key = rule.name or cls.__name__
        out[key] = cls
    return out


def python_parity_exception(cls: type[ReactionRule]) -> str | None:
    raw = getattr(cls, "rust_parity_exception", None)
    if raw is None:
        return None
    text = str(raw).strip()
    return text or None


def rust_catalog_names() -> list[str]:
    import xenosite_forest as native  # type: ignore[import-not-found]

    return list(native.RuleSet.catalog_names())


def rust_parity_exception(rule_name: str) -> str | None:
    import xenosite_forest as native  # type: ignore[import-not-found]

    raw = native.RuleSet.leaf(rule_name).parity_exception
    if raw is None:
        return None
    text = str(raw).strip()
    return text or None


def paired_rule_names() -> list[str]:
    """Same-named leaves on both sides with **no** parity exception either side."""

    py = python_leaf_classes()
    rs = set(rust_catalog_names())
    names: list[str] = []
    for name in sorted(set(py) & rs):
        if python_parity_exception(py[name]):
            continue
        if rust_parity_exception(name):
            continue
        names.append(name)
    return names


def pairing_gaps() -> list[str]:
    """Human-readable unpaired / empty-exception problems."""

    py = python_leaf_classes()
    rs = set(rust_catalog_names())
    gaps: list[str] = []
    for name, cls in sorted(py.items()):
        if name in rs:
            continue
        exc = python_parity_exception(cls)
        if not exc:
            gaps.append(
                f"Python-only {name!r} lacks rust_parity_exception "
                f"(no Rust leaf_rule twin)"
            )
    for name in sorted(rs - set(py)):
        exc = rust_parity_exception(name)
        if not exc:
            gaps.append(
                f"Rust-only {name!r} lacks parity_exception "
                f"(no Python ReactionRule twin)"
            )
    for name in sorted(set(py) & rs):
        py_exc = python_parity_exception(py[name])
        rs_exc = rust_parity_exception(name)
        if py_exc is not None and not py_exc:
            gaps.append(f"{name!r}: Python rust_parity_exception is empty")
        if rs_exc is not None and not rs_exc:
            gaps.append(f"{name!r}: Rust parity_exception is empty")
    return gaps
