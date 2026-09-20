"""Meta-test: every concrete rule's ``site_kind`` matches emitted Sites.

Uses each rule's internal ``_example_substrates`` (SMILES guaranteed to
produce metabolites). Expand those lists later so they cover all
patterns/whens — see TODO.md.
"""

from __future__ import annotations

import pytest

from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.rules import ReactionRule, TautomerRule

from .pattern_info_inventory import (
    PATTERNLESS_REACTION_RULE_BASES,
    discover_reaction_rule_classes,
    instantiate_rule,
)

# Design stubs / bases: no metabolite emit required yet.
_SKIP_EXAMPLE_REQUIREMENT: frozenset[type[ReactionRule]] = frozenset(
    {
        *PATTERNLESS_REACTION_RULE_BASES,
        TautomerRule,
    }
)


def _concrete_rules() -> list[type[ReactionRule]]:
    return [
        cls
        for cls in discover_reaction_rule_classes()
        if cls not in PATTERNLESS_REACTION_RULE_BASES
    ]


@pytest.mark.parametrize(
    "rule_cls",
    _concrete_rules(),
    ids=lambda c: c.__name__,
)
def test_site_kind_matches_emitted_sites(rule_cls: type[ReactionRule]) -> None:
    kind = getattr(rule_cls, "site_kind", None)
    assert kind in ("atom", "atom_pair"), (
        f"{rule_cls.__name__} must declare site_kind as 'atom' or 'atom_pair', got {kind!r}"
    )
    expect = 1 if kind == "atom" else 2

    examples = getattr(rule_cls, "_example_substrates", ())
    if rule_cls in _SKIP_EXAMPLE_REQUIREMENT:
        return
    assert examples, (
        f"{rule_cls.__name__} lacks _example_substrates "
        "(add short SMILES that produce metabolites; expand later for all patterns/whens)"
    )

    rule = instantiate_rule(rule_cls)
    saw_any = False
    for smiles in examples:
        mol = MolFromSmiles(smiles)
        assert mol is not None, f"{rule_cls.__name__}: bad example SMILES {smiles!r}"
        rows = list(rule.metabolites(mol))
        assert rows, (
            f"{rule_cls.__name__}: example {smiles!r} produced no metabolites"
        )
        for por in rows:
            site = por.info["site"]
            assert len(site) == expect, (
                f"{rule_cls.__name__} site_kind={kind!r} but "
                f"{smiles!r} emitted site={site!r} (len {len(site)})"
            )
            saw_any = True
    assert saw_any
