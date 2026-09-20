"""Meta-test: every concrete rule's ``site_kind`` matches emitted Sites.

Uses each rule's internal ``_example_substrates`` (SMILES guaranteed to
produce metabolites). Expand those lists later so they cover all
patterns/whens — see TODO.md.

Taxonomy:
- ``atom`` — singleton frozenset
- ``bond`` — undirected bond frozenset (sorted site ranks in unique-edit)
- ``directed_bond`` — ordered tuple of atom indexes (map order)
- ``atom_pair`` — ResonancePair ends only (frozenset)
"""

from __future__ import annotations

import pytest

from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.rules import ReactionRule, ResonancePairRule, TautomerRule

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

_SITE_KIND_LEN = {
    "atom": 1,
    "bond": 2,
    "directed_bond": 2,
    "atom_pair": 2,
}


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
    assert kind in _SITE_KIND_LEN, (
        f"{rule_cls.__name__} must declare site_kind as "
        f"'atom', 'bond', 'directed_bond', or 'atom_pair', got {kind!r}"
    )
    if kind == "atom_pair":
        assert issubclass(rule_cls, ResonancePairRule), (
            f"{rule_cls.__name__} site_kind='atom_pair' but is not a "
            f"ResonancePairRule (atom_pair is pair-rule only; use 'bond' "
            f"or 'directed_bond' for SMARTS bond sites)"
        )
    expect = _SITE_KIND_LEN[kind]

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
            if kind == "directed_bond":
                assert isinstance(site, tuple), (
                    f"{rule_cls.__name__} directed_bond must emit tuple, "
                    f"got {type(site).__name__}: {site!r}"
                )
            else:
                assert isinstance(site, frozenset), (
                    f"{rule_cls.__name__} site_kind={kind!r} must emit frozenset, "
                    f"got {type(site).__name__}: {site!r}"
                )
            saw_any = True
    assert saw_any
