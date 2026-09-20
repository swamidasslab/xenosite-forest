"""Port of ``tests/test_rules.py``: each historical example must appear as a product.

SMILES stay in ``test_rules.examples``. Tautomerization is out of scope for the
poc (see TODO.md / DROPPED.md). Conjugation examples keep full adducts
(``as_star=False``), matching the forest suite.
"""

from __future__ import annotations

import pytest

from test_rules import examples
from xenosite.refactor_poc import rules as poc_rules

from .helpers import emits_product

_CONJUGATION = frozenset(
    {"Acetylation", "Glucuronidation", "Glutathionation", "Sulfation"}
)

_CASES = [
    (rule, name, reactant, product)
    for rule, rows in examples.items()
    if rule != "Tautomerization"
    for name, reactant, product in rows
]


def _case_id(rule: str, name: str) -> str:
    return f"{rule.lower()}-{name}"


def _params():
    return [
        pytest.param(rule, name, reactant, product, id=_case_id(rule, name))
        for rule, name, reactant, product in _CASES
    ]


@pytest.mark.parametrize("rule", sorted(set(c[0] for c in _CASES)))
def test_rule_class_exists(rule):
    assert hasattr(poc_rules, rule), rule
    assert callable(getattr(poc_rules, rule))


@pytest.mark.parametrize(
    "rule, name, reactant, product",
    _params(),
)
def test_rule_emits_historical_product(rule, name, reactant, product):
    # progression: forest xfails this tetrabromo substitution under current
    # RDKit. Poc emits it.
    cls = getattr(poc_rules, rule)
    if rule in _CONJUGATION:
        instance = cls(as_star=False)
    else:
        instance = cls()

    # Historical product must appear as a yielded (single-component) mol.
    if emits_product(instance, reactant, product):
        return

    assert False, f"Failed to find {product} in {reactant} ({name})"
