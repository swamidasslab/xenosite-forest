"""FormulaDeltaMismatchWarning when product Δformula ≠ declared delta."""

from __future__ import annotations

import warnings

import pytest

from xenosite.forest.find_path import PathCounters
from xenosite.forest.rdkitutil import Smirks, mol_from_smiles
from xenosite.forest.rules import (
    FormulaDeltaMismatchWarning,
    Hydroxylation,
    SmirksReactionRule,
    _describe,
)


def test_hydroxylation_does_not_warn_on_heavy_delta():
    mol = mol_from_smiles("CC")
    rule = Hydroxylation()
    counters = PathCounters()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", FormulaDeltaMismatchWarning)
        list(rule.metabolize(mol, counters=counters))
    mismatches = [
        w for w in caught if issubclass(w.category, FormulaDeltaMismatchWarning)
    ]
    assert mismatches == []
    assert counters.formula_delta_mismatch == 0
    assert counters.formula_delta_mismatches == []


@pytest.mark.allow_formula_delta_mismatch
def test_wrong_delta_formula_warns_and_counts():
    class BogusTooMuchO(SmirksReactionRule):
        name = "bogus"
        site_kind = "atom"
        smirks = (
            (
                Smirks("[#6h3:1]>>[*:1]O"),
                _describe(adds="OO", removes="H", name="too_much_o"),
            ),
        )

    mol = mol_from_smiles("CC")
    counters = PathCounters()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", FormulaDeltaMismatchWarning)
        list(BogusTooMuchO().metabolize(mol, counters=counters))
    mismatches = [
        w for w in caught if issubclass(w.category, FormulaDeltaMismatchWarning)
    ]
    assert mismatches, "expected FormulaDeltaMismatchWarning for OO vs +O product"
    assert "too_much_o" in str(mismatches[0].message)
    assert counters.formula_delta_mismatch >= 1
    assert counters.formula_delta_mismatches
    detail = counters.formula_delta_mismatches[0]
    assert detail.pattern == "too_much_o"
    assert detail.declared_heavy.get("O") == 2
    assert detail.observed_heavy.get("O") == 1
