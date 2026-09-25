"""FormulaDeltaMismatchWarning when product Δformula ≠ declared delta."""

from __future__ import annotations

import warnings

from xenosite.forest.find_path import PathCounters
from xenosite.forest.rdkitutil import Smirks, mol_from_smiles
from xenosite.forest.rules import (
    FormulaDeltaMismatchWarning,
    Hydroxylation,
    SmirksReactionRule,
    describe,
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


def test_wrong_delta_formula_warns_and_counts():
    class BogusTooMuchO(SmirksReactionRule):
        name = "bogus"
        site_kind = "atom"
        smirks = (
            (
                Smirks("[#6h3:1]>>[*:1]O"),
                describe(adds="OO", removes="H", name="too_much_o"),
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
