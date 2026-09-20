"""Port of ``tests/test_rules.py``: each historical example must appear as a product.

SMILES stay in ``test_rules.examples``. Tautomerization is out of scope for the
forest (see TODO.md / docs/forest/DROPPED.md). Conjugation examples keep full adducts
(``as_star=False``), matching the forest suite.
"""

from __future__ import annotations

import pytest
from test_rules import examples

from xenosite.forest import rules as forest_rules

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

# Hard rules-example xfails still open after the cleavage XPASS clear
# (b2830ee / post-v0.6.1 port), emptied when kekulé parents landed. Exactly
# the 28 ``_XFAIL_IDS`` from that baseline — not the earlier quinone Reaction*
# bulk and not the early cleavage XPASS set.
_CLEARED_0_7: frozenset[str] = frozenset(
    {
        "dealkylation-Reaction94457_SRT",
        "dehydration-Reaction7389_11214_Dehydration_SRT",
        "dehydration-Reaction7389_11214_Step2",
        "dehydrogenation-Reaction3183_3757_Already_Hydroxylated1_SRT",
        "dehydrogenation-Reaction3183_3757_Already_Hydroxylated2_SRT",
        "hydrogenation-BigMolHydrogenation_SRT",
        "hydrogenation-NAPQI_Reduction_SR",
        "hydrogenation-NAPQI_Reduction_SRT",
        "hydrogenation-Reaction1224_SRT",
        "hydroxylation-Delavirdine_Hydroxylation",
        "hydroxylation-Reaction8492_99558_SRT",
        "nitrogenreduction-NitrogenReductionEx1_SRT",
        "nitrogenreduction-Reaction375_38982_Step1",
        "nitrogenreduction-Reaction85420_SRT",
        "oxygenreduction-Reaction7389_11214_Reduction_SRT",
        "quinoneformation-Epoxidation_SR_Num",
        "quinoneformation-LongRangeQuinone_SRT",
        "quinoneformation-NAPQI_Formation_SRT",
        "quinoneformation-NAPQI_QuinoneFormation_Site_Debug_SR",
        "quinoneformation-Raloxifene_SRT",
        "quinoneformation-Reaction61536_SRT",
        "quinoneformation-Reaction7751_SRT",
        "quinoneformation-Reaction9370_SRT",
        "sulfation-Reaction_162_35034_SRT",
        "sulfuroxidation-Reaction3183_3757_SulfurOxidation2_SRT",
        "sulfuroxidation-Reaction3183_3757_SulfurOxidation3_SRT",
        "sulfuroxidation-Reaction3183_3757_SulfurOxidation4_SRT",
        "sulfuroxidation-Reaction3183_3757_SulfurOxidation_SRT",
    }
)


def _case_id(rule: str, name: str) -> str:
    return f"{rule.lower()}-{name}"


def _params():
    out = []
    for rule, name, reactant, product in _CASES:
        cid = _case_id(rule, name)
        marks = [pytest.mark.cleared_0_7] if cid in _CLEARED_0_7 else []
        out.append(
            pytest.param(rule, name, reactant, product, id=cid, marks=marks)
        )
    return out


@pytest.mark.parametrize("rule", sorted(set(c[0] for c in _CASES)))
def test_rule_class_exists(rule):
    assert hasattr(forest_rules, rule), rule
    assert callable(getattr(forest_rules, rule))


@pytest.mark.parametrize(
    "rule, name, reactant, product",
    _params(),
)
def test_rule_emits_historical_product(rule, name, reactant, product):
    # progression: forest xfails this tetrabromo substitution under current
    # RDKit. Forest emits it.
    cls = getattr(forest_rules, rule)
    if rule in _CONJUGATION:
        instance = cls(as_star=False)
    else:
        instance = cls()

    # Historical product must appear as a yielded (single-component) mol.
    if emits_product(instance, reactant, product):
        return

    assert False, f"Failed to find {product} in {reactant} ({name})"
