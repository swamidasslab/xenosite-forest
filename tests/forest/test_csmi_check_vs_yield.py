"""Dedup check vs yield: warn without implying drop; yield only if unique_csmi."""

from __future__ import annotations

import logging
import warnings
from collections.abc import Generator
from typing import Any

from rdkit import Chem

from xenosite.forest.rdkitutil import Mol
from xenosite.forest.records import PatternInfo
from xenosite.forest.rules import (
    AzoSplitting,
    SiteDeduplicationWarning,
    ProductsOfReaction,
    ReactionRule,
    RuleSiteKind,
    SmirksReactionRule,
    describe,
)
from xenosite.forest.rulesets import RuleSet

_RULES_LOG = "xenosite.forest.rules"


class OverlapOhA(SmirksReactionRule):
    site_kind: RuleSiteKind = "atom"
    smirks: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h3:1]>>[*:1]O",
            describe(adds="O", removes="H", name="oh", site_map=1),
        ),
    )


class OverlapOhB(SmirksReactionRule):
    site_kind: RuleSiteKind = "atom"
    smirks: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h3:1]>>[*:1]O",
            describe(adds="O", removes="H", name="oh", site_map=1),
        ),
    )


class DoubleEmitOh(SmirksReactionRule):
    """Same hydroxylator, but metabolites re-emits the first POR (unique-edit miss)."""

    site_kind: RuleSiteKind = "atom"
    smirks: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h3:1]>>[*:1]O",
            describe(adds="O", removes="H", name="oh", site_map=1),
        ),
    )

    def metabolites(  # type: ignore[override]
        self,
        mol: Mol,
        *args: Any,
        **kwargs: Any,
    ) -> Generator[ProductsOfReaction, None, None]:
        first: ProductsOfReaction | None = None
        for por in super().metabolites(mol, *args, **kwargs):
            if first is None:
                first = por
                yield por
                yield por
            else:
                yield por


class CaptureChildUniqueCsmi(OverlapOhA):
    """Records the ``unique_csmi`` flag RuleSet passes to children."""

    seen: list[bool] = []

    def metabolize(  # type: ignore[override]
        self,
        mol: Mol,
        *args: Any,
        unique_csmi: bool = True,
        **kwargs: Any,
    ):
        type(self).seen.append(unique_csmi)
        yield from super().metabolize(
            mol, *args, unique_csmi=unique_csmi, **kwargs
        )


def _flat_csmis(rows: list) -> list[str]:
    return [p.xf.csmi for products, _ in rows for p in products]


def test_azo_cleavage_keeps_sibling_fragments_in_one_list() -> None:
    """Azo cleavage: one emission list with both aniline fragments (emission unique_csmi)."""

    mol = Chem.MolFromSmiles("c1ccc(N=Nc2ccccc2)cc1")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        rows = list(AzoSplitting().metabolize(mol, unique_csmi=True))
    assert not [
        w for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    ]
    assert len(rows) == 1
    products, info = rows[0]
    assert len(products) == 2
    assert [p.xf.csmi for p in products] == ["Nc1ccccc1", "Nc1ccccc1"]
    assert "csmi" not in info
    assert frozenset(p.xf.csmi for p in products) == frozenset({"Nc1ccccc1"})


def test_unique_csmi_false_still_packs_cleavage_siblings() -> None:
    """Yield off does not unwrap cleavage; still one list per emission."""

    mol = Chem.MolFromSmiles("c1ccc(N=Nc2ccccc2)cc1")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        rows = list(AzoSplitting().metabolize(mol, unique_csmi=False))
    assert not [
        w for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    ]
    assert len(rows) == 1
    assert _flat_csmis(rows).count("Nc1ccccc1") == 2


def test_check_warns_with_yield_off_on_unique_edit_miss() -> None:
    """Site unique-edit miss still warns when unique_csmi=False (dups yield)."""

    mol = Chem.MolFromSmiles("CC")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        rows = list(DoubleEmitOh().metabolize(mol, unique_csmi=False))
    csmi_warns = [
        w for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    ]
    assert csmi_warns, "expected SiteDeduplicationWarning from re-emitted POR"
    assert _flat_csmis(rows).count("CCO") == 2


def test_check_warns_and_yield_on_still_drops() -> None:
    """Miss warns; unique_csmi=True still suppresses the duplicate emission."""

    mol = Chem.MolFromSmiles("CC")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        rows = list(DoubleEmitOh().metabolize(mol, unique_csmi=True))
    assert [w for w in caught if issubclass(w.category, SiteDeduplicationWarning)]
    assert _flat_csmis(rows) == ["CCO"]


def test_ruleset_forces_child_yield_off() -> None:
    """RuleSet does not forward caller unique_csmi to children for yield."""

    CaptureChildUniqueCsmi.seen = []
    ruleset = RuleSet((CaptureChildUniqueCsmi,), name="Cap")
    list(ruleset.metabolize(Chem.MolFromSmiles("CC"), unique_csmi=True))
    assert CaptureChildUniqueCsmi.seen == [False]


def test_ruleset_cross_rule_bubble_then_parent_yield(caplog) -> None:
    """Both overlapping children reach parent; parent yield keeps first rule."""

    ruleset = RuleSet((OverlapOhA, OverlapOhB), name="OverlapSet")
    mol = Chem.MolFromSmiles("CC")
    with caplog.at_level(logging.INFO, logger=_RULES_LOG):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("error", SiteDeduplicationWarning)
            products = list(ruleset.metabolize(mol, unique_csmi=True))
    assert _flat_csmis(products) == ["CCO"]
    assert {type(info["rule"]).__name__ for _, info in products} == {
        "OverlapOhA"
    }
    assert not caught
    assert any(
        "kept_rule" in r.getMessage() and "dropped_rule" in r.getMessage()
        for r in caplog.records
        if r.name == _RULES_LOG and r.levelno == logging.INFO
    )


def test_ruleset_unique_csmi_false_keeps_both_rules() -> None:
    """Parent yield off: cross-rule same CSMI both emit (no parent drop)."""

    ruleset = RuleSet((OverlapOhA, OverlapOhB), name="OverlapSet")
    products = list(
        ruleset.metabolize(Chem.MolFromSmiles("CC"), unique_csmi=False)
    )
    names = [type(info["rule"]).__name__ for _, info in products]
    assert names == ["OverlapOhA", "OverlapOhB"]
    assert _flat_csmis(products) == ["CCO", "CCO"]


def test_nested_ruleset_outermost_yield_only(caplog) -> None:
    """Inner RuleSet yield-off from outer; outer uniquifies cross-child CSMI."""

    inner = RuleSet((OverlapOhA, OverlapOhB), name="Inner")
    outer = RuleSet((inner,), name="Outer")
    with caplog.at_level(logging.INFO, logger=_RULES_LOG):
        products = list(
            outer.metabolize(Chem.MolFromSmiles("CC"), unique_csmi=True)
        )
    # Inner forced unique_csmi=False by outer, so both children reach inner;
    # inner also forced False by outer, so both reach outer; outer yield keeps one.
    assert _flat_csmis(products) == ["CCO"]
    assert isinstance(products[0][1]["rule"], ReactionRule)
    assert any(
        "kept_rule" in r.getMessage() and "dropped_rule" in r.getMessage()
        for r in caplog.records
        if r.name == _RULES_LOG and r.levelno == logging.INFO
    )
