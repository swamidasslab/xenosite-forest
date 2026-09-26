"""Parametric soft-failure detector: ``unique_csmi`` yield drops (C11).

A rule should dedup from sites / patterns / ``when`` alone. When
``metabolize(..., unique_csmi=True)`` yields fewer emissions than
``unique_csmi=False``, the yield layer papered over a unique-edit or
partition miss — a **soft failure** (same family as sanitize, C10).

``SiteDeduplicationWarning`` is the check-layer unique-edit-miss signal
(same pattern + CSMI set + site ranks). It often accompanies a yield drop.

Known hits may use temporary ``pytest.xfail`` until partition / unique-edit
fixes land. An xfail here is not approval — remove the entry when green.
See ``docs/forest/RUST_PYTHON_PARITY.md`` Choices C11 / C12.
"""

from __future__ import annotations

import warnings

import pytest

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import ReactionRule, SiteDeduplicationWarning

from .pattern_info_inventory import instantiate_rule
from .rule_parity_corpus import PARITY_FUZZ_MOLS
from .rule_parity_pairs import python_leaf_classes, python_parity_exception

# Formula-delta suite gate is orthogonal; this sweep only asserts CSMI soft fails.
pytestmark = pytest.mark.allow_formula_delta_mismatch

# Temporary soft failures (C11). Drop an entry when unique-edit / pattern
# partition / C13 product-equiv handling stops relying on yield CSMI.
_XFAIL_CSMI_DEDUP: frozenset[tuple[str, str]] = frozenset(
    {
        # C13 class A: product-identical distinct directed sites (ester
        # quaternary_alcohol on both carbons of bridging O). Unequal ranks —
        # not SiteDeduplicationWarning; quiet unique_csmi drop only.
        ("Dealkylation", "CC(=O)Oc1ccccc1C(=O)O"),
        # Unique-edit miss + yield drop (halide/At geminal overlaps).
        ("OxidativeDehalogenation", "ClC(I)Cl"),
        ("OxidativeDehalogenation", "[At]C(Cl)[At]"),
        ("OxidativeDehalogenation", "ClC([At])Cl"),
    }
)


def _leaf_names() -> list[str]:
    out: list[str] = []
    for name, cls in sorted(python_leaf_classes().items()):
        if python_parity_exception(cls):
            continue
        out.append(name)
    return out


_LEAVES = _leaf_names()


def _count_emissions(
    rule: ReactionRule, smiles: str, *, unique_csmi: bool
) -> tuple[int, int]:
    """Return ``(emission_count, SiteDeduplicationWarning count)``."""

    mol = MolFromSmiles(smiles)
    assert mol is not None, smiles
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        n = len(list(rule.metabolize(mol, unique_csmi=unique_csmi)))
    n_warn = sum(
        1 for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    )
    return n, n_warn


@pytest.mark.parametrize("rule_name", _LEAVES or ["Hydroxylation"])
@pytest.mark.parametrize("smiles", PARITY_FUZZ_MOLS)
def test_no_unique_csmi_yield_drop(rule_name: str, smiles: str) -> None:
    """``unique_csmi`` must not drop emissions (sites/patterns/whens suffice)."""

    if not _LEAVES:
        pytest.skip("no Python leaf rules")

    cls = python_leaf_classes()[rule_name]
    rule = instantiate_rule(cls)
    expect_soft = (rule_name, smiles) in _XFAIL_CSMI_DEDUP

    try:
        n_off, warn_off = _count_emissions(rule, smiles, unique_csmi=False)
        n_on, warn_on = _count_emissions(rule, smiles, unique_csmi=True)

        if n_off > n_on:
            raise AssertionError(
                f"unique_csmi yield drop for {rule_name} on {smiles!r}: "
                f"unique_csmi=False → {n_off} emissions, True → {n_on} "
                f"(SiteDeduplicationWarning off={warn_off} on={warn_on}). "
                f"Dedup should come from sites/patterns/whens (C11)."
            )

        # Check layer can warn even when yield keeps both (unique_csmi=False).
        # Any warning is still a unique-edit miss soft failure.
        if warn_off or warn_on:
            raise AssertionError(
                f"SiteDeduplicationWarning for {rule_name} on {smiles!r}: "
                f"off={warn_off} on={warn_on} (unique-edit miss; C11)."
            )
    except AssertionError:
        if expect_soft:
            pytest.xfail(
                "C11 soft failure: unique_csmi / SiteDeduplicationWarning; "
                "fix unique-edit/partition then remove from _XFAIL_CSMI_DEDUP"
            )
        raise

    if expect_soft:
        pytest.fail(
            f"{rule_name} on {smiles!r} no longer hits unique_csmi soft failure — "
            f"remove from _XFAIL_CSMI_DEDUP"
        )
