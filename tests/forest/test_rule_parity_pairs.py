"""Meta-test: every Rust and Python leaf is paired, or excused on the rule.

Exception reasons are **rule attributes** (``rust_parity_exception`` /
``parity_exception``). This file does not maintain an allowlist of names.
"""

from __future__ import annotations

import pytest

from xenosite.forest.find_path_rust import native_available

from .rule_parity_pairs import (
    paired_rule_names,
    pairing_gaps,
    python_leaf_classes,
    rust_catalog_names,
)

pytestmark = pytest.mark.skipif(
    not native_available(),
    reason="xenosite-forest-native not installed",
)


def test_discover_python_and_rust_rules_nonempty():
    py = python_leaf_classes()
    rs = rust_catalog_names()
    assert py, "no Python leaf rules discovered"
    assert rs, "no Rust catalog names"


def test_every_rule_paired_or_has_attribute_exception():
    gaps = pairing_gaps()
    assert not gaps, "unpaired rules without attribute exceptions:\n" + "\n".join(
        f"  {g}" for g in gaps
    )


def test_paired_rules_exist_for_product_parity():
    paired = paired_rule_names()
    assert paired, "expected at least one paired leaf for product parity"
    # Spot-check a core Phase-I leaf is in the paired set.
    assert "Hydroxylation" in paired
