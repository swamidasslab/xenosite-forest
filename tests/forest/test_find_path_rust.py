"""Optional native ``find_path_rust`` door (requires xenosite-forest-native)."""

from __future__ import annotations

import pytest

from xenosite.forest.find_path_rust import find_path_rust, native_available

pytestmark = pytest.mark.skipif(
    not native_available(),
    reason="xenosite-forest-native not installed",
)


def test_find_path_rust_ethane_to_ethanol():
    hits, counters = find_path_rust("CC", "CCO", max_paths=1)
    assert len(hits) == 1
    assert hits[0]["steps"][0]["rule"] == "Hydroxylation"
    assert counters["billed"] >= 1
    assert counters["diversity_repush"] == 0


def test_find_path_rust_diversity_opt_in():
    hits, counters = find_path_rust(
        "COc1ccccc1",
        "Oc1ccccc1",
        max_paths=1,
        diversity=True,
    )
    assert len(hits) == 1
    assert "O" in hits[0]["smiles"] or "c" in hits[0]["smiles"]
    assert "diversity_repush" in counters
