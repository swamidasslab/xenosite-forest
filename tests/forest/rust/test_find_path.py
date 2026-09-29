"""Product-gate tests for the chematic find_path door.

Imports only ``xenosite.forest`` so missing stub re-exports fail here.
Requires the ``xenosite.forest._rust`` extension (no skip).
"""

from __future__ import annotations

from xenosite.forest import PhaseOne, available, find_path

assert available(), (
    "xenosite.forest._rust extension required for tests/forest/rust; "
    "uv sync (maturin build backend) or maturin develop "
    "--features python,extension-module"
)


def test_find_path_ethane_to_ethanol():
    hits, counters = find_path("CC", "CCO", max_paths=1)
    assert len(hits) == 1
    assert hits[0]["steps"][0]["rule"] == "Hydroxylation"
    assert counters["billed"] >= 1
    assert counters["diversity_repush"] == 0


def test_find_path_diversity_opt_in():
    hits, counters = find_path(
        "COc1ccccc1",
        "Oc1ccccc1",
        max_paths=1,
        diversity=True,
    )
    assert len(hits) == 1
    assert "O" in hits[0]["smiles"] or "c" in hits[0]["smiles"]
    assert "diversity_repush" in counters


def test_phase_one_factory_from_stub():
    rules = PhaseOne()
    assert rules is not None
    assert len(rules) >= 1
