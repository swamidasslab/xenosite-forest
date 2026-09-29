"""Product-gate tests for the chematic find_path door.

Imports only ``xenosite.forest`` so missing stub re-exports fail here.
Requires the ``xenosite.forest._rust`` extension (no skip).
"""

from __future__ import annotations

import pytest

from xenosite.forest import PhaseOne, available, find_path, normalize_tautomer

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
    assert counters["timed_out"] is False


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


def test_find_path_timeout_zero():
    hits, counters = find_path("CC", "CCO", timeout=0.0)
    assert counters["timed_out"] is True
    assert hits == []


def test_normalize_tautomer_enol_keto_agree():
    enol, changed_e = normalize_tautomer("OC=C")
    keto, changed_k = normalize_tautomer("CC=O")
    assert changed_e
    assert not changed_k
    assert enol.csmi == keto.csmi


def test_forest_mol_normalize_tautomer_method():
    keto, _ = normalize_tautomer("CC=O")
    again, changed = keto.normalize_tautomer()
    assert not changed
    assert again.csmi == keto.csmi


def test_find_path_normalize_tautomer_enol_keto():
    # Emit is normalized space — judge against normalized target CSMI.
    want, _ = normalize_tautomer("CC=O")
    hits, _ = find_path("OC=C", "CC=O", max_paths=1, max_nodes=50)
    assert len(hits) == 1
    assert hits[0]["steps"] == []
    assert hits[0]["smiles"] == want.csmi


def test_find_path_invert_target_tautomer_not_implemented():
    with pytest.raises(Exception, match="invert_target_tautomer|not implemented"):
        find_path(
            "CC=O",
            "OC=C",
            max_paths=1,
            max_nodes=50,
            invert_target_tautomer=True,
        )


def test_phase_one_factory_from_stub():
    rules = PhaseOne()
    assert rules is not None
    assert len(rules) >= 1
