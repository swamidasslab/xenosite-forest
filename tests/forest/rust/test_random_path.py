"""Product-gate tests for ``random_path``.

Requires the ``xenosite_forest`` extension (no skip).
"""

from __future__ import annotations

from xenosite.forest import PhaseOne, available, hydroxylation, random_path

assert available(), (
    "xenosite_forest extension required for tests/forest/rust; "
    "maturin develop -m crates/xenosite-forest/Cargo.toml "
    "--features python,extension-module"
)


def test_random_path_ethane_deterministic():
    a = random_path("CC", 1, max_steps=1, ruleset=hydroxylation())
    b = random_path("CC", 1, max_steps=1, ruleset=hydroxylation())
    assert a == b
    assert a["smiles"] in {"CCO", "C(C)O"}
    assert len(a["steps"]) == 1
    assert a["steps"][0]["rule"] == "Hydroxylation"
    assert a["path"][0] == "CC"
    assert a["path"][-1] == a["smiles"]
    assert len(a["patterns"]) == 1
    assert a["patterns"][0]["name"]


def test_random_path_same_seed_identical():
    rules = PhaseOne()
    a = random_path("COc1ccccc1", 42, max_steps=3, ruleset=rules)
    b = random_path("COc1ccccc1", 42, max_steps=3, ruleset=rules)
    assert a == b
    assert a["steps"], "expected a non-empty walk"


def test_random_path_different_seeds_diverge():
    rules = PhaseOne()
    a = random_path("COc1ccccc1", 1, max_steps=3, ruleset=rules)
    b = random_path("COc1ccccc1", 99, max_steps=3, ruleset=rules)
    assert a != b
    assert len(a["path"]) == len(a["steps"]) + 1
    assert len(b["path"]) == len(b["steps"]) + 1


def test_random_path_skip_filters_opt_in():
    """Shared PathwayOptions: skip multicomponent / seen only when requested."""
    from xenosite.forest import dealkylation

    open_ = random_path(
        "COc1ccccc1",
        7,
        max_steps=1,
        ruleset=dealkylation(),
    )
    filtered = random_path(
        "COc1ccccc1",
        7,
        max_steps=1,
        ruleset=dealkylation(),
        skip_multicomponent=True,
        skip_seen=True,
    )
    assert open_["steps"]
    assert filtered["steps"]
    assert "." not in filtered["smiles"]


def test_random_path_zero_steps():
    out = random_path("CC", 1, max_steps=0, ruleset=hydroxylation())
    assert out["steps"] == []
    assert out["smiles"] == "CC"
    assert out["path"] == ["CC"]
