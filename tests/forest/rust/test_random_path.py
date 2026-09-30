"""Product-gate tests for ``random_path``.

Requires the ``xenosite_forest`` extension (no skip).
"""

from __future__ import annotations

from xenosite.forest import available, hydroxylation, phase_one, random_path

assert available(), (
    "xenosite_forest extension required for tests/forest/rust; "
    "maturin develop -m crates/xenosite-forest/Cargo.toml "
    "--features python,extension-module"
)


def _steps(out):
    return out.to_dict()["steps"]


def _patterns(out):
    return out.to_dict()["patterns"]


def test_random_path_ethane_deterministic():
    a = random_path("CC", 1, max_steps=1, ruleset=hydroxylation())
    b = random_path("CC", 1, max_steps=1, ruleset=hydroxylation())
    assert a == b
    assert a.smiles in {"CCO", "C(C)O"}
    assert len(_steps(a)) == 1
    assert _steps(a)[0]["rule"] == "Hydroxylation"
    assert a.path[0] == "CC"
    assert a.path[-1] == a.smiles
    assert len(_patterns(a)) == 1
    assert _patterns(a)[0]["name"]


def test_random_path_same_seed_identical():
    rules = phase_one()
    a = random_path("COc1ccccc1", 42, max_steps=3, ruleset=rules)
    b = random_path("COc1ccccc1", 42, max_steps=3, ruleset=rules)
    assert a == b
    assert _steps(a), "expected a non-empty walk"


def test_random_path_different_seeds_diverge():
    rules = phase_one()
    a = random_path("COc1ccccc1", 1, max_steps=3, ruleset=rules)
    b = random_path("COc1ccccc1", 99, max_steps=3, ruleset=rules)
    assert a != b
    assert len(a.path) == len(_steps(a)) + 1
    assert len(b.path) == len(_steps(b)) + 1


def test_random_path_skip_filters_opt_in():
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
    assert _steps(open_)
    assert _steps(filtered)
    assert "." not in filtered.smiles


def test_random_path_zero_steps():
    out = random_path("CC", 1, max_steps=0, ruleset=hydroxylation())
    assert _steps(out) == []
    assert out.smiles == "CC"
    assert out.path == ["CC"]
