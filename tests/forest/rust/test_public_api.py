"""Product API must be importable from ``xenosite.forest`` (not ``.rust``).

These imports are the gate: if a name is missing from the public stub,
this suite fails even when the maturin module still exports it.
"""

from __future__ import annotations

import importlib

import pytest

from xenosite.forest import (  # noqa: F401 — import gap gate
    Default,
    Epoxidation,
    EpoxideOpening,
    MetabolicNetwork,
    NDealkylation,
    PhaseOne,
    QuinoneFormation,
    __version__,
    available,
    find_path,
    find_path_partial,
    normalize_tautomer,
    random_path,
)


def test_forest_rust_package_must_not_exist():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("xenosite.forest.rust")


def test_public_allowlist():
    import xenosite.forest as forest

    assert forest.__all__ == [
        "__version__",
        "available",
        "find_path",
        "find_path_partial",
        "normalize_tautomer",
        "MetabolicNetwork",
        "random_path",
        "PhaseOne",
        "Default",
        "Epoxidation",
        "QuinoneFormation",
        "EpoxideOpening",
        "NDealkylation",
        "resolve",
        "forest_xmet_sssom",
        "expand_iri",
        "to_curie",
        "BoundPattern",
        "RuleSet",
    ]
    for name in forest.__all__:
        assert hasattr(forest, name), name


def test_ruleset_factories_callable_from_stub():
    """Factories are reachable via the public stub (extension may be absent)."""

    assert callable(find_path)
    assert callable(find_path_partial)
    assert callable(normalize_tautomer)
    assert callable(random_path)
    assert callable(available)
    assert callable(PhaseOne)
    assert callable(Default)
    assert callable(Epoxidation)
    assert callable(QuinoneFormation)
    assert callable(EpoxideOpening)
    assert callable(NDealkylation)
    assert MetabolicNetwork is not None
    assert isinstance(__version__, str)


def test_sssom_exports_callable_from_stub():
    import xenosite.forest as forest

    assert callable(forest.resolve)
    assert callable(forest.forest_xmet_sssom)
    assert callable(forest.expand_iri)
    assert callable(forest.to_curie)
