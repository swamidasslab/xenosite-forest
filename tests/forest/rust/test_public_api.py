"""Product API must be importable from ``xenosite.forest`` (not ``.rust``).

These imports are the gate: if a name is missing from the public stub,
this suite fails even when the maturin module still exports it.
"""

from __future__ import annotations

import importlib

import pytest

from xenosite.forest import (  # noqa: F401 — import gap gate
    BoundPattern,
    Emission,
    ForestMol,
    MetabolicNetwork,
    PathCounters,
    PathOutcome,
    RandomPathOutcome,
    RuleSet,
    StepPlan,
    __version__,
    available,
    default_ruleset,
    epoxidation,
    epoxide_opening,
    expand_iri,
    find_path,
    find_path_partial,
    forest_xmet_sssom,
    leaf_rule,
    n_dealkylation,
    normalize_tautomer,
    phase_one,
    product_graph_bfs,
    product_graph_into,
    product_graph_ruleset,
    quinone_formation,
    random_path,
    reactivity,
    resolve,
    to_curie,
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
        "random_path",
        "product_graph_bfs",
        "product_graph_into",
        "MetabolicNetwork",
        "ForestMol",
        "PathOutcome",
        "PathCounters",
        "RandomPathOutcome",
        "Emission",
        "StepPlan",
        "RuleSet",
        "BoundPattern",
        "phase_one",
        "reactivity",
        "default_ruleset",
        "product_graph_ruleset",
        "leaf_rule",
        "epoxidation",
        "quinone_formation",
        "epoxide_opening",
        "n_dealkylation",
        "hydroxylation",
        "dehydrogenation",
        "dealkylation",
        "hydrolysis",
        "resolve",
        "forest_xmet_sssom",
        "expand_iri",
        "to_curie",
    ]
    for name in forest.__all__:
        assert hasattr(forest, name), name


def test_product_graph_ruleset_includes_qf_eh_not_tautomer():
    rs = product_graph_ruleset()
    assert rs.name == "ProductGraph"
    assert len(rs) == 3
    assert "QuinoneFormation" in rs
    assert "EpoxideHydration" in rs
    assert "PhaseOne" in rs
    assert "Tautomerization" not in rs


def test_ruleset_factories_callable_from_stub():
    assert callable(find_path)
    assert callable(find_path_partial)
    assert callable(normalize_tautomer)
    assert callable(random_path)
    assert callable(product_graph_bfs)
    assert callable(product_graph_into)
    assert callable(available)
    assert callable(phase_one)
    assert callable(default_ruleset)
    assert callable(product_graph_ruleset)
    assert callable(epoxidation)
    assert callable(quinone_formation)
    assert callable(epoxide_opening)
    assert callable(n_dealkylation)
    assert MetabolicNetwork is not None
    assert isinstance(__version__, str)
