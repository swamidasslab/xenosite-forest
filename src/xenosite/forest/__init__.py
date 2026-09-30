"""Metabolic Forest public API (Rust chematic door).

Names match ``xenosite.forest._rust`` / ``crates/xenosite-forest`` exports::

    from xenosite.forest import find_path, phase_one, MetabolicNetwork, product_graph_into

``find_path`` / ``find_path_partial`` return ``(..., network)`` — search always
records hops on a ``MetabolicNetwork`` (pass ``network=`` to extend one).

See ``docs/forest/RUST.md``. Do not extend :mod:`xenosite.forest.native` (frozen).
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("xenosite-forest")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

from ._ext import available, load
from .mol import ForestMol

__all__ = [
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
    "PartialOutcome",
    "RandomPathOutcome",
    "Emission",
    "StepPlan",
    "RuleSet",
    "BoundPattern",
    "phase_one",
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


def find_path(*args: Any, **kwargs: Any) -> Any:
    return load().find_path(*args, **kwargs)


def find_path_partial(*args: Any, **kwargs: Any) -> Any:
    return load().find_path_partial(*args, **kwargs)


def normalize_tautomer(smiles: str) -> tuple[Any, bool]:
    return load().normalize_tautomer(smiles)


def random_path(*args: Any, **kwargs: Any) -> Any:
    return load().random_path(*args, **kwargs)


def product_graph_bfs(
    start: str,
    *,
    target: str | None = None,
    max_nodes: int = 256,
    max_depth: int = 6,
    ruleset: Any | None = None,
    network: Any | None = None,
) -> Any:
    """BFS product exploration; returns a :class:`MetabolicNetwork`.

    .. warning::

       This can take a **very long time** and build a large graph even when
       ``max_nodes`` and ``max_depth`` are set. Use ``target=`` when you are
       exploring toward a single product CSMI.

    Default ``ruleset`` is :func:`product_graph_ruleset` — QuinoneFormation,
    EpoxideHydration, and Phase I core, **without** Tautomerization (unlike
    :func:`find_path`'s :func:`default_ruleset`). Pass ``network=`` to extend
    an existing graph (same object as :func:`find_path` ``network=``).
    """

    kwargs: dict[str, Any] = {
        "max_nodes": max_nodes,
        "max_depth": max_depth,
    }
    if target is not None:
        kwargs["target"] = target
    if ruleset is not None:
        kwargs["ruleset"] = ruleset
    if network is not None:
        kwargs["network"] = network
    return load().product_graph_bfs(start, **kwargs)


def product_graph_into(*args: Any, **kwargs: Any) -> None:
    load().product_graph_into(*args, **kwargs)


def phase_one() -> Any:
    return load().phase_one()


def default_ruleset() -> Any:
    return load().default_ruleset()


def product_graph_ruleset() -> Any:
    return load().product_graph_ruleset()


def leaf_rule(name: str) -> Any:
    return load().leaf_rule(name)


def epoxidation() -> Any:
    return load().epoxidation()


def quinone_formation() -> Any:
    return load().quinone_formation()


def epoxide_opening() -> Any:
    return load().epoxide_opening()


def n_dealkylation() -> Any:
    return load().n_dealkylation()


def hydroxylation() -> Any:
    return load().hydroxylation()


def dehydrogenation() -> Any:
    return load().dehydrogenation()


def dealkylation() -> Any:
    return load().dealkylation()


def hydrolysis() -> Any:
    return load().hydrolysis()


def resolve(id: str) -> Any:
    return load().resolve(id)


def forest_xmet_sssom() -> str:
    return load().forest_xmet_sssom()


def expand_iri(curie_or_iri: str) -> str:
    return load().expand_iri(curie_or_iri)


def to_curie(iri: str) -> str:
    return load().to_curie(iri)


def __getattr__(name: str) -> Any:
    """Pyclasses re-exported from ``_rust`` (same names as the extension module)."""

    if name == "MetabolicNetwork":
        return load().MetabolicNetwork
    _rust_types = (
        "BoundPattern",
        "RuleSet",
        "PathOutcome",
        "PathCounters",
        "PartialOutcome",
        "RandomPathOutcome",
        "Emission",
        "StepPlan",
        "GraphNode",
        "GraphEdge",
    )
    if name in _rust_types:
        return getattr(load(), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
