"""Metabolic Forest public API (Rust chematic door).

Names match ``xenosite.forest._rust`` / ``crates/xenosite-forest`` exports::

    from xenosite.forest import find_path, phase_one, MetabolicNetwork, product_graph_into

``find_path`` / ``find_path_partial`` both return ``(outcomes, counters, network)``
— a list of ``PathOutcome`` (``residual_cost == 0`` when exact; only
``find_path_partial`` adds closest misses with ``residual_cost > 0``). Search
always records hops on a ``MetabolicNetwork`` (pass ``network=`` to extend one).

See ``docs/forest/RUST.md``. Do not extend :mod:`xenosite.forest.native` (frozen).
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING, Any

try:
    __version__ = version("xenosite-forest")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

from ._ext import available, load
from .mol import ForestMol

if TYPE_CHECKING:
    # Signatures come from the generated ``_rust.pyi`` (``make stubs``);
    # at runtime these resolve lazily through ``__getattr__``.
    from ._rust import (
        BoundPattern,
        Emission,
        MetabolicNetwork,
        PathCounters,
        PathOutcome,
        RandomPathOutcome,
        RuleSet,
        StepPlan,
        dealkylation,
        default_ruleset,
        dehydrogenation,
        epoxidation,
        epoxide_opening,
        expand_iri,
        find_path,
        find_path_partial,
        forest_xmet_sssom,
        hydrolysis,
        hydroxylation,
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
    from ._rust import (
        GraphEdge as GraphEdge,
    )
    from ._rust import (
        GraphNode as GraphNode,
    )

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

# Public names served from the extension (everything above but the Python
# helpers), plus graph views reachable from ``MetabolicNetwork``.
_RUST_EXPORTS = frozenset(__all__) - {"__version__", "available", "ForestMol"} | {
    "GraphNode",
    "GraphEdge",
}


def __getattr__(name: str) -> Any:
    """Lazy re-export of ``_rust`` names (types come from ``_rust.pyi``)."""

    if name in _RUST_EXPORTS:
        return getattr(load(), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
