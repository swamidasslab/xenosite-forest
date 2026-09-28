"""Metabolic Forest public API (Rust chematic door).

Recommended imports::

    from xenosite.forest import find_path, PhaseOne, Epoxidation

The RDKit reference engine is :mod:`xenosite.forest.native` (optional
``[rdkit]`` extra). The frozen 0.6.x archive is :mod:`xenosite.forest.legacy``.
"""

from __future__ import annotations

from typing import Any

try:
    from ._version import __version__
except ImportError:  # pragma: no cover - missing only in incomplete checkouts
    __version__ = "0.0.0"

from ._ext import available, load

__all__ = [
    "__version__",
    "available",
    "find_path",
    "PhaseOne",
    "Epoxidation",
    "QuinoneFormation",
    "EpoxideOpening",
    "NDealkylation",
]


def find_path(
    reactant: str,
    target: str,
    *,
    max_paths: int = 1,
    max_nodes: int = 800,
    use_atom_diff: bool = True,
    lazy_closer: bool = False,
    diversity: bool = False,
    drop_skeleton_twins: bool = True,
    score: str = "log-neg-pc",
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run PhaseOne chematic ``find_path``; return ``(hits, counters)``.

    Each hit is ``{"smiles": str, "steps": [{"rule": str, "site": list[str]}]}``.
    """

    return load().find_path(
        reactant,
        target,
        max_paths=max_paths,
        max_nodes=max_nodes,
        use_atom_diff=use_atom_diff,
        lazy_closer=lazy_closer,
        diversity=diversity,
        drop_skeleton_twins=drop_skeleton_twins,
        score=score,
    )


def _ruleset(attr: str, label: str) -> Any:
    """Return a named ruleset factory result from the extension."""

    mod = load()
    factory = getattr(mod, attr, None)
    if factory is None:
        raise AttributeError(
            f"xenosite_forest has no {attr!r}; rebuild xenosite-forest-native "
            f"to use {label}"
        )
    return factory()


def PhaseOne() -> Any:
    """Composed Phase I ruleset (Rust)."""

    return _ruleset("phase_one", "PhaseOne")


def Epoxidation() -> Any:
    """Epoxidation leaf ruleset (Rust)."""

    return _ruleset("epoxidation", "Epoxidation")


def QuinoneFormation() -> Any:
    """QuinoneFormation leaf ruleset (Rust)."""

    return _ruleset("quinone_formation", "QuinoneFormation")


def EpoxideOpening() -> Any:
    """EpoxideOpening leaf ruleset (Rust)."""

    return _ruleset("epoxide_opening", "EpoxideOpening")


def NDealkylation() -> Any:
    """NDealkylation leaf ruleset (Rust)."""

    return _ruleset("n_dealkylation", "NDealkylation")


# Unexported catalog helpers (advanced / growing tests). Not in ``__all__``.
def hydroxylation() -> Any:
    return _ruleset("hydroxylation", "Hydroxylation")


def dehydrogenation() -> Any:
    return _ruleset("dehydrogenation", "Dehydrogenation")


def dealkylation() -> Any:
    return _ruleset("dealkylation", "Dealkylation")


def hydrolysis() -> Any:
    return _ruleset("hydrolysis", "Hydrolysis")


def default_ruleset() -> Any:
    return _ruleset("default_ruleset", "Default")
