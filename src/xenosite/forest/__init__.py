"""Metabolic Forest public API (Rust chematic door).

Recommended imports::

    from xenosite.forest import find_path, random_path, PhaseOne, Epoxidation

**New features** are implemented in ``crates/xenosite-forest`` and exposed
here. Do **not** add them to :mod:`xenosite.forest.native` (frozen RDKit
reference) or :mod:`xenosite.forest.legacy` (frozen 0.6.x archive). See
``docs/forest/NATIVE.md``.

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
    "random_path",
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
    timeout: float | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run PhaseOne chematic ``find_path``; return ``(hits, counters)``.

    Each hit is ``{"smiles": str, "steps": [{"rule": str, "site": list[str]}]}``.
    ``timeout`` is an optional wall-clock budget in seconds; counters include
    ``timed_out``.

    This is the Rust product door. New search features belong in the Rust crate
    and this wrapper — not in :mod:`xenosite.forest.native`.
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
        timeout=timeout,
    )


def random_path(
    reactant: str,
    seed: int,
    *,
    max_steps: int = 1,
    ruleset: Any | None = None,
    skip_multicomponent: bool = False,
    skip_seen: bool = False,
) -> dict[str, Any]:
    """Seeded random walk over a ruleset (default PhaseOne).

    Returns ``{"smiles", "path", "steps", "patterns"}``. ``path`` is reactant
    CSMI then each chosen product. Each step is
    ``{"rule", "pattern", "site", "products", "chosen"}``. Same ``seed`` is
    deterministic; different seeds diverge on a rich ruleset.

    ``skip_multicomponent`` / ``skip_seen`` are shared pathway filters (off by
    default); the same options will apply to StepSequence.apply.

    Implemented in Rust; not available on :mod:`xenosite.forest.native`.
    """

    kwargs: dict[str, Any] = {
        "max_steps": max_steps,
        "skip_multicomponent": skip_multicomponent,
        "skip_seen": skip_seen,
    }
    if ruleset is not None:
        kwargs["ruleset"] = ruleset
    return load().random_path(reactant, seed, **kwargs)


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
