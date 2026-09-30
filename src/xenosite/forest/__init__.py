"""Metabolic Forest public API (Rust chematic door).

Recommended imports::

    from xenosite.forest import find_path, random_path, PhaseOne, Default, Epoxidation

**New features** are implemented in ``crates/xenosite-forest`` and exposed
here. Do **not** add them to :mod:`xenosite.forest.native` (frozen RDKit
reference) or :mod:`xenosite.forest.legacy` (frozen 0.6.x archive). See
``docs/forest/NATIVE.md``.

The RDKit reference engine is :mod:`xenosite.forest.native` (optional
``[rdkit]`` extra). The frozen 0.6.x archive is :mod:`xenosite.forest.legacy``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("xenosite-forest")
except PackageNotFoundError:  # pragma: no cover - incomplete checkout / editable edge
    __version__ = "0.0.0"

from ._ext import available, load

__all__ = [
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


def normalize_tautomer(smiles: str) -> tuple[Any, bool]:
    """Chematic tautomer pick adopted as a tagged ForestMol.

    Returns ``(ForestMol, changed)``.
    """

    return load().normalize_tautomer(smiles)


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
    network: Any | None = None,
    normalize_tautomer: bool = False,
    invert_target_tautomer: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run default-ruleset chematic ``find_path``; return ``(hits, counters)``.

    Default ruleset is QuinoneFormation + EpoxideHydration + Tautomerization +
    PhaseOne core. Each hit is
    ``{"smiles": str, "steps": [{"rule": str, "site": list[str]}]}``.
    ``timeout`` is an optional wall-clock budget in seconds; counters include
    ``timed_out``. Pass ``network=`` a :class:`MetabolicNetwork` to record hops.
    ``normalize_tautomer`` (default False) runs chematic zwitterion → remove
    explicit H → canonical tautomer on reactant and target once before search.
    Opt in when both ends should share a chematic preferred form. Emit stays
    in that normalized target form unless ``invert_target_tautomer`` (not
    implemented yet when the target changes).

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
        network=network,
        normalize_tautomer=normalize_tautomer,
        invert_target_tautomer=invert_target_tautomer,
    )


def find_path_partial(
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
    network: Any | None = None,
    normalize_tautomer: bool = False,
    invert_target_tautomer: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Exact hits plus end-of-search closest reaches when the target is missed.

    Returns ``(exact, partials, counters)``. Each partial includes ``residual``
    with ``cost`` / ``categories``. Prefer exact when ``max_paths`` is filled.
    """

    return load().find_path_partial(
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
        network=network,
        normalize_tautomer=normalize_tautomer,
        invert_target_tautomer=invert_target_tautomer,
    )


def MetabolicNetwork() -> Any:
    """Explored reactant+metabolite graph recorded by search ``network=``."""

    return load().MetabolicNetwork()



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
            f"xenosite.forest._rust has no {attr!r}; rebuild the extension "
            f"(maturin develop --features python,extension-module) to use {label}"
        )
    return factory()


def leaf_rule(name: str) -> Any:
    """Sealed leaf ``RuleSet`` by catalog name (Rust ``LEAF_CTORS``)."""

    return load().leaf_rule(name)


def PhaseOne() -> Any:
    """Composed Phase I ruleset (Rust). No Tautomerization; see :func:`Default`."""

    return _ruleset("phase_one", "PhaseOne")


def Default() -> Any:
    """Default find_path ruleset: QF + EpoxideHydration + Tautomer + PhaseOne."""

    return _ruleset("default_ruleset", "Default")


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


def resolve(id: str) -> Any:
    """Resolve an ``xf:`` CURIE / Forest IRI to a ``RuleSet`` or ``BoundPattern``."""

    return load().resolve(id)


def forest_xmet_sssom() -> str:
    """Decompressed Forest↔XMET SSSOM TSV text (embedded at build time)."""

    return load().forest_xmet_sssom()


def expand_iri(curie_or_iri: str) -> str:
    """Expand ``xf:`` / ``xmet:`` CURIEs to absolute IRIs."""

    return load().expand_iri(curie_or_iri)


def to_curie(iri: str) -> str:
    """Compact an absolute ``xf`` / ``xmet`` IRI to a CURIE when possible."""

    return load().to_curie(iri)


def __getattr__(name: str) -> Any:
    """Lazy class exports from the Rust extension."""

    if name in ("BoundPattern", "RuleSet"):
        return getattr(load(), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")