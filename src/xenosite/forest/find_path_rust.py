"""Rust chematic doors (native extension).

Separate from the live RDKit walk in :mod:`xenosite.forest.find_path`.
Requires the ``xenosite-forest-native`` wheel (``import xenosite_forest``).

Exposes rulesets (``RuleSet.leaf`` / ``phase_one`` / …), ``find_path``, and
``enumerate`` / ``bfs`` / ``dfs``. Callers override catalogs via ``ruleset=``.
"""

from __future__ import annotations

from typing import Any


def _native() -> Any:
    try:
        import xenosite_forest as native  # type: ignore[import-not-found]
    except ImportError as e:  # pragma: no cover - optional native wheel
        raise ImportError(
            "find_path_rust requires the xenosite-forest-native extension "
            "(maturin develop -m crates/xenosite-forest/Cargo.toml "
            "--features python,extension-module)"
        ) from e
    return native


def native_available() -> bool:
    """True when ``xenosite_forest`` can be imported."""
    try:
        import xenosite_forest  # type: ignore[import-not-found]  # noqa: F401

        return True
    except ImportError:
        return False


def find_path_rust(
    reactant: str,
    target: str,
    *,
    ruleset: Any | None = None,
    max_paths: int = 1,
    max_nodes: int = 800,
    use_atom_diff: bool = True,
    lazy_closer: bool = False,
    diversity: bool = False,
    drop_skeleton_twins: bool = True,
    score: str = "log-neg-pc",
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run chematic ``find_path``; return ``(hits, counters)``.

    Default catalog is PhaseOne; pass ``ruleset`` (native ``RuleSet``) to
    override. Each hit is ``{"smiles": str, "steps": [{"rule": str, "site": list}]}``.
    """
    return _native().find_path(
        reactant,
        target,
        ruleset=ruleset,
        max_paths=max_paths,
        max_nodes=max_nodes,
        use_atom_diff=use_atom_diff,
        lazy_closer=lazy_closer,
        diversity=diversity,
        drop_skeleton_twins=drop_skeleton_twins,
        score=score,
    )


def bfs_rust(
    reactant: str,
    max_depth: int,
    *,
    ruleset: Any | None = None,
) -> list[dict[str, Any]]:
    """Native breadth-first enumerate (default: ``default_ruleset``)."""

    return _native().bfs(reactant, max_depth, ruleset=ruleset)


def dfs_rust(
    reactant: str,
    max_depth: int,
    *,
    ruleset: Any | None = None,
) -> list[dict[str, Any]]:
    """Native depth-first enumerate (default: ``default_ruleset``)."""

    return _native().dfs(reactant, max_depth, ruleset=ruleset)


def enumerate_rust(
    reactant: str,
    *,
    ruleset: Any | None = None,
    max_depth: int = 1,
    order: str = "bfs",
    max_nodes: int = 0,
    unique_csmi: bool = True,
) -> list[dict[str, Any]]:
    """Native metabolite enumerate with depth / order / dedup knobs."""

    return _native().enumerate(
        reactant,
        ruleset=ruleset,
        max_depth=max_depth,
        order=order,
        max_nodes=max_nodes,
        unique_csmi=unique_csmi,
    )
