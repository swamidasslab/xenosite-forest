"""Rust chematic ``find_path`` door (native extension).

Separate from the live RDKit walk in :mod:`xenosite.forest.find_path`.
Requires the ``xenosite-forest-native`` wheel (``import xenosite_forest``).

Example::

    from xenosite.forest.find_path_rust import find_path_rust

    hits, counters = find_path_rust("CC", "CCO", max_paths=1)
    assert hits[0]["smiles"]
    assert counters["billed"] >= 0
"""

from __future__ import annotations

from typing import Any


def _native():
    try:
        import xenosite_forest as native
    except ImportError as e:  # pragma: no cover - optional native wheel
        raise ImportError(
            "find_path_rust requires the xenosite-forest-native extension "
            "(maturin develop -m crates/xenosite-forest/Cargo.toml "
            "--features python,extension-module)"
        ) from e
    return native


def find_path_rust(
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
    Does not produce RDKit mols or archive ``Deps`` — door API only.
    """
    return _native().find_path(
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


def native_available() -> bool:
    """True when ``xenosite_forest`` can be imported."""
    try:
        import xenosite_forest  # noqa: F401

        return True
    except ImportError:
        return False
