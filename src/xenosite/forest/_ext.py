"""Load the maturin extension ``xenosite_forest`` (Rust product door)."""

from __future__ import annotations

from typing import Any


def load() -> Any:
    """Import ``xenosite_forest`` or raise a clear ``ImportError``."""

    try:
        import xenosite_forest as ext  # type: ignore[import-not-found]
    except ImportError as e:  # pragma: no cover - optional native wheel
        raise ImportError(
            "xenosite.forest requires the xenosite-forest-native extension "
            "(maturin develop -m crates/xenosite-forest/Cargo.toml "
            "--features python,extension-module)"
        ) from e
    return ext


def available() -> bool:
    """True when ``xenosite_forest`` can be imported."""

    try:
        import xenosite_forest  # type: ignore[import-not-found]  # noqa: F401

        return True
    except ImportError:
        return False
