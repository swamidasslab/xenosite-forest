"""Load the maturin extension ``xenosite.forest._rust``."""

from __future__ import annotations

from typing import Any


def load() -> Any:
    """Import ``xenosite.forest._rust`` or raise a clear ``ImportError``."""

    try:
        from . import _rust as ext  # type: ignore[attr-defined]
    except ImportError as e:  # pragma: no cover - missing wheel / editable build
        raise ImportError(
            "xenosite.forest requires the compiled Rust extension "
            "(install a platform wheel from PyPI, or: maturin develop "
            "--features python,extension-module)"
        ) from e
    return ext


def available() -> bool:
    """True when ``xenosite.forest._rust`` can be imported."""

    try:
        from . import _rust  # type: ignore[attr-defined]  # noqa: F401

        return True
    except ImportError:
        return False
