"""Load the maturin extension ``xenosite.forest._rust``."""

from __future__ import annotations

from typing import Any

_ext: Any | None = None


def load() -> Any:
    """Import ``xenosite.forest._rust`` or raise a clear ``ImportError``."""

    global _ext
    if _ext is not None:
        return _ext
    try:
        from . import _rust as ext  # type: ignore[attr-defined]
    except ImportError as e:  # pragma: no cover - missing wheel / editable build
        raise ImportError(
            "xenosite.forest requires the compiled Rust extension "
            "(install a platform wheel from PyPI, or: maturin develop "
            "--features python,extension-module)"
        ) from e
    from . import notebook

    notebook.install(ext)
    _ext = ext
    return ext


def available() -> bool:
    """True when ``xenosite.forest._rust`` can be imported."""

    try:
        from . import _rust  # type: ignore[attr-defined]  # noqa: F401

        return True
    except ImportError:
        return False
