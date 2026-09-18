"""Mark new public APIs whose surface may still change."""

from __future__ import annotations

import functools
import inspect
import warnings
from typing import Callable, Optional, TypeVar

F = TypeVar("F", bound=Callable)
C = TypeVar("C", bound=type)

_warned: set[str] = set()


class UnstableWarning(UserWarning):
    """New feature; public API may change without a deprecation cycle."""


def _emit(name: str, detail: Optional[str] = None) -> None:
    if name in _warned:
        return
    _warned.add(name)
    msg = (
        "%s is a new feature; its public API is unstable and may change "
        "without a deprecation cycle." % (name,)
    )
    if detail:
        msg = "%s %s" % (msg, detail)
    warnings.warn(msg, UnstableWarning, stacklevel=3)


def unstable(obj=None, *, name: Optional[str] = None, detail: Optional[str] = None):
    """Decorator for functions / methods / classes with an unstable public API.

    Emits :class:`UnstableWarning` once per decorated object per process.
    Filter with ``warnings.filterwarnings("ignore", category=UnstableWarning)``.
    """

    def decorate(target):
        label = name or getattr(target, "__qualname__", None) or getattr(
            target, "__name__", repr(target)
        )
        if inspect.isclass(target):
            return _decorate_class(target, label, detail)
        if not callable(target):
            raise TypeError("@unstable applies to callables and classes")
        return _decorate_callable(target, label, detail)

    if obj is None:
        return decorate
    return decorate(obj)


def _decorate_callable(fn: F, label: str, detail: Optional[str]) -> F:
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        _emit(label, detail)
        return fn(*args, **kwargs)

    wrapper.__unstable__ = True  # type: ignore[attr-defined]
    return wrapper  # type: ignore[return-value]


def _decorate_class(cls: C, label: str, detail: Optional[str]) -> C:
    init = cls.__init__

    @functools.wraps(init)
    def __init__(self, *args, **kwargs):
        _emit(label, detail)
        return init(self, *args, **kwargs)

    cls.__init__ = __init__  # type: ignore[method-assign]
    cls.__unstable__ = True  # type: ignore[attr-defined]
    return cls
