"""Copy ``_forest`` without duplicating immutable rules / pattern dicts.

Approach (see LOG.md): an explicit :func:`forest_copy` walker — not
``ReactionRule.__deepcopy__`` / pickle hooks. Magic on rules would hide the
cost model and risk surprising picklers; the forest schema already names the
three layers (``immutable`` / ``cache`` / mutable), so the copy policy lives
next to that schema.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any, cast

from xenosite.refactor_poc.records import Forest, ImmutableForest, Structure

# Top-level keys with special copy policy (not walked as mutable tree).
_IMMUTABLE_KEY = "immutable"
_CACHE_KEY = "cache"


def freeze_immutable(data: Mapping[str, Any] | ImmutableForest) -> ImmutableForest:
    """Wrap ``data`` as ``MappingProxyType`` with frozen nested maps."""

    frozen: dict[str, Any] = {}
    for key, value in dict(data).items():
        if isinstance(value, Mapping) and not isinstance(value, MappingProxyType):
            frozen[key] = MappingProxyType(dict(value))
        elif isinstance(value, MappingProxyType):
            frozen[key] = value
        else:
            # Only immutable leaves belong here; callers must not pass mutables.
            frozen[key] = value
    return cast(ImmutableForest, MappingProxyType(frozen))


def shallow_immutable(imm: ImmutableForest | Mapping[str, Any]) -> ImmutableForest:
    """New proxy shell; nested values shared by identity."""

    return cast(ImmutableForest, MappingProxyType(dict(imm)))


def start_labels_of(forest: Forest | Mapping[str, Any]) -> Mapping[int, str] | None:
    imm = forest.get(_IMMUTABLE_KEY)  # type: ignore[arg-type]
    if imm is None:
        return None
    labels = imm.get("start_labels")
    if labels is None:
        return None
    return cast(Mapping[int, str], labels)


def set_start_labels(forest: Forest, labels: dict[int, str]) -> None:
    """Replace ``immutable`` with a frozen map holding ``start_labels``."""

    forest[_IMMUTABLE_KEY] = freeze_immutable({"start_labels": labels})


def _is_reaction_rule(obj: object) -> bool:
    """True for ``ReactionRule`` (and subclasses) without importing ``rules``."""

    for cls in type(obj).__mro__:
        if cls.__name__ == "ReactionRule" and "rules" in getattr(cls, "__module__", ""):
            return True
    return False


def _is_pattern_info(obj: object) -> bool:
    """``PatternInfo`` dicts stay shared (identity used in dedup lookups)."""

    return (
        isinstance(obj, dict)
        and "possibilities" in obj
        and "span" in obj
    )


def copy_mutable(obj: Any, memo: dict[int, Any] | None = None) -> Any:
    """Deep-copy Python containers; share rules, PatternInfo, and immutables."""

    if memo is None:
        memo = {}

    if obj is None or isinstance(obj, (bool, int, float, complex, str, bytes, type)):
        return obj
    if isinstance(obj, MappingProxyType):
        return obj
    if isinstance(obj, frozenset):
        return obj
    if _is_reaction_rule(obj) or _is_pattern_info(obj):
        return obj

    obj_id = id(obj)
    if obj_id in memo:
        return memo[obj_id]

    if isinstance(obj, tuple):
        # NamedTuple / plain tuple: rebuild so nested rules stay shared.
        if hasattr(obj, "_fields"):
            copied = type(obj)(*(_copy_item(x, memo) for x in obj))
        else:
            copied = tuple(_copy_item(x, memo) for x in obj)
        memo[obj_id] = copied
        return copied

    if isinstance(obj, list):
        out_list: list[Any] = []
        memo[obj_id] = out_list
        out_list.extend(_copy_item(x, memo) for x in obj)
        return out_list

    if isinstance(obj, dict):
        out_dict: dict[Any, Any] = {}
        memo[obj_id] = out_dict
        for key, value in obj.items():
            # PatternInfo under TraceAddition["pattern"] stays shared even if
            # the shape heuristic misses (empty / partial stubs).
            if key == "pattern" and isinstance(value, dict):
                out_dict[key] = value
            else:
                out_dict[_copy_item(key, memo)] = _copy_item(value, memo)
        return out_dict

    if isinstance(obj, set):
        out_set: set[Any] = set()
        memo[obj_id] = out_set
        out_set.update(_copy_item(x, memo) for x in obj)
        return out_set

    # Unknown objects (e.g. RDKit Mol if ever reached): share by identity.
    return obj


def _copy_item(obj: Any, memo: dict[int, Any]) -> Any:
    return copy_mutable(obj, memo)


def forest_copy(forest: Forest, *, same_structure: bool = False) -> Forest:
    """Copy a forest with the three-layer policy.

    - **Mutable** tree (``atom_trace``, ``is_terminal_product``, …): deep-copied
      via :func:`copy_mutable` (rules / PatternInfo shared by identity).
    - **``immutable``**: shallow-copied (new ``MappingProxyType`` shell; nested
      values shared).
    - **``cache``** (structure-dependent ephemeral): dropped when
      ``same_structure`` is false; kept by identity when true.
    """

    out: dict[str, Any] = {}
    for key, value in forest.items():
        if key == _IMMUTABLE_KEY or key == _CACHE_KEY:
            continue
        out[key] = copy_mutable(value)

    imm = forest.get(_IMMUTABLE_KEY)
    if imm is not None:
        out[_IMMUTABLE_KEY] = shallow_immutable(imm)

    if same_structure:
        cache = forest.get(_CACHE_KEY)
        if cache is not None:
            out[_CACHE_KEY] = cache

    return cast(Forest, out)


def empty_forest() -> Forest:
    """Fresh forest shell: empty mutable cache, no immutable yet."""

    return {_CACHE_KEY: cast(Structure, {})}
