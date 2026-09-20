"""Copy ``_forest`` without duplicating immutable rules / pattern dicts.

Approach (see LOG.md): an explicit :func:`forest_copy` walker — not
``ReactionRule.__deepcopy__`` / pickle hooks. Magic on rules would hide the
cost model and risk surprising picklers; the forest schema already names the
three layers (``immutable`` / ``cache`` / mutable), so the copy policy lives
next to that schema.
"""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from types import MappingProxyType
from typing import TypeVar, cast

from xenosite.refactor_poc.records import Forest, ImmutableForest, Structure

# Top-level keys with special copy policy (not walked as mutable tree).
_IMMUTABLE_KEY = "immutable"
_CACHE_KEY = "cache"

_T = TypeVar("_T")
_ImmutableLeaf = str | int | float | bool | None | tuple[object, ...] | frozenset[object]
_ImmutableValue = _ImmutableLeaf | Mapping[object, object]


def freeze_immutable(data: Mapping[str, object] | ImmutableForest) -> ImmutableForest:
    """Wrap ``data`` as ``MappingProxyType`` with frozen nested maps."""

    frozen: dict[str, _ImmutableValue] = {}
    for key, value in dict(data).items():
        if isinstance(value, Mapping) and not isinstance(value, MappingProxyType):
            frozen[key] = MappingProxyType(dict(cast(Mapping[object, object], value)))
        elif isinstance(value, MappingProxyType):
            frozen[key] = cast(Mapping[object, object], value)
        else:
            # Only immutable leaves belong here; callers must not pass mutables.
            frozen[key] = cast(_ImmutableLeaf, value)
    return cast(ImmutableForest, MappingProxyType(frozen))


def shallow_immutable(imm: ImmutableForest | Mapping[str, object]) -> ImmutableForest:
    """New proxy shell; nested values shared by identity."""

    return cast(ImmutableForest, MappingProxyType(dict(imm)))


def start_labels_of(forest: Forest | Mapping[str, object]) -> Mapping[int, str] | None:
    imm_obj = forest.get(_IMMUTABLE_KEY)  # type: ignore[arg-type]
    if imm_obj is None:
        return None
    imm = cast(ImmutableForest, imm_obj)
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

    return isinstance(obj, dict) and "possibilities" in obj and "span" in obj


def copy_mutable(obj: _T, memo: dict[int, object] | None = None) -> _T:
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
        return cast(_T, memo[obj_id])

    if isinstance(obj, tuple):
        # NamedTuple / plain tuple: rebuild so nested rules stay shared.
        items = tuple(_copy_item(x, memo) for x in obj)
        if hasattr(obj, "_fields"):
            # Dynamic NamedTuple rebuild; field types are not known statically.
            copied = cast(object, type(obj)(*items))  # type: ignore[arg-type]
        else:
            copied = items
        memo[obj_id] = copied
        return cast(_T, copied)

    if isinstance(obj, list):
        out_list: list[object] = []
        memo[obj_id] = out_list
        out_list.extend(_copy_item(x, memo) for x in obj)
        return cast(_T, out_list)

    if isinstance(obj, dict):
        out_dict: dict[object, object] = {}
        memo[obj_id] = out_dict
        for key, value in cast(MutableMapping[object, object], obj).items():
            # PatternInfo under TraceAddition["pattern"] stays shared even if
            # the shape heuristic misses (empty / partial stubs).
            if key == "pattern" and isinstance(value, dict):
                out_dict[key] = value
            else:
                out_dict[_copy_item(key, memo)] = _copy_item(value, memo)
        return cast(_T, out_dict)

    if isinstance(obj, set):
        out_set: set[object] = set()
        memo[obj_id] = out_set
        out_set.update(_copy_item(x, memo) for x in obj)
        return cast(_T, out_set)

    # Unknown objects (e.g. RDKit Mol if ever reached): share by identity.
    return obj


def _copy_item(obj: object, memo: dict[int, object]) -> object:
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

    out: dict[str, object] = {}
    for key, value in cast(Mapping[str, object], forest).items():
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
