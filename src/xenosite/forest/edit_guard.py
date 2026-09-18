"""Guard RDKit molecules against structural edits unless explicitly opened.

``edit_guard`` installs process-wide wrappers on mutating ``Atom`` / ``Bond`` /
``RWMol`` methods and on ``Kekulize``. While a guard is active, those calls
raise :class:`MolClosedError` unless the molecule was opened with
:func:`edit_mol`. Property writes (``SetProp``, atom maps, property-cache
refresh) stay allowed.

``edit_mol`` unlocks one molecule and drops its resonance cache if a guarded
mutator actually ran (Kekulize counts). Guards are off in production until
:func:`edit_guard` is entered; tests enter one for the whole session.

.. code-block:: python

    from xenosite.forest import edit_guard, edit_mol

    with edit_guard():
        atom.SetProp("ok", "1")          # allowed
        bond.SetBondType(...)            # MolClosedError
        with edit_mol(mol):
            bond.SetBondType(...)        # allowed; resonance cleared if it ran
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from typing import Iterator

from rdkit.Chem import rdchem
from rdkit.Chem import Kekulize as _Kekulize
from rdkit.Chem.rdmolops import Kekulize as _OpsKekulize

_STACK: ContextVar[tuple] = ContextVar("xenosite_edit_guard_stack", default=())
_INSTALLED = False
_OPEN_PROP = "_forestEditOpen"
_TOKEN = 0


class MolClosedError(RuntimeError):
    """Structural edit attempted on a molecule that was not opened."""


class _Frame:
    __slots__ = ("kind", "token", "py_id", "dirty")

    def __init__(self, kind: str, token: str | None = None, py_id: int | None = None):
        self.kind = kind
        self.token = token
        self.py_id = py_id
        self.dirty = False


def _next_token() -> str:
    global _TOKEN
    _TOKEN += 1
    return str(_TOKEN)


def _mol_of(obj):
    if isinstance(obj, rdchem.Mol):
        return obj
    owning = getattr(obj, "GetOwningMol", None)
    if owning is None:
        return None
    try:
        return owning()
    except Exception:
        return None


def _tokens(mol) -> set[str]:
    """Open-tokens live on the C++ mol so GetOwningMol() wrappers see them.

    ``GetOwningMol`` returns a different Python object than the mol we opened;
    ``id()`` does not match. Props do.
    """
    if mol is None:
        return set()
    try:
        if not mol.HasProp(_OPEN_PROP):
            return set()
        return {t for t in mol.GetProp(_OPEN_PROP).split(",") if t}
    except Exception:
        return set()


def _set_tokens(mol, tokens: set[str]) -> None:
    if not hasattr(mol, "SetProp"):
        return
    if tokens:
        mol.SetProp(_OPEN_PROP, ",".join(sorted(tokens)))
    elif hasattr(mol, "HasProp") and mol.HasProp(_OPEN_PROP):
        mol.ClearProp(_OPEN_PROP)


def _guard_active(stack) -> bool:
    return any(frame.kind == "guard" for frame in stack)


def _is_open(stack, obj) -> bool:
    if any(frame.kind == "edit" and frame.py_id == id(obj) for frame in stack):
        return True
    mol = obj if isinstance(obj, rdchem.Mol) else _mol_of(obj)
    tokens = _tokens(mol)
    if not tokens:
        return False
    return any(
        frame.kind == "edit" and frame.token and frame.token in tokens for frame in stack
    )


def _mark_dirty(stack, obj) -> None:
    mol = obj if isinstance(obj, rdchem.Mol) else _mol_of(obj)
    tokens = _tokens(mol)
    for frame in stack:
        if frame.kind != "edit":
            continue
        if frame.py_id == id(obj) or (frame.token and frame.token in tokens):
            frame.dirty = True


def _refuse(obj, name: str) -> None:
    stack = _STACK.get()
    if not _guard_active(stack):
        return
    if _is_open(stack, obj):
        return
    owner = type(obj).__name__
    raise MolClosedError(
        "%s.%s refused: mol is closed. Open it with edit_mol()." % (owner, name)
    )


def _after(obj) -> None:
    _mark_dirty(_STACK.get(), obj)


def _wrap_method(orig, name: str):
    def wrapper(self, *args, **kwargs):
        _refuse(self, name)
        try:
            return orig(self, *args, **kwargs)
        finally:
            _after(self)

    wrapper._edit_guard_wrapped = True  # type: ignore[attr-defined]
    wrapper.__name__ = name
    wrapper.__wrapped__ = orig  # type: ignore[attr-defined]
    return wrapper


def _wrap_kekulize(orig):
    def wrapper(mol, *args, **kwargs):
        _refuse(mol, "Kekulize")
        _mark_dirty(_STACK.get(), mol)
        return orig(mol, *args, **kwargs)

    wrapper._edit_guard_wrapped = True  # type: ignore[attr-defined]
    wrapper.__name__ = "Kekulize"
    wrapper.__wrapped__ = orig  # type: ignore[attr-defined]
    return wrapper


_ATOM_EDITS = (
    "SetAtomicNum",
    "SetChiralTag",
    "SetFormalCharge",
    "SetHybridization",
    "SetIsAromatic",
    "SetIsotope",
    "SetMonomerInfo",
    "SetNoImplicit",
    "SetNumExplicitHs",
    "SetNumRadicalElectrons",
    "SetPDBResidueInfo",
)
_BOND_EDITS = (
    "SetBondDir",
    "SetBondType",
    "SetIsAromatic",
    "SetIsConjugated",
    "SetStereo",
    "SetStereoAtoms",
)
_RWMOL_EDITS = (
    "AddAtom",
    "AddBond",
    "InsertMol",
    "RemoveAtom",
    "RemoveBond",
    "ReplaceAtom",
    "ReplaceBond",
    "SetStereoGroups",
)


def _install() -> None:
    """Install wrappers once. Cheap no-op when the guard stack is empty."""
    global _INSTALLED
    if _INSTALLED:
        return
    pairs = (
        (rdchem.Atom, _ATOM_EDITS),
        (rdchem.Bond, _BOND_EDITS),
        (rdchem.RWMol, _RWMOL_EDITS),
        (
            rdchem.EditableMol,
            ("AddAtom", "AddBond", "RemoveAtom", "RemoveBond"),
        ),
    )
    for cls, names in pairs:
        for name in names:
            current = getattr(cls, name, None)
            if current is None or getattr(current, "_edit_guard_wrapped", False):
                continue
            setattr(cls, name, _wrap_method(current, name))

    wrapped = _wrap_kekulize(_OpsKekulize)
    import sys

    import rdkit.Chem as Chem
    import rdkit.Chem.rdmolops as rdmolops

    rdmolops.Kekulize = wrapped
    Chem.Kekulize = wrapped
    rdkit_ids = {id(_Kekulize), id(_OpsKekulize), id(Chem.Kekulize)}
    # `from rdkit... import Kekulize` binds a separate global per module.
    for mod in list(sys.modules.values()):
        if mod is None:
            continue
        name = getattr(mod, "__name__", "")
        if not (
            name.startswith("xenosite.")
            or name.startswith("rdkit.")
            or name == "rdkit"
        ):
            continue
        current = getattr(mod, "Kekulize", None)
        if current is None or getattr(current, "_edit_guard_wrapped", False):
            continue
        if id(current) in rdkit_ids or getattr(current, "__module__", "").startswith(
            "rdkit"
        ):
            try:
                mod.Kekulize = wrapped
            except Exception:
                continue
    _INSTALLED = True


@contextmanager
def edit_guard() -> Iterator[None]:
    """Reject structural edits on mols that are not inside :func:`edit_mol`.

    Nestable. Property changes are not structural and are not blocked.
    """
    _install()
    frame = _Frame("guard")
    token = _STACK.set(_STACK.get() + (frame,))
    try:
        yield
    finally:
        _STACK.reset(token)


@contextmanager
def edit_mol(mol, *, clear: str = "if_mutated") -> Iterator:
    """Open ``mol`` for structural edits.

    On exit, drop ``mol._forest['resonance']`` (and a cached standardize view)
    when a guarded mutator ran. ``clear="always"`` drops it even if every edit
    went through an unwrapped C++ path. ``clear="never"`` leaves the cache.
    """
    if clear not in ("if_mutated", "always", "never"):
        raise ValueError("clear must be 'if_mutated', 'always', or 'never'")
    _install()
    token = _next_token() if hasattr(mol, "SetProp") else None
    if token is not None:
        _set_tokens(mol, _tokens(mol) | {token})
    frame = _Frame("edit", token, id(mol))
    stack_token = _STACK.set(_STACK.get() + (frame,))
    try:
        yield mol
    finally:
        _STACK.reset(stack_token)
        if token is not None:
            _set_tokens(mol, _tokens(mol) - {token})
        if clear == "always" or (clear == "if_mutated" and frame.dirty):
            from .base import _clear_resonance

            _clear_resonance(mol)


def guard_installed() -> bool:
    return _INSTALLED
