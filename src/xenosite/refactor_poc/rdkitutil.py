"""Molecule questions for the proof of concept.

RDKit itself is imported in :mod:`xenosite.refactor_poc.rdkit_api`. This
module calls those names.

Answers about a molecule are cached on ``mol.xf.forest["cache"]``
(structure-dependent ephemeral; formerly ``"structure"``).
Accessing ``mol.xf`` mints an ephemeral :class:`Xf` (strong parent ref).
Forest attaches lazily via :func:`_require_forest` when a method needs it.
``mol.xf.has_forest`` reports wipe/absence without installing. Prefer
``mol.xf.*`` for queries and finishing; no prior ``ensure_forest`` ritual.

Typing brands (stubs; runtime is RDKit ``Mol``):

- :class:`ForestMol` — ``_forest`` present. Wipe APIs return ``NoForestMol`` /
  ``Mol``, not ``ForestMol``. Re-enter via ``mol.xf.forestmol``.
- :class:`TracingMol` — initialized atom-trace (extends ForestMol).
- :class:`NoForestMol` — no forest (constructors, reaction pieces, wipe).

Atom-trace plumbing lives under ``mol.xf.tracing``. :func:`ensure_forest` is
a thin alias of :func:`_require_forest`.

The key is one string. Process-wide data, such as parsed SMARTS reactions,
stays a module dict. Do not cache a result on a molecule this function
then edits. Edits belong on a copy from :func:`copy_mol` or :func:`rw_copy`.

:func:`ensure_kekule_parents`, :func:`parent_for_bond`, and
:func:`parents_for_ends` take a dict and do not read or write ``mol._forest``.
The rule that calls them stores that dict.

``Formula`` is a dict. Callers use keys. :class:`~xenosite.refactor_poc.records.McsResult`
and :class:`~xenosite.refactor_poc.records.FragmentSplit` are NamedTuples.
Callers use attributes.

``records.Structure`` does not list the MCS cache keys. This module writes
them anyway, because the match list has to live on the reactant:

- ``mcs_matches``: ``dict[str, McsResult]``, reactant-side embeddings
- ``mcs_targets``: ``dict[str, McsResult]``, target-side embeddings

Both are keyed by the target's canonical SMILES. ``records.py`` is not
edited here.
"""

from __future__ import annotations

import ast
from collections import defaultdict, deque
from collections.abc import Iterable, Sequence
from typing import Any, Literal, TypeGuard, TypeVar, cast, overload

from xenosite.refactor_poc.forest_copy import (
    copy_mutable,
    empty_forest,
    forest_copy,
    set_start_labels,
    shallow_immutable,
    start_labels_of,
)
from xenosite.refactor_poc.rdkit_api import (
    KEKULE_ALL,
    Atom,
    AtomCompare,
    Bond,
    BondCompare,
    BondType,
    CanonicalRankAtoms,
    ChemicalReaction,
    DisableLog,
    FindMCS,
    ForestMol,
    GetMolFrags,
    Mol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmiles,
    NoForestMol,
    ReactionFromSmarts,
    RenumberAtoms,
    ResonanceMolSupplier,
    RWMol,
    SanitizeFlags,
    SanitizeMol,
    TracingMol,
)
from xenosite.refactor_poc.records import (
    AtomPairOrbitSignature,
    BondAtomOrbitSignature,
    BondPairOrbitSignature,
    EditCounters,
    EndParents,
    Forest,
    Formula,
    FragmentSplit,
    InitializedAtomTrace,
    KekuleParents,
    McsResult,
    SiteInfo,
    SitePairOrbitTables,
    Structure,
)

DisableLog("rdApp.*")

_REACTION_CACHE: dict[str, ChemicalReaction] = {}
_MolT = TypeVar("_MolT", bound=Mol)

# Keys install_forest writes. is_tracing is true only when all of them are present.
_TRACE_KEYS = (
    "records",
    "deletes",
    "transforms",
    "additions",
    "formula",
    "delta_formula",
    "depth",
    "last_tag",
    "next_transform",
)


def _read_forest(mol: Mol) -> Forest | None:
    try:
        forest = mol._forest
    except AttributeError:
        return None
    return forest


class XfTracing:
    """Atom-trace queries and plumbing for a parent :class:`~xenosite.refactor_poc.rdkit_api.Mol`.

    Nested under ``mol.xf.tracing``. Holds a strong reference to the same
    parent as :class:`Xf`. Trace state lives on the forest.

    Read-only queries (no underscore): ``active``, ``depth``, ``atom_origin``,
    ``atom_indices``, ``atom_depths``, ``atom_root``, ``atom_added_by``,
    ``removed_roots``. These answer questions about the installed atom trace
    without exposing ``forestLabel`` tags or transform ids.

    ``atom_origin`` is the earliest index in the record (an added atom's
    birth index). ``atom_root`` is the depth-0 index, or None when the atom
    was created later. ``removed_roots`` are depth-0 indexes that left.

    Plumbing mutators (underscore): ``_ensure``, ``_stamp``, ``_install``,
    ``_trace``. Downstream code should prefer ``mol.xf.of_products`` over
    calling these directly.
    """

    __slots__ = ("_mol",)

    def __init__(self, mol: Mol) -> None:
        self._mol = mol

    @property
    def mol(self) -> Mol:
        return self._mol

    @property
    def active(self) -> bool:
        """True when the parent has an initialized atom trace."""

        return is_tracing(self.mol)

    @property
    def depth(self) -> int | None:
        """Transform depth from the root, or None when untraced.

        Trace-derived (from ``atom_trace["depth"]``), so it lives here rather
        than on :class:`Xf` beside structure queries like ``csmi`` / rings.
        """

        if not is_tracing(self.mol):
            return None
        return int(self.mol._forest["atom_trace"]["depth"])

    def atom_indices(self, idx: int) -> tuple[int, ...] | None:
        """Index frames for atom ``idx`` from origin to now.

        None when the parent is untraced or the atom has no record. Does not
        expose internal label tags.
        """

        record = self._atom_record(idx)
        if record is None:
            return None
        idxs = record.get("idx")
        if not idxs:
            return None
        return tuple(int(i) for i in idxs)

    def atom_depths(self, idx: int) -> tuple[int, ...] | None:
        """Depth frames parallel to :meth:`atom_indices`.

        None when the parent is untraced or the atom has no record.
        """

        record = self._atom_record(idx)
        if record is None:
            return None
        depths = record.get("depth")
        if not depths:
            return None
        return tuple(int(d) for d in depths)

    def atom_root(self, idx: int) -> int | None:
        """Depth-0 index of atom ``idx``, or None if it was created later.

        Distinct from :meth:`atom_origin`, which is the earliest recorded
        index even for an atom that did not exist on the reactant.
        """

        depths = self.atom_depths(idx)
        indices = self.atom_indices(idx)
        if depths is None or indices is None or 0 not in depths:
            return None
        return indices[depths.index(0)]

    def removed_roots(self) -> frozenset[int]:
        """Depth-0 indexes that are no longer in the molecule."""

        mol = self.mol
        if not is_tracing(mol):
            return frozenset()
        out: set[int] = set()
        for record in mol._forest["atom_trace"]["deletes"].values():
            depths = record.get("depth") or []
            idxs = record.get("idx") or []
            if not depths or not idxs or 0 not in depths:
                continue
            out.add(int(idxs[list(depths).index(0)]))
        return frozenset(out)

    def atom_origin(self, idx: int) -> int | None:
        """Earliest recorded index for atom ``idx``, or None if unknown."""

        history = self.atom_indices(idx)
        if history is None:
            return None
        return history[0]

    def atom_added_by(self, idx: int) -> tuple[str, frozenset[int]] | None:
        """Rule name and site for an atom created by a prior transform.

        None when the parent is untraced, the atom has no record, or the atom
        was not added by a transform. Does not expose transform ids or
        ``forestLabel`` tags — only the rule name and site the addition stored.
        """

        record = self._atom_record(idx)
        if record is None:
            return None
        added = record.get("added_by")
        if not added:
            return None
        mol = self.mol
        if not is_tracing(mol):
            return None
        trace = mol._forest["atom_trace"]
        if isinstance(added, str):
            detail = trace["additions"].get(added)
            if detail is None:
                return None
        else:
            detail = added
        name = detail.get("name")
        if name is None:
            rule = detail.get("rule")
            name = rule if isinstance(rule, str) else getattr(rule, "name", None)
        if name is None:
            return None
        from xenosite.refactor_poc.records import _flat_ints

        site = detail.get("site") or ()
        return (str(name), frozenset(_flat_ints(site)))

    def _ensure(self) -> TracingMol:
        """Initialize the atom trace on the parent if missing. Does not reset depth."""

        held = _require_forest(self.mol)
        forest = held._forest
        if "cache" not in forest:
            forest["cache"] = {}
        if "atom_trace" not in forest:
            trace: InitializedAtomTrace = {
                "records": {},
                "deletes": {},
                "transforms": [],
                "additions": {},
                "formula": molecule_formula(held),
                "delta_formula": {},
                "depth": 0,
                "last_tag": 0,
                "next_transform": 1,
            }
            forest["atom_trace"] = trace
            for atom in held.GetAtoms():
                if atom.GetAtomicNum() != 1:
                    i = atom.GetIdx()
                    trace["records"][str(i)] = {
                        "idx": [i],
                        "depth": [0],
                    }
                    trace["last_tag"] = i
        assert is_tracing(held)
        _write_forest_labels(held)
        return held

    def _stamp(self) -> TracingMol:
        """Ensure tracing and write forest labels onto atoms. Same object."""

        return self._ensure()

    def _install(self) -> TracingMol:
        """Alias of :meth:`_stamp` (search-entry wording)."""

        return self._ensure()

    def _trace(
        self, reactant: Mol, info: SiteInfo, executed: Any | None = None
    ) -> InitializedAtomTrace:
        """Record one transform from ``reactant`` onto this product mol.

        Implementation lives in :mod:`xenosite.refactor_poc.rules`. Prefer
        ``reactant.xf.of_products(...)`` over calling this directly.

        After the rule-side apply, this facade copies the reactant's
        ``immutable`` (``start_labels``) onto the product forest and restamps
        CX ``atomLabel`` so finishing stays transparent to rules / search.
        """

        from xenosite.refactor_poc.rules import _apply_forest_trace

        result = _apply_forest_trace(reactant, self.mol, info, executed=executed)
        parent = _read_forest(reactant)
        held = _read_forest(self.mol)
        if parent is not None and held is not None:
            parent_imm = parent.get("immutable")
            if parent_imm is not None:
                held["immutable"] = shallow_immutable(parent_imm)
        _restamp_start_labels(self.mol)
        return result

    def _atom_record(self, idx: int):
        """Resolve the trace record for ``idx`` via forestLabel."""

        mol = self.mol
        if not is_tracing(mol):
            return None
        atom = mol.GetAtomWithIdx(idx)
        if not atom.HasProp("forestLabel"):
            return None
        return mol._forest["atom_trace"]["records"].get(atom.GetProp("forestLabel"))


class Xf:
    """Ephemeral facade for one molecule: minted on each ``mol.xf`` read.

    Design:

    - ``Mol.xf`` is a read-only property (monkey-patched onto RDKit ``Mol``).
      Each access constructs a new :class:`Xf`; nothing is stored on the mol,
      so facades cannot be copied between molecules.
    - :class:`Xf` holds a **strong** reference to its parent. Temporary chains
      like ``MolFromSmiles(...).xf.csmi`` are safe.
    - ``has_forest`` reports whether ``_forest`` exists **without** installing.
      ``forestmol`` ensures forest and returns the parent as :class:`ForestMol`
      (typed bridge for pyright). Wipe APIs must return ``NoForestMol`` / ``Mol``.
    - :class:`TracingMol` is forest + initialized atom-trace. Prefer it for
      of_products returns, filters, and find_path walks.

    Public surface: ``has_forest``, ``forestmol``, ``csmi``, ``forest``, ``is_terminal``,
    ``clear_structure``, ring / conjugate / ``topol_equiv`` / ``formula`` /
    ``sanitize`` / ``smarts_matches``, pair-orbit
    (``atom_pair_orbit_key`` / ``bond_pair_orbit_key`` / ``bond_atom_orbit_key`` /
    ``site_pair_orbits`` / ``pair_orbit_backend``), and ``of_products``. Tracing nest:
    ``active`` / ``depth`` / ``atom_origin`` / ``atom_indices`` /
    ``atom_depths`` / ``atom_root`` / ``atom_added_by`` / ``removed_roots``,
    plus underscored ``_stamp`` / ``_ensure`` / ``_install`` / ``_trace``.
    """

    __slots__ = ("_mol",)

    def __init__(self, mol: Mol) -> None:
        # Strong ref only. Forest attaches lazily via _require_forest when needed
        # so has_forest can report a wipe without installing.
        self._mol = mol

    @property
    def mol(self) -> Mol:
        return self._mol

    @property
    def has_forest(self) -> bool:
        """True when ``_forest`` is already present. Does not install one."""

        return is_forest(self.mol)

    @property
    def forestmol(self) -> ForestMol:
        """Parent as :class:`ForestMol` after lazy attach.

        Typed bridge: plain ``Mol`` → ``ForestMol`` for call sites that need
        the brand. Same object as ``mol``; installs via :func:`_require_forest`.
        """

        return _require_forest(self.mol)

    @property
    def forest(self) -> Forest:
        """The forest on the parent; installs an empty one if missing."""

        return _require_forest(self.mol)._forest

    @property
    def tracing(self) -> XfTracing:
        """Atom-trace facade for this parent (strong ref; minted each access)."""

        return XfTracing(_require_forest(self.mol))

    @property
    def csmi(self) -> str:
        """Canonical SMILES, cached on ``structure["csmi"]`` on first read."""

        mol = _require_forest(self.mol)
        forest = mol._forest
        if "cache" not in forest:
            forest["cache"] = {}
        structure = forest["cache"]
        csmi = structure.get("csmi")
        if not csmi:
            csmi = MolToSmiles(mol, isomericSmiles=False)
            structure["csmi"] = csmi
        return csmi

    @property
    def is_terminal(self) -> bool:
        """True when this mol must not be expanded further (forest flag)."""

        return bool(self.forest.get("is_terminal_product"))

    def _mark_terminal(self, value: bool = True) -> None:
        """Set or clear the forest-level terminal flag. Plumbing for finishers."""

        self.forest["is_terminal_product"] = value

    def clear_structure(self) -> None:
        """Drop cached structure answers. Labels and the rest of the forest stay."""

        _require_forest(self.mol)._forest["cache"] = {}

    @property
    def rings(self) -> dict[int, tuple[tuple[int, ...], ...]]:
        """Per-atom ring membership (cached on structure)."""

        return _ring_membership(self.mol)

    @property
    def conjugated_systems(self) -> tuple[frozenset[int], ...]:
        """Conjugated atom sets (cached; loads resonance on first read)."""

        return _conjugated_systems(self.mol)

    @property
    def aromatic_systems(self) -> tuple[frozenset[int], ...]:
        """Aromatic connected components (cached on structure)."""

        return _aromatic_systems(self.mol)

    @property
    def topol_equiv(self) -> dict[int, int]:
        """Map each atom index to its topological equivalence class."""

        return _topol_equiv(self.mol)

    @property
    def formula(self) -> Formula:
        """Heavy-atom counts, total hydrogens, and formal charge (cached)."""

        return molecule_formula(self.mol)

    @property
    def pair_orbit_backend(self) -> Literal["nauty", "smiles", "none"]:
        """Resolved pair-orbit backend (override / env / auto). Process-wide."""

        from xenosite.refactor_poc.graph_isomorphism import get_pair_orbit_backend

        return get_pair_orbit_backend()

    @classmethod
    def set_pair_orbit_backend(
        cls, backend: Literal["nauty", "smiles", "none"] | None
    ) -> None:
        """Set process-wide pair-orbit backend, or ``None`` to clear override."""

        from xenosite.refactor_poc.graph_isomorphism import set_pair_orbit_backend

        set_pair_orbit_backend(backend)

    def atom_pair_orbit_key(
        self, site: frozenset[int]
    ) -> AtomPairOrbitSignature | None:
        """Unique-edit atom-pair signature; caches tables on forest structure."""

        from xenosite.refactor_poc.graph_isomorphism import atom_pair_orbit_key

        return atom_pair_orbit_key(_require_forest(self.mol), site)

    def bond_pair_orbit_key(
        self, bonds: frozenset[int]
    ) -> BondPairOrbitSignature | None:
        """Unique-edit bond-pair signature; caches tables on forest structure."""

        from xenosite.refactor_poc.graph_isomorphism import bond_pair_orbit_key

        return bond_pair_orbit_key(_require_forest(self.mol), bonds)

    def bond_atom_orbit_key(
        self, bond_idx: int, atom_idx: int
    ) -> BondAtomOrbitSignature:
        """Unique-edit bond–atom signature; caches tables on forest structure."""

        from xenosite.refactor_poc.graph_isomorphism import bond_atom_orbit_key

        return bond_atom_orbit_key(_require_forest(self.mol), bond_idx, atom_idx)

    def site_pair_orbits(
        self, backend: Literal["nauty", "smiles", "none"] | None = None
    ) -> SitePairOrbitTables | None:
        """Nested pair-orbit tables on ``structure``, or ``None`` for ``none``.

        Default backend is the resolved process-wide setting. Pass
        ``backend="smiles"`` to opt into RDKit isotope tables; ``"nauty"`` for
        pynauty. Forest keys: ``site_pair_orbits_nauty`` /
        ``site_pair_orbits_smiles``.
        """

        from xenosite.refactor_poc.graph_isomorphism import (
            ensure_site_pair_orbit_tables,
        )

        return ensure_site_pair_orbit_tables(_require_forest(self.mol), backend)

    def sanitize(self) -> int:
        """Sanitize a copy and cache the RDKit status code. Does not edit parent."""

        return sanitize_mol(self.mol)

    def smarts_matches(self, smarts: str) -> tuple[dict[int, int], ...]:
        """Cached substructure matches for ``smarts``, keyed by atom map."""

        return _smarts_matches(self.mol, smarts)

    def of_products(
        self,
        product_or_product_list: Mol | Sequence[Mol],
        site_info: SiteInfo,
        executed: Any | None = None,
    ) -> list[TracingMol]:
        """Stamp and trace products of this reactant; return the finished list.

        Readable finishing path: ``reactant.xf.of_products(products, site_info)``.
        Accepts one product or a sequence. Each product is taken through
        :func:`copy_mol` (keeps lifted ``forestLabel`` props and any forest;
        bare ``Chem.Mol`` drops ``_forest`` only), then
        ``product.xf.tracing._trace`` from this reactant, then
        :meth:`clear_structure`. When the executed rule (or ``site_info``'s
        rule) has ``is_terminal_rule``, each finished product is marked
        terminal. The reactant is stamped first via ``tracing._stamp``. Does
        not sanitize or split. Inputs are not edited.
        """

        reactant = self.tracing._stamp()
        if isinstance(product_or_product_list, Mol):
            products: Sequence[Mol] = (product_or_product_list,)
        else:
            products = product_or_product_list
        rule = executed if executed is not None else site_info.get("rule")
        mark_terminal = bool(getattr(rule, "is_terminal_rule", False))
        finished: list[TracingMol] = []
        for product in products:
            held = copy_mol(product)
            held.xf.tracing._trace(reactant, site_info, executed=executed)
            if mark_terminal:
                held.xf._mark_terminal(True)
            held.xf.clear_structure()
            assert is_tracing(held)
            finished.append(held)
        return finished


def _capture_start_labels(mol: ForestMol) -> None:
    """Snapshot CX ``atomLabel`` once when the forest is first established.

    An empty map is stored so a later clear cannot be mistaken for "not
    yet captured" and re-read from a stripped work copy.
    """

    forest = mol._forest
    if start_labels_of(forest) is not None:
        return
    labels: dict[int, str] = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("atomLabel"):
            labels[int(atom.GetIdx())] = atom.GetProp("atomLabel")
    set_start_labels(forest, labels)


def _restamp_start_labels(mol: Mol) -> None:
    """Rewrite captured start ``atomLabel`` props onto current start atoms."""

    forest = _read_forest(mol)
    if forest is None:
        return
    start_labels = start_labels_of(forest)
    if not start_labels:
        return
    if is_tracing(mol):
        for record in mol._forest["atom_trace"]["records"].values():
            depths = record.get("depth") or []
            idxs = record.get("idx") or []
            if not depths or not idxs or 0 not in depths:
                continue
            root = int(idxs[list(depths).index(0)])
            label = start_labels.get(root)
            if label is None:
                continue
            mol.GetAtomWithIdx(int(idxs[-1])).SetProp("atomLabel", label)
        return
    for root, label in start_labels.items():
        if root < 0 or root >= mol.GetNumAtoms():
            continue
        mol.GetAtomWithIdx(int(root)).SetProp("atomLabel", label)


def _write_forest_labels(mol: TracingMol) -> None:
    _capture_start_labels(mol)
    for tag, record in mol._forest["atom_trace"]["records"].items():
        idxs = record.get("idx")
        if not idxs:
            raise KeyError("idx")
        mol.GetAtomWithIdx(idxs[-1]).SetProp("forestLabel", tag)
    _restamp_start_labels(mol)


def _require_forest(mol: Mol) -> ForestMol:
    """Install an empty forest when missing. Same object as :class:`ForestMol`.

    The only place forest presence is coerced. ``mol.xf.forestmol`` and
    :func:`ensure_forest` call this. Captures start-atom ``atomLabel`` once.
    """

    if is_forest(mol):
        return mol
    mol._forest = empty_forest()
    assert is_forest(mol)
    _capture_start_labels(mol)
    return mol


def _mol_xf_get(self: Mol) -> Xf:
    """Mint a fresh :class:`Xf` for ``self``. Forest attaches lazily on need."""

    return Xf(self)


def _mol_xf_set(self: Mol, value: object) -> None:
    raise AttributeError("xf is read-only; each access mints a new facade")


def _mol_xf_del(self: Mol) -> None:
    raise AttributeError("xf is read-only")


# Read-only property on RDKit Mol: never stored, never copied between mols.
Mol.xf = property(_mol_xf_get, _mol_xf_set, _mol_xf_del)  # type: ignore[misc]


def is_forest(mol: Mol) -> TypeGuard[ForestMol]:
    """True when ``_forest`` is present. Does not install one."""

    return _read_forest(mol) is not None


@overload
def is_tracing(mol: TracingMol) -> TypeGuard[TracingMol]: ...
@overload
def is_tracing(mol: Mol) -> TypeGuard[TracingMol]: ...
def is_tracing(mol: Mol) -> bool:
    """True when ``atom_trace`` exists and has the keys the trace writer fills.

    Does not install a forest or a trace.
    """

    forest = _read_forest(mol)
    if forest is None:
        return False
    trace = forest.get("atom_trace")
    if trace is None:
        return False
    return all(key in trace for key in _TRACE_KEYS)


@overload
def ensure_forest(mol: TracingMol) -> TracingMol: ...
@overload
def ensure_forest(mol: ForestMol) -> ForestMol: ...
@overload
def ensure_forest(mol: Mol) -> ForestMol: ...
def ensure_forest(mol: Mol) -> ForestMol:
    """Explicit coerce to :class:`ForestMol`. Prefer ``mol.xf.forestmol``.

    Thin alias of :func:`_require_forest`. Does not copy.
    """

    return _require_forest(mol)


def wipe_forest(mol: Mol) -> NoForestMol:
    """Remove ``_forest`` if present. Same object, typed as :class:`NoForestMol`.

    Honest wipe: callers must not treat the result as ``ForestMol`` until
    ``mol.xf.forestmol`` (or ``ensure_forest``) re-attaches.
    """

    if getattr(mol, "_forest", None) is not None:
        mol._forest = None  # type: ignore[assignment]
    return mol  # type: ignore[return-value]


def _place_forest(mol: Mol, forest: Forest) -> ForestMol:
    """Write ``forest`` onto ``mol`` and return that same object as ForestMol.

    Does not store ``xf`` (minted on access). The molecule is not copied.
    """

    mol._forest = forest
    assert is_forest(mol)
    return mol


def _structure(mol: Mol) -> Structure:
    """Return ``forest["cache"]``, creating an empty dict if needed.

    Name kept for call-site stability; the forest key is ``cache``.
    """

    forest = _require_forest(mol)._forest
    if "cache" not in forest:
        forest["cache"] = {}
    return forest["cache"]


def sanitize_mol(mol: Mol) -> int:
    """Sanitize a copy and cache the RDKit status code. Does not edit ``mol``."""

    structure = _structure(mol)
    if "sanitized" in structure:
        return structure["sanitized"]

    sanitized = int(SanitizeMol(Mol(mol), catchErrors=True))
    structure["sanitized"] = sanitized
    return sanitized


def sanitize_catch(mol: Mol) -> int:
    """Sanitize ``mol`` in place. Not cached: this call edits the molecule."""

    return int(SanitizeMol(mol, catchErrors=True))


def cip_ids(mol: Mol, *, include_stereo: bool = False) -> tuple[int, ...]:
    """Per-atom CIP / topological ranks. Not uniquified (``breakTies=False``).

    Same signal as ``topol_equiv`` when ``include_stereo`` is false. With
    stereo, ``includeChirality=True``. Index ``i`` is the rank of atom ``i``.
    """

    structure = _structure(mol)
    cache_key = "cip_ids_stereo" if include_stereo else "cip_ids"
    cached = structure.get(cache_key)
    if cached is not None:
        return cast(tuple[int, ...], cached)

    sanitized = Mol(mol)
    sanitize_mol(sanitized)
    ranks = tuple(
        int(r)
        for r in CanonicalRankAtoms(
            sanitized,
            includeChirality=include_stereo,
            breakTies=False,
        )
    )
    if include_stereo:
        structure["cip_ids_stereo"] = ranks
    else:
        structure["cip_ids"] = ranks
    return ranks


def _topol_equiv(mol: Mol) -> dict[int, int]:
    """Map each atom index to its topological class. The same dict on a hit."""

    structure = _structure(mol)
    if "topol_equiv" in structure:
        return structure["topol_equiv"]

    ranks = cip_ids(mol, include_stereo=False)
    classes = {atom.GetIdx(): ranks[atom.GetIdx()] for atom in mol.GetAtoms()}
    structure["topol_equiv"] = classes
    return classes


def molecule_formula(mol: Mol) -> Formula:
    """Heavy-atom counts, total hydrogens, and formal charge.

    Explicit hydrogens are counted through ``GetTotalNumHs`` on the heavy
    atom they belong to, not as a second copy. Cached as ``cache["formula"]``.
    The trace keeps its own ``formula`` and ``delta_formula``.
    """

    structure = _structure(mol)
    cached = structure.get("formula")
    if cached is not None:
        return cached

    counts: dict[str, int] = {}
    charge = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
        hydrogens = atom.GetTotalNumHs()
        if hydrogens:
            counts["H"] = counts.get("H", 0) + hydrogens
        charge += atom.GetFormalCharge()
    formula: Formula = {"counts": counts, "charge": charge}
    structure["formula"] = formula
    return formula


@overload
def copy_mol(mol: TracingMol) -> TracingMol: ...
@overload
def copy_mol(mol: ForestMol) -> ForestMol: ...
@overload
def copy_mol(mol: Mol) -> Mol: ...
def copy_mol(mol: Mol) -> Mol:
    """``Chem.Mol`` copy that also carries a deep-copied forest.

    The constructor itself does not keep ``_forest``. This function does,
    when the source has one, so the result is not :class:`NoForestMol`.
    The new molecule is the copy. The source is not edited.
    """

    out = Mol(mol)
    if is_forest(mol):
        # Same molecular structure: keep ``cache`` by identity.
        return _place_forest(out, forest_copy(mol._forest, same_structure=True))
    return out


def rw_copy(mol: Mol) -> RWMol:
    """Editable chemistry copy. The source is not edited.

    ``_forest`` is not carried. This copy is about to be edited, and a
    copied structure cache would answer questions about the molecule
    before those edits. :func:`copy_mol` is the copy that keeps the forest.
    """

    return RWMol(Mol(mol))


def reaction_from_smarts(smarts: str) -> ChemicalReaction:
    """Parse a SMARTS reaction once. The cache is process-wide, not per mol."""

    rxn = _REACTION_CACHE.get(smarts)
    if rxn is None:
        rxn = ReactionFromSmarts(smarts)
        rxn._setImplicitPropertiesFlag(False)
        _REACTION_CACHE[smarts] = rxn
    return rxn


def run_reactants(smarts: str, mol: Mol) -> tuple[tuple[NoForestMol, ...], ...]:
    """Run one cached reaction on ``mol``. Empty when RDKit refuses the run."""

    reaction = reaction_from_smarts(smarts)
    try:
        product_sets = reaction.RunReactants((mol,))
    except (RuntimeError, ValueError):
        return ()
    if not product_sets:
        return ()
    return tuple(product_sets)


@overload
def cannonicalize_order(
    mol: TracingMol, tracing_reset: bool = True
) -> tuple[TracingMol, str]: ...
@overload
def cannonicalize_order(mol: Mol, tracing_reset: bool = True) -> tuple[Mol, str]: ...
def cannonicalize_order(mol: Mol, tracing_reset: bool = True) -> tuple[Mol, str]:
    """Renumber into canonical SMILES order. Returns the new mol and that SMILES.

    ``RenumberAtoms`` builds a new molecule and the forest is copied onto it.
    The input is not edited. A traced input comes back traced.
    """

    csmi = MolToSmiles(mol, isomericSmiles=False)
    smiles_order = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))

    renumber_map = [0] * mol.GetNumAtoms()
    for new_pos, old_idx in enumerate(smiles_order):
        renumber_map[old_idx] = new_pos

    source = ensure_forest(mol)
    # Atom order changed: drop structure-dependent cache, then seed csmi.
    renumbered = _place_forest(
        RenumberAtoms(mol, renumber_map),
        forest_copy(source._forest, same_structure=False),
    )
    renumbered._forest["cache"] = {"csmi": csmi}

    if tracing_reset:
        _reordered_forest_labels(renumbered)

    return renumbered, csmi

# TODO: this function should ensure Mol atoms exactly matches the forest labels,
# or throw error. It's only a helper that's meant to work if labels were dropped from mol.
def _reordered_forest_labels(mol: Mol) -> None:
    if not is_tracing(mol):
        raise KeyError("atom_trace")
    trace = mol._forest["atom_trace"]
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        if atom.GetAtomicNum() != 1:
            tag = atom.GetProp("forestLabel")
            record = trace["records"][tag]
            if "idx" not in record:
                raise KeyError("idx")
            record["idx"][-1] = index
    _restamp_start_labels(mol)


def restamp_product_forest_last_layer(mol: Mol) -> None:
    """Rewrite the product's last ``atom_trace`` idx layer and restamp props.

    Product-only: does not edit a parent mol. Use after opt-in canonical
    emission remaps chemistry onto the lex representative (or after
    ``RenumberAtoms`` on a product). Updates ``idx[-1]`` from current atom
    positions via ``forestLabel``, then restamps CX ``atomLabel``.
    """

    if not is_tracing(mol):
        return
    _reordered_forest_labels(mol)
    _write_forest_labels(mol)




def mol_from_smiles(smiles: str) -> NoForestMol:
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("could not parse %r" % (smiles,))
    return mol


@overload
def as_mol(value: str) -> NoForestMol: ...
@overload
def as_mol(value: _MolT) -> _MolT: ...
def as_mol(value: Mol | str) -> Mol:
    if value is None:
        raise ValueError("mol is required")
    if isinstance(value, str):
        return mol_from_smiles(value)
    return value


def canon_smiles(value: Mol | str) -> str:
    """Canonical SMILES with atom-map numbers cleared on a copy."""

    mol = as_mol(value)
    copied = Mol(mol)
    for atom in copied.GetAtoms():
        atom.SetAtomMapNum(0)
    return MolToSmiles(copied, isomericSmiles=False)


def _bond_key(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _current_bond_map(mol: Mol) -> dict[tuple[int, int], float]:
    bonds = {}
    for bond in mol.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[_bond_key(left, right)] = bond.GetBondTypeAsDouble()
    return bonds


def _connected_components(mol: Mol, atoms: Iterable[int]) -> list[frozenset[int]]:
    atoms = set(atoms)
    seen: set[int] = set()
    systems = []
    for start in atoms:
        if start in seen:
            continue
        comp: set[int] = set()
        queue = deque([start])
        while queue:
            index = queue.popleft()
            if index in comp:
                continue
            comp.add(index)
            for neighbor in mol.GetAtomWithIdx(index).GetNeighbors():
                other = neighbor.GetIdx()
                if other in atoms and other not in comp:
                    queue.append(other)
        seen |= comp
        if len(comp) >= 2:
            systems.append(frozenset(comp))
    return systems


def _load_resonance(mol: Mol) -> None:
    """Cache kekulé bond maps and conjugated-atom sets on structure."""

    structure = _structure(mol)
    if "resonance_bonds" in structure:
        return

    base = Mol(mol)
    maps = []
    groups = []
    try:
        if SanitizeMol(base, catchErrors=True):
            raise ValueError("unsanitizable")
        supplier = ResonanceMolSupplier(base, KEKULE_ALL)
        n_groups = supplier.GetNumConjGrps()
        grouped = defaultdict(set)
        for atom in base.GetAtoms():
            group = supplier.GetAtomConjGrpIdx(atom.GetIdx())
            if 0 <= group < n_groups:
                grouped[group].add(atom.GetIdx())
        groups = [frozenset(values) for values in grouped.values() if len(values) >= 2]
        seen = set()
        for res in supplier:
            if res is None:
                continue
            bond_map = _current_bond_map(res)
            key = tuple(sorted(bond_map.items()))
            if key in seen:
                continue
            seen.add(key)
            maps.append(bond_map)
    except (ValueError, RuntimeError):
        maps = []
        groups = []

    if not groups:
        aromatic = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
        groups = _connected_components(mol, aromatic)
        if not groups:
            conjugated: set[int] = set()
            for bond in base.GetBonds():
                if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() >= 1.5:
                    conjugated.add(bond.GetBeginAtomIdx())
                    conjugated.add(bond.GetEndAtomIdx())
            groups = _connected_components(base, conjugated)

    if not maps:
        maps = [_current_bond_map(base)]

    structure["resonance_bonds"] = tuple(maps)
    structure["conjugated_systems"] = tuple(groups)


def resonance_bond_maps(mol: Mol) -> tuple[dict[tuple[int, int], float], ...]:
    _load_resonance(mol)
    structure = _structure(mol)
    if "resonance_bonds" not in structure:
        raise KeyError("resonance_bonds")
    return structure["resonance_bonds"]


def _conjugated_systems(mol: Mol) -> tuple[frozenset[int], ...]:
    _load_resonance(mol)
    structure = _structure(mol)
    if "conjugated_systems" not in structure:
        raise KeyError("conjugated_systems")
    return structure["conjugated_systems"]


def _aromatic_systems(mol: Mol) -> tuple[frozenset[int], ...]:
    structure = _structure(mol)
    if "aromatic_systems" not in structure:
        aromatic = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
        structure["aromatic_systems"] = tuple(_connected_components(mol, aromatic))
    return structure["aromatic_systems"]


def _pi_center(mol: Mol, index: int) -> bool:
    atom = mol.GetAtomWithIdx(index)
    if atom.GetIsAromatic():
        return True
    for bond in atom.GetBonds():
        if bond.GetIsAromatic() or bond.GetBondType() in (BondType.DOUBLE, BondType.TRIPLE):
            return True
    return False


def _conjugated_bond(mol: Mol, bond: Bond) -> bool:
    """True when the bond belongs to one conjugated system.

    Aromatic bonds cross only inside a ring, so a biaryl linker does not
    join two rings. A single bond crosses only from carbon to N, O, or S
    on a pi center. Fused aromatic bonds stay in one system.
    """

    order = bond.GetBondType()
    if order in (BondType.DOUBLE, BondType.TRIPLE):
        return True
    if bond.GetIsAromatic() and bond.IsInRing():
        return True
    if order != BondType.SINGLE:
        return False
    left = bond.GetBeginAtomIdx()
    right = bond.GetEndAtomIdx()
    elements = {
        mol.GetAtomWithIdx(left).GetAtomicNum(),
        mol.GetAtomWithIdx(right).GetAtomicNum(),
    }
    if 6 not in elements or not (elements & {7, 8, 16}):
        return False
    return _pi_center(mol, left) or _pi_center(mol, right)


def _conjugated_component(
    mol: Mol, start: int
) -> tuple[frozenset[int], frozenset[tuple[int, int]]]:
    atoms = {start}
    bonds: set[tuple[int, int]] = set()
    stack = [start]
    while stack:
        index = stack.pop()
        for bond in mol.GetAtomWithIdx(index).GetBonds():
            if not _conjugated_bond(mol, bond):
                continue
            other = bond.GetOtherAtomIdx(index)
            bonds.add(_bond_key(index, other))
            if other not in atoms:
                atoms.add(other)
                stack.append(other)
    return frozenset(atoms), frozenset(bonds)


def _kekule_slots(
    cache: KekuleParents,
) -> tuple[
    list[Mol],
    list[dict[tuple[int, int], float]],
    dict[frozenset[int], tuple[int, ...]],
    dict[tuple[tuple[int, int], float], int],
]:
    parents = cache.get("parents")
    orders = cache.get("orders")
    systems = cache.get("systems")
    by_order = cache.get("by_order")
    if parents is None:
        parents = []
        cache["parents"] = parents
    if orders is None:
        orders = []
        cache["orders"] = orders
    if systems is None:
        systems = {}
        cache["systems"] = systems
    if by_order is None:
        by_order = {}
        cache["by_order"] = by_order
    return parents, orders, systems, by_order


def _bond_order_sums(mol: Mol) -> dict[int, float]:
    return {
        atom.GetIdx(): sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        for atom in mol.GetAtoms()
    }


def move_charge_with_bonds(
    mol: Mol,
    before: dict[int, float],
    aromatic: set[int] | None = None,
) -> None:
    """Move formal charge when a bond-order flip would leave it behind.

    The oxygen whose bond order rose by one loses a negative charge. The
    oxygen whose bond order fell gains it. A neutral carbon keeps charge 0
    and moves hydrogen instead, because that hydrogen has to travel with
    the bond.

    A neutral aromatic atom is already the right charge. Kekulizing its
    1.5-order bonds is not a flip that should mint ``[n-]`` or ``[n+]``.
    Charged atoms still follow the bond, including a nitro oxygen.
    """

    aromatic = aromatic or set()
    for atom in mol.GetAtoms():
        old = before.get(atom.GetIdx())
        if old is None:
            continue
        new = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        delta = int(round(new - old))
        if delta == 0:
            continue
        neutral = atom.GetFormalCharge() == 0
        if neutral and atom.GetAtomicNum() == 6:
            _shift_hydrogens(atom, -delta)
        elif neutral and atom.GetIdx() in aromatic:
            continue
        else:
            atom.SetFormalCharge(atom.GetFormalCharge() + delta)


def _shift_hydrogens(atom: Atom, change: int) -> None:
    try:
        implicit = atom.GetNumImplicitHs()
    except RuntimeError:
        atom.UpdatePropertyCache(strict=False)
        implicit = atom.GetNumImplicitHs()
    total = atom.GetNumExplicitHs() + implicit
    updated = total + change
    if updated < 0:
        return
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(updated)


def _write_assignment(
    mol: Mol,
    atoms: frozenset[int],
    bonds: frozenset[tuple[int, int]],
    seed: tuple[int, int],
) -> tuple[NoForestMol, dict[tuple[int, int], float]] | None:
    """One kekulé assignment of ``atoms``. Other bonds stay as they were."""

    if seed[0] not in atoms or seed[1] not in atoms:
        return None
    adj: dict[int, list[int]] = {atom: [] for atom in atoms}
    for left, right in bonds:
        adj[left].append(right)
        adj[right].append(left)
    for nbrs in adj.values():
        nbrs.sort()
    carbons = sorted(atom for atom in atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6)
    doubles: dict[int, int] = {seed[0]: seed[1], seed[1]: seed[0]}

    def place(idx: int) -> bool:
        if idx == len(carbons):
            return True
        carbon = carbons[idx]
        if carbon in doubles:
            return place(idx + 1)
        for nbr in adj[carbon]:
            if nbr in doubles:
                continue
            doubles[carbon] = nbr
            doubles[nbr] = carbon
            if place(idx + 1):
                return True
            del doubles[carbon]
            del doubles[nbr]
        return False

    if not place(0):
        return None
    before = _bond_order_sums(mol)
    aromatic = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    rw = RWMol(Mol(mol))
    written: dict[tuple[int, int], float] = {}
    for left, right in bonds:
        bond = rw.GetBondBetweenAtoms(left, right)
        if bond is None:
            continue
        is_double = doubles.get(left) == right
        bond.SetBondType(BondType.DOUBLE if is_double else BondType.SINGLE)
        bond.SetIsAromatic(False)
        written[_bond_key(left, right)] = 2.0 if is_double else 1.0
    for atom in atoms:
        rw.GetAtomWithIdx(atom).SetIsAromatic(False)
    move_charge_with_bonds(rw, before, aromatic)
    return rw.GetMol(), written


def aromatic_parent_atoms(mol: Mol, start: int, end: int) -> frozenset[int] | None:
    """Aromatic atoms of the conjugated component that holds both ends.

    A biaryl single bond does not join two rings, and an exocyclic amide is
    not aromatic, so neither is kekulized with this ring. None when the ends
    do not share that component.
    """

    atoms, _bonds = _conjugated_component(mol, start)
    if end not in atoms:
        return None
    aromatic = frozenset(
        index for index in atoms if mol.GetAtomWithIdx(index).GetIsAromatic()
    )
    if start not in aromatic or end not in aromatic or len(aromatic) < 2:
        return None
    return aromatic


def _bonds_within(mol: Mol, atoms: frozenset[int]) -> frozenset[tuple[int, int]]:
    """Bonds whose ends are both in ``atoms``. The system's own edges."""

    bonds: set[tuple[int, int]] = set()
    for bond in mol.GetBonds():
        left = bond.GetBeginAtomIdx()
        right = bond.GetEndAtomIdx()
        if left in atoms and right in atoms:
            bonds.add(_bond_key(left, right))
    return frozenset(bonds)


def _store_assignments(
    mol: Mol,
    atoms: frozenset[int],
    bonds: frozenset[tuple[int, int]],
    cache: KekuleParents,
) -> tuple[int, ...]:
    """One parent per kekulé assignment of ``atoms``. Cached on that atom set."""

    parents, orders, systems, by_order = _kekule_slots(cache)
    held = systems.get(atoms)
    if held is not None:
        return held
    indexes: list[int] = []
    seen: set[tuple[tuple[tuple[int, int], float], ...]] = set()
    for bond in sorted(bonds):
        assigned = _write_assignment(mol, atoms, bonds, bond)
        if assigned is None:
            continue
        parent, bond_orders = assigned
        signature = tuple(sorted(bond_orders.items()))
        if signature in seen:
            continue
        seen.add(signature)
        index = len(parents)
        parents.append(parent)
        orders.append(bond_orders)
        indexes.append(index)
        for key, order in bond_orders.items():
            slot = (key, order)
            if slot not in by_order:
                by_order[slot] = index
    found = tuple(indexes)
    systems[atoms] = found
    return found


def ensure_kekule_parents(
    mol: Mol,
    left: int,
    right: int,
    cache: KekuleParents,
) -> tuple[int, ...]:
    """One parent per assignment of the system that contains ``(left, right)``.

    Other conjugated systems stay aromatic. Writes ``cache``. Does not read
    or write ``mol._forest``. Returns the parent indexes of that system.
    """

    atoms, bonds = _conjugated_component(mol, left)
    seed = _bond_key(left, right)
    if seed not in bonds and mol.GetBondBetweenAtoms(left, right) is not None:
        bonds = frozenset((*bonds, seed))
        if right not in atoms:
            atoms = frozenset((*atoms, right))
    return _store_assignments(mol, atoms, bonds, cache)


def parent_for_bond(
    cache: KekuleParents, left: int, right: int, order: float = 2.0
) -> Mol | None:
    """The cached parent in which ``(left, right)`` has ``order``.

    Does not search. Does not read or write ``mol._forest``.
    """

    by_order = cache.get("by_order")
    parents = cache.get("parents")
    if by_order is None or parents is None:
        return None
    index = by_order.get((_bond_key(left, right), order))
    if index is None:
        return None
    return parents[index]


def _ensure_atoms(mol: Mol, atom: int, cache: KekuleParents) -> frozenset[int]:
    atoms, bonds = _conjugated_component(mol, atom)
    _parents, _orders, systems, _by_order = _kekule_slots(cache)
    if atoms in systems:
        return atoms
    if not bonds:
        systems[atoms] = ()
        return atoms
    left, right = min(bonds)
    ensure_kekule_parents(mol, left, right, cache)
    return atoms


def parents_for_ends(
    mol: Mol,
    start: int,
    end: int,
    cache: KekuleParents,
    atoms: frozenset[int] | None = None,
) -> EndParents:
    """Parents covering ``start`` and ``end``.

    ``atoms``, when given, is the system the caller already chose (an
    aromatic component, or a conjugated one). Assignments stay inside that
    set, so an exocyclic amide is not rewritten just because a ring carbon
    touches it. Without ``atoms``, each end's conjugated component is used.

    Same system: that system's assignments. Different systems: each system's
    assignments, not a product of every system. Does not read or write
    ``mol._forest``.
    """

    if atoms is not None:
        bonds = _bonds_within(mol, atoms)
        if not bonds:
            return EndParents(parents=(), same_system=True)
        indexes = _store_assignments(mol, atoms, bonds, cache)
        parents = cache.get("parents") or []
        return EndParents(
            parents=tuple(parents[index] for index in indexes),
            same_system=True,
        )

    start_atoms = _ensure_atoms(mol, start, cache)
    end_atoms = _ensure_atoms(mol, end, cache)
    _parents, _orders, systems, _by_order = _kekule_slots(cache)
    parents = cache.get("parents") or []
    if start_atoms == end_atoms:
        indexes = systems.get(start_atoms, ())
        return EndParents(
            parents=tuple(parents[index] for index in indexes),
            same_system=True,
        )
    indexes = systems.get(start_atoms, ()) + systems.get(end_atoms, ())
    return EndParents(
        parents=tuple(parents[index] for index in indexes),
        same_system=False,
    )


def _ring_membership(mol: Mol) -> dict[int, tuple[tuple[int, ...], ...]]:
    structure = _structure(mol)
    if "rings" not in structure:
        work = Mol(mol)
        SanitizeMol(work, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        atom_rings = work.GetRingInfo().AtomRings()
        structure["rings"] = {
            idx: tuple(ring for ring in atom_rings if idx in ring)
            for idx in range(work.GetNumAtoms())
        }
    return structure["rings"]


def _smarts_matches(mol: Mol, smarts: str) -> tuple[dict[int, int], ...]:
    cache = _structure(mol).setdefault("smarts_matches", {})
    if smarts not in cache:
        query = MolFromSmarts(smarts)
        hits = []
        if query is not None:
            mapnos = [atom.GetAtomMapNum() for atom in query.GetAtoms()]
            for match in mol.GetSubstructMatches(query):
                mapped = {
                    mapno: idx for idx, mapno in zip(match, mapnos) if mapno
                }
                if 1 in mapped:
                    hits.append(mapped)
        cache[smarts] = tuple(hits)
    return cache[smarts]


def _bump(counters: EditCounters | None, name: str, amount: int = 1) -> None:
    if counters is None:
        return
    current = getattr(counters, name)
    if not isinstance(current, int):
        raise TypeError(name)
    setattr(counters, name, current + amount)


def _sanitize_piece(frag: Mol) -> bool:
    """True when ``frag`` sanitizes. The second try drops explicit H.

    Pair edits can leave ``[CH2]`` on a carbon whose new bonds already
    use those hydrogens. ``O=[CH2][CH2]=O`` is that draft of glyoxal.
    The bond orders are the product. The explicit count is not.
    """

    if not SanitizeMol(frag, catchErrors=True):
        return True
    for atom in frag.GetAtoms():
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    return not SanitizeMol(frag, catchErrors=True)


def carry_forest(src: Mol, dst: Mol) -> Mol:
    """Deep-copy ``src._forest`` onto ``dst``, remapped by ``forestLabel``.

    Prefer react → split → :func:`~xenosite.refactor_poc.rules.forest_trace`,
    so fragments are usually untraced when split and tracing is installed per
    piece. When ``src`` is already traced (copy / rare re-split), each call
    still installs a new forest dict: live ``records`` keep only labels on
    ``dst`` with the current index rewritten; cleavage siblings are dropped
    from live records (not left with stale indices). Prior ``deletes`` are
    copied only for same-size mols. Structure caches are cleared.
    """

    if not is_forest(src):
        return dst

    src_forest = src._forest
    child: Forest = empty_forest()
    src_imm = src_forest.get("immutable")
    if src_imm is not None:
        child["immutable"] = shallow_immutable(src_imm)

    try:
        same_size = src.GetNumAtoms() == dst.GetNumAtoms()
    except Exception:
        same_size = False

    if same_size and "is_terminal_product" in src_forest:
        child["is_terminal_product"] = src_forest["is_terminal_product"]

    trace = src_forest.get("atom_trace")
    if trace is not None and "records" in trace:
        label_to_idx: dict[str, int] = {}
        for atom in dst.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            if atom.HasProp("forestLabel"):
                label_to_idx[atom.GetProp("forestLabel")] = atom.GetIdx()

        records: dict = {}
        for tag, rec in trace["records"].items():
            label = str(tag)
            if label not in label_to_idx:
                # Cleavage sibling / foreign — absent on this fragment.
                continue
            rec_copy = copy_mutable(rec)
            idxs = rec_copy.get("idx")
            if not idxs:
                raise KeyError("idx")
            depths = rec_copy.get("depth")
            frame = int(trace["depth"])
            if depths is None or frame not in depths:
                raise KeyError(
                    "atom_trace label %s has no idx at depth %s" % (tag, frame)
                )
            rec_copy["idx"] = list(idxs)
            rec_copy["idx"][-1] = label_to_idx[label]
            records[label] = rec_copy

        carried = {
            "records": records,
            "deletes": (
                copy_mutable(trace.get("deletes") or {}) if same_size else {}
            ),
            "transforms": list(trace.get("transforms") or []),
            "additions": copy_mutable(trace.get("additions") or {}),
            "delta_formula": copy_mutable(trace.get("delta_formula") or {}),
            "depth": int(trace["depth"]),
            "last_tag": int(trace["last_tag"]),
            "next_transform": int(trace.get("next_transform") or 1),
            "formula": {"counts": {}, "charge": 0},
        }
        child["atom_trace"] = carried  # type: ignore[typeddict-item]
        held = _place_forest(dst, child)
        carried["formula"] = molecule_formula(held)
        child["cache"] = {}
        _restamp_start_labels(held)
        return held

    held = _place_forest(dst, child)
    _restamp_start_labels(held)
    return held


def sanitized_fragments(mol: Mol, counters: EditCounters | None = None) -> FragmentSplit:
    """Split, drop the dealkylation leaving group, sanitize, carry forest.

    Each piece is a connected mol. Empty pieces when any fragment fails
    sanitize. Callers read ``pieces``.
    """

    if isinstance(mol, RWMol):
        mol = mol.GetMol()
    frags = list(GetMolFrags(mol, asMols=True, sanitizeFrags=False)) or [mol]
    out = []
    for frag in frags:
        if any(atom.HasProp("dealk-noncarbon") for atom in frag.GetAtoms()):
            continue
        if not _sanitize_piece(frag):
            _bump(counters, "sanitize_dropped")
            return FragmentSplit(pieces=())
        out.append(carry_forest(mol, frag))
    return FragmentSplit(pieces=tuple(out))


def split_fragments(raw: Mol) -> FragmentSplit:
    """One mol, or the connected components of a disconnected product.

    Carries forest onto each piece. Structure caches are cleared.
    """

    try:
        groups = GetMolFrags(raw)
    except ValueError:
        return FragmentSplit(pieces=(raw,))
    if len(groups) <= 1:
        return FragmentSplit(pieces=(raw,))
    frags = list(GetMolFrags(raw, asMols=True, sanitizeFrags=False)) or [raw]
    return FragmentSplit(pieces=tuple(carry_forest(raw, frag) for frag in frags))


def _without_atoms(mol: Mol, drop: set[int]) -> tuple[NoForestMol, dict[int, int]]:
    """Copy with ``drop`` removed. The map sends each new index to the old one."""

    survivors = [index for index in range(mol.GetNumAtoms()) if index not in drop]
    new_to_old = {new: old for new, old in enumerate(survivors)}
    editable = RWMol(Mol(mol))
    for index in sorted(drop, reverse=True):
        editable.RemoveAtom(index)
    remainder = editable.GetMol()
    SanitizeMol(remainder, catchErrors=True)
    return remainder, new_to_old


def _remember(
    hits: tuple[tuple[int, ...], ...],
    into: list[tuple[int, ...]],
    seen: set[tuple[int, ...]],
) -> None:
    for hit in hits:
        if hit in seen:
            continue
        seen.add(hit)
        into.append(hit)


def _placements(
    reactant: Mol, target: Mol
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """Full-size embeddings, then smaller matches on the uncovered remainder.

    A later round deletes atoms a previous embedding already used and runs
    the same MCS again, so a smaller leftover is not truncated away.
    Orientations of one atom set stay; this does not pick a winner.
    """

    reactant_hits: list[tuple[int, ...]] = []
    target_hits: list[tuple[int, ...]] = []
    reactant_seen: set[tuple[int, ...]] = set()
    target_seen: set[tuple[int, ...]] = set()
    placed: set[frozenset[int]] = set()

    def absorb(piece: Mol, index_of: dict[int, int]) -> bool:
        query = _mcs_query(piece, target)
        if query is None or query.GetNumAtoms() < 2:
            return False
        raw = _full_matches(piece, query)
        if not raw or not _full_matches(target, query):
            return False
        pending: list[tuple[int, ...]] = []
        fresh: set[frozenset[int]] = set()
        for hit in raw:
            orig = tuple(index_of[index] for index in hit)
            atoms = frozenset(orig)
            if len(atoms) < 2 or orig in reactant_seen:
                continue
            if atoms in placed and atoms not in fresh:
                continue
            pending.append(orig)
            fresh.add(atoms)
        if not fresh:
            return False
        _remember(tuple(pending), reactant_hits, reactant_seen)
        _remember(_full_matches(target, query), target_hits, target_seen)
        placed.update(fresh)
        return True

    identity = {index: index for index in range(reactant.GetNumAtoms())}
    absorb(reactant, identity)
    covered = {index for hit in reactant_hits for index in hit}
    for _round in range(3):
        if len(covered) >= reactant.GetNumAtoms() - 1:
            break
        remainder, new_to_old = _without_atoms(reactant, covered)
        if remainder.GetNumHeavyAtoms() < 2:
            break
        if not absorb(remainder, new_to_old):
            break
        covered = {index for hit in reactant_hits for index in hit}
    return tuple(reactant_hits), tuple(target_hits)


def mcs_matches(reactant: Mol, target: Mol) -> McsResult:
    """Every placement of ``target`` on ``reactant``, not only the best.

    The first round is every full-size embedding. Later rounds match the
    uncovered remainder, so a smaller placement is not truncated away. The
    reactant structure holds both sides, keyed by the target's canonical
    SMILES. This function does not score them.
    """

    structure = _structure(reactant)
    cache: dict[str, McsResult] = structure.setdefault("mcs_matches", {})
    targets: dict[str, McsResult] = structure.setdefault("mcs_targets", {})
    held_target = ensure_forest(target)
    key = held_target.xf.csmi
    cached = cache.get(key)
    if cached is not None:
        return cached

    reactant_hits, target_hits = _placements(reactant, target)
    found = McsResult(embeddings=reactant_hits)
    cache[key] = found
    targets[key] = McsResult(embeddings=target_hits)
    return found


def mcs_target_matches(reactant: Mol, target: Mol) -> McsResult:
    """Target-side embeddings for the same MCS query. Filled with :func:`mcs_matches`."""

    mcs_matches(reactant, target)
    held_target = ensure_forest(target)
    key = held_target.xf.csmi
    structure = _structure(reactant)
    if "mcs_targets" not in structure:
        raise KeyError("mcs_targets")
    return structure["mcs_targets"][key]


def _mcs_query(reactant: Mol, target: Mol) -> NoForestMol | None:
    mcs = FindMCS(
        [reactant, target],
        atomCompare=AtomCompare.CompareElements,
        bondCompare=BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        timeout=2,
    )
    if mcs.numAtoms <= 0 or mcs.canceled:
        return None
    return MolFromSmarts(mcs.smartsString)


def _full_matches(mol: Mol, query: Mol | None) -> tuple[tuple[int, ...], ...]:
    if query is None:
        return ()
    size = query.GetNumAtoms()
    return tuple(
        match
        for match in mol.GetSubstructMatches(query, uniquify=False)
        if len(match) == size
    )

