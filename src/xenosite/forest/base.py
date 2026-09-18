"""Defines reaction archetypes to be parameterized or customized in rules.py"""

# Standard Library
import ast
import collections
import itertools
import logging
import re
from collections import defaultdict, deque
from contextlib import contextmanager
from copy import deepcopy

# Third Party
from .edit_guard import _OPEN_PROP, edit_mol
from .utils import (
    clean,
    merge,
    load,
    unmapped_smiles,
    canon_smi,
    is_rdkit_valid,
    note_sanitize_drop,
    _mol_smiles,
    refresh_mol,
)
from rdkit import Chem, rdBase
from rdkit.Chem.rdmolfiles import (
    MolToSmiles,
    MolFromSmiles,
    MolFromSmarts,
    CanonicalRankAtoms,
)
from rdkit.Chem.rdchem import (
    Mol,
    BondType,
    Atom,
    RWMol,
    ResonanceMolSupplier,
    KEKULE_ALL,
)
from rdkit.Chem.rdmolops import (
    SanitizeMol,
    Kekulize,
    SanitizeFlags,
    CombineMols,
    FragmentOnBonds,
    RenumberAtoms,
)
from rdkit.Chem import AllChem

# Prevents spammy rdkit messages
rdBase.DisableLog("rdApp.*")

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Private mol._forest bag (not a public API).
#
# Schema (reserved keys):
#   resonance   -> ResonanceCache (lazy joined forms / bfs_all_pairs)
#   atom_refs   -> AtomRefsIndex (creation lookup for AtomRef.added_by)
#                 Keys are (rule, site) as metabolize saw them on the parent;
#                 AtomRef.resolve projects site through atom_trace frames.
#   atom_trace  -> {
#                    "records": {label: record},   # live labels only
#                    "removed": [event, ...],      # chemical-removal history
#                    "last_tag": int,             # next label allocator
#                  }
#                 label matches the stable ``_forestLabel`` atom prop (not CX
#                 ``atomLabel``, which conjugates use for display names)
#                 stamped once; never rewritten. Remap via carry_forest after
#                 Mol()/GetMolFrags.
#
#                 record (live):
#                   idx, depth  — parallel GetIdx history at each tag depth
#                   added_by?   — optional (rule, site) when created by a step
#
#                 removed event:
#                   depth       — tag depth at which labels left this structure
#                   removed_by  — rule name / (rule, site) / True
#                   records     — {label: record} snapshot of those labels
#
#                 Cleavage: do **not** record sibling labels in ``removed``.
#                 carry_forest / tag drop labels absent from this mol; they
#                 belong on the other fragment. ``removed`` is only for atoms
#                 chemically deleted from this structure.
#
#                 AtomRef.resolve (origin): find the live label whose idx at
#                 the ref depth equals origin, then read that label's idx at
#                 the mol's latest tagged depth.
#
# Share with _copy_forest (same dict object) when caches should be shared.
# Product handoff / edited mols must **not** keep a parent's resonance cache
# (same atom count ≠ same bonding). Prefer carry_forest (no resonance share)
# or _clear_resonance after in-place edits. install_product_forest always
# drops resonance on reaction products.
# Do not store public phase1_steps JSON here (mol prop only).
# ---------------------------------------------------------------------------

_RESONANCE_CACHE_ENABLED = True


def _set_resonance_cache_enabled(enabled):
    """Private: enable/disable mol._forest resonance caching (tests / parity)."""
    global _RESONANCE_CACHE_ENABLED
    _RESONANCE_CACHE_ENABLED = bool(enabled)


@contextmanager
def _resonance_cache_disabled():
    """Private: temporarily disable resonance caching for parity comparisons."""
    prev = _RESONANCE_CACHE_ENABLED
    _set_resonance_cache_enabled(False)
    try:
        yield
    finally:
        _set_resonance_cache_enabled(prev)


def _forest_state(mol):
    """Get-or-create the private ``mol._forest`` dict."""
    forest = getattr(mol, "_forest", None)
    if forest is None:
        forest = {}
        mol._forest = forest
    return forest


def _clear_resonance(mol):
    """Drop cached resonance after any structural / kekule edit to ``mol``.

    Also drops a cached standardized (kekulized) view — it is only valid for
    the pre-edit structure.
    """
    forest = getattr(mol, "_forest", None)
    if forest is not None:
        forest.pop("resonance", None)
        forest.pop("standardized_mol", None)
    return mol


def _copy_forest(src, dst):
    """Share ``src._forest`` onto ``dst`` (same object) when present.

    Only for unedited views of the same structure. Prefer :func:`carry_forest`
    (independent dict, no resonance share) when ``dst`` may be sanitized,
    kekulized, or otherwise edited.
    """
    forest = getattr(src, "_forest", None)
    if forest is not None:
        dst._forest = forest


class AtomRefsIndex:
    """Maps AtomRef.added_by ``(rule, site)`` -> current GetIdx on a mol.

    ``site`` is whatever GetIdx frozenset metabolize saw on the parent at
    recording time (often the current frame of that step, not depth-0).
    ``AtomRef.resolve`` projects the caller's site through ``atom_trace``
    starting at ``AtomRef.depth`` (the frame those site idxs were written in).
    """

    __slots__ = ("_entries",)

    def __init__(self, entries=None):
        # (rule_name, frozenset[int]) -> int
        self._entries = dict(entries or ())

    def copy(self):
        return AtomRefsIndex(self._entries)

    def set(self, added_by, idx: int) -> None:
        rule, at = added_by
        self._entries[(rule, frozenset(at))] = int(idx)

    def lookup(self, added_by):
        rule, at = added_by
        return self._entries.get((rule, frozenset(at)))

    def items(self):
        return self._entries.items()

    def remap_from_parent(self, parent_mol, child_mol):
        """Rebuild entries with child GetIdx via react_atom_idx from parent."""
        out = AtomRefsIndex()
        for key, parent_idx in self._entries.items():
            for atom in child_mol.GetAtoms():
                if (
                    atom.HasProp("react_atom_idx")
                    and int(atom.GetProp("react_atom_idx")) == parent_idx
                ):
                    out._entries[key] = atom.GetIdx()
                    break
                if (
                    atom.HasProp("current_idx")
                    and int(atom.GetProp("current_idx")) == parent_idx
                ):
                    out._entries[key] = atom.GetIdx()
                    break
        return out


def _atom_refs_index(mol) -> AtomRefsIndex:
    forest = _forest_state(mol)
    index = forest.get("atom_refs")
    if index is None:
        index = AtomRefsIndex()
        forest["atom_refs"] = index
    return index


def _record_atom_creations(rule_name, origin_at, before, after, index: AtomRefsIndex):
    """Record atoms created by ``rule_name`` at ``origin_at`` onto ``index``."""
    origin_at = frozenset(origin_at or ())
    if not origin_at:
        return

    old_on_after = set()
    for atom in after.GetAtoms():
        if atom.GetAtomMapNum() > 0 or atom.HasProp("react_atom_idx"):
            old_on_after.add(atom.GetIdx())
    new_atoms = [
        atom for atom in after.GetAtoms() if atom.GetIdx() not in old_on_after
    ]
    if not new_atoms:
        return

    chosen = None
    if len(new_atoms) == 1:
        chosen = new_atoms[0]
    else:
        # ``origin_at`` is the formation site as metabolize / resolve saw it
        # (current-frame idxs on ``before``). Prefer a new atom bonded to one.
        targets = {int(x) for x in origin_at}
        for atom in new_atoms:
            for t in targets:
                if t < after.GetNumAtoms() and atom.GetIdx() != t and after.GetBondBetweenAtoms(
                    atom.GetIdx(), t
                ):
                    chosen = atom
                    break
            if chosen is not None:
                break
        if chosen is None:
            chosen = new_atoms[0]

    index.set((rule_name, origin_at), chosen.GetIdx())


def install_product_forest(
    parent,
    product,
    *,
    rule_name=None,
    origin_site=None,
    record_creation=True,
):
    """Install/remap ``_forest`` on ``product``; optionally record creations.

    Does **not** share ``resonance`` with ``parent`` (reaction products need a
    fresh cache). Remaps ``atom_refs`` from parent when the product does not
    already carry an index. Preserves ``atom_trace`` from tagging.
    """
    parent_forest = getattr(parent, "_forest", None) or {}
    child = getattr(product, "_forest", None)
    if child is None:
        child = {}
        product._forest = child

    # Reaction products must not reuse the reactant's resonance forms.
    child.pop("resonance", None)

    # Keep product atom_trace from tagging; do not pull parent's live dict.
    if "atom_refs" not in child or child.get("atom_refs") is None:
        prev = parent_forest.get("atom_refs")
        if prev is not None:
            child["atom_refs"] = prev.remap_from_parent(parent, product)
        else:
            child["atom_refs"] = AtomRefsIndex()

    if (
        record_creation
        and rule_name is not None
        and origin_site is not None
    ):
        _record_atom_creations(
            rule_name, origin_site, parent, product, child["atom_refs"]
        )
    return child["atom_refs"]


def carry_forest(src, dst, share_resonance=False):
    """Copy ``src._forest`` onto ``dst``, remapping idxs via stable ``_forestLabel``.

    RDKit ``Mol()`` / ``GetMolFrags`` drop Python attrs but keep atom props.
    Each heavy atom keeps a unique ``_forestLabel`` stamped once; this rebuilds
    ``atom_trace`` / ``atom_refs`` current idxs from those labels.

    Live ``atom_trace["records"]`` are filtered to labels on ``dst`` (cleavage
    siblings are dropped, not written to ``removed``). The ``removed`` event
    list is copied for same-size mols only so chemical-removal history survives
    ``copy_mol``; fragments start without the parent's removal log.

    Resonance is **not** carried by default. Same atom count does not imply the
    same bonding / resonance forms. Pass ``share_resonance=True`` only for an
    unedited identity view of the same structure; clear via
    :func:`_clear_resonance` after any subsequent edit.
    """
    src_forest = getattr(src, "_forest", None)
    if src_forest is None:
        return dst

    try:
        same_size = src.GetNumAtoms() == dst.GetNumAtoms()
    except Exception:
        same_size = False

    label_to_idx = {}
    for atom in dst.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.HasProp(AtomTracker.atom_tag_prop_name):
            label_to_idx[atom.GetProp(AtomTracker.atom_tag_prop_name)] = atom.GetIdx()

    child = {}
    if share_resonance and "resonance" in src_forest:
        child["resonance"] = src_forest["resonance"]

    trace = src_forest.get("atom_trace")
    if trace and "records" in trace:
        records = {}
        for tag, rec in trace["records"].items():
            label = str(tag)
            if label not in label_to_idx:
                # Cleavage sibling / foreign — drop (do not record as removed).
                continue
            rec_copy = deepcopy(rec)
            if rec_copy["depth"]:
                last_i = len(rec_copy["depth"]) - 1
                rec_copy["idx"][last_i] = label_to_idx[label]
            records[tag] = rec_copy

        removed = []
        if same_size:
            removed = [
                AtomTracker._copy_removed_event(event)
                for event in trace.get("removed") or ()
            ]

        child["atom_trace"] = {
            "records": records,
            "removed": removed,
            "last_tag": trace.get("last_tag", max(records) if records else 0),
        }

    prev_refs = src_forest.get("atom_refs")
    if prev_refs is not None:
        # Remap creation idxs by ``_forestLabel`` when the created atom still carries it.
        remapped = AtomRefsIndex()
        # Build idx->label on src for entries, then label->idx on dst.
        src_idx_to_label = {}
        for atom in src.GetAtoms():
            if atom.HasProp(AtomTracker.atom_tag_prop_name):
                src_idx_to_label[atom.GetIdx()] = atom.GetProp(
                    AtomTracker.atom_tag_prop_name
                )
        for key, src_idx in prev_refs.items():
            label = src_idx_to_label.get(src_idx)
            if label is not None and label in label_to_idx:
                remapped.set(key, label_to_idx[label])
            else:
                # Fall back to react_atom_idx remap.
                for atom in dst.GetAtoms():
                    if (
                        atom.HasProp("react_atom_idx")
                        and int(atom.GetProp("react_atom_idx")) == src_idx
                    ):
                        remapped.set(key, atom.GetIdx())
                        break
        child["atom_refs"] = remapped

    dst._forest = child
    return dst


def copy_mol(mol):
    """``Chem.Mol(mol)`` that also carries ``_forest`` (via :func:`carry_forest`).

    Use this instead of bare ``Chem.Mol`` / ``Mol`` whenever the copy must keep
    atom-trace or atom_refs history. Resonance is not shared (same size ≠ same
    bonding); a later edit must not see the source's forms. Throwaway copies
    for SMILES / sanitize probes can still use ``Chem.Mol`` directly.
    """
    if mol is None:
        return None
    out = Chem.Mol(mol)
    carry_forest(mol, out)  # share_resonance=False
    # An open-for-edit token must not leak onto copies.
    if out.HasProp(_OPEN_PROP):
        out.ClearProp(_OPEN_PROP)
    # Mol-level props (LAST_TAG, phase1_steps, …) are copied by Chem.Mol;
    # Python ``_forest`` is not — carry_forest restores it.
    return out


def _resonance_cache(mol):
    """Return the mol's ResonanceCache when caching is enabled, else None."""
    if not _RESONANCE_CACHE_ENABLED:
        return None
    forest = _forest_state(mol)
    cache = forest.get("resonance")
    if cache is None:
        cache = ResonanceCache()
        forest["resonance"] = cache
    return cache


class _ModeResonanceCache(object):
    """Lazily fills joined resonance forms for one systems mode (conjugated/aromatic)."""

    def __init__(self, mode):
        self.mode = mode
        self._items = []  # (joined_mol, system)
        self._done = False
        self._gen = None
        self.compute_count = 0

    @property
    def forms_materialized(self):
        return len(self._items)

    @property
    def exhausted(self):
        return self._done

    def iter_joined(self, resonate, mol):
        """Yield ``(Mol(copy), system)`` pull-through; fill cache on demand."""
        i = 0
        while True:
            if i < len(self._items):
                joined, system = self._items[i]
                i += 1
                yield copy_mol(joined), system
                continue
            if self._done:
                return
            if self._gen is None:
                self.compute_count += 1
                self._gen = self._produce(resonate, mol)
            try:
                item = next(self._gen)
            except StopIteration:
                self._done = True
                self._gen = None
                return
            self._items.append(item)

    def _produce(self, resonate, mol):
        for res_frag, rest, bonds, system in resonate._resfrags(
            mol, output_systems=True
        ):
            joined = resonate.join_fragments([res_frag] + list(rest), bonds)
            yield joined, system


class ResonanceCache(object):
    """Per-substrate resonance state stored on ``mol._forest['resonance']``."""

    def __init__(self):
        self._modes = {}
        self.bfs_paths = None
        self.bfs_compute_count = 0

    def mode(self, flag):
        entry = self._modes.get(flag)
        if entry is None:
            entry = _ModeResonanceCache(flag)
            self._modes[flag] = entry
        return entry


def can_smi(line="", rdmol=None):
    """Converts rdmol or line (in SMILES) to canonical SMILES."""

    if isinstance(rdmol, (list, tuple)):
        return [can_smi(rdmol=x) for x in rdmol]

    if rdmol:
        SanitizeMol(rdmol, catchErrors=True)
        line = unmapped_smiles(rdmol)

    if "." in line:
        return list(itertools.chain(*[can_smi(line=x) for x in line.split(".")]))

    split = line.split()
    smi = split[0]

    rdmol = MolFromSmiles(smi)
    if rdmol:
        out = unmapped_smiles(rdmol)
    else:
        out = smi

    out = re.sub(r"\[CH*\]", "C", out)
    out = re.sub(r"\[H\]", "", out)
    out = re.sub(r"\[C\]", "C", out)

    return [out]


def can_smi_set(rdmols):
    return frozenset(itertools.chain(*[can_smi(rdmol=x) for x in rdmols]))


class AtomTracker(object):
    tag_name = "ATOM_INDEX_PATHS"
    last_tag_name = "LAST_TAG"
    previous_index_prop_name = "current_idx"
    atom_tag_prop_name = "_forestLabel"

    def __init__(self, *args, **kwargs):
        super(AtomTracker, self).__init__()

    @staticmethod
    def topol_equiv(mol):
        copy = Mol(mol)
        SanitizeMol(copy, catchErrors=True)
        """ Get dict mapping atom indexes to topological IDs.

        >>> from rdkit import Chem
        >>> mol = MolFromSmiles('CCC')
        >>> AtomTracker.topol_equiv(mol)
        {0: 0, 1: 2, 2: 0}

        """

        return {
            a.GetIdx(): i
            for a, i in zip(
                mol.GetAtoms(),
                CanonicalRankAtoms(copy, includeChirality=False, breakTies=False),
            )
        }

    @staticmethod
    def site_to_topol_site(site, topol_equiv):
        """Map a ``(rule_name, atom_idxs)`` site to a topology identity.

        Rule suffixes (``_SmartsReactionRuleRxn0``, ``_12``, …) are stripped so
        matches share a pathway label. Atom indices become a sorted multiset of
        topological ranks (same signal as predict ``atoms.cipRank`` /
        XenoSite UI dedup).
        """
        name = site[0]
        if isinstance(name, str):
            name = name.split("_", 1)[0]
        ranks = tuple(sorted(topol_equiv[x] for x in site[1]))
        return name, ranks

    @staticmethod
    def add_current_idx_as_atom_prop(mol, propname="idx"):
        [
            a.SetProp(propname, str(a.GetIdx()))
            for a in mol.GetAtoms()
            if a.GetAtomicNum() != 1
        ]

    @classmethod
    def _copy_tag_records(cls, records):
        return deepcopy(records)

    @classmethod
    def _copy_removed_event(cls, event):
        """Deep-copy a ``removed`` event, filling default keys."""
        return {"depth": 0, "removed_by": True, "records": {}, **deepcopy(event)}

    @classmethod
    def _normalize_atom_trace(cls, trace):
        """Return ``atom_trace`` with default keys filled in."""
        out = {"records": {}, "removed": [], "last_tag": 0, **dict(trace or {})}
        out.setdefault("records", {})
        out.setdefault("removed", [])
        if out.get("last_tag") is None:
            records = out["records"] or {}
            out["last_tag"] = max(records) if records else 0
        return out

    @classmethod
    def _load_atom_trace(cls, mol, strict=True):
        """Return full ``atom_trace`` dict from ``_forest`` or legacy prop."""
        forest = getattr(mol, "_forest", None) or {}
        trace = forest.get("atom_trace")
        if trace is not None and "records" in trace:
            normalized = cls._normalize_atom_trace(trace)
            forest["atom_trace"] = normalized
            return normalized

        try:
            raw = mol.GetProp(cls.tag_name)
        except KeyError as err:
            if strict:
                raise err
            return {"records": {}, "removed": [], "last_tag": 0}

        try:
            records = ast.literal_eval(raw)
        except (SyntaxError, ValueError) as err:
            if strict:
                raise err
            return {"records": {}, "removed": [], "last_tag": 0}

        try:
            last_tag = int(mol.GetProp(cls.last_tag_name))
        except (KeyError, ValueError):
            last_tag = max(records) if records else 0

        trace = {"records": records, "removed": [], "last_tag": last_tag}
        # Migrate once onto _forest so later reads skip literal_eval.
        _forest_state(mol)["atom_trace"] = trace
        return trace

    @classmethod
    def tags(cls, record, depth=None, idx=None, strict=True, compact=False, **kwargs):
        # Tag records use 0-based RDKit indices (GetIdx()). compact_tags()
        # converts those to 1-based atom numbers (SMILES :N) when compact=True.

        if isinstance(record, list):
            return itertools.chain(
                *[
                    cls.tags(x, depth=depth, idx=idx, strict=strict, compact=compact)
                    for x in record
                ]
            )
        if isinstance(record, Mol):
            trace = cls._load_atom_trace(record, strict=strict)
            record = trace["records"]
        elif not isinstance(record, dict):
            raise ValueError("Must submit RDKit Mol or dict")

        if depth is not None:
            record = {
                tag: data
                for tag, data in list(record.items())
                if depth in data["depth"]
            }

        if idx is not None:
            record = {
                tag: data for tag, data in list(record.items()) if idx in data["idx"]
            }

        if idx is None and depth is None:
            record = cls._copy_tag_records(record)

        if compact:
            return cls.compact_tags(record)
        else:
            return record

    @staticmethod
    def compact_tags(record, adjust_root_by=1):
        # 1-based atom numbers (SMILES / map convention). Internal idx lists stay
        # 0-based RDKit GetIdx(); this is the bridge used by AtomTrace.
        return {
            k: {d: i + adjust_root_by for d, i in zip(v["depth"], v["idx"])}
            for k, v in list(record.items())
        }

    @classmethod
    def depths(cls, record, strict=True):
        if isinstance(record, Mol):
            record = cls.tags(record, strict=strict)
        elif not isinstance(record, dict):
            raise ValueError("Must submit RDKit Mol or dict")

        return sorted(
            set(itertools.chain(*[x["depth"] for x in list(record.values())]))
        )

    @classmethod
    def next_depth(cls, previous_tags):
        depth = 1

        previous = cls.depths(previous_tags)
        if previous:
            depth += max(previous)

        return depth

    @classmethod
    def metabolite_index_to_reversed_index_record(
        cls, metabolite, exact_depth=2, strict=True
    ):

        depths = cls.depths(metabolite, strict=strict)
        if len(depths) < exact_depth:
            return None

        reversed_depth = list(reversed(depths[-exact_depth:]))

        idx_record = cls.tags(metabolite, compact=True)

        metabolite_index_to_reversed_index_record = defaultdict(list)

        for depth_to_idx in list(idx_record.values()):
            if set(reversed_depth) != set(depth_to_idx):
                continue

            metabolite_idx = depth_to_idx[reversed_depth[-1]]

            for depth in reversed_depth:
                metabolite_index_to_reversed_index_record[metabolite_idx].append(
                    depth_to_idx[depth]
                )

        return metabolite_index_to_reversed_index_record

    def initialize_tags(self, mol):
        # Tags use 0-based RDKit indices (GetIdx()). Public AtomTrace converts to
        # 1-based atom numbers. Each heavy atom gets a unique ``_forestLabel`` once;
        # ``_forest["atom_trace"]`` maps those labels to idxs across depths.

        forest = getattr(mol, "_forest", None) or {}
        if forest.get("atom_trace") is not None:
            return
        # Labels without forest: ``Chem.Mol()`` / RDKit copy dropped ``_forest``
        # but kept ``_forestLabel``. Rebuild a live snapshot from those labels
        # so metabolize/tag can continue (pre-copy history is unavailable;
        # current idxs become depth 0 of the restored trace).
        labeled = [
            a
            for a in mol.GetAtoms()
            if a.GetAtomicNum() != 1 and a.HasProp(self.atom_tag_prop_name)
        ]
        if labeled:
            records = {}
            for atom in labeled:
                label = atom.GetProp(self.atom_tag_prop_name)
                try:
                    tag = int(label)
                except ValueError:
                    tag = label
                records[tag] = {"idx": [atom.GetIdx()], "depth": [0]}
            self._save_tags(mol, records)
            return

        if mol.HasProp(self.tag_name):
            # Legacy string prop — migrate via tags().
            self.tags(mol, strict=False)
            return

        initial_tags = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue

            atom.SetProp("initial", "1")
            idx = atom.GetIdx()
            initial_tags[idx] = {"idx": [idx], "depth": [0]}
            atom.SetProp(self.atom_tag_prop_name, str(idx))

        self._save_tags(mol, initial_tags)

    def tag(self, product, reactant, strict=True, **kwargs):
        """Copy live atom tags from reactant onto product; advance depth.

        Live ``records`` hold only labels present on ``product``. Labels on the
        reactant frontier that are absent here are dropped (cleavage siblings
        belong on the other fragment — not written to ``removed``). Prior
        chemical-removal events on ``reactant`` are carried forward in
        ``removed``. Pass ``record_removals=True`` to append absent frontier
        labels as a new ``removed`` event (in-place chemical deletion).
        """

        prev_trace = self._load_atom_trace(reactant, strict=strict)
        previous_tags = prev_trace.get("records") or {}
        removed = [
            self._copy_removed_event(ev) for ev in (prev_trace.get("removed") or ())
        ]
        next_depth = self.next_depth(previous_tags)

        frontier = next_depth - 1
        new_tags = {}
        absent = {}
        for tag, rec in previous_tags.items():
            if frontier not in rec["depth"]:
                continue
            record = deepcopy(rec)
            new_tags[tag] = record

        old_to_new_atom_indexes = self._old_to_new_atom_indexes(product)

        new_atom_indexes = [
            a.GetIdx() for a in product.GetAtoms() if a.GetAtomicNum() != 1
        ]

        for atom_unique_tag, record in list(new_tags.items()):
            last_index = record["idx"][-1]
            if last_index in old_to_new_atom_indexes:
                new_index = old_to_new_atom_indexes[last_index]
                record["idx"].append(new_index)
                record["depth"].append(record["depth"][-1] + 1)
                new_atom_indexes.remove(new_index)
            else:
                # Absent on this product — cleavage sibling or chemical loss.
                absent[atom_unique_tag] = new_tags.pop(atom_unique_tag)

        if absent and kwargs.get("record_removals"):
            removed.append(
                self._copy_removed_event(
                    {
                        "depth": next_depth,
                        "removed_by": kwargs.get(
                            "removed_by", getattr(self, "name", True)
                        ),
                        "records": absent,
                    }
                )
            )
        # else: drop absences (default) — do not record cleavage.

        self._save_tags(
            product,
            self._tag_new_atoms(
                new_tags, new_atom_indexes, next_depth, self._next_tag(reactant)
            ),
            removed=removed,
        )

    @classmethod
    def _next_tag(cls, mol):
        forest = getattr(mol, "_forest", None) or {}
        trace = forest.get("atom_trace")
        if trace is not None and "last_tag" in trace:
            try:
                return int(trace["last_tag"]) + 1
            except (TypeError, ValueError):
                pass
        try:
            last_tag = mol.GetProp(cls.last_tag_name)
        except KeyError:
            last_tag = "0"

        try:
            last_tag_depth = int(last_tag)
        except ValueError:
            last_tag_depth = 0

        return last_tag_depth + 1

    def _tag_new_atoms(self, new_tags, untagged_atom_indexes, next_depth, next_tag):

        for idx in untagged_atom_indexes:
            new_tags[next_tag] = {"idx": [idx], "depth": [next_depth]}
            next_tag += 1

        return new_tags

    def _old_to_new_atom_indexes(self, mol):
        # 0-based RDKit indices. current_idx / react_atom_idx are GetIdx() on
        # the reactant; values are GetIdx() on this product.

        mapping = {}
        for a in mol.GetAtoms():
            if a.HasProp(self.previous_index_prop_name):
                mapping[int(a.GetProp(self.previous_index_prop_name))] = a.GetIdx()
            elif a.HasProp("react_atom_idx"):
                mapping[int(a.GetProp("react_atom_idx"))] = a.GetIdx()
        return mapping

    def _clear_atom_maps(self, mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

    def _current_to_origin(self, tags, last_depth, level=0):
        """Map current 0-based idx -> origin GetIdx() at ``level``, or None if added."""
        origin_of = {}
        for rec in tags.values():
            if last_depth not in rec["depth"]:
                continue
            cur = rec["idx"][rec["depth"].index(last_depth)]
            if level in rec["depth"]:
                origin_of[cur] = rec["idx"][rec["depth"].index(level)]
            else:
                origin_of[cur] = None
        return origin_of

    def _stamp_origin_maps(self, mol, tags=None, level=0):
        """Set molAtomMapNumber to the 1-based origin atom number at ``level``.

        New atoms stay 0 (unmapped). Hydrogens are never mapped.
        """
        from .trace import atom_no

        if tags is None:
            try:
                tags = self.tags(mol)
            except (KeyError, SyntaxError, ValueError):
                return

        depths = self.depths(tags)
        if not depths:
            return
        last_depth = max(depths)
        origin_of = self._current_to_origin(tags, last_depth, level=level)

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                atom.SetAtomMapNum(0)
                continue
            orig = origin_of.get(atom.GetIdx())
            if orig is None:
                atom.SetAtomMapNum(0)
            else:
                atom.SetAtomMapNum(atom_no(orig))

    def _reactant_aligned_order(self, mol, origin_of):
        """Atom order following level-0 reactant order, with additions interleaved.

        Walk surviving atoms in increasing origin index. After each survivor,
        emit new atoms attached to it (BFS over the new-atom subgraph). If a
        new atom bridges two survivors, attach it after the lower origin.
        Each disconnected fragment is ordered on its own.
        """
        n_atoms = mol.GetNumAtoms()
        order = []
        emitted = set()

        def neighbors(idx):
            return [a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors()]

        def is_heavy_new(idx):
            return origin_of.get(idx) is None and mol.GetAtomWithIdx(idx).GetAtomicNum() != 1

        def component_min_survivor_origin(start):
            seen = set()
            stack = [start]
            min_orig = None
            while stack:
                a = stack.pop()
                if a in seen:
                    continue
                seen.add(a)
                for ni in neighbors(a):
                    if is_heavy_new(ni):
                        stack.append(ni)
                    elif origin_of.get(ni) is not None:
                        o = origin_of[ni]
                        if min_orig is None or o < min_orig:
                            min_orig = o
            return min_orig

        for frag in Chem.GetMolFrags(mol, asMols=False):
            survivors = sorted(
                (i for i in frag if origin_of.get(i) is not None),
                key=lambda i: origin_of[i],
            )
            for s in survivors:
                if s not in emitted:
                    order.append(s)
                    emitted.add(s)
                queue = deque()
                for ni in neighbors(s):
                    if ni in emitted or not is_heavy_new(ni):
                        continue
                    if component_min_survivor_origin(ni) == origin_of[s]:
                        queue.append(ni)
                while queue:
                    nidx = queue.popleft()
                    if nidx in emitted:
                        continue
                    order.append(nidx)
                    emitted.add(nidx)
                    for ni in neighbors(nidx):
                        if ni not in emitted and is_heavy_new(ni):
                            queue.append(ni)
            for i in frag:
                if i not in emitted and mol.GetAtomWithIdx(i).GetAtomicNum() != 1:
                    order.append(i)
                    emitted.add(i)
            for i in frag:
                if i not in emitted:
                    order.append(i)
                    emitted.add(i)

        for i in range(n_atoms):
            if i not in emitted:
                order.append(i)
        return order

    def _align_and_stamp(self, product):
        """Renumber to reactant-aligned order, rewrite last-depth tags, stamp maps."""
        try:
            tags = self.tags(product)
        except (KeyError, SyntaxError, ValueError):
            return product

        depths = self.depths(tags)
        if not depths:
            return product
        last_depth = max(depths)
        origin_of = self._current_to_origin(tags, last_depth, level=0)
        new_order = self._reactant_aligned_order(product, origin_of)

        if new_order != list(range(product.GetNumAtoms())):
            # RDKit RenumberAtoms drops molecule-level props (e.g. phase1_steps).
            mol_props = {
                name: product.GetProp(name) for name in product.GetPropNames()
            }
            # Preserve private _forest across renumber (RDKit drops Python attrs).
            forest = getattr(product, "_forest", None)
            product = RenumberAtoms(product, new_order)
            if forest is not None:
                product._forest = forest
            for name, value in mol_props.items():
                if not product.HasProp(name):
                    product.SetProp(name, value)
            old_to_new = {old: new for new, old in enumerate(new_order)}
            for rec in tags.values():
                if last_depth not in rec["depth"]:
                    continue
                i = rec["depth"].index(last_depth)
                rec["idx"][i] = old_to_new[rec["idx"][i]]
            self._save_tags(product, tags)

        # Keep current_idx as the reactant GetIdx() so a wrapping RuleSet.metabolize
        # can tag() again. Next-step metabolize resets it via add_current_idx_as_atom_prop.
        self._stamp_origin_maps(product, tags, level=0)
        return product

    def _save_tags(self, mol, tags, removed=None):
        if isinstance(tags, defaultdict):
            tags = {x: y for x, y in list(tags.items())}

        last_tag = max(tags) if tags else 0
        forest = _forest_state(mol)
        prev = forest.get("atom_trace") or {}
        if removed is None:
            removed = list(prev.get("removed") or [])
        forest["atom_trace"] = self._normalize_atom_trace(
            {"records": tags, "removed": removed, "last_tag": last_tag}
        )

        # One stable ``_forestLabel`` per heavy atom (never rewritten). History
        # lives only in ``_forest["atom_trace"]`` keyed by that label.
        # (CX ``atomLabel`` is reserved for conjugate display names.)
        if tags:
            depths = self.depths(tags)
            last_depth = max(depths) if depths else None
            for tag, rec in tags.items():
                if last_depth is None or last_depth not in rec["depth"]:
                    continue
                idx = rec["idx"][rec["depth"].index(last_depth)]
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 1:
                    continue
                label = str(tag)
                if not atom.HasProp(self.atom_tag_prop_name):
                    atom.SetProp(self.atom_tag_prop_name, label)

        # Do not write the growing ATOM_INDEX_PATHS string. Keep LAST_TAG for
        # callers / legacy _next_tag fallbacks that still read the mol prop.
        mol.SetProp(self.last_tag_name, str(last_tag))

    @staticmethod
    def all_atom_prop_names(mol):
        return frozenset.union(*[frozenset(a.GetPropNames()) for a in mol.GetAtoms()])

    @staticmethod
    def set_atom_prop(mol, idx, prop, val):
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp(prop, str(val))


class ConjugatedSystems(object):
    """Fragments molecules on the basis of the their aromatic or conjugated systemes.
    Currently, for most purposes, such as quinone formation, what is really of interest is
    AROMATIC systems, so this is the default.
    However, for other purposes, such as general hydrogenation reactions, of interest
    are the more broadly defined CONJUGATED systems.
    """

    flag = "conjugated"

    def fragments(self, mol):
        """Returns the aromatic or conjugated systems of mol.

        >>> CS = ConjugatedSystems()
        >>> mol = MolFromSmiles("C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1")
        >>> fragment, other_fragments, system, bond_types = next(CS.fragments(mol))
        >>> system
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
        >>> canon_smi(fragment) == canon_smi('[10*]C1=CC2=CC=CC=C2C=C1')
        True
        >>> canon_smi(other_fragments) == canon_smi(['[7*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1'])
        True
        """

        # This makes a copy to prevent the input molecule being unexpectedly modified
        mol = copy_mol(mol)

        # The input molecule will be fragmented and rejoined.
        # There is no guarantee that the atom ordering will remain the same.
        # Consequently, it is important to use the SetProp function on rdkit atoms to track atoms.

        AtomTracker.add_current_idx_as_atom_prop(mol)
        AtomTracker.add_current_idx_as_atom_prop(mol, propname="reactant_idx2")

        SanitizeMol(mol, catchErrors=True)

        systems = self.systems(mol)

        if len(systems) == 1:
            yield mol, [], systems[0], {}

        else:
            for system in systems:
                fragments, bond_types = self._separate_system(mol, system)

                yield fragments[0], fragments[1:], system, bond_types

    def _separate_system(self, mol, system):
        """Fragments mol based on atom indexes specified in system.

        >>> mol = MolFromSmiles("C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1")
        >>> CS = ConjugatedSystems()
        >>> fragments,bond_types = CS._separate_system(mol,set([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        >>> canon_smi(fragments) == canon_smi(['[10*]C1=CC2=CC=CC=C2C=C1', '[7*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1'])
        True

        """
        mol = copy_mol(mol)

        bonds_to_break = []
        for idx in system:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("conjugated_system", "1")

            for bond in atom.GetBonds():
                bond_atom_indexes = set([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

                if not bond_atom_indexes.issubset(system):
                    bonds_to_break.append(bond_atom_indexes)

        bonds_broken, bond_types = self._fragment_on_bonds(mol, bonds_to_break)

        fragments = Chem.GetMolFrags(bonds_broken, asMols=True, sanitizeFrags=False)
        ordered_fragments = sorted(
            fragments,
            key=lambda f: 1
            not in [a.HasProp("conjugated_system") for a in f.GetAtoms()],
        )
        return ordered_fragments, bond_types

    def _add_orphan_double_bond_systems(self, sys_id_to_atom_id, mol):
        """Adds any double bonds in mol to sys_id_to_atom_id that are not already specified."""

        if sys_id_to_atom_id:
            assigned = set.union(*list(sys_id_to_atom_id.values()))
        else:
            assigned = set([])

        unassigned = [
            set([b.GetBeginAtomIdx(), b.GetEndAtomIdx()]) - assigned
            for b in mol.GetBonds()
            if b.GetBondTypeAsDouble() in [1.5, 2]
        ]

        for group in unassigned:
            if not group:
                continue
            if sys_id_to_atom_id:
                sys_id_to_atom_id[max(sys_id_to_atom_id.keys()) + 1] = group
            else:
                sys_id_to_atom_id[0] = group

    def systems(self, mol):
        """Returns disjoint conjugated systems.

        >>> mol = MolFromSmiles("C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C")
        >>> ConjugatedSystems().systems(mol)
        [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {11, 12, 13, 14, 15, 16, 17, 18}]

        """
        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)

        with edit_mol(mol):
            try:
                Kekulize(mol, clearAromaticFlags=True)
            except ValueError:
                pass

        try:
            supplier = ResonanceMolSupplier(mol, KEKULE_ALL)
        except ValueError:
            try:
                supplier = ResonanceMolSupplier(mol)
            except ValueError:
                return [set()]

        atoms_to_groups = {
            a.GetIdx(): supplier.GetAtomConjGrpIdx(a.GetIdx()) for a in mol.GetAtoms()
        }

        groups = list(range(supplier.GetNumConjGrps()))

        sys_id_to_atom_id = collections.defaultdict(set)

        for idx, group in atoms_to_groups.items():
            if group in groups:
                sys_id_to_atom_id[group].add(idx)

        self._add_orphan_double_bond_systems(sys_id_to_atom_id, mol)

        return list(sys_id_to_atom_id.values())

    def _fragment_on_bonds(self, mol, bonds_to_break_by_atom_idx):
        """Breaks mol on all atom pairs specified,
        after first setting atom properties to enable later reconstruction of the input molecule.

        >>> AS = ConjugatedSystems()
        >>> mol = MolFromSmiles('C=CC(=CCC(Cl)C(F))C')
        >>> fragments,bond_types = AS._fragment_on_bonds(mol,[(0,1),(7,8)])
        >>> canon_smi(fragments) == canon_smi('*=CC(C)=CCC(Cl)C[8*].[1*]=C.[7*]F')
        True
        >>> list(bond_types.keys())
        [frozenset({0, 1}), frozenset({8, 7})]
        >>> list(bond_types.values())
        [rdkit.Chem.rdchem.BondType.DOUBLE, rdkit.Chem.rdchem.BondType.SINGLE]

        """
        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        with edit_mol(mol):
            try:
                Kekulize(mol, clearAromaticFlags=True)
            except ValueError:
                pass

        # The input molecule will be fragmented and rejoined,
        # and there is no guarantee that the atom ordering will remain the same.
        # Consequently, it is important to use the SetProp function on rdkit atoms to track atoms.
        # [a.SetProp('idx', str(a.GetIdx())) for a in mol.GetAtoms()]
        AtomTracker.add_current_idx_as_atom_prop(mol)

        if not bonds_to_break_by_atom_idx:
            return mol, {}

        bonds_to_break = []
        bond_types = {}

        for idx1, idx2 in bonds_to_break_by_atom_idx:
            bond_to_break = mol.GetBondBetweenAtoms(idx1, idx2)
            bond_types[frozenset([idx1, idx2])] = bond_to_break.GetBondType()
            bonds_to_break.append(bond_to_break.GetIdx())

        return FragmentOnBonds(mol, bonds_to_break, addDummies=True), bond_types

    def join_fragments(self, fragments, bond_types):
        """Calls _combine_fragments, _remove_dummy_atoms, and _restore_bonds in order to cleanly
        reproduce the original molecule prior to _fragment_on_bonds."""
        mol = self._combine_fragments(fragments)
        mol = self._remove_dummy_atoms(mol)
        return self._restore_bonds(mol, bond_types)

    def _combine_fragments(self, fragments):
        """Combines rdkit molecules into a single molecule and reorders by their original atom
        indexes. Note that this method does not add bonds, this is performed by the
        _restore_bonds method."""
        fragments = list(fragments)
        combined = fragments.pop()
        while fragments:
            combined = CombineMols(combined, fragments.pop())
        return self._recover_original_atom_order(combined)

    def _restore_bonds(self, mol, bond_types):
        """Adds bonds to mol. bond_types is a dict mapping from tuples of atom indexes to RDKit
        bond types."""

        with edit_mol(mol):
            for idxs, bond in bond_types.items():
                if mol.GetBondBetweenAtoms(*idxs):
                    mol.GetBondBetweenAtoms(*idxs).SetBondType(bond)
                else:
                    emol = Chem.rdchem.EditableMol(mol)
                    with edit_mol(emol):
                        emol.AddBond(*idxs, order=bond)
                    mol = emol.GetMol()
        return mol

    def _remove_dummy_atoms(self, mol):
        """Delete FragmentOnBonds attachment dummies only.

        Conjugation products may already carry a bare ``*`` adduct. Removing those
        shifts atom indices and breaks ``_restore_bonds`` (Full BFS depth≥2).
        """

        def _fragment_dummy_idxs(m):
            return [
                atom.GetIdx()
                for atom in m.GetAtoms()
                if atom.GetSymbol() == "*"
                and atom.HasProp("_forest_fragment_dummy")
                and atom.GetBoolProp("_forest_fragment_dummy")
            ]

        dummies = _fragment_dummy_idxs(mol)
        while dummies:
            emol = Chem.rdchem.EditableMol(mol)
            with edit_mol(emol):
                emol.RemoveAtom(dummies.pop())
            mol = emol.GetMol()
            dummies = _fragment_dummy_idxs(mol)
        return mol

    def _add_idx_prop_to_dummy_atoms(self, combined):
        """Detect FragmentOnBonds dummies and assign idx props."""
        dummies = [x for x in combined.GetAtoms() if not x.HasProp("idx")]
        originals = [
            int(x.GetProp("idx")) for x in combined.GetAtoms() if x.HasProp("idx")
        ]
        start = max(originals) + 1 if originals else 0
        for idx, dummy in enumerate(dummies, start=start):
            dummy.SetProp("idx", str(idx))
            # Distinguish from conjugation star adducts (also atomic num 0).
            dummy.SetBoolProp("_forest_fragment_dummy", True)

    def _recover_original_atom_order(self, combined):
        """Uses the idx atom property to reorder atoms by their indexes prior to fragmentation."""
        self._add_idx_prop_to_dummy_atoms(combined)

        original_order = [
            a.GetIdx()
            for a in sorted(combined.GetAtoms(), key=lambda x: int(x.GetProp("idx")))
        ]

        SanitizeMol(combined, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)

        return RenumberAtoms(combined, original_order)


class AromaticSystems(ConjugatedSystems):
    flag = "aromatic"

    def systems(self, mol):
        """Returns disjoint aromatic systems.

        >>> mol = MolFromSmiles("C1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=C")
        >>> AromaticSystems().systems(mol)
        [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {11, 12, 13, 14, 15, 16}]

        """

        atoms2neighbors = {
            a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()
        }

        valid = [a.GetIdx() for a in mol.GetAromaticAtoms()]

        groups = []
        for idx in valid:
            groups.append([idx] + [x for x in atoms2neighbors[idx] if x in valid])

        return merge(groups)


class QueryMol(object):
    """Tools for querying a molecule with SMARTS
    or for intrinsic propertys inferrable from an rdmol."""

    def __init__(self, query_smarts=None, *args, **kwargs):
        super(QueryMol, self).__init__(*args, **kwargs)

        if query_smarts is None:
            self.queries = []
            self.query_smarts = []

        else:
            self.query_smarts = query_smarts
            self.queries = []
            for entry in query_smarts:
                name, smarts, options = self._parse_query_smarts_entry(entry)
                query, queryidx2mapid = self._prepare_query(smarts)
                self.queries.append((name, query, queryidx2mapid, options))

        self.valid_modifications = set(
            itertools.chain(
                *[[x[0]] if isinstance(x[0], str) else x[0] for x in self.queries]
            )
        )

    @staticmethod
    def _parse_query_smarts_entry(entry):
        """Normalize ``(mods, smarts)`` or ``(mods, smarts, options)`` entries.

        ``options`` may include:

        - ``pathway``: optional uncommon chemotype name (off unless enabled)
        - ``one_sided`` / ``pair_limit``: max ends in a pair that may use this pathway
        - ``recommend``: target SMARTS; guided search opts the pathway in when it
          matches and formula hints agree
        - ``formula_hint``: :class:`~xenosite.forest.path_context.FormulaHint`,
          :class:`~xenosite.forest.path_context.FormulaAny`,
          :class:`~xenosite.forest.path_context.FormulaMatch`, a sequence of
          those, or ``callable(mol, match)`` — effect may depend on the hit
        """
        if not isinstance(entry, (tuple, list)) or len(entry) not in (2, 3):
            raise TypeError(
                "query_smarts entry must be (mods, smarts) or "
                "(mods, smarts, options), got %r" % (entry,)
            )
        name, smarts = entry[0], entry[1]
        options = dict(entry[2]) if len(entry) == 3 and entry[2] else {}
        return name, smarts, options

    @staticmethod
    def match_options(match) -> dict:
        """Options dict from a :meth:`match_queries` hit (empty if legacy 2-tuple)."""
        if isinstance(match, tuple) and len(match) >= 3:
            return match[2] or {}
        return {}

    @staticmethod
    def match_modifications(match):
        return match[1]

    def pair_query_options_allowed(self, match1, match2, active_pathways=()) -> bool:
        """True if a pair of query hits respects pathway / one-sided options.

        Used by resonance-pair rules so optional SMARTS stay configurable via
        the query list rather than hard-coded chemotype branches.
        """
        active = frozenset(active_pathways or ())
        opts1 = self.match_options(match1)
        opts2 = self.match_options(match2)
        pathway_counts = {}
        for opts in (opts1, opts2):
            pw = opts.get("pathway")
            if not pw:
                continue
            if pw not in active:
                return False
            pathway_counts[pw] = pathway_counts.get(pw, 0) + 1
        for opts in (opts1, opts2):
            pw = opts.get("pathway")
            if not pw:
                continue
            limit = opts.get("pair_limit")
            if limit is None and opts.get("one_sided"):
                limit = 1
            if limit is not None and pathway_counts.get(pw, 0) > int(limit):
                return False
        return True

    def pathways_recommended_by_queries(self, mol, target) -> frozenset:
        """Pathway names whose ``recommend`` SMARTS match ``target`` and formula hint allows.

        Intended for guided search: uncommon SMARTS stay off in broad enumeration
        but opt in when the target (+ formula) asks for them.
        """
        from .path_context import (
            FormulaHint,
            NEUTRAL,
            formula_compatible,
            heavy_formula_equal,
        )

        out = set()
        if target is None:
            return frozenset()
        for entry in self.query_smarts:
            _name, _smarts, options = self._parse_query_smarts_entry(entry)
            pw = options.get("pathway")
            if not pw:
                continue
            hint = options.get("formula_hint", NEUTRAL)
            if hint is None:
                hint = NEUTRAL
            # Pre-match: possibility set must be compatible (sound).
            if not formula_compatible(hint, mol, target, match=None):
                continue
            rec = options.get("recommend")
            if not rec:
                continue
            try:
                q = MolFromSmarts(rec) if isinstance(rec, str) else rec
                if q is not None and target.HasSubstructMatch(q):
                    out.add(pw)
            except Exception:
                continue
        return frozenset(out)

    def neighbors(self, mol, idx):
        """Get the indexes neighboring idx in mol."""
        return set([n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()])

    def _bfs_atom_path(self, mol, start, end, alternate_bonds=2):
        """Returns a list of all paths in mol between start and end atom.

        >>> mol = MolFromSmiles('c1ccccc1')
        >>> Resonate()._bfs_atom_path(mol,0,1,alternate_bonds=None)
        [[0, 1], [0, 5, 4, 3, 2, 1]]

        Double- vs single-bond alternating paths are a Kekulé pairing. Which
        edge is double depends on RDKit; together they are always both routes.

        >>> doubles = Resonate()._bfs_atom_path(MolFromSmiles('c1ccccc1'), 0, 1)
        >>> singles = Resonate()._bfs_atom_path(MolFromSmiles('c1ccccc1'), 0, 1, alternate_bonds=1)
        >>> sorted(doubles + singles, key=len)
        [[0, 1], [0, 5, 4, 3, 2, 1]]
        >>> len(doubles) == len(singles) == 1
        True

        """

        # This will be modified by the search function defined below.
        valid_paths = []

        def search(mol, end, path):
            """Searches for end idx by extending path."""

            if end == path[-1]:
                # valid_paths is modified in place, as it was initialized
                # outside of this search function.
                valid_paths.append(path)

            else:
                current_endpoint = mol.GetAtomWithIdx(path[-1])
                neighbors = [n.GetIdx() for n in current_endpoint.GetNeighbors()]

                if alternate_bonds:

                    if len(path) % 2:
                        valid_neighbors = [
                            n
                            for n in neighbors
                            if int(
                                mol.GetBondBetweenAtoms(
                                    path[-1], n
                                ).GetBondTypeAsDouble()
                            )
                            == alternate_bonds
                        ]

                    else:
                        valid_neighbors = [
                            n
                            for n in neighbors
                            if int(
                                mol.GetBondBetweenAtoms(
                                    path[-1], n
                                ).GetBondTypeAsDouble()
                            )
                            != alternate_bonds
                        ]

                else:
                    valid_neighbors = neighbors

                valid = set(valid_neighbors) - set(path)

                for idx in valid:
                    search(mol, end, path=path + [idx])

        if alternate_bonds:
            SanitizeMol(mol, SanitizeFlags.SANITIZE_CLEANUP)
            with edit_mol(mol):
                try:
                    Kekulize(mol, clearAromaticFlags=True)
                except ValueError:
                    pass

        # Execute recursive search that modifies valid_paths
        search(mol, end, path=[start])

        return sorted(valid_paths, key=len)

    def bfs_all_pairs(self, mol, **kwargs):
        """
        >>> mol = MolFromSmiles('c1ccccc1')
        >>> all_paths = QueryMol().bfs_all_pairs(mol,alternate_bonds=None)
        >>> all_paths[frozenset([2, 4])]
        [[2, 3, 4], [2, 1, 0, 5, 4]]
        """

        cache = _resonance_cache(mol)
        # Only memoize the default alternating-bond BFS used by resonance rules.
        if cache is not None and not kwargs:
            if cache.bfs_paths is None:
                cache.bfs_compute_count += 1
                cache.bfs_paths = self._bfs_all_pairs_uncached(mol, **kwargs)
            return cache.bfs_paths

        return self._bfs_all_pairs_uncached(mol, **kwargs)

    def _bfs_all_pairs_uncached(self, mol, **kwargs):
        AtomTracker.add_current_idx_as_atom_prop(mol)

        all_paths = {}

        current2originalidx = {
            a.GetIdx(): int(a.GetProp("idx"))
            for a in mol.GetAtoms()
            if a.GetAtomicNum() != 0
        }

        for pair in itertools.combinations(current2originalidx, 2):

            paths = self._bfs_atom_path(mol, *pair, **kwargs)

            original_idxs_path = [[current2originalidx[a] for a in p] for p in paths]

            all_paths[frozenset([current2originalidx[x] for x in pair])] = (
                original_idxs_path
            )

        return all_paths

    def rings(self, mol):
        """Construct mapping from each atom index to a list of tuples specifiying the atom indexes
        of all rings containing that atom.

        >>> mol = MolFromSmiles('c1ccc2ccccc2c1Cc1ccccc1')
        >>> rings = QueryMol().rings(mol)
        >>> rings[3]
        [(0, 9, 8, 3, 2, 1), (4, 5, 6, 7, 8, 3)]
        >>> rings[10]
        []

        """

        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        ring_info = mol.GetRingInfo()
        return {
            a.GetIdx(): [r for r in ring_info.AtomRings() if a.GetIdx() in r]
            for a in mol.GetAtoms()
        }

    def in_ring_size(self, rings, idx, size):
        """Returns True if idx in any rings of the specified size, otherwise returns False."""
        return size in [
            len(y)
            for y in itertools.chain(*[x for x in list(rings.values())])
            if idx in y
        ]

    def _queryidx2mapid(self, query):
        """Returns dict mapping from atom indexes to map IDs."""
        return {
            a.GetIdx(): int(a.GetProp("molAtomMapNumber"))
            for a in query.GetAtoms()
            if a.HasProp("molAtomMapNumber")
        }

    def match(self, mol, query=None, smarts=None, queryidx2mapid=None):
        """Returns dict mapping from mapids in query to atom indexes in mol.

        >>> mol = MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
        >>> next(EditMol().match(mol,smarts='[#6;R:1][#7H1:2]'))
        {1: 4, 2: 3}

        """
        if not query:
            if smarts:
                query, queryidx2mapid = self._prepare_query(smarts)
            else:
                raise NotImplementedError(
                    "Must submit either query mol or smarts pattern"
                )

        if not queryidx2mapid:
            queryidx2mapid = self._queryidx2mapid(query)

        # The ordering of the indices in each match corresponds to the atom ordering in the query.
        for match in mol.GetSubstructMatches(query):
            matchidx2queryidx = {y: x for x, y in enumerate(match)}
            mapid2queryidx = {
                queryidx2mapid[q]: i
                for i, q in list(matchidx2queryidx.items())
                if q in queryidx2mapid
            }

            yield mapid2queryidx

    @staticmethod
    def _process_smarts_match_dict(match_dict):
        out = {}

        for k, v in list(match_dict.items()):
            vs = [x for x in v if x]
            if vs:
                if isinstance(vs[0][0], (list, tuple)):
                    out[k] = frozenset(
                        itertools.chain(*[itertools.chain(*z) for z in vs])
                    )
                else:
                    out[k] = frozenset(itertools.chain(*vs))

        return out

    @classmethod
    def smarts_match(cls, mol, smarts, index=0):
        """
        >>> QueryMol.smarts_match('C=CC','[#6:1]~[#6:2]')
        ((0, 1), (1, 2))

        >>> QueryMol.smarts_match('C=CC','[#6:1]=[#6:2]')
        ((0, 1),)

        >>> QueryMol.smarts_match('C=CC',['[#6]=[#6]','[#6]-[#6]'])
        [((0, 1),), ((1, 2),)]

        >>> QueryMol.smarts_match('C=CC',{'S1':'[#6]=[#6]','S2':'[#6]-[#6]'})
        {'S1': frozenset({0, 1}), 'S2': frozenset({1, 2})}

        >>> QueryMol.smarts_match('C=CC',{'S1':['C~C~C','[#6]=[#6]'],'S2':'[#6]-[#6]'})
        {'S1': frozenset({0, 1, 2}), 'S2': frozenset({1, 2})}
        
        """

        if isinstance(mol, str):
            mol = MolFromSmiles(mol)

        if not mol.HasProp("standardized"):
            mol = cls.standardize(mol)
            if not mol:
                return None

        if isinstance(smarts, list):
            return [cls.smarts_match(mol, x, index=index) for x in smarts]

        if isinstance(smarts, dict):
            return cls._process_smarts_match_dict(
                {
                    k: cls.smarts_match(mol, v, index=index)
                    for k, v in list(smarts.items())
                }
            )

        query = MolFromSmarts(smarts)

        try:
            matches = mol.GetSubstructMatches(query)
        except:
            raise ValueError("Problem preparing query from %s" % smarts)

        return tuple([tuple([x + index for x in y]) for y in matches])

    def _prepare_query(self, smarts):
        """Preprocess smarts query.

        >>> query_mol, queryidx2mapid = EditMol()._prepare_query('[#6;R:1][#7H1:2]')
        >>> type(query_mol)
        <class 'rdkit.Chem.rdchem.Mol'>
        >>> queryidx2mapid
        {0: 1, 1: 2}

        """

        query = MolFromSmarts(smarts)

        try:
            assert query
        except AssertionError:
            raise ValueError("Problem with this pattern: %s" % smarts)

        return query, self._queryidx2mapid(query)

    @staticmethod
    def standardize(inmol):
        """Copies inmol to avoid modifying the input, and makes sure the structure is sanitized.

        >>> bad_input_mol = MolFromSmiles('C#C(C)C',sanitize=False)
        >>> good_input_mol = MolFromSmiles('C#C(C)',sanitize=False)
        >>> EditMol().standardize(bad_input_mol)
        False
        >>> mol = EditMol().standardize(good_input_mol)
        >>> isinstance(mol,Mol)
        True

        """

        # Reuse one kekulized prototype per substrate. Return a Chem.Mol copy
        # so callers cannot corrupt the prototype. Do not share resonance with
        # the aromatic parent or across edited working copies.
        forest = _forest_state(inmol)
        cached = forest.get("standardized_mol")
        if cached is not None:
            out = Chem.Mol(cached)
            carry_forest(cached, out)
            return out

        mol = Mol(inmol)
        carry_forest(inmol, mol)
        try:
            with edit_mol(mol):
                SanitizeMol(mol)
                Kekulize(mol, clearAromaticFlags=True)
        except ValueError:
            return False

        _clear_resonance(mol)
        # Ensure RingInfo is valid after kekulize (RDKit 2026).
        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        refresh_mol(mol)
        mol.SetProp("standardized", "1")
        forest["standardized_mol"] = mol
        out = Chem.Mol(mol)
        carry_forest(mol, out)
        return out

    def match_queries(self, mol):
        """Returns dict mapping from each atom index to a list of matches to self.queries.

        Each element of the list is a tuple ``(mapids, modifications, options)``.
        ``options`` comes from an optional third field on ``query_smarts`` entries.

        >>> query_smarts=[('name1', '[#6R:1][#7,#8h1:2]'),('name2','[#6R:1][#7:2][#6:3]')]
        >>> RRR = EditMol(query_smarts=query_smarts)
        >>> mol = MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
        >>> RRR.match_queries(mol)
        {4: [({1: 4, 2: 3}, 'name1', {}), ({1: 4, 2: 3, 3: 1}, 'name2', {})], 7: [({1: 7, 2: 8}, 'name1', {})]}



        """
        if isinstance(mol, str):
            mol = MolFromSmiles(mol)

        if not mol.HasProp("standardized"):
            mol = self.standardize(mol)
            if not mol:
                return dict()

        mapped_matches = collections.defaultdict(list)

        for modifications, query, queryidx2mapid, options in self.queries:
            for map2queryidx in self.match(
                mol, query=query, queryidx2mapid=queryidx2mapid
            ):

                mapped_matches[map2queryidx[min(map2queryidx)]].append(
                    (map2queryidx, modifications, options)
                )

        return dict(mapped_matches)


class EditMol(QueryMol):
    """Collection of standarized molecule modification functions."""

    bonds = {
        1: BondType.SINGLE,
        1.5: BondType.AROMATIC,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
    }

    def break_bond(self, emol, idx1, idx2):
        """Breaks bond between idx1 and idx2 in emol."""
        with edit_mol(emol):
            bond = emol.GetBondBetweenAtoms(idx1, idx2)
            if bond:
                emol.RemoveBond(idx1, idx2)
                self.adjust_hydrogen_count(emol, idx1, 1)
                self.adjust_hydrogen_count(emol, idx2, 1)

    def change_bond(self, emol, idx1, idx2, new_bond_type=1):
        """Change bond between idx1 and idx2 in emol to new_bond_type."""
        with edit_mol(emol):
            bond = emol.GetBondBetweenAtoms(idx1, idx2)
            if bond:
                bond.SetBondType(self.bonds[new_bond_type])

    def add_atom(self, emol, idx, new_atomic_num=1, new_bond_type=1):
        """Add atom to emol with new_atomic_num to atom with idx with new_bond_type."""
        with edit_mol(emol):
            i = emol.AddAtom(Atom(new_atomic_num))
            emol.AddBond(idx, i, self.bonds[new_bond_type])

    def replace_atom(self, emol, idx, new_atomic_num=1):
        """Changes atom with idx in emol to atom with new_atomic_num."""
        with edit_mol(emol):
            emol.GetAtomWithIdx(idx).SetAtomicNum(new_atomic_num)

    def set_charge(self, emol, idx, charge):
        """Sets idx charge in emol."""
        with edit_mol(emol):
            emol.GetAtomWithIdx(idx).SetFormalCharge(charge)

    def adjust_hydrogen_count(self, mol, atom, change):
        """Adjusts the number of explicit hydrogens on atom."""

        if isinstance(atom, int):
            atom = mol.GetAtomWithIdx(atom)

        try:
            implicit = atom.GetNumImplicitHs()
        except RuntimeError:
            atom.UpdatePropertyCache(strict=False)
            implicit = atom.GetNumImplicitHs()
        total_hydrogens = atom.GetNumExplicitHs() + implicit
        if total_hydrogens > 0:
            with edit_mol(mol):
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(total_hydrogens + change)

    def _bonds_from_atom_path(self, mol, atoms):
        """Converts list of atom indexes to a list of the bond objects between each atom pair.
        >>> bonds = EditMol()._bonds_from_atom_path(MolFromSmiles('C=CC=CCCl'),[0,1,2,3,4])
        >>> [(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in bonds]
        [(0, 1), (1, 2), (2, 3), (3, 4)]
        """

        for pos in range(len(atoms) - 1):
            try:
                bond = mol.GetBondBetweenAtoms(*atoms[pos : pos + 2])
            except RuntimeError:
                return
            if bond:
                yield bond

    def _correct_hydrogens_of_endpoint(self, mol, atomidx, newbond):
        """If newbond is single, add one hydrogen to atomidx.
        If newbond is double, remove one hydrogen from atomidx."""
        if newbond.GetBondType() == BondType.SINGLE:
            self.adjust_hydrogen_count(mol, atomidx, 1)
        elif newbond.GetBondType() == BondType.DOUBLE:
            self.adjust_hydrogen_count(mol, atomidx, -1)

    def swap_bonds_along_path(self, mol, atoms, attach_path=False, **kwargs):
        """Exchange double and single bonds in mol along the bond path specified by atoms.

        >>> mol = MolFromSmiles('C=CC=CCCl')
        >>> EditMol().swap_bonds_along_path(mol,atoms=[0,1,2,3,4])
        >>> canon_smi(mol) == canon_smi('CC=CC=CCl')
        True

        >>> mol = MolFromSmiles('C=CC=CCCl')
        >>> bonds = list(mol.GetBonds())[:-1]
        >>> EditMol().swap_bonds_along_path(mol,[0,1,2,3,4])
        >>> canon_smi(mol) == canon_smi('CC=CC=CCl')
        True

        """
        with edit_mol(mol):
            bond = False
            for i, bond in enumerate(self._bonds_from_atom_path(mol, atoms)):

                if bond.GetBondType() == BondType.DOUBLE:
                    bond.SetBondType(BondType.SINGLE)

                elif bond.GetBondType() == BondType.SINGLE:
                    bond.SetBondType(BondType.DOUBLE)

                if i == 0:  # fix hydrogens at begininng
                    self._correct_hydrogens_of_endpoint(mol, atoms[0], bond)

            # fix hydrogens at end
            if bond:
                self._correct_hydrogens_of_endpoint(mol, atoms[-1], bond)

            if attach_path:
                mol.SetProp("path", str(atoms))

        # Bond edits invalidate any cached resonance forms. edit_mol only
        # sees the mutation when the guard is installed, which it is not yet.
        _clear_resonance(mol)

    def apply_modifications(self, mol, modifications, **kwargs):
        """Makes an editable copy of mol and successively applies each submitted modification."""

        emol = RWMol(Mol(mol))
        with edit_mol(emol):
            for item in modifications:
                map2queryidx, modification = item[0], item[1]
                valid = self.modify(emol, map2queryidx, modification, **kwargs)
                if not valid:
                    return False
        return emol.GetMol()

    def modify(self, emol, mapid2atomidx, modifications):
        """Should return True if all modifications were successful, otherwise return False."""
        raise NotImplementedError


class Resonate(ConjugatedSystems, EditMol):
    """Generates effecient resonance structures by resonating each atom system independently."""

    def __call__(self, mol, first_yield_input=True):
        """Generates effecient resonance structures by resonating each atom system independently.

        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> smiles = [MolToSmiles(x, kekuleSmiles=True) for x in Resonate()(mol)]
        >>> all(MolFromSmiles(s) for s in smiles)
        True
        >>> len(smiles) > 1
        True

        """
        return self.resonance_structures(mol, first_yield_input=first_yield_input)

    def resonance_structures(self, mol, first_yield_input=True):
        """Generates effecient resonance structures by resonating each atom system independently.

        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> smiles = [MolToSmiles(x, kekuleSmiles=True) for x in Resonate().resonance_structures(mol)]
        >>> all(MolFromSmiles(s) for s in smiles)
        True

        """

        if first_yield_input:
            yield mol

        AtomTracker.add_current_idx_as_atom_prop(mol)

        cache = _resonance_cache(mol)
        if cache is not None:
            for joined, _system in cache.mode(self.flag).iter_joined(self, mol):
                yield joined
            return

        for res_frag, rest_of_molecule, bond_types in self._resfrags(mol):

            yield self.join_fragments([res_frag] + rest_of_molecule, bond_types)

    def resonate_with_pair_paths(self, mol, valid_atoms=None):
        """
        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> res_struct, pair, path = next(Resonate().resonate_with_pair_paths(mol))
        >>> canon_smi(res_struct) == canon_smi('NCCCCC1=C2C=C(CC3=CC=CC(C=CC4=CC=CC=C4)=C3)C=CC2=CC=C1')
        True
        >>> pair, path[0], path[-1]
        ((5, 6), 5, 6)
        """

        all_system_paths = self.bfs_all_pairs(mol)

        cache = _resonance_cache(mol)
        if cache is not None:
            joined_iter = cache.mode(self.flag).iter_joined(self, mol)
            for res_struct, system in joined_iter:
                for pair in itertools.combinations(system, 2):
                    if valid_atoms and not set(pair).issubset(valid_atoms):
                        continue
                    for path in all_system_paths[frozenset(pair)]:
                        yield copy_mol(res_struct), pair, path
            return

        for system_fragment, the_rest, bonds, system in self._resfrags(
            mol, output_systems=True
        ):

            res_struct = self.join_fragments([system_fragment] + the_rest, bonds)

            for pair in itertools.combinations(system, 2):

                if valid_atoms and not set(pair).issubset(valid_atoms):
                    continue

                for path in all_system_paths[frozenset(pair)]:
                    yield copy_mol(res_struct), pair, path

    def _resfrags(self, mol, output_systems=False):
        """
        >>> RES = Resonate()
        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> res_frag, other_fragments, bonds, system = next(RES._resfrags(mol,output_systems=True))
        >>> system
        {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}
        >>> canon_smi(res_frag) == canon_smi('[4*]C1=C2C=C([15*])C=CC2=CC=C1')
        True
        >>> canon_smi(other_fragments) == canon_smi(['[5*]CCCCN', '[12*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1'])
        True
        """

        for fragment_to_resonate, other_fragments, system, bond_types in self.fragments(
            mol
        ):

            # Sanitization is required to make the ResonanceMolSupplier work correctly.
            SanitizeMol(fragment_to_resonate, catchErrors=True)

            try:
                supplier = ResonanceMolSupplier(fragment_to_resonate, KEKULE_ALL)

            except ValueError:
                try:
                    supplier = ResonanceMolSupplier(mol)
                except ValueError:
                    supplier = [fragment_to_resonate]

            for resmol in supplier:

                # I tried this but it hampered speed and accuracy:
                # Kekulize(resmol, clearAromaticFlags=True)

                if not resmol:
                    continue

                outputs = [resmol, other_fragments, bond_types]

                if output_systems:
                    outputs.append(system)

                yield tuple(outputs)


def _site_atom_set(site):
    """Normalize a metabolize site to a frozenset of atom indices."""
    if isinstance(site, tuple) and len(site) >= 2 and not isinstance(site[0], int):
        return frozenset(site[1])
    if isinstance(site, (list, tuple, set, frozenset)):
        return frozenset(site)
    return frozenset([site])


def _site_matches_filters(site, include_sites, exclude_sites) -> bool:
    atoms = _site_atom_set(site)
    if exclude_sites is not None:
        excluded = {_site_atom_set(s) for s in exclude_sites}
        if atoms in excluded:
            return False
    if include_sites is not None:
        included = {_site_atom_set(s) for s in include_sites}
        if atoms not in included:
            return False
    return True


def _include_atom_sets(include_sites):
    """frozenset of atom-index frozensets, or ``None`` if unrestricted."""
    if include_sites is None:
        return None
    return frozenset(_site_atom_set(s) for s in include_sites)


class ReactionRule(AtomTracker):
    """Template for all rules."""

    sites_on = "atoms"
    phase1_sites_on = "bonds"
    # When True, phase1_steps returns a degenerate singleton plan for valid sites.
    phase1_equivalent = False

    def __init__(
        self,
        name=None,
        sites_on=None,
        phase1_sites_on=None,
        longname=None,
        *args,
        **kwargs
    ):
        """sites_on is used by report.
        phase1_sites_on is used by the --phase1 option. If no phase1_sites_on, default to sites_on.
        """
        super(ReactionRule, self).__init__(*args, **kwargs)

        if not name:
            self.name = str(self).split(".")[-1].split()[0]
        else:
            self.name = name

        if longname is None:
            self.longname = self.name
        else:
            self.longname = longname

        try:
            assert "_" not in self.name
        except AssertionError as err:
            err.args = ("Cannot have '_' in rule name",)
            raise

        if sites_on is not None:
            self.sites_on = sites_on

        if phase1_sites_on is not None:
            self.phase1_sites_on = phase1_sites_on

    def __call__(self, mol, **kwargs):
        for site, metabolites in self.metabolize(mol, **kwargs):
            yield site, metabolites

    def __iter__(self):
        return iter([self])

    def phase1_steps(self, mol, site, **kwargs):
        """Return one Phase1-equivalent :class:`~xenosite.forest.step_plan.StepPlan`.

        Public signature is ``(mol, site)``. Subclasses that are Phase I set
        ``phase1_equivalent = True`` for a degenerate singleton plan. Others
        raise ``NotImplementedError``. Private underscore kwargs may be used
        by subclasses for efficiency.
        """
        from .step_plan import StepPlan

        if not self.phase1_equivalent:
            raise NotImplementedError(
                "%s does not define Phase1-equivalent steps" % self.name
            )
        from .step_plan import StepPlan, _origin_refs

        want = self._cast_sites(site)[0][1]
        for outsite, _products in self.metabolize(
            mol,
            tag_atoms=False,
            only_emit_topologically_distinct_sites=False,
            attach_phase1_steps=False,
        ):
            emitted = outsite[1] if isinstance(outsite, tuple) else outsite
            if frozenset(emitted) == frozenset(want):
                # Site idxs are GetIdx on ``mol`` — stamp that frame's depth.
                return StepPlan.singleton(self.name, _origin_refs(want, mol))
        return StepPlan.empty()

    def format_site(self, site, just_rule_name=False):

        if just_rule_name and isinstance(site, (tuple, list)):
            return site[0]

        if isinstance(site, tuple):
            return tuple(
                [self.format_site(x, just_rule_name=just_rule_name) for x in site]
            )

        if isinstance(site, list):
            return [self.format_site(x, just_rule_name=just_rule_name) for x in site]

        if isinstance(site, str):
            return site.split("_")[0]

        return site

    def metabolize(
        self,
        mol,
        only_emit_topologically_distinct_sites=True,
        tag_atoms=True,
        format_output_site=True,
        only_largest_fragment=False,
        do_not_tag_atoms=False,
        only_unique=False,
        strict=True,
        attach_phase1_steps=False,
        include_sites=None,
        exclude_sites=None,
        **kwargs
    ):
        if only_unique:
            only_emit_topologically_distinct_sites = False
            unique = False
            unique_smi = []

        # Maps must not be present for CanonicalRankAtoms or SMARTS matching.
        self._clear_atom_maps(mol)
        topol_equiv = self.topol_equiv(mol)
        if only_emit_topologically_distinct_sites:
            try:
                topol_equiv = self.topol_equiv(mol)
            except:
                # SanitizeMol(mol, SanitizeFlags.SANITIZE_CLEANUP, catchErrors=True)
                topol_equiv = self.topol_equiv(mol)

        tagging = tag_atoms and not do_not_tag_atoms
        if tagging:
            self.initialize_tags(mol)
            AtomTracker.add_current_idx_as_atom_prop(
                mol, propname=AtomTracker.previous_index_prop_name
            )

        seen = []
        skipped = []
        products_iter = self.metabolites(
            mol,
            format_output_site=format_output_site,
            do_not_tag_atoms=do_not_tag_atoms,
            strict=strict,
            attach_phase1_steps=attach_phase1_steps,
            include_sites=include_sites,
            exclude_sites=exclude_sites,
            **kwargs
        )
        while True:
            if tagging:
                self._clear_atom_maps(mol)
            try:
                site, metabolites = next(products_iter)
            except StopIteration:
                if tagging:
                    self._stamp_origin_maps(mol)
                break

            metabolites = [x for x in metabolites if x]
            if not metabolites:
                continue

            valid = []
            for metabolite in metabolites:
                if is_rdkit_valid(metabolite):
                    valid.append(metabolite)
                    continue
                note_sanitize_drop(1)
                _log.debug(
                    "Dropping RDKit-invalid %s metabolite %s",
                    self.name,
                    _mol_smiles(metabolite),
                )
            if not valid:
                continue
            metabolites = valid

            if only_largest_fragment:
                metabolites = [
                    sorted(
                        metabolites, key=lambda x: x.GetNumAtoms(), reverse=True
                    )[0]
                ]

            if only_emit_topologically_distinct_sites:
                # Match XenoSite UI: pathway + topo ranks + product SMILES.
                # Do not collapse distinct products that share a topo site class.
                topsite = self.site_to_topol_site(site, topol_equiv)
                identity = (topsite[0], topsite[1], can_smi_set(metabolites))
                if identity in seen:
                    continue
                seen.append(identity)

            if tagging:
                aligned = []
                # Formation site atoms are reactant-frame idxs when unformatted.
                origin_site = site[1] if isinstance(site, tuple) and len(site) == 2 else site
                if not isinstance(origin_site, frozenset):
                    try:
                        origin_site = frozenset(origin_site)
                    except TypeError:
                        origin_site = frozenset()
                for metabolite in metabolites:
                    self.tag(
                        metabolite,
                        reactant=mol,
                        strict=strict,
                        # Single product → chemical deletion goes to removed[].
                        # Multi-fragment cleavage → siblings are on other
                        # products; drop absences (do not record).
                        record_removals=len(metabolites) == 1,
                    )
                    metabolite = self._align_and_stamp(metabolite)
                    install_product_forest(
                        mol,
                        metabolite,
                        rule_name=self.name,
                        origin_site=origin_site,
                        record_creation=True,
                    )
                    aligned.append(metabolite)
                metabolites = aligned
                self._stamp_origin_maps(mol)

            if format_output_site:
                outsite = self.format_site(site)
            else:
                outsite = site

            if only_unique:
                unique_metabolites = can_smi_set(metabolites)

                if unique_metabolites in unique_smi:
                    continue
                else:
                    unique_smi.append(unique_metabolites)

            if include_sites is not None or exclude_sites is not None:
                if not _site_matches_filters(outsite, include_sites, exclude_sites):
                    continue

            if attach_phase1_steps:
                self._attach_phase1_steps_to_products(mol, outsite, metabolites)

            yield outsite, metabolites

    def _attach_phase1_steps_to_products(self, mol, outsite, metabolites):
        """Stamp phase1_steps on products that do not already carry a plan."""
        from .step_plan import StepPlan

        missing = [m for m in metabolites if StepPlan.try_from_mol(m) is None]
        if not missing:
            return
        # Subclasses (e.g. QuinoneFormation) may stamp only selected fragments.
        if any(StepPlan.try_from_mol(m) is not None for m in metabolites):
            return
        if not self.phase1_equivalent:
            raise NotImplementedError(
                "attach_phase1_steps is not supported for %s" % self.name
            )
        from .step_plan import _origin_refs

        site = outsite[1] if isinstance(outsite, tuple) else outsite
        # Site idxs are GetIdx on ``mol`` (may be mid-depth; created atoms
        # often have no depth-0 entry).
        plan = StepPlan.singleton(self.name, _origin_refs(site, mol))
        for metabolite in missing:
            plan.attach_to_mol(metabolite)

    def is_redundant(self, peer_rules) -> bool:
        """True if this rule should be dropped given the rest of the ruleset."""
        return False

    def is_terminal_product(self, mol) -> bool:
        """True if ``mol`` must not be expanded further in guided path search."""
        return False

    def formula_hints(self):
        """Per-reaction / per-query :class:`~xenosite.forest.path_context.FormulaHint`\\ s.

        Prefer ``formula_hint`` on each SMARTS entry's options dict (reaction
        ``smarts`` / ``rxns`` or ``query_smarts``). Legacy class attribute
        ``formula_effects`` is still read when options omit hints.
        ``None`` entries mean “unknown — do not filter”. Empty / missing ⇒ no
        formula-based SMARTS skipping (``could_help`` stays permissive unless
        overridden).
        """
        opts = getattr(self, "smarts_options", None)
        if opts and any("formula_hint" in o for o in opts):
            return [o.get("formula_hint") for o in opts]
        # query_smarts options (ResonancePairRule / EditMol)
        qsmarts = getattr(self, "query_smarts", None) or []
        if qsmarts:
            from_query = []
            for entry in qsmarts:
                try:
                    _n, _s, o = QueryMol._parse_query_smarts_entry(entry)
                except Exception:
                    from_query.append(None)
                    continue
                from_query.append(o.get("formula_hint") if o else None)
            if from_query and any(h is not None for h in from_query):
                return from_query
        return getattr(self, "formula_effects", None) or getattr(
            type(self), "formula_effects", None
        )

    def smarts_compatible(self, rxn_index: int, mol, target, match=None) -> bool:
        """False ⇒ skip reaction SMARTS ``rxn_index`` for paths toward ``target``.

        ``match`` optional: when provided, match-polymorphic ``formula_hint``
        resolvers can refine the effect; otherwise the pre-match possibility
        set is used (sound skip).
        """
        from .path_context import formula_compatible

        hints = self.formula_hints()
        if not hints:
            return True
        if rxn_index < 0 or rxn_index >= len(hints):
            return True
        hint = hints[rxn_index]
        if hint is None:
            return True
        return formula_compatible(hint, mol, target, match=match)

    def elements_may_add(self):
        """Element symbols this rule may introduce, or ``None`` if unknown.

        Default: derive from :meth:`formula_hints` positive ``delta`` entries.
        Rules without hints return ``None`` (cannot prove impossibility).
        """
        from .path_context import elements_may_add_from_hints

        return elements_may_add_from_hints(self.formula_hints())

    def could_help(self, mol, target, ctx) -> bool:
        """False ⇒ this rule cannot appear on any path from ``mol`` to ``target``.

        Default: if :meth:`formula_hints` are declared, True when any hint is
        compatible with ``mol``→``target``; otherwise True (permissive).
        """
        from .path_context import any_hint_compatible

        hints = self.formula_hints()
        if not hints:
            return True
        return any_hint_compatible(hints, mol, target)

    def is_cleavage(self) -> bool:
        """True if any formula hint may shrink / fragment the molecule."""
        from .path_context import hint_includes_cleave

        hints = self.formula_hints()
        if not hints:
            return False
        return any(h is not None and hint_includes_cleave(h) for h in hints)

    def cleave_alone(self) -> bool:
        """True ⇒ join the cleavage peer pass when the target is smaller.

        The guided cleavage-only pass expands every ``is_cleavage()`` rule as
        peers (primary and secondary CLEAVE-tagged chemistry together). Override
        this when a rule should join that pass even without a CLEAVE formula
        hint (e.g. oxidative dehalogenation: C–X cleavage with ADD_O hints).
        """
        return False

    def sites_toward(self, mol, target, ctx):
        """Required sites for a direct path, or ``None`` for unrestricted."""
        return None

    def sites_maybe(self, mol, target, ctx):
        """Optional (Maybe) sites that do not themselves advance MCS/delta."""
        return ()

    def child_may_reach(self, parent, child, target, ctx) -> bool:
        """False ⇒ reject ``child`` before enqueue (Impossible product).

        Default: formula distance to ``target`` must not increase.
        """
        if ctx is None:
            return True
        try:
            return ctx.distance_proxy(child) <= ctx.distance_proxy(parent)
        except Exception:
            return True

    def enumerate_for_path(
        self,
        mol,
        ctx,
        expand_phase1_plans=True,
        include_sites=None,
        exclude_sites=None,
        only_emit_topologically_distinct_sites=True,
        **kwargs,
    ):
        """Yield ``("plan", expr, site)`` or ``("hop", site, products)`` for guided search.

        Batches metabolize once per rule/mol. Never loops per-site generation.
        Topologically equivalent site orbits are pruned by default.
        Passes ``toward_target`` so per-SMARTS formula hints can skip reactions.

        When :meth:`sites_toward` returns a list and phase1 plans are enabled,
        emit plans for those sites from reactant SMARTS + MCS pruning **without**
        calling ``metabolize`` / ``RunReactants`` first.
        """
        from .step_plan import pathway_alternatives

        if not self.could_help(mol, ctx.target, ctx):
            return

        toward = self.sites_toward(mol, ctx.target, ctx)
        maybe = self.sites_maybe(mol, ctx.target, ctx)
        if include_sites is None and toward is not None:
            include_sites = list(toward) + list(maybe or ())

        supports_phase1 = False
        if expand_phase1_plans:
            try:
                supports_phase1 = bool(getattr(self, "phase1_equivalent", False))
                if not supports_phase1:
                    supports_phase1 = (
                        type(self).phase1_steps is not ReactionRule.phase1_steps
                    )
            except Exception:
                supports_phase1 = False

        # Cleavage toward-sites: plan without materializing every SMARTS hit.
        if (
            include_sites is not None
            and supports_phase1
            and expand_phase1_plans
            and toward is not None
        ):
            for site in include_sites:
                atoms = _site_atom_set(site)
                try:
                    expr = self.phase1_steps(mol, atoms)
                except NotImplementedError:
                    continue
                if pathway_alternatives(expr):
                    yield ("plan", expr, site)
            return

        topol = None  # metabolize topo+product identity is the orbit prune
        for site, products in self.metabolize(
            mol,
            include_sites=include_sites,
            exclude_sites=exclude_sites,
            tag_atoms=False,
            only_emit_topologically_distinct_sites=only_emit_topologically_distinct_sites,
            attach_phase1_steps=False,
            toward_target=ctx.target,
            **kwargs,
        ):
            atoms = _site_atom_set(site)
            if supports_phase1:
                try:
                    expr = self.phase1_steps(mol, atoms)
                except NotImplementedError:
                    yield ("hop", site, products)
                    continue
                if pathway_alternatives(expr):
                    yield ("plan", expr, site)
                    continue
            yield ("hop", site, products)

    def metabolites(self, mol, **kwargs):
        """Should return a tuple of lists. The first element will be the site, the second element
        will be a list of metabolites. Any rule that could cause molecule fragmention should
        separate the fragments within the outputted list of metabolites. One way to do this is to
        return clean(product), which will separate all fragments and sanitize each one.
        """
        raise NotImplementedError

    def metabolites_from_sites(
        self,
        mol,
        sites,
        just_smiles=False,
        deplete_sites=False,
        tag_atoms=True,
        **kwargs
    ):
        """Generates all metabolites matching sites.
        If deplete_sites=True, then stop iteration once a metabolite for each site has been found.
        Otherwise, return all metabolites with matching sites."""

        sites = self._cast_sites(sites)
        # Restrict generation before clean/RunReactants when possible — full
        # metabolize + post-filter still builds every resonance product.
        include = [_site_atom_set(s) for s in sites]

        for site, metabolite in self.metabolize(
            mol, tag_atoms=tag_atoms, include_sites=include, **kwargs
        ):
            if tuple(site) in sites:

                if just_smiles:
                    yield list(map(unmapped_smiles, metabolite))

                else:
                    yield tuple(site), metabolite

                if deplete_sites:
                    sites.remove(tuple(site))
                    if not sites:
                        break

    def _cast_sites(self, sites):
        """Converts sites in various formats to a list of tuples, with the first element of
        tuple specifing the rule name, and the second element being a frozenset of atom indexes.
        """

        if isinstance(sites, int):
            return [(self.name, frozenset([sites]))]

        if isinstance(sites, tuple):
            return [(self.name, frozenset(sites))]

        elif isinstance(sites, frozenset):
            return [(self.name, sites)]

        elif isinstance(sites, list):

            if isinstance(sites[0], frozenset):
                return [(self.name, x) for x in sites]

            elif isinstance(sites[0], list):
                return [tuple(x) for x in sites]

            elif isinstance(sites[0], int):
                return [(self.name, frozenset(sites))]

            return sites

        else:
            raise NotImplementedError(
                "sites must be an integer, a frozenset, or a list"
            )


class SmartsReactionRule(ReactionRule):
    """Performs reactions specified by SMARTS.

    Each ``smarts`` / ``rxns`` entry may be a plain SMARTS string or
    ``(smarts, options)`` where ``options`` may include ``formula_hint``
    (and the same keys as ``query_smarts`` options: ``pathway``, …).
    """

    parameters = ["Reaction SMARTS"]
    smarts = []
    mapid_site = []

    @staticmethod
    def _parse_smarts_entry(entry):
        """Return ``(smarts_str, options_dict)`` for a reaction SMARTS list item."""
        if isinstance(entry, str):
            return entry, {}
        if isinstance(entry, (tuple, list)) and len(entry) >= 1:
            smarts = entry[0]
            options = dict(entry[1]) if len(entry) > 1 and entry[1] else {}
            if not isinstance(smarts, str):
                raise TypeError(
                    "reaction SMARTS entry must start with a string, got %r" % (entry,)
                )
            return smarts, options
        raise TypeError(
            "reaction SMARTS entry must be a string or (smarts, options), got %r"
            % (entry,)
        )

    def __init__(self, rxns=None, mapid_site=None, *args, **kwargs):

        super(SmartsReactionRule, self).__init__(*args, **kwargs)

        if mapid_site:
            self.mapid_site = mapid_site

        if rxns:
            if isinstance(rxns, str):
                rxns = [rxns]
            entries = list(rxns)
        elif self.smarts:
            if isinstance(self.smarts, str):
                entries = [self.smarts]
            else:
                entries = list(self.smarts)
        else:
            entries = []

        parsed = [self._parse_smarts_entry(e) for e in entries]
        self.smarts = [s for s, _o in parsed]
        self.smarts_options = [o for _s, o in parsed]
        self.rxns = [self._smarts2rxns(s, **kwargs) for s in self.smarts]

    def formula_hints(self):
        """Hints from each reaction SMARTS entry's ``formula_hint`` option.

        Falls back to legacy class ``formula_effects`` when options omit hints.
        """
        opts = getattr(self, "smarts_options", None) or []
        if opts and any("formula_hint" in o for o in opts):
            return [o.get("formula_hint") for o in opts]
        return super(SmartsReactionRule, self).formula_hints()

    def iter_reactant_site_matches(self, mol):
        """Yield formation sites from reactant-side SMARTS only (no ``RunReactants``)."""
        mapids = list(self.mapid_site or [])
        if not mapids:
            return
        seen = set()
        for smarts in self.smarts or ():
            lhs = smarts.split(">>", 1)[0].strip()
            # Drop reaction-product grouping parens sometimes left on LHS.
            if lhs.startswith("(") and lhs.endswith(")"):
                lhs = lhs[1:-1]
            query = Chem.MolFromSmarts(lhs)
            if query is None:
                continue
            for match in mol.GetSubstructMatches(query):
                by_map = {}
                for qi, atom in enumerate(query.GetAtoms()):
                    mid = int(atom.GetAtomMapNum() or 0)
                    if mid and qi < len(match):
                        by_map[mid] = int(match[qi])
                try:
                    site = frozenset(by_map[m] for m in mapids)
                except KeyError:
                    continue
                if len(site) < 1 or site in seen:
                    continue
                seen.add(site)
                yield site

    def sites_toward(self, mol, target, ctx):
        """When T is smaller, keep cleavage sites whose MCS side can hold T.

        Candidate bonds come from reactant SMARTS matches; pruning uses the
        reactant↔target MCS and a bridge split on the SMILES graph — cleavage
        products are not materialized.

        Sites are **prioritized** by MCS frontier / embedding disagreement.
        Only sites passing :func:`~xenosite.forest.path_context.cleavage_site_safe_to_drop_required`
        thresholds are hard-dropped from Required (deep interior; deep exterior
        when MCS already matches T's size). Other non-frontier sites stay, later
        in the list.
        """
        if not self.is_cleavage():
            return super(SmartsReactionRule, self).sites_toward(mol, target, ctx)
        if ctx is None or target is None:
            return None
        try:
            if target.GetNumHeavyAtoms() >= mol.GetNumHeavyAtoms():
                return None
        except Exception:
            return None
        from .path_context import (
            cleavage_site_may_reach,
            cleavage_site_priority,
            cleavage_site_safe_to_drop_required,
        )

        out = []
        seen = set()
        for site in self.iter_reactant_site_matches(mol):
            key = frozenset(site)
            if key in seen:
                continue
            seen.add(key)
            site_fs = frozenset(site)
            if not cleavage_site_may_reach(mol, site_fs, target, ctx):
                continue
            if cleavage_site_safe_to_drop_required(mol, site_fs, ctx, target):
                continue
            out.append(site_fs)
        out.sort(key=lambda s: cleavage_site_priority(mol, s, ctx))
        return out

    def metabolites(self, mol, kekulize=True, toward_target=None, **kwargs):
        """By default, mol will be kekulized.

        ``toward_target``: optional product mol; reaction SMARTS whose
        :meth:`smarts_compatible` is False are skipped (formula hints).

        Kekulize runs on a copy so the caller's bonding / resonance cache are
        not mutated or cleared.
        """
        if kekulize:
            mol = copy_mol(mol)
            self._kekulize(mol)

        self._remove_props(mol)
        self._clear_atom_maps(mol)
        refresh_mol(mol)

        for rxn_num, rxn in enumerate(self.rxns):
            if toward_target is not None and not self.smarts_compatible(
                rxn_num, mol, toward_target
            ):
                continue
            self._clear_atom_maps(mol)
            try:
                reactant_products = rxn.RunReactants((mol,))
            except RuntimeError:
                _log.debug(
                    "Skipping %s rxn %d on unsanitizable reactant %s",
                    self.name,
                    rxn_num,
                    _mol_smiles(mol),
                )
                continue

            for prod_num, prod in enumerate(reactant_products):

                products = list(prod)
                site = self._get_site(products, mapid_site=self.mapid_site)
                fullsite = self._get_site(products)

                # RunReactants replaces atoms, so the props needed to be copied back
                self._copy_props(
                    products, {"react_atom_idx": AtomTracker.previous_index_prop_name}
                )

                yield (self.name + "_SmartsReactionRuleRxn%d" % (rxn_num), site), clean(
                    products
                )

    def _smarts2rxns(self, smarts, use_implicit_properties=False, **kwargs):
        """Converts list of SMARTS reactions to RDKit reactions."""
        rxn = AllChem.ReactionFromSmarts(smarts)
        if not use_implicit_properties:
            rxn._setImplicitPropertiesFlag(False)
        return rxn

    def _kekulize(self, mol):
        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        with edit_mol(mol):
            try:
                Kekulize(mol, clearAromaticFlags=True)
            except ValueError:
                pass
        refresh_mol(mol)
        _clear_resonance(mol)

    def _copy_props(
        self,
        mol,
        props_to_copy=None,
    ):

        if props_to_copy is None:
            props_to_copy = {
                "react_atom_idx": AtomTracker.previous_index_prop_name
            }

        if isinstance(mol, list):
            return [self._copy_props(x, props_to_copy) for x in mol]

        for atom in mol.GetAtoms():
            for old, new in list(props_to_copy.items()):
                if atom.HasProp(old):
                    atom.SetProp(new, atom.GetProp(old))

    def _remove_props(self, mol, props_to_remove=["old_mapno", "react_atom_idx"]):

        if isinstance(mol, list):
            return [self._remove_props(x) for x in mol]

        for atom in mol.GetAtoms():
            for prop in props_to_remove:
                if atom.HasProp(prop):
                    atom.ClearProp(prop)

    def _get_site(self, product, mapid_site=None):
        """Return the sites in product based on the reactant_idx property assigned by
        rxns.RunReactants in self.metabolites."""

        if isinstance(product, list):
            return frozenset.union(
                *[self._get_site(p, mapid_site=mapid_site) for p in product]
            )

        mapno2idx = {
            int(a.GetProp("react_atom_idx")): int(a.GetProp("old_mapno"))
            for a in product.GetAtoms()
            if a.HasProp("old_mapno") and a.HasProp("react_atom_idx")
        }
        if mapid_site:
            return frozenset([r for r, m in mapno2idx.items() if m in mapid_site])
        return frozenset(mapno2idx)


class ResonanceRule(Resonate, SmartsReactionRule):
    """Base for rules that apply SMARTS reactions to resonance structures."""

    parameters = ["Reaction SMARTS", "Endpoint Modifications"]

    def metabolites(self, mol, **kwargs):
        for num, res_struct in enumerate(self.resonance_structures(mol)):
            for site, products in SmartsReactionRule.metabolites(self,
                res_struct, kekulize=False, **kwargs
            ):
                yield (site[0] + "_ResonanceRule%d " % num, site[1]), products


class ResonancePairRule(ResonanceRule):
    parameters = ["Reaction SMARTS", "Endpoint Modifications", "System Type"]
    """Base for rules that will apply EditMol modifications to pairs of atoms within
    resonance structures."""


if __name__ == "__main__":
    import doctest

    doctest.testmod()
