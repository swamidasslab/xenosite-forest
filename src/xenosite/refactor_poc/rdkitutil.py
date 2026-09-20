"""Molecule questions for the proof of concept.

RDKit itself is imported in :mod:`xenosite.refactor_poc.rdkit_api`. This
module calls those names.

Answers about a molecule are cached on ``get_forest(ensure_forest(mol))["structure"]``.
The caller rebinds ``mol = ensure_forest(mol)`` when the molecule had no forest.
``get_forest`` only reads a molecule that already has one.
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
import copy
from collections import defaultdict, deque
from collections.abc import Iterable
from typing import TypeGuard, TypeVar, overload

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
    ForestNoTracingMol,
    ForestTracingMol,
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
)
from xenosite.refactor_poc.records import (
    EditCounters,
    EndParents,
    Forest,
    Formula,
    FragmentSplit,
    KekuleParents,
    McsResult,
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


@overload
def is_forest(mol: ForestTracingMol) -> TypeGuard[ForestTracingMol]: ...
@overload
def is_forest(mol: ForestNoTracingMol) -> TypeGuard[ForestNoTracingMol]: ...
@overload
def is_forest(mol: ForestMol) -> TypeGuard[ForestMol]: ...
@overload
def is_forest(mol: Mol) -> TypeGuard[ForestMol]: ...
def is_forest(mol: Mol) -> bool:
    """True when ``_forest`` is present. Does not install one."""

    return _read_forest(mol) is not None


@overload
def is_tracing(mol: ForestTracingMol) -> TypeGuard[ForestTracingMol]: ...
@overload
def is_tracing(mol: Mol) -> TypeGuard[ForestTracingMol]: ...
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
def ensure_forest(mol: ForestTracingMol) -> ForestTracingMol: ...
@overload
def ensure_forest(mol: ForestNoTracingMol) -> ForestNoTracingMol: ...
@overload
def ensure_forest(mol: ForestMol) -> ForestMol: ...
@overload
def ensure_forest(mol: Mol) -> ForestMol: ...
def ensure_forest(mol: Mol) -> ForestMol:
    """Install a missing forest on ``mol`` and return that same object.

    An existing forest is left alone, including its depth. A traced molecule
    stays traced. This does not copy and does not raise. Callers rebind:
    ``mol = ensure_forest(mol)``.
    """

    if is_forest(mol):
        return mol
    structure: Structure = {}
    mol._forest = {"structure": structure}
    assert is_forest(mol)
    return mol


def get_forest(mol: ForestMol) -> Forest:
    """The forest already on ``mol``. Does not install and does not copy."""

    return mol._forest


def _place_forest(mol: Mol, forest: Forest) -> ForestMol:
    """Write ``forest`` onto ``mol`` and return that same object.

    The parameter is :class:`Mol`, not :class:`NoForestMol`, so the write is
    the base attribute. The molecule is not copied.
    """

    mol._forest = forest
    assert is_forest(mol)
    return mol


def _structure(mol: Mol) -> Structure:
    forest = get_forest(ensure_forest(mol))
    if "structure" not in forest:
        raise KeyError("structure")
    return forest["structure"]


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


def topol_equiv(mol: Mol) -> dict[int, int]:
    """Map each atom index to its topological class. The same dict on a hit."""

    structure = _structure(mol)
    if "topol_equiv" in structure:
        return structure["topol_equiv"]

    sanitized = Mol(mol)
    sanitize_mol(sanitized)

    classes = {
        atom.GetIdx(): rank
        for atom, rank in zip(
            mol.GetAtoms(),
            CanonicalRankAtoms(sanitized, includeChirality=False, breakTies=False),
        )
    }
    structure["topol_equiv"] = classes
    return classes


def molecule_formula(mol: Mol) -> Formula:
    """Heavy-atom counts, total hydrogens, and formal charge.

    Explicit hydrogens are counted through ``GetTotalNumHs`` on the heavy
    atom they belong to, not as a second copy. Cached as ``structure["formula"]``.
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
def copy_mol(mol: ForestTracingMol) -> ForestTracingMol: ...
@overload
def copy_mol(mol: ForestNoTracingMol) -> ForestNoTracingMol: ...
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
        return _place_forest(out, copy.deepcopy(get_forest(mol)))
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
    mol: ForestTracingMol, tracing_reset: bool = True
) -> tuple[ForestTracingMol, str]: ...
@overload
def cannonicalize_order(mol: Mol, tracing_reset: bool = True) -> tuple[ForestMol, str]: ...
def cannonicalize_order(mol: Mol, tracing_reset: bool = True) -> tuple[ForestMol, str]:
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
    renumbered = _place_forest(
        RenumberAtoms(mol, renumber_map), copy.deepcopy(get_forest(source))
    )
    get_forest(renumbered)["structure"] = {"csmi": csmi}

    if tracing_reset:
        _reordered_forest_labels(renumbered)

    return renumbered, csmi

# TODO: this function should ensure Mol atoms exactly matches the forest labels,
# or throw error. It's only a helper that's meant to work if labels were dropped from mol.
def _reordered_forest_labels(mol: ForestMol) -> None:
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


def get_csmi(mol: Mol) -> str:
    structure = _structure(mol)
    csmi = structure.get("csmi")
    if not csmi:
        csmi = MolToSmiles(mol, isomericSmiles=False)
        structure["csmi"] = csmi
    return csmi


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


def conjugated_systems(mol: Mol) -> tuple[frozenset[int], ...]:
    _load_resonance(mol)
    structure = _structure(mol)
    if "conjugated_systems" not in structure:
        raise KeyError("conjugated_systems")
    return structure["conjugated_systems"]


def aromatic_systems(mol: Mol) -> tuple[frozenset[int], ...]:
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


def move_charge_with_bonds(mol: Mol, before: dict[int, float]) -> None:
    """Move formal charge when a bond-order flip would leave it behind.

    The oxygen whose bond order rose by one loses a negative charge. The
    oxygen whose bond order fell gains it. A neutral carbon keeps charge 0
    and moves hydrogen instead, because that hydrogen has to travel with
    the bond.
    """

    for atom in mol.GetAtoms():
        old = before.get(atom.GetIdx())
        if old is None:
            continue
        new = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        delta = int(round(new - old))
        if delta == 0:
            continue
        if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0:
            _shift_hydrogens(atom, -delta)
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
    move_charge_with_bonds(rw, before)
    return rw.GetMol(), written


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

    parents, orders, systems, by_order = _kekule_slots(cache)
    atoms, bonds = _conjugated_component(mol, left)
    seed = _bond_key(left, right)
    if seed not in bonds and mol.GetBondBetweenAtoms(left, right) is not None:
        bonds = frozenset((*bonds, seed))
        if right not in atoms:
            atoms = frozenset((*atoms, right))
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


def parents_for_ends(mol: Mol, start: int, end: int, cache: KekuleParents) -> EndParents:
    """Parents covering ``start`` and ``end``.

    Same conjugated system: that system's assignments. Different systems:
    each system's assignments, not a product of every system. Does not read
    or write ``mol._forest``.
    """

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


def ring_membership(mol: Mol) -> dict[int, tuple[tuple[int, ...], ...]]:
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


def smarts_matches(mol: Mol, smarts: str) -> tuple[dict[int, int], ...]:
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


def sanitized_fragments(mol: Mol, counters: EditCounters | None = None) -> FragmentSplit:
    """Split, drop the dealkylation leaving group, sanitize.

    Empty pieces when any fragment fails. Callers read ``pieces``.
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
        out.append(frag)
    return FragmentSplit(pieces=tuple(out))


def split_fragments(raw: Mol) -> FragmentSplit:
    """One mol, or the fragments of a disconnected reaction product."""

    try:
        groups = GetMolFrags(raw)
    except ValueError:
        return FragmentSplit(pieces=(raw,))
    if len(groups) <= 1:
        return FragmentSplit(pieces=(raw,))
    frags = list(GetMolFrags(raw, asMols=True, sanitizeFrags=False)) or [raw]
    return FragmentSplit(pieces=tuple(frags))


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
    key = get_csmi(target)
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
    key = get_csmi(target)
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

