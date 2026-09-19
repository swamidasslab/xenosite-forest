"""Molecule questions for the proof of concept.

RDKit itself is imported in :mod:`xenosite.refactor_poc.rdkit_api`. This
module calls those names.

Answers about a molecule are cached on ``get_forest(mol)["structure"]``.
The key is one string. Process-wide data, such as parsed SMARTS reactions,
stays a module dict. Do not cache a result on a molecule this function
then edits. Edits belong on a copy from :func:`copy_mol` or :func:`rw_copy`.

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
from typing import Any

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
    GetMolFrags,
    Mol,
    MolFromSmarts,
    MolFromSmiles,
    MolToSmiles,
    RWMol,
    ReactionFromSmarts,
    RenumberAtoms,
    ResonanceMolSupplier,
    SanitizeFlags,
    SanitizeMol,
)

from xenosite.refactor_poc.records import (
    Forest,
    Formula,
    FragmentSplit,
    McsResult,
    Structure,
)

DisableLog("rdApp.*")

_REACTION_CACHE: dict[str, ChemicalReaction] = {}


def get_forest(mol: Mol, new_structure: bool = False) -> Forest:
    # Declared on the type, but a fresh RDKit mol has not been stamped yet.
    try:
        forest: Forest | None = mol._forest
    except AttributeError:
        forest = None
    if forest is None:
        structure: Structure = {}
        forest = {"structure": structure}
        mol._forest = forest

    if new_structure:
        forest = copy.deepcopy(forest)
        del forest["structure"]

    return forest


def _structure(mol: Mol) -> Structure:
    forest = get_forest(mol)
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


def copy_mol(mol: Mol) -> Mol:
    """``Chem.Mol`` copy that also carries a deep-copied forest."""

    out = Mol(mol)
    forest = getattr(mol, "_forest", None)
    if forest is not None:
        out._forest = copy.deepcopy(forest)
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


def run_reactants(smarts: str, mol: Mol) -> tuple[tuple[Mol, ...], ...]:
    """Run one cached reaction on ``mol``. Empty when RDKit refuses the run."""

    reaction = reaction_from_smarts(smarts)
    try:
        product_sets = reaction.RunReactants((mol,))
    except (RuntimeError, ValueError):
        return ()
    if not product_sets:
        return ()
    return tuple(product_sets)


def cannonicalize_order(mol: Mol, tracing_reset: bool = True) -> tuple[Mol, str]:
    """Renumber into canonical SMILES order. Returns the new mol and that SMILES.

    The mol is not a record field, so this stays a tuple. The SMILES is also
    the only structure entry copied onto the renumbered molecule.
    """

    csmi = MolToSmiles(mol, isomericSmiles=False)
    smiles_order = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))

    renumber_map = [0] * mol.GetNumAtoms()
    for new_pos, old_idx in enumerate(smiles_order):
        renumber_map[old_idx] = new_pos

    renumbered = RenumberAtoms(mol, renumber_map)
    forest = get_forest(mol)
    renumbered._forest = copy.deepcopy(forest)
    renumbered._forest["structure"] = {"csmi": csmi}

    if tracing_reset:
        _reordered_forest_labels(renumbered)

    return renumbered, csmi

def _reordered_forest_labels(mol: Mol) -> None:
    forest = get_forest(mol)
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        if atom.GetAtomicNum() != 1:
            tag = atom.GetProp("forestLabel")
            if "atom_trace" not in forest:
                raise KeyError("atom_trace")
            trace = forest["atom_trace"]
            if "records" not in trace:
                raise KeyError("records")
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


def mol_from_smiles(smiles: str) -> Mol:
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("could not parse %r" % (smiles,))
    return mol


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


def _bump(counters: Any, name: str, amount: int = 1) -> None:
    if counters is None:
        return
    setattr(counters, name, getattr(counters, name) + amount)


def sanitized_fragments(mol: Mol, counters: Any = None) -> FragmentSplit:
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
        if SanitizeMol(frag, catchErrors=True):
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


def mcs_matches(reactant: Mol, target: Mol) -> McsResult:
    """Every full-size embedding of ``target`` on ``reactant``, not only the best.

    The reactant structure holds the NamedTuple, keyed by the target's
    canonical SMILES. A later search can read every embedding. This function
    does not score them.
    """

    structure = _structure(reactant)
    cache: dict[str, McsResult] = structure.setdefault("mcs_matches", {})
    targets: dict[str, McsResult] = structure.setdefault("mcs_targets", {})
    key = get_csmi(target)
    cached = cache.get(key)
    if cached is not None:
        return cached

    query = _mcs_query(reactant, target)
    found = McsResult(embeddings=_full_matches(reactant, query))
    cache[key] = found
    targets[key] = McsResult(embeddings=_full_matches(target, query))
    return found


def mcs_target_matches(reactant: Mol, target: Mol) -> McsResult:
    """Target-side embeddings for the same MCS query. Filled with :func:`mcs_matches`."""

    mcs_matches(reactant, target)
    key = get_csmi(target)
    structure = _structure(reactant)
    if "mcs_targets" not in structure:
        raise KeyError("mcs_targets")
    return structure["mcs_targets"][key]


def _mcs_query(reactant: Mol, target: Mol) -> Mol | None:
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

