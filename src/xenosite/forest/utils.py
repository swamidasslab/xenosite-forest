"""Molecule loading and cleanup helpers."""

import itertools
import logging

from rdkit import Chem, rdBase

# Prevents spammy rdkit messages during sanitization probes.
rdBase.DisableLog("rdApp.*")

_log = logging.getLogger(__name__)


def refresh_mol(mol):
    """Fill valence caches without failing on odd valences (RDKit 2026).

    Resonance copies from ``join_fragments`` can look valid but have an empty
    implicit-H cache. RDKit 2026 asserts those caches exist inside
    ``RunReactants`` and ``MolToSmiles``.
    """
    if mol is None:
        return mol
    if isinstance(mol, (list, tuple)):
        for item in mol:
            refresh_mol(item)
        return mol
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return mol


def unmapped_smiles(mol, **kwargs):
    """SMILES with atom-map numbers cleared so structure identity ignores `:N` maps."""
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    refresh_mol(mol)
    return Chem.MolToSmiles(mol, **kwargs)


def canon_smi(obj, **kwargs):
    """Canonical SMILES for structure identity, stable across RDKit versions.

    Accepts a SMILES string, an RDKit mol, or a sequence of those. Mol inputs
    have atom-map numbers cleared first. Each structure is parsed and re-emitted
    with this RDKit, so two writings of the same molecule compare equal.
    Do not pass kekuleSmiles=True: Kekulé form is not unique across RDKit versions.

    >>> canon_smi('[10*]C1=CC2=CC=CC=C2C=C1') == canon_smi('[10*]C1=CC2=C(C=CC=C2)C=C1')
    True
    """
    if isinstance(obj, (list, tuple)):
        return [canon_smi(x, **kwargs) for x in obj]
    if isinstance(obj, str):
        smi = obj
    else:
        smi = unmapped_smiles(obj, **kwargs)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    return Chem.MolToSmiles(mol, **kwargs)


def load(inp, single=False):
    """Construct RDKit molecules from SMILES, files, mol blocks, or existing mols.

    Args:
        inp: a string, RDKit mol, or a list/tuple of those.
        single: if True, return only the first molecule.

    Returns:
        a list of RDKit molecules, or the first molecule if single=True.
    """
    if isinstance(inp, (list, tuple)):
        if single:
            return load(inp[0], single=True) if inp else None
        return list(itertools.chain.from_iterable(load(x) for x in inp))

    if isinstance(inp, Chem.Mol):
        mols = [inp]
    elif isinstance(inp, str):
        mols = _load_from_string(inp)
    else:
        raise ValueError("Must submit a string, RDKit molecule, or a list/tuple of those.")

    mols = [mol for mol in mols if mol]
    if single:
        return mols[0] if mols else None
    return mols


def _load_from_string(inp):
    if inp.endswith(".sdf"):
        mols = list(Chem.ForwardSDMolSupplier(inp, strictParsing=False, removeHs=True))
        if None in mols:
            mols = list(
                Chem.ForwardSDMolSupplier(
                    inp, strictParsing=False, sanitize=False, removeHs=True
                )
            )
        return load(mols)

    if inp.endswith(".smi"):
        with open(inp) as handle:
            return [Chem.MolFromSmiles(line.split()[0]) for line in handle]

    if "$$$$" in inp:
        mol = Chem.MolFromMolBlock(inp, removeHs=True)
        if not mol:
            raise ValueError("Error constructing RDKit molecule from mol block")
        return [mol]

    if " " not in inp:
        mol = Chem.MolFromSmiles(inp)
        if not mol:
            raise ValueError("Error constructing RDKit molecule from %s" % inp)
        return [mol]

    return load(inp.split())


def _mol_smiles(mol):
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return "<unwritable>"


def is_rdkit_valid(mol):
    """True if ``mol`` sanitizes and round-trips through SMILES."""
    if mol is None:
        return False
    smi = _mol_smiles(mol)
    if not smi or smi == "<unwritable>":
        return False
    return Chem.MolFromSmiles(smi) is not None


def _sanitize_kekulize(mol, reset_hs=False):
    if reset_hs:
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_CLEANUP, catchErrors=True)
    for atom in mol.GetAtoms():
        atom.SetNoImplicit(False)
    Chem.SanitizeMol(mol)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except ValueError:
        pass
    refresh_mol(mol)
    smi = Chem.MolToSmiles(mol)
    if Chem.MolFromSmiles(smi) is None:
        raise ValueError("SMILES round-trip failed: %s" % smi)
    return mol


def sanitize_reason(mol):
    """Short RDKit sanitization error for ``mol``, or None if it is valid."""
    try:
        _sanitize_kekulize(Chem.Mol(mol), reset_hs=False)
        return None
    except Exception as err:
        msg = str(err).strip().split("\n")[0]
        return msg or err.__class__.__name__


def sanitize_metabolite(mol):
    """Return a sanitized copy of ``mol``, or None if it is RDKit-invalid.

    Prefer keeping existing hydrogens (so imidazole [nH] is not lost). If that
    fails, reset explicit hydrogens — some reaction products keep leftover Hs —
    and sanitize again.
    """
    for reset_hs in (False, True):
        try:
            return _sanitize_kekulize(Chem.Mol(mol), reset_hs=reset_hs)
        except Exception:
            continue
    return None


def clean(mol):
    """Split mol into sanitized, kekulized fragments.

    If any fragment is RDKit-invalid, the whole product set is dropped.
    Details are logged at DEBUG only so default runs stay quiet.
    Leaving leftover fragments from a failed reaction would emit chemically
    incomplete structures (for example acetaldehyde from a failed quinone
    formation).
    """
    if isinstance(mol, (list, tuple)):
        parts = [clean(x) for x in mol]
        if any(len(part) == 0 for part in parts):
            return []
        return list(itertools.chain.from_iterable(parts))

    out = []
    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
        sanitized = sanitize_metabolite(frag)
        if sanitized is None:
            _log.debug(
                "Dropping RDKit-invalid metabolite %s (%s)",
                _mol_smiles(frag),
                sanitize_reason(frag) or "RDKit sanitization failed",
            )
            return []
        out.append(sanitized)
    return out


_PARENT_ATOM_PROPS = ("current_idx", "react_atom_idx")


def _atom_from_reactant(atom):
    """True when ``atom`` was copied from the reactant (has a parent index prop)."""
    return any(atom.HasProp(p) for p in _PARENT_ATOM_PROPS)


def collapse_conjugate_to_star(product):
    """Replace newly added conjugate atoms with a dummy ``*`` at each attachment.

    Parent atoms are those carrying ``current_idx`` / ``react_atom_idx`` from
    ``SmartsReactionRule`` / AtomTracker. Atoms without those props are treated
    as the conjugate group.
    """
    if product is None:
        return product

    new_atoms = set()
    for atom in product.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if not _atom_from_reactant(atom):
            new_atoms.add(atom.GetIdx())
    if not new_atoms:
        return Chem.Mol(product)

    attach = set()
    for bond in product.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        a_new, b_new = a in new_atoms, b in new_atoms
        if a_new == b_new:
            continue
        attach.add(b if a_new else a)
    if not attach:
        return Chem.Mol(product)

    def _after_remove(idx):
        return idx - sum(1 for n in new_atoms if n < idx)

    rw = Chem.RWMol(product)
    for idx in sorted(new_atoms, reverse=True):
        rw.RemoveAtom(idx)
    for pa in sorted(_after_remove(a) for a in attach):
        dummy = rw.AddAtom(Chem.Atom(0))
        rw.AddBond(pa, dummy, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol, catchErrors=True)
    refresh_mol(mol)
    return mol


def label_star_atoms(mol, label):
    """Set CX ``atomLabel`` on dummy (atomic number 0) atoms. Mutates ``mol``."""
    if mol is None or not label:
        return mol
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetProp("atomLabel", label)
    return mol


def mol_to_cxsmiles(mol, isomericSmiles=False):
    """Non-isomeric CXSMILES when possible; else plain SMILES.

    Clears non-``atomLabel`` atom props so the CX block only carries labels.
    The token before ``|`` is valid on its own and depicts a bare ``*``.
    """
    if mol is None:
        return None
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        label = atom.GetProp("atomLabel") if atom.HasProp("atomLabel") else None
        for prop in list(atom.GetPropNames()):
            atom.ClearProp(prop)
        if label is not None:
            atom.SetProp("atomLabel", label)
    refresh_mol(mol)
    try:
        out = Chem.MolToCXSmiles(mol, isomericSmiles)
    except Exception:
        return unmapped_smiles(mol, isomericSmiles=isomericSmiles)
    return out or unmapped_smiles(mol, isomericSmiles=isomericSmiles)


def apply_star_conjugate(mol, label=None):
    """Collapse a full conjugate to a star adduct; optionally set ``atomLabel``."""
    collapsed = collapse_conjugate_to_star(mol)
    if label:
        label_star_atoms(collapsed, label)
    return collapsed


def merge(intervals):
    """Merge a list of overlapping lists into a single list of disjoint sets.

    >>> merge([[0,1,2,3,4],[0,5,6],[7,8,9]])
        [{0, 1, 2, 3, 4, 5, 6}, {8, 9, 7}]

    """
    if not intervals:
        return []

    sorted_intervals = sorted(map(set, intervals))
    merged = [sorted_intervals.pop(0)]

    while sorted_intervals:
        interval = sorted_intervals.pop(0)
        overlap = [x for x in merged if x & interval]
        if overlap:
            try:
                assert len(overlap) == 1
            except AssertionError:
                return []
            already_established = overlap.pop()
            already_established |= interval
        else:
            merged.append(interval)

    return merged
