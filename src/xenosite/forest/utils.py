"""Molecule loading and cleanup helpers."""

import itertools
import warnings

from rdkit import Chem, rdBase

# Prevents spammy rdkit messages during sanitization probes.
rdBase.DisableLog("rdApp.*")


def unmapped_smiles(mol, **kwargs):
    """SMILES with atom-map numbers cleared so structure identity ignores `:N` maps."""
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
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


def _warn_invalid_metabolite(mol, reason=None):
    smi = _mol_smiles(mol)
    reason = reason or sanitize_reason(mol) or "RDKit sanitization failed"
    warnings.warn(
        "Dropping RDKit-invalid metabolite %s (%s)" % (smi, reason),
        UserWarning,
        stacklevel=3,
    )


def clean(mol):
    """Split mol into sanitized, kekulized fragments.

    If any fragment is RDKit-invalid, the whole product set is dropped with a
    warning. Leaving leftover fragments from a failed reaction would emit
    chemically incomplete structures (for example acetaldehyde from a failed
    quinone formation).
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
            _warn_invalid_metabolite(frag)
            return []
        out.append(sanitized)
    return out


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
