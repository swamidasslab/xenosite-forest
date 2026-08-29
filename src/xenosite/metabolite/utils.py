"""Molecule loading and cleanup helpers."""

import itertools

from rdkit import Chem


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


def clean(mol):
    """Split mol into a list of sanitized, kekulized fragments."""
    if isinstance(mol, (list, tuple)):
        return list(itertools.chain.from_iterable(clean(x) for x in mol))

    out = []
    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
        for atom in frag.GetAtoms():
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        Chem.SanitizeMol(frag, Chem.SanitizeFlags.SANITIZE_CLEANUP, catchErrors=True)
        for atom in frag.GetAtoms():
            atom.SetNoImplicit(False)
        Chem.SanitizeMol(frag, catchErrors=True)
        try:
            Chem.Kekulize(frag, clearAromaticFlags=True)
        except ValueError:
            pass
        out.append(frag)
    return out


def merge(intervals):
    """Merge a list of overlapping lists into a single list of disjoint sets.

    >>> merge([[0,1,2,3,4],[0,5,6],[7,8,9]])
    [set([0, 1, 2, 3, 4, 5, 6]), set([8, 9, 7])]

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
