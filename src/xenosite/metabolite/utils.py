# Standard Library
import importlib
import itertools
import io
import pickle
import pkgutil
import re
import sys
import tempfile
import unittest
import os
import tarfile


import numpy as np

# import openbabel
import pandas as pd

# import pybel
from rdkit import Chem


def interchange_rdmols_and_pymols(inp, loaded=True, skip=None):
    if not loaded:
        inp = load(inp, single=True, skip=skip)
    if isinstance(inp, pybel.Molecule):
        return pymol_to_rdmol(inp, skip=skip)
    elif isinstance(inp, Chem.Mol):
        return rdmol_to_pymol(inp, skip=skip)


def pymol_to_rdmol(pymol, skip=None):
    if skip is None:
        skip = []

    sdf_file = tempfile.NamedTemporaryFile(suffix="sdf").name
    with open(sdf_file, "w") as F:
        F.write(pymol.write("sdf"))

    supp = Chem.SDMolSupplier(sdf_file, sanitize=False, strictParsing=False)

    try:
        rdmol = next(supp)
    except StopIteration:
        rdmol = None
    os.remove(sdf_file)

    return rdmol


def rdmol_to_pymol(rdmol, hydrogens=True, skip=None):
    import pybel

    if skip is None:
        skip = []

    data = rdmol.GetPropsAsDict()

    Chem.SanitizeMol(rdmol, Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
    try:
        Chem.Kekulize(rdmol, clearAromaticFlags=True)
    except ValueError:
        pass

    try:
        block = Chem.MolToMolBlock(rdmol)

    except ValueError:
        try:
            block = Chem.MolToMolBlock(rdmol, kekulize=False)
        except ValueError:
            return None

    pymol = pybel.readstring("sdf", block)
    if hydrogens:
        pymol.removeh()
        pymol.addh()

    for k, v in list(data.items()):
        if k not in skip:
            pymol.data[k] = v
    return pymol


def load(inp, single=False, moltype="rdmol", skip=None):
    """Constructs rdkit molecules (default) or pybel molecules from various inputs.

    Args:
        - inp: the input. Can be a string or list of strings.
        - single: return only the first molecule
        - moltype: the type of molecule to output, rdmol (default) or pymol

    Returns:
        a list of rdkit or pybel molecules OR just the first one if single=True

    """
    if moltype == "pybel":
        import pybel

    if skip is None:
        skip = []

    if isinstance(inp, (list, tuple)):
        if single:
            if inp:
                return load(inp[0], single=single, moltype=moltype, skip=skip)
            else:
                return None
        else:
            return list(
                itertools.chain(*[load(x, moltype=moltype, skip=skip) for x in inp])
            )

    elif isinstance(inp, str):
        mol = None
        mols = []

        if inp.endswith(".sdf"):
            mols = list(
                Chem.ForwardSDMolSupplier(inp, strictParsing=False, removeHs=True)
            )

            if None in mols:
                mols = list(
                    Chem.ForwardSDMolSupplier(
                        inp, strictParsing=False, sanitize=False, removeHs=True
                    )
                )

        elif inp.endswith(".smi"):
            mols = []
            for line in open(inp):
                mols.append(Chem.MolFromSmiles(line.split()[0]))

        if mols:
            if len(mols) == 1:
                mol = mols.pop()
            else:
                return load(mols, single=single, moltype=moltype, skip=skip)

        else:
            if "$$$$" in inp:
                mol = Chem.MolFromMolBlock(inp, removeHs=True)
            elif " " not in inp:
                mol = Chem.MolFromSmiles(inp)

            else:
                return load(inp.split(), moltype=moltype, skip=skip)

        if not mol:
            raise ValueError("Error constructing RDKit Molecule from %s" % inp)

        return load([mol], single=single, moltype=moltype, skip=skip)

    elif isinstance(inp, Chem.Mol):
        inp = [inp]

    elif isinstance(inp, pybel.Molecule):
        if moltype == "rdmol":
            inp = [pymol_to_rdmol(inp, skip=skip)]
        elif moltype == "pymol":
            inp = [inp]
        else:
            raise ValueError("Valid values for moltype are rdmol or pymol")

    else:
        raise ValueError("Must submit a string or a list or tuple of strings")

    outmols = []

    for outmol in inp:
        if not outmol:
            continue

        if moltype == "rdmol":
            if isinstance(outmol, Chem.Mol):
                outmols.append(outmol)
            elif isinstance(outmol, pybel.Molecule):
                outmols.append(pymol_to_rdmol(outmol, skip=skip))

        if moltype == "pymol":
            if isinstance(outmol, Chem.Mol):
                outmols.append(rdmol_to_pymol(outmol, skip=skip))
            elif isinstance(outmol, pybel.Molecule):
                outmols.append(outmol)

    successfuly_constructed_outmols = [x for x in outmols if x]

    if single:
        if successfuly_constructed_outmols:
            final = successfuly_constructed_outmols[0]
    else:
        final = successfuly_constructed_outmols

    return final


def clean(mol):
    """Splits mol into a list of sanitized, kekulized fragments."""

    if isinstance(mol, (list, tuple)):
        return list(itertools.chain(*[clean(x) for x in mol]))

    out = []

    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):

        [(a.SetNoImplicit(True), a.SetNumExplicitHs(0)) for a in frag.GetAtoms()]
        Chem.SanitizeMol(frag, Chem.SanitizeFlags.SANITIZE_CLEANUP, catchErrors=True)

        [a.SetNoImplicit(False) for a in frag.GetAtoms()]

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


def can_smi(mol):
    """Converts mol or line (in SMILES) to canonical SMILES."""

    if mol is None:
        return frozenset([])

    if isinstance(mol, (list, tuple)):
        return frozenset.union(*[can_smi(x) for x in mol])

    elif isinstance(mol, str):
        return can_smi(Chem.MolFromSmiles(mol.split()[0]))

    elif isinstance(mol, pybel.Molecule):
        return can_smi(Chem.MolFromMolBlock(mol.write("sdf")))

    elif isinstance(mol, Chem.Mol):
        Chem.SanitizeMol(mol, catchErrors=True)
        line = Chem.MolToSmiles(mol)

    else:

        raise ValueError("Must input Chem.Mol, pybel.Molecule, string, or a list/tuple")

    if "." in line:

        frags = []

        for frag in line.split("."):
            can_frag = can_smi(frag)
            frags.append(can_frag)
        largest_fragment = frags[pd.Series(frags).map(lambda x: len(str(x))).idxmax()]
        return largest_fragment

    split = line.split()
    smi = split[0]

    mol = Chem.MolFromSmiles(smi)

    if mol:
        mol = Chem.RenumberAtoms(
            mol, Chem.CanonicalRankAtoms(mol, includeChirality=False)
        )
        out = Chem.CanonSmiles(Chem.MolToSmiles(mol))

    else:
        out = smi

    out = re.sub(r"\[CH*\]", "C", out)
    out = re.sub(r"\[H\]", "", out)
    out = re.sub(r"\[C\]", "C", out)

    return frozenset([out])


def can_smi_new(mol):
    """Converts mol or line (in SMILES) to canonical SMILES."""
    if isinstance(mol, pybel.Molecule):
        mol = Chem.MolFromMolBlock(
            mol.write("sdf"), sanitize=False, strictParsing=False
        )

    if isinstance(mol, Chem.Mol):
        Chem.SanitizeMol(mol, catchErrors=True)
        try:
            mol = Chem.RemoveHs(mol)
        except:
            pass

        try:
            Chem.rdmolops.RemoveStereochemistry(mol)
        except:
            pass

        mol = Chem.RenumberAtoms(
            mol, Chem.CanonicalRankAtoms(mol, includeChirality=False)
        )
        out = Chem.CanonSmiles(Chem.MolToSmiles(mol))
        return frozenset([out])

    elif isinstance(mol, str):

        if mol.endswith(".sdf"):
            mols = [
                x
                for x in Chem.SDMolSupplier(
                    mol.encode(), sanitize=False, strictParsing=False
                )
            ]

        elif mol.endswith(".smi"):
            mols = [Chem.MolFromSmiles(x.split()[0], sanitize=False) for x in open(mol)]

        else:
            mols = [Chem.MolFromSmiles(mol.split()[0], sanitize=False)]

        return [can_smi(x) for x in mols]

    else:
        raise ValueError("Expected input type", type(mol))


def sdf2cansmi(sdf, outfile=None):
    out = {}

    for num, mol in enumerate(Chem.SDMolSupplier(sdf.encode()), start=1):
        out[can_smi([mol])] = num

    if outfile is not None:
        with open(outfile, "w") as F:
            F.write(str(out))

    return out


def sdf_to_can_smi_dict_file(sdf, outfile):
    out = {}
    for num, mol in enumerate(Chem.SDMolSupplier(sdf), start=1):
        out[can_smi([mol])] = num

    with open(outfile, "w") as F:
        F.write(str(out))


def cansmi_list(SMILES):
    all_cansmi = set([])
    for smi in SMILES:
        canonical_smiles = can_smi(smi)
        all_cansmi.add(canonical_smiles)
    return all_cansmi


class PyMolPredictorBase(object):
    def pandas_to_csv(self, df):
        S = io.StringIO()
        df.to_csv(S, sep="\t")
        S.seek(0)
        return S

    def csv_to_pandas(self, call):
        S = io.StringIO()
        call(S)
        S.seek(0)
        return pd.read_csv(S, sep="\t", index_col="ID")

    def predict(self, pymol):
        raise NotImplemented

    def structure(self, pymol, site):
        raise NotImplemented

    def all_structures(self, pymol):
        raise NotImplemented

    def scored_metabolites(self, pymol):
        raise NotImplemented

    def read(self, fmt, string):

        mol = pybel.readstring(fmt, string)

        mol.removeh()
        # mol = pybel.Molecule(mol)
        mol.convertdbonds()

        heavy = [x.idx for x in mol.atoms if not x.OBAtom.IsHydrogen()]
        light = [x.idx for x in mol.atoms if x.OBAtom.IsHydrogen()]
        mol.OBMol.RenumberAtoms(heavy + light)

        openbabel.obErrorLog.SetOutputLevel(0)
        mol = pybel.readstring("sdf", mol.write("sdf"))
        openbabel.obErrorLog.SetOutputLevel(1)

        return mol

    def predict_structures(self, pymol):
        pred = self.predict(pymol)["site"]
        poss = self.all_structures(pymol)

        out = {}
        for site, structure in poss.items():
            b1_site = frozenset([x + 1 for x in site])
            out[structure] = pred[b1_site]
            # out[poss[site]] = pred[b1_site]

        return out


class RebasingUnpickler(pickle.Unpickler):
    def __init__(self, base, *args, **kw):
        pickle.Unpickler.__init__(self, *args, **kw)
        self.base = base

    def find_class(self, module, name):
        try:
            return pickle.Unpickler.find_class(
                self, ".".join([self.base, module]), name
            )
        except ImportError:
            return pickle.Unpickler.find_class(self, module, name)


class SubstitutingRebasingUnpickler(pickle.Unpickler):
    def __init__(self, substitutions, *args, **kwargs):
        pickle.Unpickler.__init__(self, *args, **kwargs)
        self.subs = substitutions

    def find_class(self, module, name):
        try:
            nmod = self.subs[module]
            return pickle.Unpickler.find_class(self, nmod, name)
        except KeyError:
            return pickle.Unpickler.find_class(self, module, name)


class TestBase(unittest.TestCase):
    def assertAtLeastOneNotAlmostEqual(self, vals1, vals2, places=None):

        s1 = set(vals1.keys())
        s2 = set(vals2.keys())

        combined = s1 & s2

        matches = []

        for k in combined:
            v1 = vals1[k]
            v2 = vals2[k]

            if places:
                v1 = np.round(v1, places)
                v2 = np.round(v2, places)

            matches.append(v1 == v2)

        try:
            self.assertIn(False, matches)

        except AssertionError:

            compared = pd.concat([pd.Series(vals1), pd.Series(vals2)], axis=1)

            matches = pd.Series(matches, index=compared.index)
            matches.name = "matches"

            compared = pd.concat([matches, compared], axis=1)

            print("\n" * 3, file=sys.stderr)
            compared.to_csv(sys.stdout, sep="\t", header=0)
            print("\n" * 3, file=sys.stderr)

            raise AssertionError

    def number_of_values(self, vals, correct_number_of_values=1, decimals=None):
        if decimals:
            vals = np.round(vals, decimals)
        self.assertEqual(correct_number_of_values, len(set(vals)))

    def cansmi_list(self, SMILES):
        all_cansmi = set([])
        for smi in SMILES:
            canonical_smiles = can_smi(smi)
            self.assertIsNotNone(canonical_smiles)
            all_cansmi.add(canonical_smiles)
        return all_cansmi

    def equal_SMILES(self, SMILES1, SMILES2, NOT_EQUAL=False):
        cansmi1 = can_smi(SMILES1)
        cansmi2 = can_smi(SMILES2)

        self.assertIsNotNone(cansmi1)
        self.assertIsNotNone(cansmi2)

        if NOT_EQUAL:
            self.assertNotEqual(cansmi1, cansmi2)
        else:
            self.assertEqual(cansmi1, cansmi2)

    def equal_SMILES_population(self, SMILES1_list, SMILES2_list, NOT_EQUAL=False):

        cansmi1 = set(self.cansmi_list(SMILES1_list))
        cansmi2 = set(self.cansmi_list(SMILES2_list))

        if NOT_EQUAL:
            self.assertNotEqual(cansmi1, cansmi2)
        else:
            self.assertSetEqual(
                cansmi1,
                cansmi2,
                "Missed:\n%s" % "\n".join(cansmi2 - cansmi1)
                + "\n\n"
                + "Full Set:\n%s" % "\n".join(SMILES2_list),
            )

    def predictions_not_consistent(self, preds1, preds2, places=None):

        for outer_key in preds1:

            vals1 = preds1[outer_key]
            vals2 = preds2[outer_key]

            if isinstance(vals1, dict):
                self.assertAtLeastOneNotAlmostEqual(vals1, vals2, places=places)

    def prediction_consistency(
        self, preds1, preds2, places=None, infer_common_keys=False
    ):

        if infer_common_keys:
            keys = set(preds1.keys()) & set(preds2.keys())

        else:
            keys = list(preds1.keys())

        for outer_key in keys:

            vals1 = preds1[outer_key]
            vals2 = preds2[outer_key]

            if isinstance(vals1, float):
                if places:
                    self.assertAlmostEqual(vals1, vals2, places)

                else:
                    self.assertEqual(vals1, vals2)

            else:
                for inner_key in vals1:

                    v1 = vals1[inner_key]
                    v2 = vals2[inner_key]

                    if places:
                        self.assertAlmostEqual(v1, v2, places)
                    else:
                        self.assertEqual(v1, v2)
