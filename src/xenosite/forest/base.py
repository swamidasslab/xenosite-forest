"""Defines reaction archetypes to be parameterized or customized in rules.py"""

# Standard Library
import ast
import collections
import copy
import itertools
import re
from collections import defaultdict, deque

# Third Party
from .utils import clean, merge, load
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


def can_smi(line="", rdmol=None):
    """Converts rdmol or line (in SMILES) to canonical SMILES."""

    if isinstance(rdmol, (list, tuple)):
        return [can_smi(rdmol=x) for x in rdmol]

    if rdmol:
        SanitizeMol(rdmol, catchErrors=True)
        line = MolToSmiles(rdmol)

    if "." in line:
        return list(itertools.chain(*[can_smi(line=x) for x in line.split(".")]))

    split = line.split()
    smi = split[0]

    rdmol = MolFromSmiles(smi)
    if rdmol:
        out = MolToSmiles(rdmol)
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
        return site[0], frozenset([topol_equiv[x] for x in site[1]])

    @staticmethod
    def add_current_idx_as_atom_prop(mol, propname="idx"):
        [
            a.SetProp(propname, str(a.GetIdx()))
            for a in mol.GetAtoms()
            if a.GetAtomicNum() != 1
        ]

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
            try:
                record = ast.literal_eval(record.GetProp(cls.tag_name))

            except (KeyError, SyntaxError, ValueError) as err:

                if strict:
                    raise err
                else:
                    record = {}
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
            record = copy.deepcopy(record)

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
        # 1-based atom numbers.

        if not mol.HasProp(self.tag_name):

            initial_tags = {}
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    continue

                atom.SetProp("initial", "1")
                idx = atom.GetIdx()
                initial_tags[idx] = {"idx": [idx], "depth": [0]}

            self._save_tags(mol, initial_tags)
            self._stamp_origin_maps(mol, initial_tags)

    def tag(self, product, reactant, strict=True, **kwargs):
        """Copies all atom tags from reactant and updates the tags in product."""

        previous_tags = self.tags(reactant, strict=strict)
        next_depth = self.next_depth(previous_tags)
        # ext_depth = max(self.depths(previous_tags)) + 1

        new_tags = self.tags(previous_tags, depth=next_depth - 1, strict=strict)

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

        # Copy all missing tags
        for tag in set(previous_tags) - set(new_tags):
            new_tags[tag] = previous_tags[tag]

        self._save_tags(
            product,
            self._tag_new_atoms(
                new_tags, new_atom_indexes, next_depth, self._next_tag(reactant)
            ),
        )

    @classmethod
    def _next_tag(cls, mol):
        try:
            last_tag = mol.GetProp(cls.last_tag_name)
        except KeyError:
            last_tag = "0"

        try:
            last_tag_depth = int(last_tag)
        except ValueError:
            last_tag_depth = 0

        return last_tag_depth + 1
        # return int(mol.GetProp(self.last_tag_name)) + 1

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
            product = RenumberAtoms(product, new_order)
            old_to_new = {old: new for new, old in enumerate(new_order)}
            for rec in tags.values():
                if last_depth not in rec["depth"]:
                    continue
                i = rec["depth"].index(last_depth)
                rec["idx"][i] = old_to_new[rec["idx"][i]]
            self._save_tags(product, tags)

        self.add_current_idx_as_atom_prop(
            product, propname=self.previous_index_prop_name
        )
        self._stamp_origin_maps(product, tags, level=0)
        return product

    def _save_tags(self, mol, tags):
        if isinstance(tags, defaultdict):
            tags = {x: y for x, y in list(tags.items())}

        mol.SetProp(self.tag_name, str(tags))
        mol.SetProp(self.last_tag_name, str(max(tags)))

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
        >>> MolToSmiles(fragment)
        '[10*]C1=CC2=CC=CC=C2C=C1'

        >>> [MolToSmiles(x) for x in other_fragments]
        ['[7*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1']
        """

        # This makes a copy to prevent the input molecule being unexpectedly modified
        mol = Mol(mol)

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
        >>> [MolToSmiles(f) for f in fragments]
        ['[10*]C1=CC2=CC=CC=C2C=C1', '[7*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1']

        """
        mol = Mol(mol)

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
        >>> MolToSmiles(fragments)
        '*=CC(C)=CCC(Cl)C[8*].[1*]=C.[7*]F'
        >>> list(bond_types.keys())
        [frozenset({0, 1}), frozenset({8, 7})]
        >>> list(bond_types.values())
        [rdkit.Chem.rdchem.BondType.DOUBLE, rdkit.Chem.rdchem.BondType.SINGLE]

        """
        SanitizeMol(mol, SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
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

        for idxs, bond in bond_types.items():
            if mol.GetBondBetweenAtoms(*idxs):
                mol.GetBondBetweenAtoms(*idxs).SetBondType(bond)
            else:
                emol = Chem.rdchem.EditableMol(mol)
                emol.AddBond(*idxs, order=bond)
                mol = emol.GetMol()
        return mol

    def _remove_dummy_atoms(self, mol):
        """Deletes dummy atoms from mol introduced by FragmentOnBonds."""
        dummies = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "*"]
        while dummies:
            emol = Chem.rdchem.EditableMol(mol)
            emol.RemoveAtom(dummies.pop())
            mol = emol.GetMol()
            dummies = [
                atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "*"
            ]
        return mol

    def _add_idx_prop_to_dummy_atoms(self, combined):
        """Detects dummy atoms and adds idx prop."""
        dummies = [x for x in combined.GetAtoms() if not x.HasProp("idx")]
        originals = [
            int(x.GetProp("idx")) for x in combined.GetAtoms() if x.HasProp("idx")
        ]
        for idx, dummy in enumerate(dummies, start=max(originals) + 1):
            dummy.SetProp("idx", str(idx))

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
            self.queries = [
                (name,) + self._prepare_query(smarts) for name, smarts in query_smarts
            ]

        self.valid_modifications = set(
            itertools.chain(
                *[[x[0]] if isinstance(x[0], str) else x[0] for x in self.queries]
            )
        )

    def neighbors(self, mol, idx):
        """Get the indexes neighboring idx in mol."""
        return set([n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()])

    def _bfs_atom_path(self, mol, start, end, alternate_bonds=2):
        """Returns a list of all paths in mol between start and end atom.

        >>> mol = MolFromSmiles('c1ccccc1')
        >>> Resonate()._bfs_atom_path(mol,0,1,alternate_bonds=None)
        [[0, 1], [0, 5, 4, 3, 2, 1]]

        >>> Resonate()._bfs_atom_path(mol,0,1)
        [[0, 5, 4, 3, 2, 1]]

        >>> Resonate()._bfs_atom_path(mol,0,1,alternate_bonds=1)
        [[0, 1]]

        """

        # This will be modified by the search function defined below.
        valid_paths = []

        def search(mol, end, path=None):
            """Searches for end idx by extending path."""

            if path is None:
                path = []

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

        mol = Mol(inmol)
        try:
            SanitizeMol(mol)
            Kekulize(mol, clearAromaticFlags=True)
        except ValueError:
            return False

        mol.SetProp("standardized", "1")
        return mol

    def match_queries(self, mol):
        """Returns dict mapping from each atom index to a list of matches to self.queries.

        Each element of the list is a tuple. The first value in each tuple is a mapping from
        mapids to atom indexes, as produced by self.match. The second value in each tuple
        is the name of each query as listed in self.queries.

        >>> query_smarts=[('name1', '[#6R:1][#7,#8h1:2]'),('name2','[#6R:1][#7:2][#6:3]')]
        >>> RRR = EditMol(query_smarts=query_smarts)
        >>> mol = MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
        >>> RRR.match_queries(mol)
        {4: [({1: 4, 2: 3}, 'name1'), ({1: 4, 2: 3, 3: 1}, 'name2')], 7: [({1: 7, 2: 8}, 'name1')]}



        """
        if isinstance(mol, str):
            mol = MolFromSmiles(mol)

        if not mol.HasProp("standardized"):
            mol = self.standardize(mol)
            if not mol:
                return dict()

        mapped_matches = collections.defaultdict(list)

        for modifications, query, queryidx2mapid in self.queries:
            for map2queryidx in self.match(
                mol, query=query, queryidx2mapid=queryidx2mapid
            ):

                mapped_matches[map2queryidx[min(map2queryidx)]].append(
                    (map2queryidx, modifications)
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
        bond = emol.GetBondBetweenAtoms(idx1, idx2)
        if bond:
            emol.RemoveBond(idx1, idx2)
            self.adjust_hydrogen_count(emol, idx1, 1)
            self.adjust_hydrogen_count(emol, idx2, 1)

    def change_bond(self, emol, idx1, idx2, new_bond_type=1):
        """Change bond between idx1 and idx2 in emol to new_bond_type."""
        bond = emol.GetBondBetweenAtoms(idx1, idx2)
        if bond:
            bond.SetBondType(self.bonds[new_bond_type])

    def add_atom(self, emol, idx, new_atomic_num=1, new_bond_type=1):
        """Add atom to emol with new_atomic_num to atom with idx with new_bond_type."""
        i = emol.AddAtom(Atom(new_atomic_num))
        emol.AddBond(idx, i, self.bonds[new_bond_type])

    def replace_atom(self, emol, idx, new_atomic_num=1):
        """Changes atom with idx in emol to atom with new_atomic_num."""
        emol.GetAtomWithIdx(idx).SetAtomicNum(new_atomic_num)

    def set_charge(self, emol, idx, charge):
        """Sets idx charge in emol."""
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
        >>> MolToSmiles(mol)
        'CC=CC=CCl'

        >>> mol = MolFromSmiles('C=CC=CCCl')
        >>> bonds = list(mol.GetBonds())[:-1]
        >>> EditMol().swap_bonds_along_path(mol,[0,1,2,3,4])
        >>> MolToSmiles(mol)
        'CC=CC=CCl'

        """
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

    def apply_modifications(self, mol, modifications, **kwargs):
        """Makes an editable copy of mol and successively applies each submitted modification."""

        emol = RWMol(Mol(mol))
        for map2queryidx, modification in modifications:
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

        for res_frag, rest_of_molecule, bond_types in self._resfrags(mol):

            yield self.join_fragments([res_frag] + rest_of_molecule, bond_types)

    def resonate_with_pair_paths(self, mol, valid_atoms=None):
        """
        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> res_struct, pair, path = next(Resonate().resonate_with_pair_paths(mol))
        >>> MolToSmiles(res_struct, kekuleSmiles=True), pair, path[0], path[-1]
        ('NCCCCC1=C2C=C(CC3=CC=CC(C=CC4=CC=CC=C4)=C3)C=CC2=CC=C1', (5, 6), 5, 6)
        """

        all_system_paths = self.bfs_all_pairs(mol)

        for system_fragment, the_rest, bonds, system in self._resfrags(
            mol, output_systems=True
        ):

            res_struct = self.join_fragments([system_fragment] + the_rest, bonds)

            for pair in itertools.combinations(system, 2):

                if valid_atoms and not set(pair).issubset(valid_atoms):
                    continue

                for path in all_system_paths[frozenset(pair)]:
                    yield Mol(res_struct), pair, path

    def _resfrags(self, mol, output_systems=False):
        """
        >>> RES = Resonate()
        >>> mol = MolFromSmiles('NCCCCC1=CC=CC2=C1C=C(C=C2)CC3=CC=CC(=C3)C=Cc1ccccc1')
        >>> res_frag, other_fragments, bonds, system = next(RES._resfrags(mol,output_systems=True))
        >>> system
        {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}
        >>> MolToSmiles(res_frag, kekuleSmiles=True)
        '[4*]C1=C2C=C([15*])C=CC2=CC=C1'
        >>> [MolToSmiles(f, kekuleSmiles=True) for f in other_fragments]
        ['[5*]CCCCN', '[12*]CC1=CC=CC(C=CC2=CC=CC=C2)=C1']
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


class ReactionRule(AtomTracker):
    """Template for all rules."""

    sites_on = "atoms"
    phase1_sites_on = "bonds"

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

    def format_site(self, site, just_rule_name=False):

        if isinstance(site, tuple):
            return tuple(
                [self.format_site(x, just_rule_name=just_rule_name) for x in site]
            )

        if isinstance(site, list):
            return [self.format_site(x, just_rule_name=just_rule_name) for x in site]

        if isinstance(site, str):
            return site.split("_")[0]

        if just_rule_name and isinstance(site, (tuple, list)):
            return site[0]
        else:
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
        **kwargs
    ):
        if only_unique:
            only_emit_topologically_distinct_sites = False
            unique = False
            unique_smi = []

        topol_equiv = self.topol_equiv(mol)
        if only_emit_topologically_distinct_sites:
            try:
                topol_equiv = self.topol_equiv(mol)
            except:
                # SanitizeMol(mol, SanitizeFlags.SANITIZE_CLEANUP, catchErrors=True)
                topol_equiv = self.topol_equiv(mol)

        tagging = tag_atoms and not do_not_tag_atoms
        if tagging:
            self._clear_atom_maps(mol)
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

            if only_largest_fragment:
                metabolites = [x for x in metabolites if x]
                if metabolites:
                    metabolites = [
                        sorted(
                            metabolites, key=lambda x: x.GetNumAtoms(), reverse=True
                        )[0]
                    ]

            if only_emit_topologically_distinct_sites:
                topsite = self.site_to_topol_site(site, topol_equiv)

                if topsite in seen:
                    continue

                seen.append(topsite)

            if tagging:
                aligned = []
                for metabolite in metabolites:
                    self.tag(metabolite, reactant=mol, strict=strict)
                    aligned.append(self._align_and_stamp(metabolite))
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

            yield outsite, [x for x in metabolites if x]

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

        for site, metabolite in self.metabolize(mol, **kwargs):
            if tuple(site) in sites:

                if just_smiles:
                    yield list(map(MolToSmiles, metabolite))

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
    """Performs reactions specified by SMARTS."""

    parameters = ["Reaction SMARTS"]
    smarts = []
    mapid_site = []

    def __init__(self, rxns=None, mapid_site=None, *args, **kwargs):

        super(SmartsReactionRule, self).__init__(*args, **kwargs)

        if mapid_site:
            self.mapid_site = mapid_site

        if rxns:
            if isinstance(rxns, str):
                rxns = [rxns]
            self.smarts = rxns
        elif self.smarts:
            if isinstance(self.smarts, str):
                self.smarts = [self.smarts]
        else:
            self.smarts = []

        self.rxns = [self._smarts2rxns(rxn, **kwargs) for rxn in self.smarts]

    def metabolites(self, mol, kekulize=True, **kwargs):
        """By default, mol will be kekulized."""
        if kekulize:
            self._kekulize(mol)

        self._remove_props(mol)
        self._clear_atom_maps(mol)

        for rxn_num, rxn in enumerate(self.rxns):
            self._clear_atom_maps(mol)
            for prod_num, prod in enumerate(rxn.RunReactants([mol])):

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
        try:
            Kekulize(mol, clearAromaticFlags=True)
        except ValueError:
            pass

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
