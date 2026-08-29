from __future__ import annotations
import networkx as nx
from typing import NamedTuple, Iterable
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol
from typing import Generator
from xenosite.metabolite.base import AtomTracker
from xenosite.metabolite import rulesets
import ast
import tqdm


class Rxn(NamedTuple):
    reactant: str
    product: str
    type: str
    site: tuple
    map: tuple[tuple[int, ...], tuple[int, ...]]  # from reactant idx to product indx

    def __eq__(self, other):
        return self[:-1] == other[:-1]

    def __hash__(self):
        return hash(self[:-1])


def reordering(mol: Mol):
    reorder = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))
    return {v: n for n, v in enumerate(reorder)}


def metabolites(rule: rulesets.RuleSet, reactant: str) -> Generator[Rxn, None, None]:
    r = MolFromSmiles(reactant)
    assert r, reactant
    rout = MolToSmiles(r, isomericSmiles=False)
    react_reord = reordering(r)

    for product, site, mols in rule.find_path(r):
        if not MolFromSmiles(product):
            continue

        p = mols[1]
        assert type(p) == Mol

        if p.GetNumAtoms() < 2:
            continue

        record: dict = AtomTracker.tags(p)  # type: ignore
        prod_reord = reordering(p)
        mapping = []

        for _, info in record.items():
            if len(info["depth"]) != 2:
                continue
            i = info["idx"]
            mapping.append((react_reord[i[0]], prod_reord[i[1]]))

        site_type, site_idx = site[0]

        mapping: list[tuple[int, int]] = sorted(mapping)
        frm = tuple(x[0] for x in mapping)
        to = tuple(x[1] for x in mapping)

        yield Rxn(
            rout,
            product,
            site_type,
            tuple(sorted(react_reord[s] for s in site_idx)),
            (frm, to),
        )


class MetaboliteNetwork(nx.DiGraph):

    def __init__(self, mols: Mol | str | Iterable[Mol | str] | None = None):
        nx.DiGraph.__init__(self)

        if not mols:
            return

        if isinstance(mols, Mol) or isinstance(mols, str):
            mols = [mols]

        for m in mols:
            self.add_molecule(m)

    def __str__(self):
        return super(nx.DiGraph, self).__str__()

    def __repr__(self):
        return super(nx.DiGraph, self).__repr__()

    def add_molecule(self, mol: str | Mol):
        smiles = self.canonize(mol)
        self.add_node(smiles)
        return smiles

    def canonize(self, x: str | Mol):
        if isinstance(x, str):
            m = MolFromSmiles(x)
            assert m
            return MolToSmiles(m, isomericSmiles=False)

        if isinstance(x, Mol):
            return MolToSmiles(x, isomericSmiles=False)

        raise ValueError(f"Invalid type {type(x)}")

    def paths(self, reactant: str, product: str, extra_depth=1):
        r = self.canonize(reactant)
        p = self.canonize(product)

        shortest = None
        for path in nx.all_simple_paths(self, r, p):
            if shortest is None:
                shortest = path

            if len(shortest) + extra_depth <= len(path):
                break

            yield path

    def expand(self, rule):
        keys = [
            n
            for n in list(self.nodes)
            if len(n) > 1 and not self.nodes[n].get("expanded", False)
        ]

        for n in tqdm.tqdm(keys):
            try:
                for m in metabolites(rule, n):
                    e = n, self.canonize(m.product)
                    self.add_edge(*e)
                    self.edges[e]["rxn"] = self.edges[e].get("rxn", set())
                    self.edges[e]["rxn"].add(m)
            except Exception as e:
                print(f"Error expanding {n}")
            finally:
                self.nodes[n]["expanded"] = True
