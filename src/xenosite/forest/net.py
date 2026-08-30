from __future__ import annotations
import networkx as nx
from typing import NamedTuple, Iterable
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol
from typing import Generator
from xenosite.forest.base import AtomTracker
from xenosite.forest import rulesets
from xenosite.forest.utils import unmapped_smiles
import ast
import tqdm


class Rxn(NamedTuple):
    reactant: str
    product: str
    type: str
    site: tuple
    map: tuple[tuple[int, ...], tuple[int, ...]]  # SMILES output order, not AtomTrace atom numbers

    def __eq__(self, other):
        return self[:-1] == other[:-1]

    def __hash__(self):
        return hash(self[:-1])


def reordering(mol: Mol):
    """Map mol GetIdx() -> position in unmapped canonical SMILES output order."""
    copy = Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    MolToSmiles(copy, isomericSmiles=False)
    order = ast.literal_eval(copy.GetProp("_smilesAtomOutputOrder"))
    return {v: n for n, v in enumerate(order)}


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

        # Tag idx values are 0-based GetIdx(); remap to SMILES output order
        # (not AtomTrace 1-based atom numbers).
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
            return unmapped_smiles(m, isomericSmiles=False)

        if isinstance(x, Mol):
            return unmapped_smiles(x, isomericSmiles=False)

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
            except Exception:
                pass
            finally:
                self.nodes[n]["expanded"] = True

    def copy(self):
        """Deep copy of the network (nodes, edges, and reaction metadata)."""
        other = MetaboliteNetwork()
        other.add_nodes_from((n, dict(d)) for n, d in self.nodes(data=True))
        for u, v, data in self.edges(data=True):
            other.add_edge(u, v, **{k: (set(val) if k == "rxn" else val)
                                    for k, val in data.items()})
        return other

    def prune(
        self,
        *,
        nodes=None,
        max_generation=None,
        reaction_types=None,
        path=None,
        roots=None,
        inplace=False,
    ):
        """Filter the network and return a copy (or modify in place).

        Parameters
        ----------
        nodes:
            Keep only these molecules.
        max_generation:
            Keep molecules with BFS generation ``<= max_generation``.
        reaction_types:
            Keep only edges whose reaction names intersect this set.
        path:
            Keep only the molecules (and consecutive edges) on this path.
        roots:
            Roots used when computing generations.
        inplace:
            If True, replace this graph; otherwise return a new network.
        """
        from xenosite.forest.draw import filter_graph

        view = filter_graph(
            self,
            nodes=nodes,
            max_generation=max_generation,
            reaction_types=reaction_types,
            path=path,
            roots=roots,
        )
        # Preserve isolates when only capping generation.
        if nodes is None and reaction_types is None and path is None and max_generation is not None:
            from xenosite.forest.draw import generations

            dist = generations(self, roots=roots)
            for n, g in dist.items():
                if g <= max_generation and n not in view:
                    view.add_node(n, **dict(self.nodes[n]))

        other = MetaboliteNetwork()
        other.add_nodes_from((n, dict(d)) for n, d in view.nodes(data=True))
        for u, v, data in view.edges(data=True):
            other.add_edge(u, v, **{k: (set(val) if k == "rxn" else val)
                                    for k, val in data.items()})

        if inplace:
            self.clear()
            self.add_nodes_from((n, dict(d)) for n, d in other.nodes(data=True))
            for u, v, data in other.edges(data=True):
                self.add_edge(u, v, **data)
            return self
        return other

    def draw(self, highlight=None, **kwargs):
        """Layered network figure with molecule drawings (Jupyter-friendly HTML/SVG).

        Common options: ``nodes``, ``mol_size``, ``max_generation``,
        ``reaction_types``, ``background``, ``mol_background``, ``mol_border``,
        ``show_labels``, ``show_legend``, ``title``, ``fade_others``.
        See ``xenosite.forest.draw.draw_network``.
        """
        from xenosite.forest.draw import draw_network

        return draw_network(self, highlight=highlight, **kwargs)

    def grid(self, mols=None, **kwargs):
        """Grid of metabolite structures, ordered by generation."""
        from xenosite.forest.draw import draw_grid

        return draw_grid(self, mols=mols, **kwargs)

    def draw_path(self, path, **kwargs):
        """Linear reaction scheme for one path through the network."""
        from xenosite.forest.draw import draw_path

        return draw_path(self, path, **kwargs)

    def save_draw(self, path, highlight=None, **kwargs):
        """Draw the network and save to SVG/PDF/PS/PNG. See ``Drawing.save``."""
        return self.draw(highlight=highlight, **kwargs).save(path)

    def _repr_html_(self):
        return self.draw().html
