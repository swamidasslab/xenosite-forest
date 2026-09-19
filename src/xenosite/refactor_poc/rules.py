"""Define specific reaction rules."""

# Standard Library
import itertools
from collections import defaultdict, deque
import json
import copy

# Third Party
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import CanonicalRankAtoms, RenumberAtoms, MolToSmiles
from rdkit.Chem import GetMolFrags, SanitizeMol
from rdkit.Chem.rdchem import (
    Atom,
    BondType,
    KEKULE_ALL,
    Mol,
    ResonanceMolSupplier,
    RWMol,
)
from typing import Any, NamedTuple, TypedDict
from collections.abc import Callable, Generator


rdBase.DisableLog("rdApp.*")


def get_forest(mol, new_structure=False) -> dict:

    forest = getattr(mol, "_forest", None)
    if forest is None:
        forest = {"structure": {}}
        mol._forest = forest

    if new_structure:
        forest = copy.deepcopy(forest)
        del forest["structure"]

    return forest


def sanitize_mol(mol):
    forest = get_forest(mol)

    if "sanitized" in forest["structure"]:
        return forest["structure"]["sanitized"]

    sanitized = SanitizeMol(Mol(mol), catchErrors=True)

    forest["structure"]["sanitized"] = sanitized
    return sanitized


def topol_equiv(mol):
    """Get dict mapping atom indexes to topological IDs.

    >>> from rdkit import Chem
    >>> mol = MolFromSmiles('CCC')
    >>> AtomTracker.topol_equiv(mol)
    {0: 0, 1: 2, 2: 0}
    """

    forest = get_forest(mol)

    if "topol_equiv" in forest["structure"]:
        return forest["structure"]["topol_equiv"]

    sanitized = Mol(mol)
    sanitize_mol(sanitized)

    topol_equiv = {
        a.GetIdx(): i
        for a, i in zip(
            mol.GetAtoms(),
            CanonicalRankAtoms(sanitized, includeChirality=False, breakTies=False),
        )
    }
    forest["structure"]["topol_equiv"] = topol_equiv
    return topol_equiv


def set_terminal_product(mol, value=True):
    forest = get_forest(mol)
    forest["is_terminal_product"] = value
    return mol


class ProductsOfReaction(NamedTuple):
    info: dict[str, Any]
    products: list[Mol]


class ReactionRule:
    """Template for all rules."""

    is_terminal_rule = False

    def _clear_atom_maps(self, mol):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    def __init__(
        self,
        name=None,
        sites_on=None,
        longname=None,
        *args,
        **kwargs,
    ):
        # super().__init__( *args, **kwargs)

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

    def __call__(self, mol, **kwargs):
        yield from self.metabolize(mol, **kwargs)

    def __iter__(self):
        return iter([self])

    def metabolize(
        self,
        mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        unique_csmi=True,
        **kwargs,
    ) -> Generator[ProductsOfReaction, None, None]:
        assert mol is not None
        install_forest(mol)

        if self.is_terminal_product(mol):
            return

        # Maps must not be present for CanonicalRankAtoms or SMARTS matching.
        self._clear_atom_maps(mol)

        seen = set()

        for por in self.metabolites(
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            **kwargs,
        ):
            info = por.info
            products = por.products
            # print("INFO", info, len(products))
            site = info["site"]
            top_site = self._top_site(site, mol)

            sig = (top_site, tuple(info["rule"]))

            if sig in seen:
                continue

            seen.add(sig)

            for n, p in enumerate(products):
                sanitize_mol(p)

                if self.is_terminal_rule:
                    set_terminal_product(p)

                forest_trace(mol, p, info["rule"], site)

                p, csmi = cannonicalize_order(p)

                if unique_csmi:
                    if csmi in seen:
                        continue
                    seen.add(csmi)

                i = info.copy()
                i["product_index"] = n
                i["product_count"] = len(products)
                i["csmi"] = csmi

                yield p, i

    def _top_site(self, site, mol):
        te = topol_equiv(mol)
        if type(site) == int:
            return te[site]
        else:
            return frozenset(int(te[s]) for s in site)

    def is_terminal_product(self, mol) -> bool:
        """True if ``mol`` must not be expanded further in guided path search."""

        forest = get_forest(mol)
        if "is_terminal_product" in forest["structure"]:
            return forest["structure"]["is_terminal_product"]

        return False

    def metabolites(
        self,
        mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        **kwargs,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Should return a tuple of lists. The first element will be the site, the second element
        will be a list of metabolites. Any rule that could cause molecule fragmention should
        separate the fragments within the outputted list of metabolites. One way to do this is to
        return clean(product), which will separate all fragments and sanitize each one.
        """
        raise NotImplementedError


def install_forest(mol):
    """Install the forest data structure in the molecule."""
    forest = mol._forest = get_forest(mol)

    forest["structure"] = forest.get("structure", {})

    if "atom_trace" not in forest:
        forest["atom_trace"] = {
            "records": {},
            "deletes": {},
            "transforms": [],
            "depth": 0,
            "last_tag": 0,
        }

        for a in mol.GetAtoms():
            if a.GetAtomicNum() != 1:
                i = a.GetIdx()
                forest["atom_trace"]["records"][str(i)] = {
                    "idx": [i],
                    "depth": [0],
                }

                forest["atom_trace"]["last_tag"] = i

    stamp_forest_labels(mol)

    return forest


def stamp_forest_labels(mol):
    forest = get_forest(mol)
    if "atom_trace" not in forest:
        install_forest(mol)

    for tag, record in forest["atom_trace"]["records"].items():
        i = record["idx"][-1]
        atom = mol.GetAtomWithIdx(i)
        atom.SetProp("forestLabel", tag)

    return mol


import ast


def cannonicalize_order(mol, tracing_reset=True):
    csmi = MolToSmiles(mol, isomericSmiles=False)

    smiles_order = ast.literal_eval(mol.GetProp("_smilesAtomOutputOrder"))

    renumber_map = [0] * mol.GetNumAtoms()
    for new_pos, old_idx in enumerate(smiles_order):
        renumber_map[old_idx] = new_pos

    r = RenumberAtoms(mol, renumber_map)

    forest = get_forest(mol)

    r._forest = copy.deepcopy(forest)
    r._forest["structure"] = {"csmi": csmi}

    if tracing_reset:
        reordered_forest_labels(r)

    return r, csmi


def get_csmi(mol):
    forest = get_forest(mol)
    csmi = forest["structure"].get("csmi", None)

    if not csmi:
        csmi = MolToSmiles(mol, isomericSmiles=False)
        forest["structure"]["csmi"] = csmi
    return csmi


def reordered_forest_labels(mol):
    forest = get_forest(mol)
    # if "atom_trace" not in forest:
    #     install_forest(mol)

    for atom in mol.GetAtoms():
        i = atom.GetIdx()

        if atom.GetAtomicNum() != 1:
            assert atom.HasProp("forestLabel")
            tag = atom.GetProp("forestLabel")
            record = forest["atom_trace"]["records"][tag]
            record["idx"][-1] = i


def forest_trace(reactant, product, rule, site):
    """Trace the forest of the molecule."""

    stamp_forest_labels(reactant)
    forest = get_forest(product, new_structure=True)
    parent_forest = get_forest(reactant)

    if "atom_trace" not in parent_forest:
        install_forest(reactant)

    trace = forest["parent_atom_trace"] = copy.deepcopy(parent_forest["atom_trace"])

    depth = trace["depth"] = trace["depth"] + 1

    records = trace["records"]
    new_records = {}

    for product_atom in product.GetAtoms():
        product_idx = product_atom.GetIdx()
        if product_atom.HasProp("forestLabel"):
            tag = product_atom.GetProp("forestLabel")
            new_records[tag] = records.pop(tag)
            new_records[tag]["idx"].append(product_idx)
            new_records[tag]["depth"].append(depth)

        else:
            tag = trace["last_tag"] = trace["last_tag"] + 1
            new_records[str(tag)] = {
                "idx": [product_idx],
                "depth": [depth],
                "added_by": {"rule": rule, "site": site, "depth": depth - 1},
            }
            product_atom.SetProp("forestLabel", str(tag))

    for tag, record in records.items():
        record["removed_by"] = (rule, site)
        trace["deletes"][tag] = record

    trace["records"] = new_records

    product._forest["atom_trace"] = trace

    return trace


class When(TypedDict, total=False):
    """Constraint that picks one branch of a SMARTS OR once atoms are known."""

    map: int
    z: int
    h: int


class Effect(TypedDict, total=False):
    """One concrete outcome.

    Declared on a possibility, then filled in from the matched atoms:
    ``symbol`` / ``h`` / ``site_aromatic`` are the site atom;
    ``partner`` / ``partner_h`` are the atom the SMARTS OR was ambiguous about.
    """

    adds: str
    removes: str
    cleaves: bool
    dearomatizes: bool
    methide: bool
    needs: str
    partner: str
    partner_h: int
    symbol: str
    h: int
    site_aromatic: bool
    when: When


class PatternInfo(TypedDict, total=False):
    """What a SMARTS pattern can do, and how to apply it.

    ``filter_rules(rule, info)`` sees this whole object, before any match.
    ``possibilities`` is every branch the SMARTS OR allows. ``span`` collapses
    that list: a bare value means every branch agrees, a tuple means the
    site has not chosen yet.

    ``filter_sites(site, info)`` sees one resolved :class:`Effect` in
    ``info["options"]``. Pair rules also set ``info["ends"]`` to the effect
    at each end, so a merged ``adds`` / ``removes`` string is not the only
    record of which end did what.

    ``edit`` / ``site_map`` / ``skip_same_rings`` are mechanism, not chemistry.
    """

    possibilities: tuple[Effect, ...]
    span: dict[str, Any]
    edit: str
    site_map: int
    skip_same_rings: bool


_EFFECT_DEFAULTS = {
    "adds": "",
    "removes": "",
    "cleaves": False,
    "dearomatizes": False,
    "methide": False,
    "needs": "",
}

_SYMBOL = {
    1: "H",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    16: "S",
    17: "Cl",
    35: "Br",
    53: "I",
}


def _span(possibilities):
    """Certain values stay bare. Disagreeing values become a tuple."""

    keys = []
    for possibility in possibilities:
        for key in possibility:
            if key != "when" and key not in keys:
                keys.append(key)
    span = {}
    for key in keys:
        values = []
        for possibility in possibilities:
            value = possibility.get(key, _EFFECT_DEFAULTS.get(key))
            if value not in values:
                values.append(value)
        span[key] = values[0] if len(values) == 1 else tuple(values)
    return span


def branches(whens, site_map=1, removes_partner=False, **effect):
    """Copy ``effect`` once per ``when``. The site atom and the OR atom differ.

    A ``when`` on ``site_map`` records ``h`` / ``symbol`` (the site was the
    ambiguous atom). A ``when`` on any other map records ``partner`` /
    ``partner_h``.
    """

    out = []
    for when in whens:
        when = dict(when)
        item = dict(effect)
        item["when"] = when
        symbol = _SYMBOL.get(when.get("z", -1))
        if when.get("map") == site_map:
            if "h" in when:
                item["h"] = when["h"]
            if symbol:
                item["symbol"] = symbol
        else:
            if symbol:
                item["partner"] = symbol
                if removes_partner:
                    item["removes"] = symbol
            if "h" in when:
                item["partner_h"] = when["h"]
        out.append(item)
    return tuple(out)


def describe(*possibilities, edit=None, site_map=1, skip_same_rings=False, **single):
    """Build a :class:`PatternInfo`.

    One outcome: ``describe(adds="O", removes="H")``.
    Several: ``describe(*branches(...), edit="dealkylate")``.
    """

    if single and possibilities:
        raise TypeError("pass one effect as keywords, or several dicts, not both")
    if single:
        possibilities = (single,)
    if not possibilities:
        raise TypeError("a pattern needs at least one possibility")
    normalized = []
    for possibility in possibilities:
        effect = dict(_EFFECT_DEFAULTS)
        effect.update(possibility)
        normalized.append(effect)
    info = {
        "possibilities": tuple(normalized),
        "span": _span(normalized),
        "site_map": site_map,
    }
    if edit is not None:
        info["edit"] = edit
    if skip_same_rings:
        info["skip_same_rings"] = True
    return info


def may(info, key, value=True):
    """True if any possibility has this outcome. For ``adds`` / ``removes`` / ``needs``, ``value`` may be a substring."""

    for possibility in info["possibilities"]:
        have = possibility.get(key, _EFFECT_DEFAULTS.get(key))
        if isinstance(value, str) and isinstance(have, str) and key in (
            "adds",
            "removes",
            "needs",
        ):
            if value in have:
                return True
        elif have == value:
            return True
    return False


def must(info, key, value=True):
    """True if every possibility has this outcome."""

    possibilities = info["possibilities"]
    if not possibilities:
        return False
    return all(
        (
            isinstance(value, str)
            and isinstance(possibility.get(key, ""), str)
            and key in ("adds", "removes", "needs")
            and value in possibility.get(key, "")
        )
        or possibility.get(key, _EFFECT_DEFAULTS.get(key)) == value
        for possibility in possibilities
    )


def _when_matches(mol, mapped, when):
    idx = mapped.get(when.get("map"))
    if idx is None:
        return False
    atom = mol.GetAtomWithIdx(idx)
    if "z" in when and atom.GetAtomicNum() != when["z"]:
        return False
    if "h" in when and atom.GetTotalNumHs() != when["h"]:
        return False
    return True


def resolve_effect(mol, mapped, info):
    """Narrow ``info`` to the one effect this match actually is.

    Uses the caller's mol, not a kekulé copy: aromatic flags must still be set.
    """

    possibilities = info.get("possibilities") or (info,)
    chosen = None
    for possibility in possibilities:
        when = possibility.get("when")
        if when is None or _when_matches(mol, mapped, when):
            chosen = possibility
            break
    if chosen is None:
        chosen = possibilities[0]
    effect = dict(_EFFECT_DEFAULTS)
    effect.update(chosen)
    effect.pop("when", None)
    site_map = info.get("site_map", 1)
    if site_map in mapped:
        atom = mol.GetAtomWithIdx(mapped[site_map])
        effect["symbol"] = atom.GetSymbol()
        effect["h"] = atom.GetTotalNumHs()
        effect["site_aromatic"] = atom.GetIsAromatic()
    when = chosen.get("when") or {}
    branch = when.get("map")
    if branch in mapped and branch != site_map:
        partner = mol.GetAtomWithIdx(mapped[branch])
        effect["partner"] = partner.GetSymbol()
        effect["partner_h"] = partner.GetTotalNumHs()
    return effect


def merge_effects(left, right, both_aromatic):
    """One effect for a pair. Per-end detail stays on ``info["ends"]``."""

    needs = (left.get("needs") or "") + (right.get("needs") or "")
    can = left.get("dearomatizes") or right.get("dearomatizes")
    return {
        "adds": left.get("adds", "") + right.get("adds", ""),
        "removes": left.get("removes", "") + right.get("removes", ""),
        "cleaves": bool(left.get("cleaves") or right.get("cleaves")),
        "dearomatizes": bool(can and both_aromatic),
        "methide": bool(left.get("methide")) ^ bool(right.get("methide")),
        "needs": needs,
        "partners": (left.get("partner", ""), right.get("partner", "")),
        "aromatic": bool(both_aromatic),
    }


class SmartsReactionRule(ReactionRule):
    """Performs reactions specified by SMARTS.

    Each entry is ``(smarts, options)`` with :class:`PatternInfo` options.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = ()

    rxns: list[tuple[Any, dict[str, Any]]]

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.rxns = [(self._smarts2rxns(s, **kwargs), opt) for s, opt in self.smarts]

    def metabolites(
        self,
        mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        context_mol=None,
        **kwargs,
    ):
        """Run this rule's SMARTS reactions.

        ``context_mol`` is the unsubstituted parent when ``mol`` is a kekulé
        copy. Resolution reads aromatic flags from it.
        """

        context = mol if context_mol is None else context_mol

        for rxn_num, (rxn, pattern) in enumerate(self.rxns):
            if not filter_rules(self, pattern):
                continue

            try:
                reactant_products = rxn.RunReactants((mol,))  # pyright: ignore[reportAttributeAccessIssue]
            except RuntimeError:
                continue

            for prod_num, prod in enumerate(reactant_products):
                products = [self._lift_forest_labels(mol, p) for p in prod]
                site = self._get_site(products)
                mapped = {}
                for product in products:
                    mapno2idx, _ = self._get_product_mappings(product)
                    mapped.update(mapno2idx)
                info = {
                    "site": site,
                    "rule": self,
                    "options": resolve_effect(context, mapped, pattern),
                    "rxn_num": rxn_num,
                }
                if filter_sites(site, info):
                    yield ProductsOfReaction(
                        info=info,
                        products=products,
                    )

    def _smarts2rxns(self, smarts, use_implicit_properties=False, **kwargs):
        """Converts SMARTS reactions to RDKit reactions."""
        rxn = AllChem.ReactionFromSmarts(smarts)  # pyright: ignore[reportAttributeAccessIssue]
        if not use_implicit_properties:
            rxn._setImplicitPropertiesFlag(False)
        return rxn

    def _lift_forest_labels(self, reactant, product):
        """Carry the forest labels from the reactant to the product."""

        for a in product.GetAtoms():
            if a.HasProp("react_atom_idx"):
                react_idx = int(a.GetProp("react_atom_idx"))
                reactant_atom = reactant.GetAtomWithIdx(react_idx)
                if reactant_atom.HasProp("forestLabel"):
                    a.SetProp("forestLabel", reactant_atom.GetProp("forestLabel"))

        return product

    def _get_product_mappings(self, product):
        mapno2idx = {}
        reactant2idx = {}
        for a in product.GetAtoms():
            prod_idx = a.GetIdx()
            react_idx = (
                int(a.GetProp("react_atom_idx"))
                if a.HasProp("react_atom_idx")
                else None
            )
            mapno = int(a.GetProp("old_mapno")) if a.HasProp("old_mapno") else None

            if react_idx is not None and mapno is not None:
                mapno2idx[mapno] = react_idx

            if react_idx is not None:
                reactant2idx[react_idx] = prod_idx

        return mapno2idx, reactant2idx

    def _get_site(self, products):
        """Return the sites in product based on the reactant_idx property assigned by
        rxns.RunReactants in self.metabolites."""

        site = set()
        for product in products:
            mapno2idx, _ = self._get_product_mappings(product)
            site = site | set(mapno2idx.values())
        return frozenset(site)


# Bond order stored in the structure cache. Not RDKit mols: those do not
# belong in a dict that is deep-copied onto every product.
_BOND = {
    1: BondType.SINGLE,
    1.0: BondType.SINGLE,
    2: BondType.DOUBLE,
    2.0: BondType.DOUBLE,
    3: BondType.TRIPLE,
    3.0: BondType.TRIPLE,
    1.5: BondType.AROMATIC,
}


def _bond_key(i, j):
    return (i, j) if i < j else (j, i)


def _current_bond_map(mol):
    bonds = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[_bond_key(i, j)] = bond.GetBondTypeAsDouble()
    return bonds


def _connected_components(mol, atoms):
    atoms = set(atoms)
    seen = set()
    systems = []
    for start in atoms:
        if start in seen:
            continue
        comp = set()
        queue = deque([start])
        while queue:
            i = queue.popleft()
            if i in comp:
                continue
            comp.add(i)
            for nbr in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nbr.GetIdx()
                if j in atoms and j not in comp:
                    queue.append(j)
        seen |= comp
        if len(comp) >= 2:
            systems.append(frozenset(comp))
    return systems


def _load_resonance(mol):
    """Cache kekulé bond maps and conjugated-atom sets on ``_forest['structure']``."""

    structure = get_forest(mol)["structure"]
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
        groups = [frozenset(v) for v in grouped.values() if len(v) >= 2]
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
        aromatic = {a.GetIdx() for a in mol.GetAromaticAtoms()}
        groups = _connected_components(mol, aromatic)
        if not groups:
            conjugated = set()
            for bond in base.GetBonds():
                if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() >= 1.5:
                    conjugated.add(bond.GetBeginAtomIdx())
                    conjugated.add(bond.GetEndAtomIdx())
            groups = _connected_components(base, conjugated)

    if not maps:
        maps = [_current_bond_map(base)]

    structure["resonance_bonds"] = tuple(maps)
    structure["conjugated_systems"] = tuple(groups)


def resonance_bond_maps(mol):
    _load_resonance(mol)
    return get_forest(mol)["structure"]["resonance_bonds"]


def conjugated_systems(mol):
    _load_resonance(mol)
    return get_forest(mol)["structure"]["conjugated_systems"]


def aromatic_systems(mol):
    structure = get_forest(mol)["structure"]
    if "aromatic_systems" not in structure:
        aromatic = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
        structure["aromatic_systems"] = tuple(_connected_components(mol, aromatic))
    return structure["aromatic_systems"]


def system_neighbors(mol, system):
    system = set(system)
    neighbors = {i: [] for i in system}
    for i in system:
        for nbr in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nbr.GetIdx()
            if j in system and j > i:
                neighbors[i].append(j)
                neighbors[j].append(i)
    return neighbors


def odd_anchor_pairs(anchors, neighbors):
    """Pairs of anchors separated by an odd number of bonds.

    Ortho (1) and para (3) pass. Meta (2) does not. The walk may cross
    non-anchor atoms; ``neighbors`` is the whole system.
    """

    anchors = [a for a in anchors if a in neighbors]
    pairs = []
    for i, start in enumerate(anchors):
        dist = {start: 0}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for nbr in neighbors[node]:
                if nbr not in dist:
                    dist[nbr] = dist[node] + 1
                    queue.append(nbr)
        for end in anchors[i + 1 :]:
            d = dist.get(end)
            if d is not None and d % 2 == 1:
                pairs.append((start, end))
    return pairs


def alternating_path(bond_map, start, end, neighbors):
    """Shortest path whose first bond is double and whose bonds then alternate."""

    if start == end or start not in neighbors or end not in neighbors:
        return None
    queue = deque([(start, 2.0, (start,))])
    seen = {(start, 2)}
    while queue:
        node, want, path = queue.popleft()
        next_want = 1.0 if int(want) == 2 else 2.0
        for nbr in neighbors[node]:
            if nbr in path:
                continue
            order = bond_map.get(_bond_key(node, nbr))
            if order is None or int(order) != int(want):
                continue
            nxt = path + (nbr,)
            if nbr == end:
                return list(nxt)
            state = (nbr, int(next_want))
            if state in seen:
                continue
            seen.add(state)
            queue.append((nbr, next_want, nxt))
    return None


def overlay_kekule(mol, bond_map):
    """Copy ``mol`` and set kekulé bond orders. Atom indexes stay put."""

    rw = RWMol(Mol(mol))
    for bond in rw.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = bond_map.get(_bond_key(i, j))
        if order is None:
            continue
        bond.SetBondType(_BOND.get(order, BondType.SINGLE))
        bond.SetIsAromatic(False)
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)
    return rw


def adjust_hydrogens(mol, idx, change):
    atom = mol.GetAtomWithIdx(idx)
    try:
        implicit = atom.GetNumImplicitHs()
    except RuntimeError:
        atom.UpdatePropertyCache(strict=False)
        implicit = atom.GetNumImplicitHs()
    total = atom.GetNumExplicitHs() + implicit
    updated = total + change
    if total > 0 and updated >= 0:
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(updated)


def swap_bonds_along_path(mol, atoms):
    """Flip single and double bonds along ``atoms``. Adjust H at the ends."""

    bond = None
    for i in range(len(atoms) - 1):
        bond = mol.GetBondBetweenAtoms(atoms[i], atoms[i + 1])
        if bond is None:
            return False
        if bond.GetBondType() == BondType.DOUBLE:
            bond.SetBondType(BondType.SINGLE)
        elif bond.GetBondType() == BondType.SINGLE:
            bond.SetBondType(BondType.DOUBLE)
        if i == 0:
            _correct_end_hydrogens(mol, atoms[0], bond)
    if bond is not None:
        _correct_end_hydrogens(mol, atoms[-1], bond)
    return True


def _correct_end_hydrogens(mol, idx, bond):
    if bond.GetBondType() == BondType.SINGLE:
        adjust_hydrogens(mol, idx, 1)
    elif bond.GetBondType() == BondType.DOUBLE:
        adjust_hydrogens(mol, idx, -1)


def ring_membership(mol):
    structure = get_forest(mol)["structure"]
    if "rings" not in structure:
        work = Mol(mol)
        SanitizeMol(work, Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
        atom_rings = work.GetRingInfo().AtomRings()
        structure["rings"] = {
            idx: tuple(ring for ring in atom_rings if idx in ring)
            for idx in range(work.GetNumAtoms())
        }
    return structure["rings"]


def _same_rings(rings, i, j):
    return rings.get(i, ()) == rings.get(j, ())


def edit_single_to_double(rw, mapped, info, rings):
    a, b = mapped.get(1), mapped.get(2)
    if a is None or b is None:
        return False
    if info.get("skip_same_rings") and _same_rings(rings, a, b):
        return False
    bond = rw.GetBondBetweenAtoms(a, b)
    if bond is None:
        return False
    bond.SetBondType(BondType.DOUBLE)
    return True


def edit_add_carbonyl_o(rw, mapped, info, rings):
    carbon = mapped.get(1)
    if carbon is None:
        return False
    oxygen = rw.AddAtom(Atom(8))
    rw.AddBond(carbon, oxygen, BondType.DOUBLE)
    return True


def edit_replace_halogen(rw, mapped, info, rings):
    carbon, halogen = mapped.get(1), mapped.get(2)
    if carbon is None or halogen is None:
        return False
    atom = rw.GetAtomWithIdx(halogen)
    atom.SetAtomicNum(8)
    atom.SetFormalCharge(0)
    bond = rw.GetBondBetweenAtoms(carbon, halogen)
    if bond is None:
        return False
    bond.SetBondType(BondType.DOUBLE)
    return True


def edit_iminium(rw, mapped, info, rings):
    if not edit_single_to_double(rw, mapped, {}, rings):
        return False
    rw.GetAtomWithIdx(mapped[2]).SetFormalCharge(1)
    return True


def edit_dealkylate(rw, mapped, info, rings):
    hetero, alkyl = mapped.get(2), mapped.get(3)
    if hetero is None or alkyl is None:
        return False
    if not edit_single_to_double(rw, mapped, info, rings):
        return False
    rw.GetAtomWithIdx(alkyl).SetProp("dealk-noncarbon", "1")
    if rw.GetBondBetweenAtoms(hetero, alkyl) is None:
        return False
    rw.RemoveBond(hetero, alkyl)
    adjust_hydrogens(rw, hetero, 1)
    adjust_hydrogens(rw, alkyl, 1)
    return True


EDITS = {
    "single_to_double": edit_single_to_double,
    "add_carbonyl_o": edit_add_carbonyl_o,
    "replace_halogen": edit_replace_halogen,
    "iminium": edit_iminium,
    "dealkylate": edit_dealkylate,
}


def smarts_matches(mol, smarts):
    cache = get_forest(mol)["structure"].setdefault("smarts_matches", {})
    if smarts not in cache:
        query = Chem.MolFromSmarts(smarts)
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


def sanitized_fragments(mol):
    """Split, drop the dealkylation leaving group, sanitize. Empty if any frag fails."""

    if isinstance(mol, RWMol):
        mol = mol.GetMol()
    frags = list(GetMolFrags(mol, asMols=True, sanitizeFrags=False)) or [mol]
    out = []
    for frag in frags:
        if any(atom.HasProp("dealk-noncarbon") for atom in frag.GetAtoms()):
            continue
        if SanitizeMol(frag, catchErrors=True):
            return []
        out.append(frag)
    return out


def _merge_options(left, right, dearomatizes):
    return merge_effects(left, right, dearomatizes)


def _site_atoms(mapped, info):
    key = info.get("site_map", 1)
    if key not in mapped:
        return None
    return mapped[key]


def _paths_for_pair(bond_maps, start, end, neighbors):
    found = []
    for bond_map in bond_maps:
        path = alternating_path(bond_map, start, end, neighbors)
        if path:
            found.append((bond_map, path))
    found.sort(key=lambda item: len(item[1]))
    return found


def pair_metabolites(rule, mol, filter_rules, filter_sites):
    """Shared ResonancePairRule loop.

    ``filter_rules`` drops endpoint patterns before they are matched.
    ``filter_sites`` drops a pair before any kekulé overlay or bond edit.
    """

    active = [
        (smarts, info)
        for smarts, info in rule.endpoints
        if filter_rules(rule, info)
    ]
    if not active:
        return

    if rule.systems == "aromatic":
        systems = aromatic_systems(mol)
    else:
        systems = conjugated_systems(mol)
    if not systems:
        return

    hits = defaultdict(list)
    for smarts, info in active:
        for mapped in smarts_matches(mol, smarts):
            hits[mapped[1]].append((mapped, info))
    if len(hits) < 2:
        return

    bond_maps = None
    rings = None
    for system in systems:
        anchors = [atom for atom in hits if atom in system]
        neighbors = system_neighbors(mol, system)
        for start, end in odd_anchor_pairs(anchors, neighbors):
            combos = []
            both_aromatic = (
                mol.GetAtomWithIdx(start).GetIsAromatic()
                and mol.GetAtomWithIdx(end).GetIsAromatic()
            )
            for (map1, info1), (map2, info2) in itertools.product(
                hits[start], hits[end]
            ):
                site_a = _site_atoms(map1, info1)
                site_b = _site_atoms(map2, info2)
                if site_a is None or site_b is None or site_a == site_b:
                    continue
                site = frozenset((site_a, site_b))
                end1 = resolve_effect(mol, map1, info1)
                end2 = resolve_effect(mol, map2, info2)
                preview = {
                    "site": site,
                    "rule": rule,
                    "options": merge_effects(end1, end2, both_aromatic),
                    "ends": (end1, end2),
                    "path_ends": frozenset((start, end)),
                }
                if filter_sites(site, preview):
                    combos.append((map1, info1, map2, info2, preview))
            if not combos:
                continue
            if bond_maps is None:
                bond_maps = resonance_bond_maps(mol)
            paths = _paths_for_pair(bond_maps, start, end, neighbors)
            if not paths:
                continue
            for map1, info1, map2, info2, preview in combos:
                if rings is None and (
                    info1.get("skip_same_rings") or info2.get("skip_same_rings")
                ):
                    rings = ring_membership(mol)
                ring_table = rings or {}
                for bond_map, path in paths:
                    rw = overlay_kekule(mol, bond_map)
                    edit1 = EDITS.get(info1.get("edit", ""))
                    edit2 = EDITS.get(info2.get("edit", ""))
                    if edit1 is None or edit2 is None:
                        break
                    if not edit1(rw, map1, info1, ring_table):
                        break
                    if not edit2(rw, map2, info2, ring_table):
                        break
                    if not swap_bonds_along_path(rw, path):
                        continue
                    products = sanitized_fragments(rw)
                    if not products:
                        continue
                    info = dict(preview)
                    info["path"] = tuple(path)
                    yield ProductsOfReaction(info=info, products=products)
                    break


class ResonanceRule(SmartsReactionRule):
    """Run this rule's SMARTS on each kekulé form.

    Epoxidation is the rule this exists for: a double bond moves between
    kekulé forms, so the SMARTS has to see each of them. The forms are bond
    maps in ``_forest['structure']``, not copies kept on the class.
    """

    def metabolites(
        self,
        mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        context_mol=None,
        **kwargs,
    ):
        if not self.rxns:
            return
        context = mol if context_mol is None else context_mol
        for bond_map in resonance_bond_maps(mol):
            form = overlay_kekule(mol, bond_map).GetMol()
            yield from SmartsReactionRule.metabolites(
                self,
                form,
                filter_rules=filter_rules,
                filter_sites=filter_sites,
                context_mol=context,
                **kwargs,
            )


class ResonancePairRule(ResonanceRule):
    """Two endpoint matches joined by an alternating path.

    Subclasses set ``endpoints`` to ``(smarts, PatternInfo)`` pairs and
    ``systems`` to ``"conjugated"`` or ``"aromatic"``. ``info["edit"]`` names
    a module-level edit applied at that end; the path's bonds are then flipped.

    Dehydrogenation and quinone formation both use this. Hydroxylation does
    not: it is a one-atom SMARTS reaction.
    """

    endpoints: tuple[tuple[str, PatternInfo], ...] = ()
    systems = "conjugated"

    def metabolites(
        self,
        mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        **kwargs,
    ):
        yield from ResonanceRule.metabolites(
            self,
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            **kwargs,
        )
        yield from pair_metabolites(self, mol, filter_rules, filter_sites)


class Hydroxylation(SmartsReactionRule):
    """Adds a hydroxyl to carbon.

    Both patterns add OH and remove one H. ``[#6h]`` is h=1, 2, or 3;
    ``[#6h2]`` is the subset with at least two hydrogens. The match records
    which of those the atom actually is.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h:1]>>[*:1]O",
            describe(
                *branches(
                    (
                        {"map": 1, "z": 6, "h": 1},
                        {"map": 1, "z": 6, "h": 2},
                        {"map": 1, "z": 6, "h": 3},
                    ),
                    adds="O",
                    removes="H",
                )
            ),
        ),
        (
            "[#6h2:1]>>[*:1]O",
            describe(
                *branches(
                    ({"map": 1, "z": 6, "h": 2}, {"map": 1, "z": 6, "h": 3}),
                    adds="O",
                    removes="H",
                )
            ),
        ),
    )


class Dehydrogenation(ResonancePairRule):
    """Drop H2 from one bond, or across a conjugated path.

    The path case turns a hydroquinone into a quinone: each end is already
    an OH or NH, and the bonds between them flip. Two hydroxylations of
    benzene, then this rule, is benzoquinone.

    ``filter_rules`` sees ``span`` (what is still ambiguous). ``filter_sites``
    sees the branch the two atoms selected, including whether the path
    actually dearomatizes.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6h:1]-[#8H1:2]>>[*:1]=[*:2]",
            describe(removes="HH", partner="O"),
        ),
        (
            "[#6h:1]-[#7D1H2,#7D2H1:2]>>[*:1]=[*:2]",
            describe(
                *branches(
                    ({"map": 2, "z": 7, "h": 2}, {"map": 2, "z": 7, "h": 1}),
                    removes="HH",
                )
            ),
        ),
        (
            "[#6h:1]-[#6D1H3,#6D2H2,#6D3H1:2]>>[*:1]=[*:2]",
            describe(
                *branches(
                    (
                        {"map": 2, "z": 6, "h": 3},
                        {"map": 2, "z": 6, "h": 2},
                        {"map": 2, "z": 6, "h": 1},
                    ),
                    removes="HH",
                )
            ),
        ),
    )

    endpoints: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]-[#8H:2]",
            describe(
                removes="H",
                partner="O",
                dearomatizes=True,
                edit="single_to_double",
                site_map=2,
            ),
        ),
        (
            "[#6:1]-[#7D1H2,#7D2H1:2]",
            describe(
                *branches(
                    ({"map": 2, "z": 7, "h": 2}, {"map": 2, "z": 7, "h": 1}),
                    site_map=2,
                    removes="H",
                    dearomatizes=True,
                ),
                edit="single_to_double",
                site_map=2,
            ),
        ),
    )


class QuinoneFormation(ResonancePairRule):
    """Dearomatize an aromatic pair into a quinone, imine, or methide.

    A bare aromatic CH can gain its carbonyl oxygen here (``needs`` ``"O"``).
    The same quinone is hydroxylation of that carbon, then dehydrogenation.
    Use that split when a search should hydroxylate only atoms that still
    need oxygen, and dearomatize only the ring.

    The exocyclic single-to-double SMARTS is one pattern with several
    partners (O, N, alkyl C). ``span['partner']`` is that tuple until the
    match; the resolved effect names the partner this site actually has.
    """

    systems = "aromatic"

    endpoints: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6R:1]-[#8H,#7D1H2,#7D2H1,#6D1H3,#6D2H2,#6D3H1:2]",
            describe(
                *branches(
                    (
                        {"map": 2, "z": 8},
                        {"map": 2, "z": 7, "h": 2},
                        {"map": 2, "z": 7, "h": 1},
                    ),
                    removes="H",
                    dearomatizes=True,
                ),
                *branches(
                    (
                        {"map": 2, "z": 6, "h": 3},
                        {"map": 2, "z": 6, "h": 2},
                        {"map": 2, "z": 6, "h": 1},
                    ),
                    removes="H",
                    dearomatizes=True,
                    methide=True,
                ),
                edit="single_to_double",
                site_map=1,
                skip_same_rings=True,
            ),
        ),
        (
            "[#6D2H1;R:1]",
            describe(
                adds="O",
                removes="H",
                dearomatizes=True,
                needs="O",
                edit="add_carbonyl_o",
                site_map=1,
            ),
        ),
        (
            "[#6H0R:1]-[F,Cl,Br,I:2]",
            describe(
                *branches(
                    (
                        {"map": 2, "z": 9},
                        {"map": 2, "z": 17},
                        {"map": 2, "z": 35},
                        {"map": 2, "z": 53},
                    ),
                    adds="O",
                    dearomatizes=True,
                    removes_partner=True,
                ),
                edit="replace_halogen",
                site_map=1,
            ),
        ),
        (
            "[#6H0R:1]-[#7D3:2]",
            describe(
                partner="N",
                dearomatizes=True,
                edit="iminium",
                site_map=1,
                skip_same_rings=True,
            ),
        ),
        (
            "[#6R:1]-[#7,#8:2]-[#6:3]",
            describe(
                *branches(
                    ({"map": 2, "z": 7}, {"map": 2, "z": 8}),
                    cleaves=True,
                    dearomatizes=True,
                ),
                edit="dealkylate",
                site_map=1,
                skip_same_rings=True,
            ),
        ),
    )

