"""Reactant-to-target search over the proof-of-concept rules.

One atom diff drives ``filter_rules`` and ``filter_sites``. A rule is a
SMARTS pattern plus ``describe()``; it does not grow this search. The yield
is a phase-I :class:`~xenosite.forest.step_plan.Deps` plus the cleavage
fragments that were not expanded.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import GetMolFrags, SanitizeMol, rdFMCS

from xenosite.forest.step_plan import AtomRef, Deps, Step
from xenosite.refactor_poc.rules import (
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
    RuleSet,
    cannonicalize_order,
    forest_trace,
    install_forest,
)


# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------


class PathCounters:
    """Billed work for one search.

    ``mol_edits`` counts a ``RunReactants`` or a kekulé overlay that a filter
    has already accepted. ``billed`` is edits plus nodes visited. It is a
    measurement, not a second search mode.
    """

    def __init__(self):
        self.rule_expansions = 0
        self.sites_considered = 0
        self.sites_skipped = 0
        self.mol_edits = 0
        self.sanitize_dropped = 0
        self.nodes = 0

    @property
    def billed(self):
        return self.mol_edits + self.nodes


# ---------------------------------------------------------------------------
# Outcome
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CleavageSide:
    """Fragment discarded at a cleavage. Not searched.

    ``site`` is the cleavage site on the parent. ``opens`` are earlier
    ring-open sites on the same walk, oldest first. ``side`` is the
    canonical SMILES of the discarded fragment.
    """

    site: frozenset
    side: str
    opens: tuple = ()

    def span_sites(self):
        return tuple(self.opens) + (self.site,)


def _site_keys(site):
    if site is None:
        return frozenset()
    keys = set()
    for item in site:
        origin = getattr(item, "origin", None)
        if origin is not None:
            keys.add(int(origin))
            continue
        try:
            keys.add(int(item))
        except (TypeError, ValueError):
            keys.add(str(item))
    return frozenset(keys)


@dataclass(frozen=True)
class Maybe:
    """Uncleared cleavage fragments. Not a step, and not searched."""

    entries: tuple = ()

    def __bool__(self):
        return bool(self.entries)

    def sides(self):
        return tuple(entry.side for entry in self.entries)

    def allows(self, rule_name=None, site=None, *, side=None):
        """True when ``side`` is a discarded fragment, or ``site`` overlaps one.

        The bifurcating cleavage site itself does not pass: that step is
        already in the plan. A different reaction on those atoms does.
        """

        if not self.entries:
            return False
        if side is not None:
            want = canon_smiles(side)
            return any(entry.side == want for entry in self.entries)
        if site is None:
            return True
        keys = _site_keys(site)
        if not keys:
            return False
        for entry in self.entries:
            if keys == _site_keys(entry.site):
                continue
            for span in entry.span_sites():
                if keys & _site_keys(span):
                    return True
        return False


@dataclass(frozen=True)
class PathOutcome:
    """Required phase-I plan, plus cleavage fragments left on :class:`Maybe`."""

    plan: Deps
    maybe: Maybe
    smiles: str

    def allows(self, rule_name=None, site=None, *, side=None):
        return self.maybe.allows(rule_name, site, side=side)


# ---------------------------------------------------------------------------
# Atom diff
# ---------------------------------------------------------------------------


def _as_mol(value):
    if isinstance(value, str):
        mol = Chem.MolFromSmiles(value)
        if mol is None:
            raise ValueError("could not parse %r" % (value,))
        return mol
    return value


def _hydrogens(atom):
    try:
        return atom.GetTotalNumHs()
    except RuntimeError:
        return 0


def canon_smiles(value):
    mol = _as_mol(value)
    copy = Chem.Mol(mol)
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy, isomericSmiles=False)


class AtomDiff:
    """Local reactant-to-target differences. This is the whole heuristic.

    A mapped atom can need an element, lose aromaticity, or lose hydrogens.
    An unmapped reactant atom, or a mapped bond that is absent in the target,
    is cleavage. ``filter_rules`` sees the molecule-level summary. ``filter_sites``
    sees one site.
    """

    def __init__(
        self,
        reactant,
        target,
        mapping,
        needs_oxygen,
        needs_carbonyl,
        needs_alcohol,
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta,
        n_extra,
        bond_order_mismatches,
        bond_raises,
    ):
        self.reactant = reactant
        self.target = target
        self.mapping = mapping
        self.needs_oxygen = frozenset(needs_oxygen)
        self.needs_carbonyl = frozenset(needs_carbonyl)
        self.needs_alcohol = frozenset(needs_alcohol)
        self.cleaved = frozenset(cleaved)
        self.cleavage_bonds = set(cleavage_bonds)
        self.loses_aromaticity = frozenset(loses_aromaticity)
        self.h_delta = dict(h_delta)
        self.n_extra = n_extra
        self.bond_order_mismatches = bond_order_mismatches
        self.bond_raises = set(bond_raises)
        self.reactant_heavy = reactant.GetNumHeavyAtoms()
        self.target_heavy = target.GetNumHeavyAtoms()

    @property
    def target_smaller(self):
        return self.target_heavy < self.reactant_heavy

    @property
    def has_cleavage(self):
        return bool(self.cleaved or self.cleavage_bonds)

    @property
    def h_loss(self):
        return any(delta < 0 for delta in self.h_delta.values())

    def site_is_cleavage(self, atoms):
        """True when ``atoms`` is the bond that separates kept from gone.

        Touching a leaving-group atom is not enough. The site has to be
        the bridge itself, or one atom of that bridge.
        """

        atoms = set(atoms)
        if not atoms:
            return False
        for bond in self.cleavage_bonds:
            if bond <= atoms or atoms <= bond:
                return True
        return False

    def cost(self):
        """How much of the reactant still disagrees with the target.

        Used only to refuse a child that moved away. Not a ranker.
        """

        h_off = sum(1 for delta in self.h_delta.values() if delta)
        return (
            3 * len(self.cleaved)
            + 3 * self.n_extra
            + 2 * len(self.needs_oxygen)
            + len(self.loses_aromaticity)
            + h_off
            + 3 * len(self.cleavage_bonds)
            + self.bond_order_mismatches
        )


def _mapping_score(reactant, target, r_match, t_match):
    """Prefer a pairing that keeps rings on rings and bond orders put."""

    score = 0
    paired = {}
    for r_idx, t_idx in zip(r_match, t_match):
        paired[r_idx] = t_idx
        ra = reactant.GetAtomWithIdx(r_idx)
        ta = target.GetAtomWithIdx(t_idx)
        if ra.GetAtomicNum() != ta.GetAtomicNum():
            score -= 50
        if ra.GetIsAromatic() == ta.GetIsAromatic():
            score += 4
        else:
            score -= 12
        if ra.IsInRing() == ta.IsInRing():
            score += 3
        else:
            score -= 8
        if _hydrogens(ra) == _hydrogens(ta):
            score += 1
    for bond in reactant.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in paired or j not in paired:
            continue
        other = target.GetBondBetweenAtoms(paired[i], paired[j])
        if other is None:
            score -= 4
            continue
        if abs(other.GetBondTypeAsDouble() - bond.GetBondTypeAsDouble()) < 0.2:
            score += 1
        else:
            score -= 1
    return score


def _best_mapping(reactant, target):
    mcs = rdFMCS.FindMCS(
        [reactant, target],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        matchValences=False,
        ringMatchesRingOnly=False,
        completeRingsOnly=False,
        timeout=2,
    )
    if mcs.numAtoms <= 0 or mcs.canceled:
        return {}
    query = Chem.MolFromSmarts(mcs.smartsString)
    if query is None:
        return {}
    r_matches = reactant.GetSubstructMatches(query, uniquify=False)
    t_matches = target.GetSubstructMatches(query, uniquify=False)
    if not r_matches or not t_matches:
        return {}

    # Symmetric molecules can have dozens of embeddings. Score a bounded set.
    r_matches = r_matches[:24]
    t_matches = t_matches[:24]
    best = None
    best_score = None
    for r_match in r_matches:
        for t_match in t_matches:
            score = _mapping_score(reactant, target, r_match, t_match)
            if best_score is None or score > best_score:
                best_score = score
                best = (r_match, t_match)
    r_match, t_match = best
    return {r: t for r, t in zip(r_match, t_match)}


def atom_diff(reactant, target):
    """Pair reactant atoms with target atoms and record the local change.

    String inputs are parsed with ``MolFromSmiles``. Indexes then refer to
    that parse, which is stable for a given SMILES.
    """

    reactant = _as_mol(reactant)
    target = _as_mol(target)
    mapping = _best_mapping(reactant, target)
    image = set(mapping.values())

    needs_oxygen = set()
    needs_carbonyl = set()
    needs_alcohol = set()
    for atom in target.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in image:
            continue
        for neighbor in atom.GetNeighbors():
            bond = target.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond is None:
                continue
            for r_idx, t_idx in mapping.items():
                if t_idx != neighbor.GetIdx():
                    continue
                needs_oxygen.add(r_idx)
                if bond.GetBondTypeAsDouble() >= 1.5:
                    needs_carbonyl.add(r_idx)
                else:
                    needs_alcohol.add(r_idx)

    cleaved = set()
    for atom in reactant.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.GetIdx() not in mapping:
            cleaved.add(atom.GetIdx())

    cleavage_bonds = set()
    loses_aromaticity = set()
    h_delta = {}
    bond_order_mismatches = 0
    bond_raises = set()
    for r_idx, t_idx in mapping.items():
        ra = reactant.GetAtomWithIdx(r_idx)
        ta = target.GetAtomWithIdx(t_idx)
        if ra.GetIsAromatic() and not ta.GetIsAromatic():
            loses_aromaticity.add(r_idx)
        h_delta[r_idx] = _hydrogens(ta) - _hydrogens(ra)

    for bond in reactant.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in mapping or j not in mapping:
            # A cut that can reach the target joins a kept atom to a leaving
            # atom. Bonds inside the leaving group are not that cut.
            if (i in mapping) != (j in mapping):
                cleavage_bonds.add(frozenset((i, j)))
            continue
        other = target.GetBondBetweenAtoms(mapping[i], mapping[j])
        if other is None:
            cleavage_bonds.add(frozenset((i, j)))
            continue
        delta = other.GetBondTypeAsDouble() - bond.GetBondTypeAsDouble()
        if abs(delta) >= 0.2:
            bond_order_mismatches += 1
        if delta >= 0.2:
            bond_raises.add(frozenset((i, j)))

    n_extra = sum(
        1
        for atom in target.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetIdx() not in image
    )
    return AtomDiff(
        reactant,
        target,
        mapping,
        needs_oxygen,
        needs_carbonyl,
        needs_alcohol,
        cleaved,
        cleavage_bonds,
        loses_aromaticity,
        h_delta,
        n_extra,
        bond_order_mismatches,
        bond_raises,
    )


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------


def _span_values(span, key, default):
    if key not in span:
        return (default,)
    value = span[key]
    if isinstance(value, tuple):
        return value
    return (value,)


def _any_span(span, key, pred, default):
    return any(pred(value) for value in _span_values(span, key, default))


def _all_span(span, key, pred, default):
    values = _span_values(span, key, default)
    return bool(values) and all(pred(value) for value in values)


def _effect_adds_oxygen(effect):
    return "O" in (effect.get("adds") or "") or "O" in (effect.get("needs") or "")


def _pattern_could_help(info, diff):
    """``filter_rules`` sees ``span`` before any match."""

    span = info.get("span") or {}
    can_cleave = _any_span(span, "cleaves", bool, False)
    # The kept piece is smaller. A non-cleaving pattern cannot get there,
    # and editing the intact parent is the cost this filter exists to avoid.
    if diff.target_smaller and not can_cleave:
        return False
    if (
        _all_span(span, "adds", lambda value: "O" in (value or ""), "")
        and not can_cleave
        and not diff.needs_oxygen
    ):
        return False
    if _all_span(span, "cleaves", bool, False) and not diff.has_cleavage:
        return False
    if (
        _all_span(span, "dearomatizes", bool, False)
        and not diff.loses_aromaticity
    ):
        return False
    drops_h_only = (
        not can_cleave
        and not _any_span(span, "adds", lambda value: bool(value), "")
        and _any_span(span, "removes", lambda value: bool(value) and "H" in value, "")
    )
    if drops_h_only and not diff.h_loss and not diff.loses_aromaticity:
        return False
    return True


def _alkyl_bond_raises(mol, atom_idx, diff):
    """True when an exocyclic C-C bond at ``atom_idx`` is higher in the target."""

    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 6 or neighbor.GetIsAromatic():
            continue
        if frozenset((atom_idx, neighbor.GetIdx())) in diff.bond_raises:
            return True
    return False


def _site_could_help(site, info, diff, mol):
    """``filter_sites`` sees one resolved effect and the local atom diff."""

    effect = info.get("options") or {}
    atoms = set(site)
    if effect.get("cleaves"):
        return diff.site_is_cleavage(atoms)

    ends = info.get("ends")
    end_atoms = info.get("end_atoms")
    if ends and end_atoms and len(tuple(ends)) == len(tuple(end_atoms)):
        for atom, end in zip(end_atoms, ends):
            if _effect_adds_oxygen(end) and atom not in diff.needs_oxygen:
                return False
            # An alkyl partner turns the ring bond into an exocyclic double
            # bond (methide). Skip it unless that C-C bond is higher in the target.
            if (end.get("partner") or "") == "C" and not _alkyl_bond_raises(
                mol, atom, diff
            ):
                return False
    elif _effect_adds_oxygen(effect) and not effect.get("dearomatizes"):
        oxygen_sites = [atom for atom in atoms if atom in diff.needs_oxygen]
        if not oxygen_sites:
            return False
        # The local change is a carbonyl on a ring that stops being aromatic.
        # A bare hydroxylation does not do that; the dearomatizing edit does.
        if all(
            atom in diff.needs_carbonyl and atom in diff.loses_aromaticity
            for atom in oxygen_sites
        ):
            return False

    if effect.get("dearomatizes"):
        scope = atoms | set(info.get("path_ends") or ())
        if not (scope & set(diff.loses_aromaticity)):
            return False

    removes = effect.get("removes") or ""
    if "H" in removes and not _effect_adds_oxygen(effect) and not effect.get("cleaves"):
        scope = atoms | set(info.get("path_ends") or ())
        loses_h = any(diff.h_delta.get(atom, 0) < 0 for atom in scope)
        if not loses_h and not (scope & set(diff.loses_aromaticity)):
            return False
    return True


def _filters(diff, enabled, mol):
    if not enabled:
        return (lambda rule, info: True), (lambda site, info: True)

    def filter_rules(rule, info):
        return _pattern_could_help(info, diff)

    def filter_sites(site, info):
        return _site_could_help(site, info, diff, mol)

    return filter_rules, filter_sites


def _rule_can_cleave(rule):
    patterns = list(getattr(rule, "smarts", ()) or ())
    patterns = patterns + list(getattr(rule, "endpoints", ()) or ())
    for _smarts, info in patterns:
        if _any_span(info.get("span") or {}, "cleaves", bool, False):
            return True
    return False


# ---------------------------------------------------------------------------
# Steps and dependency
# ---------------------------------------------------------------------------


def _atom_ref(mol, idx):
    """Origin index, or ``added_by`` when this atom was created by a step."""

    atom = mol.GetAtomWithIdx(idx)
    forest = getattr(mol, "_forest", None) or {}
    records = (forest.get("atom_trace") or {}).get("records") or {}
    if atom.HasProp("forestLabel"):
        record = records.get(atom.GetProp("forestLabel")) or {}
        added = record.get("added_by")
        if added:
            rule = added["rule"]
            name = rule if isinstance(rule, str) else rule.name
            return AtomRef(added_by=(name, frozenset(added["site"])))
    return idx


def _step(mol, rule_name, site):
    return Step(rule_name, [_atom_ref(mol, idx) for idx in site])


def _steps_for(mol, info):
    """Phase-I steps for one accepted edit.

    QuinoneFormation is not itself a step. The hop stands in for the
    hydroxylations that supply each oxygen and the dehydrogenation that
    follows them.
    """

    if info["rule"].name == "QuinoneFormation":
        return _quinone_phase1(mol, info)
    return (_step(mol, info["rule"].name, info["site"]),)


def _bonded(mol, idx, atomic_num):
    """Neighbor of ``idx`` with this atomic number, if the mol already has one."""

    atom = mol.GetAtomWithIdx(idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == atomic_num:
            return neighbor.GetIdx()
    return None


def _quinone_phase1(mol, info):
    """Hydroxylations that supply missing oxygens, then one dehydrogenation.

    An end that already carries oxygen keeps that atom. An end that
    ``needs`` oxygen becomes a Hydroxylation; the dehydrogenation site
    points at the oxygen that step would add. No ``end_maps`` field is
    required: the partner atom is the neighbor already on ``mol``.
    """

    ends = tuple(info.get("ends") or ())
    end_atoms = tuple(info.get("end_atoms") or ())
    hydroxylations = []
    dh_refs = []
    for end, atom in zip(ends, end_atoms):
        partner = end.get("partner") or ""
        adds_oxygen = "O" in (end.get("needs") or "") or (
            "O" in (end.get("adds") or "") and partner != "O"
        )
        if adds_oxygen:
            step = Step("Hydroxylation", [_atom_ref(mol, atom)])
            hydroxylations.append(step)
            origin = next(iter(step.site)).origin
            dh_refs.append(AtomRef(added_by=("Hydroxylation", frozenset({origin}))))
            continue
        atomic_num = {"O": 8, "N": 7, "C": 6, "S": 16}.get(partner)
        hetero = _bonded(mol, atom, atomic_num) if atomic_num else None
        if hetero is not None:
            dh_refs.append(_atom_ref(mol, hetero))
    if not dh_refs:
        return (_step(mol, "Dehydrogenation", info["site"]),)
    return tuple(hydroxylations) + (Step("Dehydrogenation", dh_refs),)


def _deps(steps):
    """A later step depends on an earlier one when its site was ``added_by`` it."""

    edges = []
    for later, step in enumerate(steps):
        for ref in step.site:
            if ref.added_by is None:
                continue
            rule_name, site = ref.added_by
            wanted = frozenset(site)
            for earlier, previous in enumerate(steps):
                if previous.rule != rule_name:
                    continue
                origins = frozenset(
                    item.origin for item in previous.site if item.origin is not None
                )
                if origins == wanted:
                    edges.append((earlier, later))
    return Deps(steps, edges)


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------


def default_ruleset():
    """The poc catalog as one rule. Phase I is the ``Deps`` search yields."""

    return RuleSet(
        (Dealkylation, QuinoneFormation, Hydroxylation, Dehydrogenation),
        name="Poc",
    )


def _pieces(raw):
    """One mol, or the fragments of a disconnected reaction product."""

    try:
        groups = Chem.GetMolFrags(raw)
    except ValueError:
        return [raw]
    if len(groups) <= 1:
        return [raw]
    return list(GetMolFrags(raw, asMols=True, sanitizeFrags=False)) or [raw]


def _finish(parent, raw_products, info, counters):
    """Sanitize and trace each fragment. Failed sanitizes are dropped."""

    finished = []
    for raw in raw_products:
        for piece in _pieces(raw):
            if SanitizeMol(piece, catchErrors=True):
                counters.sanitize_dropped += 1
                continue
            for atom in piece.GetAtoms():
                atom.SetAtomMapNum(0)
            forest_trace(parent, piece, info["rule"], info["site"])
            ordered, smiles = cannonicalize_order(piece)
            finished.append((ordered, smiles))
    return finished


def _keep_fragment(finished, target):
    """The fragment closest to ``target``. The rest were cleaved off."""

    best = None
    best_cost = None
    for mol, smiles in finished:
        if smiles == canon_smiles(target):
            cost = -1
        else:
            cost = atom_diff(mol, target).cost()
        if best is None or cost < best_cost:
            best = (mol, smiles)
            best_cost = cost
    discarded = [item for item in finished if item[0] is not best[0]]
    return best, discarded


@dataclass
class _Walk:
    mol: object
    steps: tuple
    sides: tuple
    opens: tuple


def find_path(
    reactant,
    target,
    ruleset=None,
    counters=None,
    *,
    use_filters=True,
    max_paths=1,
    max_nodes=800,
):
    """Yield phase-I plans that turn ``reactant`` into ``target``.

    ``ruleset`` is one :class:`~xenosite.refactor_poc.rules.RuleSet`. Filters
    come from :func:`atom_diff` and are applied to each child pattern; the
    set does not hide them. Pass ``use_filters=False`` to bill the same pair
    with every site edited. ``counters`` is optional; pass one in when the
    test needs ``billed``.
    """

    reactant = _as_mol(reactant)
    target_mol = _as_mol(target)
    reactant = Chem.Mol(reactant)
    install_forest(reactant)
    target_smiles = canon_smiles(target_mol)
    if ruleset is None:
        ruleset = default_ruleset()
    if counters is None:
        counters = PathCounters()

    queue = deque([_Walk(reactant, (), (), ())])
    seen = {canon_smiles(reactant)}
    found = 0

    while queue and found < max_paths and counters.nodes < max_nodes:
        walk = queue.popleft()
        counters.nodes += 1
        here = canon_smiles(walk.mol)
        if here == target_smiles:
            yield PathOutcome(
                plan=_deps(walk.steps),
                maybe=Maybe(walk.sides),
                smiles=here,
            )
            found += 1
            continue

        diff = atom_diff(walk.mol, target_mol)
        parent_cost = diff.cost()
        filter_rules, filter_sites = _filters(diff, use_filters, walk.mol)
        # Cleavage children run first. The set still passes each pattern to
        # the filters; this only picks an order.
        order_key = None
        if diff.target_smaller or diff.has_cleavage:
            order_key = lambda rule: (0 if _rule_can_cleave(rule) else 1, rule.name)

        for por in ruleset.metabolites(
            walk.mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            counters=counters,
            order_key=order_key,
        ):
            finished = _finish(walk.mol, por.products, por.info, counters)
            if not finished:
                continue
            kept, discarded = _keep_fragment(finished, target_mol)
            if kept is None:
                continue
            child, child_smiles = kept
            closer = (
                child_smiles == target_smiles
                or atom_diff(child, target_mol).cost() < parent_cost
            )
            if not closer:
                continue
            if child_smiles in seen and child_smiles != target_smiles:
                continue
            seen.add(child_smiles)

            cleaves = bool((por.info.get("options") or {}).get("cleaves"))
            if len(finished) == 1 and cleaves:
                opens = walk.opens + (por.info["site"],)
                sides = walk.sides
            elif len(finished) > 1:
                opens = walk.opens
                sides = walk.sides + tuple(
                    CleavageSide(
                        site=por.info["site"],
                        side=smiles,
                        opens=walk.opens,
                    )
                    for _mol, smiles in discarded
                )
            else:
                opens = walk.opens
                sides = walk.sides

            steps = walk.steps + _steps_for(walk.mol, por.info)
            queue.append(_Walk(child, steps, sides, opens))
