"""Define specific reaction rules."""

from __future__ import annotations

# Standard Library
import copy
import itertools
from collections import defaultdict, deque
from collections.abc import Callable, Generator, Iterable, Mapping, Sequence
from typing import NamedTuple

from xenosite.refactor_poc.rdkit_api import (
    ChemicalReaction,
    MolFromSmarts,
    MolToSmarts,
    RenumberAtoms,
)
from xenosite.refactor_poc.rdkitutil import (
    Atom,
    Bond,
    BondType,
    ForestMol,
    ForestTracingMol,
    Mol,
    RWMol,
    aromatic_systems,
    cannonicalize_order,
    conjugated_systems,
    copy_mol,
    ensure_forest,
    ensure_kekule_parents,
    is_tracing,
    molecule_formula,
    move_charge_with_bonds,
    parent_for_bond,
    parents_for_ends,
    reaction_from_smarts,
    resonance_bond_maps,
    ring_membership,
    run_reactants,
    rw_copy,
    sanitize_catch,
    sanitize_mol,
    sanitized_fragments,
    smarts_matches,
    topol_equiv,
)
from xenosite.refactor_poc.records import (
    AtomTrace,
    Effect,
    Forest,
    Formula,
    InitializedAtomTrace,
    KekuleParents,
    PatternInfo,
    Site,
    TraceAddition,
    When,
)


def set_terminal_product(mol: Mol, value: bool = True) -> ForestMol:
    held = ensure_forest(mol)
    held._forest["is_terminal_product"] = value
    return held


class ProductsOfReaction(NamedTuple):
    """One edit from :meth:`ReactionRule.metabolites`, before tracing.

    ``info`` describes the edit. ``products`` are the mols it made.
    :meth:`ReactionRule.metabolize` turns each of those mols into a
    ``(product, info)`` pair and is what callers should use.
    """

    info: dict[str, object]
    products: list[Mol]


class ReactionRule:
    """A reaction, and the contract every reaction keeps.

    Call :meth:`metabolize`. Subclasses implement :meth:`metabolites` only.
    The invariants below are the library's. A new rule does not restate them
    and does not bypass them.
    """

    is_terminal_rule: bool = False
    name: str | None
    longname: str | None
    sites_on: object | None = None

    def _clear_atom_maps(self, mol: Mol) -> Mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    def __init__(
        self,
        name: str | None = None,
        sites_on: object | None = None,
        longname: str | None = None,
        *args: object,
        **kwargs: object,
    ) -> None:
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
            if self.name is not None:
                assert "_" not in self.name
        except AssertionError as err:
            err.args = ("Cannot have '_' in rule name",)
            raise

        if sites_on is not None:
            self.sites_on = sites_on

    def __call__(
        self, mol: Mol, **kwargs: object
    ) -> Generator[tuple[ForestTracingMol, dict[str, object]], None, None]:
        yield from self.metabolize(mol, **kwargs)

    def __iter__(self) -> Iterable[ReactionRule]:
        return iter([self])

    def metabolize(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda rule, info: True,
        filter_sites: FilterSites = lambda site, info: True,
        unique_csmi: bool = True,
        **kwargs: object,
    ) -> Generator[tuple[ForestTracingMol, dict[str, object]], None, None]:
        """Apply this rule and yield ``(product, info)`` pairs.

        This is the method callers use. It keeps the invariants below.
        :meth:`metabolites` supplies the chemistry and must not try to.

        The mol that was passed in:

        - Its atoms, bonds, charges, hydrogens, and atom-map numbers are
          unchanged.
        - If it has no ``_forest`` trace, one is installed and left there.
        - If it already has a trace, that trace's depth is not changed.

        Every product:

        - Is a sanitized mol with its own ``_forest``.
        - Has ``atom_trace["depth"]`` one greater than the parent.
        - Has ``atom_trace["formula"]`` equal to the atom counts and formal
          charge of that product, hydrogens included.
        - Records each new transform once, as ``R1``, ``R2``, and so on.
          A new atom's ``added_by`` is that id. The site, the rule
          hierarchy, the resolved effect, the name, ``phase1``, and the
          depth of the site's index frame live under
          ``atom_trace["additions"][id]``. The change in formula lives
          under ``atom_trace["delta_formula"][id]``.
        - Carries ``info["csmi"]``, the canonical SMILES of that product.
          The same SMILES is not yielded twice.

        ``filter_rules(rule, pattern_info)`` sees the pattern before a
        match. ``filter_sites(site, info)`` sees the resolved effect
        before an edit. Either may refuse. A refusal edits nothing.
        """
        assert mol is not None
        # The caller's chemistry is not edited. A mol with no trace gets one,
        # at its current depth, so products can sit one step below it.
        ensure_tracing(mol)
        # Matching, map clearing, and forest stamps happen on a copy.
        mol = ensure_tracing(_work_copy(mol))

        if self.is_terminal_product(mol):
            return

        # Maps must not be present for CanonicalRankAtoms or SMARTS matching.
        self._clear_atom_maps(mol)

        seen: set[str] = set()

        for por in self.metabolites(
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            **kwargs,
        ):
            info = por.info
            products = por.products
            # print("INFO", info, len(products))

            # Same canonical SMILES is one outcome. Two sites in one atom
            # class can still be different molecules (ortho quinone and para).
            for n, p in enumerate(products):
                sanitize_mol(p)

                if self.is_terminal_rule:
                    set_terminal_product(p)

                forest_trace(mol, p, info, executed=self)

                ordered, csmi = cannonicalize_order(p)
                assert is_tracing(ordered)
                # canonicalize deep-copies the forest; keep the rule's
                # PatternInfo object as the addition's pattern link.
                pattern = info.get("pattern")
                if pattern is not None:
                    tid = ordered._forest["atom_trace"]["transforms"][-1]
                    ordered._forest["atom_trace"]["additions"][tid]["pattern"] = (
                        pattern
                    )

                if unique_csmi:
                    if csmi in seen:
                        continue
                    seen.add(csmi)

                i = dict(info)
                i["product_index"] = n
                i["product_count"] = len(products)
                i["csmi"] = csmi

                yield ordered, i

    def _top_site(self, site: object, mol: Mol) -> int | frozenset[int]:
        te = topol_equiv(mol)
        if type(site) is int:
            return te[site]
        else:
            return frozenset(int(te[s]) for s in site)

    def is_terminal_product(self, mol: Mol) -> bool:
        """True if ``mol`` must not be expanded further in guided path search."""

        forest = ensure_forest(mol)._forest
        structure = forest.get("structure")
        if not structure or "is_terminal_product" not in structure:
            return False
        return structure["is_terminal_product"]

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda rule, info: True,
        filter_sites: FilterSites = lambda site, info: True,
        **kwargs: object,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Yield one :class:`ProductsOfReaction` per edit. Subclasses override this.

        Do not override :meth:`metabolize`. It is what keeps the caller-facing
        invariants: the input chemistry is unchanged, the parent gains a
        ``_forest`` if it lacked one, and every product is one depth below
        that parent with a single ``additions`` record per transform.

        Each yielded item:

        - ``info["site"]`` is the atom or atom pair this edit is about.
        - ``info["rule"]`` is this rule. A containing ruleset is not
          substituted for it.
        - ``info["options"]`` is the one resolved effect at that site.
        - ``products`` is a list of mols. A cleavage puts each fragment in
          that list. It does not return one mol that is several pieces.

        ``filter_rules(rule, pattern_info)`` is called before a match and
        can see the pattern's ``span``. False means that pattern is skipped.
        ``filter_sites(site, info)`` is called after the effect is resolved
        and before any edit. False means that site is skipped.

        This method must not edit the mol it is given. ``metabolize`` has
        already handed it a copy. It must not attach ``_forest`` to the
        caller's mol, and it must not set product depth. ``metabolize``
        does both after this method returns.
        """
        raise NotImplementedError


FilterRules = Callable[[ReactionRule, PatternInfo], bool]
FilterSites = Callable[[object, dict[str, object]], bool]


def __getattr__(name: str) -> object:
    """``RuleSet`` lives in ``rulesets``. Keep ``from .rules import RuleSet`` working."""

    if name == "RuleSet":
        from .rulesets import RuleSet

        return RuleSet
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def formula_delta(before: Formula, after: Formula) -> Formula:
    """Change in atom counts and formal charge from ``before`` to ``after``."""

    keys = set(before.get("counts") or {}) | set(after.get("counts") or {})
    counts: dict[str, int] = {}
    for key in keys:
        delta = (after.get("counts") or {}).get(key, 0) - (
            before.get("counts") or {}
        ).get(key, 0)
        if delta:
            counts[key] = delta
    return {
        "counts": counts,
        "charge": after.get("charge", 0) - before.get("charge", 0),
    }


def _rule_name(rule: object) -> str | None:
    if rule is None:
        return None
    if isinstance(rule, str):
        return rule
    return getattr(rule, "name", type(rule).__name__)


def _work_copy(mol: Mol) -> Mol:
    """A mol the rule may stamp. The caller's object is left alone."""

    return copy_mol(mol)


def ensure_tracing(mol: Mol) -> ForestTracingMol:
    """Forest is present and the trace is initialized. Does not reset depth.

    A missing trace is started the way :func:`install_forest` always has.
    An existing trace, including its depth, is left in place.
    """

    held = ensure_forest(mol)
    forest = held._forest
    if "structure" not in forest:
        forest["structure"] = {}
    if "atom_trace" not in forest:
        trace: InitializedAtomTrace = {
            "records": {},
            "deletes": {},
            "transforms": [],
            "additions": {},
            "formula": molecule_formula(held),
            "delta_formula": {},
            "depth": 0,
            "last_tag": 0,
            "next_transform": 1,
        }
        forest["atom_trace"] = trace
        for a in held.GetAtoms():
            if a.GetAtomicNum() != 1:
                i = a.GetIdx()
                trace["records"][str(i)] = {
                    "idx": [i],
                    "depth": [0],
                }
                trace["last_tag"] = i
    assert is_tracing(held)
    _write_forest_labels(held)
    return held


def install_forest(mol: Mol) -> Forest:
    """Install the forest data structure in the molecule."""

    return ensure_tracing(mol)._forest


def _write_forest_labels(mol: ForestTracingMol) -> None:
    for tag, record in mol._forest["atom_trace"]["records"].items():
        idx = record.get("idx")
        if not idx:
            raise KeyError("idx")
        mol.GetAtomWithIdx(idx[-1]).SetProp("forestLabel", tag)


def stamp_forest_labels(mol: Mol) -> ForestTracingMol:
    return ensure_tracing(mol)


def reordered_forest_labels(mol: ForestTracingMol) -> None:
    trace = mol._forest["atom_trace"]
    # if "atom_trace" not in forest:
    #     install_forest(mol)

    for atom in mol.GetAtoms():
        i = atom.GetIdx()

        if atom.GetAtomicNum() != 1:
            assert atom.HasProp("forestLabel")
            tag = atom.GetProp("forestLabel")
            record = trace["records"][tag]
            idx = record.get("idx")
            if not idx:
                raise KeyError("idx")
            idx[-1] = i


def _site_tuple(site: object) -> Site:
    if isinstance(site, int):
        return (site,)
    if isinstance(site, (frozenset, set, list, tuple)):
        return tuple(sorted(site))
    raise TypeError(site)


def _trace_info(info: Mapping[str, object]) -> dict[str, object]:
    """Pattern fields worth keeping, without a second copy of the rule object."""

    kept: dict[str, object] = {}
    for key, value in info.items():
        if key == "rule":
            kept[key] = _rule_name(value)
        elif key == "rule_chain":
            if isinstance(value, (tuple, list)):
                kept[key] = tuple(_rule_name(item) for item in value)
            else:
                kept[key] = value
        elif key in ("options", "pattern"):
            continue
        else:
            kept[key] = value
    return kept


def _as_effect(value: object) -> Effect:
    """Copy known effect fields from an open dict."""

    effect: Effect = {
        "adds": "",
        "removes": "",
        "cleaves": False,
        "leave_count": None,
        "breaks_ring": False,
        "dearomatizes": False,
        "methide": False,
        "needs": "",
    }
    if not isinstance(value, dict):
        return effect
    adds = value.get("adds")
    if isinstance(adds, str):
        effect["adds"] = adds
    removes = value.get("removes")
    if isinstance(removes, str):
        effect["removes"] = removes
    cleaves = value.get("cleaves")
    if isinstance(cleaves, bool):
        effect["cleaves"] = cleaves
    leave_count = value.get("leave_count")
    if leave_count is None or isinstance(leave_count, int):
        effect["leave_count"] = leave_count
    breaks_ring = value.get("breaks_ring")
    if isinstance(breaks_ring, bool):
        effect["breaks_ring"] = breaks_ring
    dearomatizes = value.get("dearomatizes")
    if isinstance(dearomatizes, bool):
        effect["dearomatizes"] = dearomatizes
    methide = value.get("methide")
    if isinstance(methide, bool):
        effect["methide"] = methide
    needs = value.get("needs")
    if isinstance(needs, str):
        effect["needs"] = needs
    partner = value.get("partner")
    if isinstance(partner, str):
        effect["partner"] = partner
    partner_h = value.get("partner_h")
    if isinstance(partner_h, int):
        effect["partner_h"] = partner_h
    symbol = value.get("symbol")
    if isinstance(symbol, str):
        effect["symbol"] = symbol
    h = value.get("h")
    if isinstance(h, int):
        effect["h"] = h
    site_aromatic = value.get("site_aromatic")
    if isinstance(site_aromatic, bool):
        effect["site_aromatic"] = site_aromatic
    when = value.get("when")
    if isinstance(when, dict):
        effect["when"] = when
    return effect


def forest_trace(
    reactant: Mol,
    product: Mol,
    info: Mapping[str, object],
    executed: ReactionRule | None = None,
) -> AtomTrace:
    """Record one transform on the product's atom trace.

    A new atom's ``added_by`` is an id such as ``R1``. The site, the rule
    hierarchy, the resolved effect, and the formula change live once, under
    ``atom_trace["additions"][id]``. ``depth`` is the reactant index frame
    the site is written in.
    """

    rule = info.get("rule")
    site = info.get("site")
    parent = stamp_forest_labels(reactant)
    held = ensure_forest(product)
    trace = copy.deepcopy(parent._forest["atom_trace"])
    trace.setdefault("additions", {})
    trace.setdefault("delta_formula", {})
    trace.setdefault("transforms", [])
    trace.setdefault("next_transform", 1)

    number = trace["next_transform"]
    trace["next_transform"] = number + 1
    transform_id = "R%d" % number

    depth = trace["depth"] = trace["depth"] + 1
    frame = depth - 1
    chain: list[ReactionRule] = []
    rule_chain = info.get("rule_chain") or ()
    extras: tuple[object, ...]
    if isinstance(rule_chain, tuple):
        extras = rule_chain
    elif isinstance(rule_chain, list):
        extras = tuple(rule_chain)
    else:
        extras = ()
    for item in extras + ((executed,) if executed is not None else ()):
        if item is rule or not isinstance(item, ReactionRule):
            continue
        if any(item is seen for seen in chain):
            continue
        chain.append(item)
    if isinstance(rule, ReactionRule) and all(rule is not seen for seen in chain):
        chain.append(rule)

    before = trace.get("formula") or molecule_formula(reactant)
    after = molecule_formula(held)
    trace["formula"] = after
    trace["delta_formula"][transform_id] = formula_delta(before, after)
    addition: TraceAddition = {
        "site": _site_tuple(site),
        "rules": tuple(chain),
        "info": _trace_info(info),
        "effect": _as_effect(info.get("options")),
        "name": _rule_name(rule),
        "phase1": info.get("phase1"),
        "depth": frame,
        "pattern": info.get("pattern") if isinstance(info.get("pattern"), dict) else None,
    }
    trace["additions"][transform_id] = addition
    trace["transforms"].append(transform_id)

    records = trace["records"]
    new_records = {}

    for product_atom in held.GetAtoms():
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
                "added_by": transform_id,
            }
            product_atom.SetProp("forestLabel", str(tag))

    for tag, record in records.items():
        record["removed_by"] = transform_id
        trace["deletes"][tag] = record

    trace["records"] = new_records
    held._forest["atom_trace"] = trace

    return trace


_EFFECT_DEFAULTS = {
    "adds": "",
    "removes": "",
    "cleaves": False,
    "leave_count": None,
    "breaks_ring": False,
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
    15: "P",
    16: "S",
    17: "Cl",
    35: "Br",
    53: "I",
    85: "At",
}


def _span(possibilities: Sequence[Mapping[str, object]]) -> dict[str, object]:
    """Certain values stay bare. Disagreeing values become a tuple."""

    keys: list[str] = []
    for possibility in possibilities:
        for key in possibility:
            if key != "when" and key not in keys:
                keys.append(key)
    span: dict[str, object] = {}
    for key in keys:
        values: list[object] = []
        for possibility in possibilities:
            value = possibility.get(key, _EFFECT_DEFAULTS.get(key))
            if value not in values:
                values.append(value)
        span[key] = values[0] if len(values) == 1 else tuple(values)
    return span


def branches(
    whens: Sequence[When],
    site_map: int = 1,
    removes_partner: bool = False,
    **effect: object,
) -> tuple[Effect, ...]:
    """Copy ``effect`` once per ``when``. The site atom and the OR atom differ.

    A ``when`` on ``site_map`` records ``h`` / ``symbol`` (the site was the
    ambiguous atom). A ``when`` on any other map records ``partner`` /
    ``partner_h``.
    """

    out: list[Effect] = []
    for when in whens:
        when = dict(when)
        item = _as_effect(effect)
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


def describe(
    *possibilities: Effect,
    edit: str | None = None,
    site_map: int | tuple[int, ...] = 1,
    pin: tuple[int, ...] | None = None,
    skip_same_rings: bool = False,
    name: str | None = None,
    **single: object,
) -> PatternInfo:
    """Build a :class:`PatternInfo`.

    One outcome: ``describe(adds="O", removes="H")``.
    Several: ``describe(*branches(...), edit="dealkylate")``.
    ``name`` distinguishes this pattern from the others on the same rule.
    """

    if single and possibilities:
        raise TypeError("pass one effect as keywords, or several dicts, not both")
    normalized: list[Effect]
    if single:
        normalized = [_as_effect(single)]
    else:
        if not possibilities:
            raise TypeError("a pattern needs at least one possibility")
        normalized = []
        for possibility in possibilities:
            effect = _as_effect(_EFFECT_DEFAULTS)
            effect.update(possibility)
            normalized.append(effect)
    info: PatternInfo = {
        "possibilities": tuple(normalized),
        "span": _span(normalized),
        "site_map": site_map,
    }
    if edit is not None:
        info["edit"] = edit
    if pin:
        info["pin"] = tuple(pin)
    if skip_same_rings:
        info["skip_same_rings"] = True
    if name is not None:
        info["name"] = name
    return info


def _assign_pattern_names(patterns: Iterable[PatternInfo]) -> None:
    """Give each unnamed pattern a numbered name on that pattern object."""

    next_num = 1
    for pattern in patterns:
        if pattern.get("name"):
            continue
        pattern["name"] = str(next_num)
        next_num += 1


def may(info: PatternInfo, key: str, value: object = True) -> bool:
    """True if any possibility has this outcome.

    For ``adds`` / ``removes`` / ``needs``, ``value`` may be a substring.
    """

    for possibility in info["possibilities"]:
        have = possibility.get(key, _EFFECT_DEFAULTS.get(key))
        if (
            isinstance(value, str)
            and isinstance(have, str)
            and key
            in (
                "adds",
                "removes",
                "needs",
            )
        ):
            if value in have:
                return True
        elif have == value:
            return True
    return False


def must(info: PatternInfo, key: str, value: object = True) -> bool:
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


def _when_matches(mol: Mol, mapped: Mapping[int, int], when: When) -> bool:
    idx = mapped.get(when.get("map", -1))
    if idx is None:
        return False
    atom = mol.GetAtomWithIdx(idx)
    if "z" in when and atom.GetAtomicNum() != when["z"]:
        return False
    if "h" in when and atom.GetTotalNumHs() != when["h"]:
        return False
    return True


def resolve_effect(mol: Mol, mapped: Mapping[int, int], info: PatternInfo) -> Effect:
    """Narrow ``info`` to the one effect this match actually is.

    Uses the caller's mol, not a kekulé copy: aromatic flags must still be set.
    """

    possibilities = info.get("possibilities") or (_as_effect(info),)
    chosen: Effect | None = None
    for possibility in possibilities:
        when = possibility.get("when")
        if when is None or _when_matches(mol, mapped, when):
            chosen = possibility
            break
    if chosen is None:
        chosen = possibilities[0]
    effect = _as_effect(_EFFECT_DEFAULTS)
    effect.update(chosen)
    effect.pop("when", None)
    site_map = info.get("site_map", 1)
    if isinstance(site_map, int) and site_map in mapped:
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
    if effect.get("cleaves"):
        effect["breaks_ring"] = _cleavage_breaks_ring(mol, mapped, site_map)
    return effect


def _cleavage_breaks_ring(
    mol: Mol, mapped: Mapping[int, int], site_map: int | tuple[int, ...]
) -> bool:
    """True when the two site atoms share a ring, so the cleaved bond is in it."""

    if not isinstance(site_map, tuple) or len(site_map) != 2:
        return False
    left = mapped.get(site_map[0])
    right = mapped.get(site_map[1])
    if not isinstance(left, int) or not isinstance(right, int):
        return False
    rings = ring_membership(mol)
    return bool(set(rings.get(left, ())) & set(rings.get(right, ())))


def merge_effects(
    left: Effect, right: Effect, both_aromatic: bool
) -> dict[str, object]:
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


def _bump(counters: object | None, name: str, amount: int = 1) -> None:
    if counters is None:
        return
    current = getattr(counters, name)
    if not isinstance(current, int):
        raise TypeError(name)
    setattr(counters, name, current + amount)


def _isotope_smarts(smarts: str, pin: tuple[int, ...]) -> str:
    """Restrict a SMARTS reaction to the pinned maps.

    The matched atoms are stamped with isotope ``8000 + map number``. The
    isotope is applied to the whole atom query, so an OR list does not keep
    the label on only its first alternative. Those atoms are written first:
    ``SubstructMatch`` starts at query atom 0 and does not reorder.
    """

    import re

    wanted = set(pin)

    def repl(match):
        mapno = int(match.group(2))
        if mapno not in wanted:
            return match.group(0)
        isotope = 8000 + mapno
        return "[" + match.group(1) + "&" + str(isotope) + "*:" + match.group(2) + "]"

    rewritten = re.sub(r"\[([^\[\]]*):(\d+)\]", repl, smarts)
    return _isotope_atoms_first(rewritten, wanted)


def _isotope_atoms_first(smarts: str, pin: set[int]) -> str:
    """Put each isotope-bearing reactant atom before the others."""

    if ">>" not in smarts or not pin:
        return smarts
    reactant, product = smarts.split(">>", 1)
    query = MolFromSmarts(reactant)
    if query is None:
        return smarts
    pinned = [
        idx
        for _mapno, idx in sorted(
            (atom.GetAtomMapNum(), atom.GetIdx())
            for atom in query.GetAtoms()
            if atom.GetAtomMapNum() in pin
        )
    ]
    if not pinned:
        return smarts
    rest = [idx for idx in range(query.GetNumAtoms()) if idx not in set(pinned)]
    order = pinned + rest
    if order == list(range(query.GetNumAtoms())):
        return smarts
    rewritten = MolToSmarts(RenumberAtoms(query, order))
    if MolFromSmarts(rewritten) is None:
        return smarts
    return rewritten + ">>" + product


def _site_indexes(mapped: Mapping[int, int], pattern: PatternInfo) -> frozenset[int]:
    """Atom indexes the pattern calls the site. Defaults to map 1."""

    key = pattern.get("site_map", 1)
    if isinstance(key, (list, tuple)):
        idxs = [mapped[k] for k in key if k in mapped]
    elif key in mapped:
        idxs = [mapped[key]]
    else:
        idxs = list(mapped.values())
    return frozenset(idxs)


def react_at(
    rule: SmartsReactionRule,
    smarts: str,
    mol: Mol,
    mapped: Mapping[int, int],
    counters: object | None = None,
    pin: tuple[int, ...] | None = None,
) -> list[Mol]:
    """Run ``smarts`` on one match. One call is one ``mol_edits``.

    ``pin`` is the pattern's map numbers. Each of those atoms is stamped.
    Absent ``pin`` stamps every mapped atom, which is the same list.
    """

    _bump(counters, "mol_edits")
    chosen = tuple(
        mapno
        for mapno in (pin if pin is not None else tuple(mapped))
        if mapno in mapped
    )
    if not chosen:
        chosen = tuple(mapped)
    stamped = rw_copy(mol)
    for mapno in chosen:
        stamped.GetAtomWithIdx(mapped[mapno]).SetIsotope(8000 + mapno)
    try:
        product_sets = run_reactants(_isotope_smarts(smarts, chosen), stamped)
    except (RuntimeError, ValueError):
        return []
    if not product_sets:
        return []

    products: list[Mol] = []
    for prod in product_sets[0]:
        for atom in prod.GetAtoms():
            if atom.GetIsotope() >= 8000:
                atom.SetIsotope(0)
            atom.SetAtomMapNum(0)
        products.append(rule._lift_forest_labels(mol, prod))
    return products


class SmartsReactionRule(ReactionRule):
    """Performs reactions specified by SMARTS.

    Each entry is ``(smarts, options)`` with :class:`PatternInfo` options.
    Matches are filtered before ``RunReactants``. One SMARTS is applied once
    per topological site; a later SMARTS with the same formula effect still runs.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = ()

    rxns: list[tuple[str, ChemicalReaction, PatternInfo]]

    def __init__(self, *args: object, **kwargs: object) -> None:

        super().__init__(*args, **kwargs)

        self.rxns = [
            (smarts, self._smarts2rxns(smarts, **kwargs), opt)
            for smarts, opt in self.smarts
        ]
        patterns = [opt for _smarts, _rxn, opt in self.rxns]
        endpoints = getattr(self, "endpoints", ()) or ()
        patterns.extend(opt for _smarts, opt in endpoints)
        _assign_pattern_names(patterns)

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda rule, info: True,
        filter_sites: FilterSites = lambda site, info: True,
        context_mol: Mol | None = None,
        **kwargs: object,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`ReactionRule.metabolites`.

        ``filter_rules`` sees this rule and the pattern, including ``span``,
        before SMARTS runs. ``filter_sites`` sees the resolved effect in
        ``info["options"]`` before ``RunReactants``. ``context_mol`` is the
        unsubstituted parent when ``mol`` is a kekulé copy. Aromatic flags
        are read from it. ``counters``, when passed, records one
        ``rule_expansions`` per call, one ``sites_considered`` per match,
        ``sites_skipped`` when a filter or a topological duplicate refuses
        the site, and one ``mol_edits`` inside :func:`react_at`.
        """

        counters = kwargs.get("counters")
        context = mol if context_mol is None else context_mol
        _bump(counters, "rule_expansions")
        seen: set[object] = set()
        ranks = topol_equiv(context)

        for work in _kekule_forms(mol):
            for rxn_num, (smarts, _rxn, pattern) in enumerate(self.rxns):
                if not filter_rules(self, pattern):
                    continue

                reactant = smarts.split(">>", 1)[0]
                for mapped in smarts_matches(work, reactant):
                    site = _site_indexes(mapped, pattern)
                    if not site:
                        continue
                    effect = resolve_effect(context, mapped, pattern)
                    info: dict[str, object] = {
                        "site": site,
                        "rule": self,
                        "options": effect,
                        "rxn_num": rxn_num,
                        "pattern": pattern,
                    }
                    _bump(counters, "sites_considered")
                    if not filter_sites(site, info):
                        _bump(counters, "sites_skipped")
                        continue
                    # Same map roles and the same incident bond orders are one
                    # edit. Equivalent carbons share a rank. Another Kekulé
                    # writing, or swapping which atom is map 1, is not.
                    signature = (
                        tuple(
                            (mapno, ranks[idx])
                            for mapno, idx in sorted(mapped.items())
                        ),
                        _incident_orders(work, ranks, mapped),
                        rxn_num,
                        effect.get("adds"),
                        effect.get("removes"),
                        bool(effect.get("cleaves")),
                        bool(effect.get("dearomatizes")),
                    )
                    if signature in seen:
                        _bump(counters, "sites_skipped")
                        continue
                    products = react_at(
                        self, smarts, work, mapped, counters, pattern.get("pin")
                    )
                    if not products:
                        continue
                    seen.add(signature)
                    yield ProductsOfReaction(info=info, products=products)

    def _smarts2rxns(
        self,
        smarts: str,
        use_implicit_properties: bool = False,
        **kwargs: object,
    ) -> ChemicalReaction:
        """Converts SMARTS reactions to RDKit reactions."""
        return reaction_from_smarts(smarts)

    # TODO: should this method be moved to library? Needed by this class,
    # but isn't specific to its quirks.
    def _lift_forest_labels(self, reactant: Mol, product: Mol) -> Mol:
        """Carry the forest labels from the reactant to the product."""

        for a in product.GetAtoms():
            if a.HasProp("react_atom_idx"):
                react_idx = int(a.GetProp("react_atom_idx"))
                reactant_atom = reactant.GetAtomWithIdx(react_idx)
                if reactant_atom.HasProp("forestLabel"):
                    a.SetProp("forestLabel", reactant_atom.GetProp("forestLabel"))

        return product

    # this method should probably stay here, because it has to do with how
    # reactions are process.
    # TODO: A reaction helper in rdkitutil could apply reactions and return
    # products while maintaining the forest instead of maintaing that logic here.
    def _get_product_mappings(
        self, product: Mol
    ) -> tuple[dict[int, int], dict[int, int]]:
        mapno2idx: dict[int, int] = {}
        reactant2idx: dict[int, int] = {}
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

    def _get_site(self, products: Sequence[Mol]) -> frozenset[int]:
        """Return the sites in product based on the reactant_idx property assigned by
        rxns.RunReactants in self.metabolites."""

        site: set[int] = set()
        for product in products:
            mapno2idx, _ = self._get_product_mappings(product)
            site = site | set(mapno2idx.values())
        return frozenset(site)


# Bond order stored in the structure cache. Not RDKit mols: those do not
# belong in a dict that is deep-copied onto every product.
# Int keys: 1.0, 2.0, and 3.0 hash the same as 1, 2, and 3, so a
# GetBondTypeAsDouble() lookup still hits these entries.
_BOND = {
    1: BondType.SINGLE,
    2: BondType.DOUBLE,
    3: BondType.TRIPLE,
    1.5: BondType.AROMATIC,
}


def _bond_key(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def _current_bond_map(mol: Mol) -> dict[tuple[int, int], float]:
    bonds: dict[tuple[int, int], float] = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[_bond_key(i, j)] = bond.GetBondTypeAsDouble()
    return bonds


def _connected_components(mol: Mol, atoms: Iterable[int]) -> list[frozenset[int]]:
    atoms = set(atoms)
    seen: set[int] = set()
    systems: list[frozenset[int]] = []
    for start in atoms:
        if start in seen:
            continue
        comp: set[int] = set()
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


def system_neighbors(mol: Mol, system: Iterable[int]) -> dict[int, list[int]]:
    system = set(system)
    neighbors: dict[int, list[int]] = {i: [] for i in system}
    for i in system:
        for nbr in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nbr.GetIdx()
            if j in system and j > i:
                neighbors[i].append(j)
                neighbors[j].append(i)
    return neighbors


def odd_anchor_pairs(
    anchors: Sequence[int], neighbors: Mapping[int, Sequence[int]]
) -> list[tuple[int, int]]:
    """Pairs of anchors separated by an odd number of bonds.

    Ortho (1) and para (3) pass. Meta (2) does not. The walk may cross
    non-anchor atoms; ``neighbors`` is the whole system.
    """

    anchors = [a for a in anchors if a in neighbors]
    pairs: list[tuple[int, int]] = []
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


def alternating_path(
    bond_map: Mapping[tuple[int, int], float],
    start: int,
    end: int,
    neighbors: Mapping[int, Sequence[int]],
) -> list[int] | None:
    """Shortest path whose first bond is double and whose bonds then alternate."""

    if start == end or start not in neighbors or end not in neighbors:
        return None
    queue: deque[tuple[int, float, tuple[int, ...]]] = deque([(start, 2.0, (start,))])
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


def _kekule_forms(mol: Mol) -> tuple[Mol, ...]:
    """Kekulé copies when any atom is aromatic. Indexes stay put.

    One form is not enough: a ring-bond cleavage on the other Kekulé writing
    is a different product, and it will not sanitize while atoms stay
    aromatic. A mol that is already kekulé is returned as given. Forest
    labels are copied so the product can be traced.
    """

    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return (mol,)
    maps = resonance_bond_maps(mol)
    if not maps:
        return (mol,)
    forms = []
    for bond_map in maps:
        work = overlay_kekule(mol, bond_map).GetMol()
        for atom in mol.GetAtoms():
            if not atom.HasProp("forestLabel"):
                continue
            work.GetAtomWithIdx(atom.GetIdx()).SetProp(
                "forestLabel", atom.GetProp("forestLabel")
            )
        forms.append(work)
    return tuple(forms)


def _incident_orders(
    mol: Mol, ranks: dict[int, int], mapped: Mapping[int, int]
) -> tuple[tuple[int, int, float], ...]:
    """Bond orders touching the matched atoms, in rank space.

    Two Kekulé forms of the same site differ here. Two equivalent carbons
    do not, so they stay one edit.
    """

    idxs = set(mapped.values())
    bonds: list[tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i not in idxs and j not in idxs:
            continue
        a, b = sorted((ranks[i], ranks[j]))
        bonds.append((a, b, bond.GetBondTypeAsDouble()))
    return tuple(sorted(bonds))


def overlay_kekule(mol: Mol, bond_map: Mapping[tuple[int, int], float]) -> RWMol:
    """Copy ``mol`` and set kekulé bond orders. Atom indexes stay put.

    Bond orders move without their charge unless that charge is put back.
    A neutral carbon moves hydrogen instead.
    """

    before = {
        atom.GetIdx(): sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
        for atom in mol.GetAtoms()
    }
    rw = rw_copy(mol)
    for bond in rw.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = bond_map.get(_bond_key(i, j))
        if order is None:
            continue
        bond.SetBondType(_BOND.get(order, BondType.SINGLE))
        bond.SetIsAromatic(False)
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)
    move_charge_with_bonds(rw, before)
    return rw


def adjust_hydrogens(mol: Mol, idx: int, change: int) -> None:
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


def swap_bonds_along_path(mol: Mol, atoms: Sequence[int]) -> bool:
    """Flip single and double bonds along ``atoms``. Adjust H at the ends.

    False when a bond is missing or none of them flipped. A triple bond
    stays a triple bond. That copy is not a product.
    """

    bond = None
    flipped = False
    for i in range(len(atoms) - 1):
        bond = mol.GetBondBetweenAtoms(atoms[i], atoms[i + 1])
        if bond is None:
            return False
        if bond.GetBondType() == BondType.DOUBLE:
            bond.SetBondType(BondType.SINGLE)
            flipped = True
        elif bond.GetBondType() == BondType.SINGLE:
            bond.SetBondType(BondType.DOUBLE)
            flipped = True
        if i == 0:
            _correct_end_hydrogens(mol, atoms[0], bond)
    if bond is not None:
        _correct_end_hydrogens(mol, atoms[-1], bond)
    return flipped


def _correct_end_hydrogens(mol: Mol, idx: int, bond: Bond) -> None:
    if bond.GetBondType() == BondType.SINGLE:
        adjust_hydrogens(mol, idx, 1)
    elif bond.GetBondType() == BondType.DOUBLE:
        adjust_hydrogens(mol, idx, -1)


def _same_rings(
    rings: Mapping[int, tuple[tuple[int, ...], ...]], i: int, j: int
) -> bool:
    return rings.get(i, ()) == rings.get(j, ())


def edit_single_to_double(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
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


def edit_add_carbonyl_o(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
    carbon = mapped.get(1)
    if carbon is None:
        return False
    oxygen = rw.AddAtom(Atom(8))
    rw.AddBond(carbon, oxygen, BondType.DOUBLE)
    return True


def edit_replace_halogen(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
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


def edit_iminium(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
    if not edit_single_to_double(rw, mapped, {}, rings):
        return False
    rw.GetAtomWithIdx(mapped[2]).SetFormalCharge(1)
    return True


def edit_keep(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
    """Leave the endpoint bond alone. The path flip is the reaction."""

    return 1 in mapped


def edit_dealkylate(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: Mapping[str, object],
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
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


EDITS: dict[
    str,
    Callable[
        [
            RWMol,
            Mapping[int, int],
            Mapping[str, object],
            Mapping[int, tuple[tuple[int, ...], ...]],
        ],
        bool,
    ],
] = {
    "single_to_double": edit_single_to_double,
    "add_carbonyl_o": edit_add_carbonyl_o,
    "replace_halogen": edit_replace_halogen,
    "iminium": edit_iminium,
    "dealkylate": edit_dealkylate,
    "keep": edit_keep,
}


def _merge_options(
    left: Effect, right: Effect, dearomatizes: bool
) -> dict[str, object]:
    return merge_effects(left, right, dearomatizes)


def _site_atoms(mapped: Mapping[int, int], info: PatternInfo) -> int | None:
    key = info.get("site_map", 1)
    if not isinstance(key, int) or key not in mapped:
        return None
    return mapped[key]


def _kekule_cache(mol: Mol) -> KekuleParents:
    """The dict the resonance rules store. Helpers never touch ``_forest``."""

    forest = ensure_forest(mol)._forest
    if "structure" not in forest:
        raise KeyError("structure")
    structure = forest["structure"]
    cache = structure.get("kekule_parents")
    if cache is None:
        fresh: KekuleParents = {
            "parents": [],
            "orders": [],
            "systems": {},
            "by_order": {},
        }
        structure["kekule_parents"] = fresh
        return fresh
    return cache


def _reactant_parent(mol: Mol, mapped: dict[int, int], cache: KekuleParents) -> Mol | None:
    """Kekulé parent whose bond orders match an aromatic hit. Else ``mol``."""

    left = mapped.get(1)
    right = mapped.get(2)
    if left is None or right is None:
        return mol
    bond = mol.GetBondBetweenAtoms(left, right)
    if bond is None:
        return mol
    begin = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
    end = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
    if bond.GetIsAromatic():
        order = 2.0
    elif begin.GetIsAromatic() or end.GetIsAromatic():
        order = bond.GetBondTypeAsDouble()
    else:
        return mol
    ensure_kekule_parents(mol, left, right, cache)
    return parent_for_bond(cache, left, right, order)


class ResonanceRule(SmartsReactionRule):
    """Match once on the aromatic parent, then react on a cached kekulé parent.

    The pattern's ``=,:`` bond matches aromatic bonds. One parent is cached
    per assignment of the conjugated system that contains the match. Other
    systems stay aromatic. The dict is stored on the molecule's forest.
    The helpers that fill it do not read ``_forest``.
    """

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda rule, info: True,
        filter_sites: FilterSites = lambda site, info: True,
        context_mol: Mol | None = None,
        **kwargs: object,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`SmartsReactionRule.metabolites`.

        SMARTS runs once, on this mol. A hit picks the cached parent by bond
        order. The reaction runs on that copy. ``context_mol`` keeps the
        original aromatic flags for :func:`resolve_effect`.
        """

        if not self.rxns:
            return
        counters = kwargs.get("counters")
        context = mol if context_mol is None else context_mol
        _bump(counters, "rule_expansions")
        cache = _kekule_cache(mol)
        seen: set[object] = set()
        ranks = topol_equiv(context)

        for rxn_num, (smarts, _rxn, pattern) in enumerate(self.rxns):
            if not filter_rules(self, pattern):
                continue
            reactant = smarts.split(">>", 1)[0]
            for mapped in smarts_matches(mol, reactant):
                site = _site_indexes(mapped, pattern)
                if not site:
                    continue
                effect = resolve_effect(context, mapped, pattern)
                info: dict[str, object] = {
                    "site": site,
                    "rule": self,
                    "options": effect,
                    "rxn_num": rxn_num,
                    "pattern": pattern,
                }
                _bump(counters, "sites_considered")
                if not filter_sites(site, info):
                    _bump(counters, "sites_skipped")
                    continue
                work = _reactant_parent(mol, mapped, cache)
                if work is None:
                    _bump(counters, "sites_skipped")
                    continue
                signature = (
                    tuple(
                        (mapno, ranks[idx])
                        for mapno, idx in sorted(mapped.items())
                    ),
                    _incident_orders(work, ranks, mapped),
                    rxn_num,
                    effect.get("adds"),
                    effect.get("removes"),
                    bool(effect.get("cleaves")),
                    bool(effect.get("dearomatizes")),
                )
                if signature in seen:
                    _bump(counters, "sites_skipped")
                    continue
                products = react_at(
                    self, smarts, work, mapped, counters, pattern.get("pin")
                )
                if not products:
                    continue
                seen.add(signature)
                yield ProductsOfReaction(info=info, products=products)


class ResonancePairRule(ResonanceRule):
    """Two endpoint matches joined by an alternating path.

    Subclasses set ``endpoints`` to ``(smarts, PatternInfo)`` pairs and
    ``systems`` to ``"conjugated"`` or ``"aromatic"``. ``info["edit"]`` names
    a module-level edit applied at that end; the path's bonds are then flipped.

    Dehydrogenation and quinone formation edit the ends, then flip the path.
    Hydrogenation names ``keep``: the flip is the whole reaction.
    Hydroxylation does not use this. It is a one-atom SMARTS reaction.
    """

    endpoints: tuple[tuple[str, PatternInfo], ...] = ()
    systems: str = "conjugated"

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda rule, info: True,
        filter_sites: FilterSites = lambda site, info: True,
        context_mol: Mol | None = None,
        **kwargs: object,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`ReactionRule.metabolites`.

        ``filter_rules`` runs before any pair is built. ``filter_sites``
        sees the merged effect and can refuse the pair before a kekulé
        form is overlaid.
        """
        yield from ResonanceRule.metabolites(
            self,
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            context_mol=context_mol,
            **kwargs,
        )
        yield from self.pair_metabolites(
            mol, filter_rules, filter_sites, counters=kwargs.get("counters")
        )

    def pair_metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules,
        filter_sites: FilterSites,
        counters: object | None = None,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Endpoint loop.

        ``filter_rules`` drops endpoint patterns before they are matched.
        ``filter_sites`` drops a pair before a cached parent is copied.
        """

        _bump(counters, "rule_expansions")

        active = [
            (smarts, info)
            for smarts, info in self.endpoints
            if filter_rules(self, info)
        ]
        if not active:
            return

        if self.systems == "aromatic":
            systems = aromatic_systems(mol)
        else:
            systems = conjugated_systems(mol)
        if not systems:
            return

        hits: dict[int, list[tuple[dict[int, int], PatternInfo]]] = defaultdict(list)
        for smarts, info in active:
            for mapped in smarts_matches(mol, smarts):
                hits[mapped[1]].append((mapped, info))
        if len(hits) < 2:
            return

        rings: dict[int, tuple[tuple[int, ...], ...]] | None = None
        cache: KekuleParents | None = None
        for system in systems:
            anchors = [atom for atom in hits if atom in system]
            neighbors = system_neighbors(mol, system)
            for start, end in odd_anchor_pairs(anchors, neighbors):
                combos: list[
                    tuple[
                        dict[int, int],
                        PatternInfo,
                        dict[int, int],
                        PatternInfo,
                        dict[str, object],
                    ]
                ] = []
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
                    preview: dict[str, object] = {
                        "site": site,
                        "rule": self,
                        "options": merge_effects(end1, end2, both_aromatic),
                        "ends": (end1, end2),
                        "end_atoms": (site_a, site_b),
                        "end_maps": (map1, map2),
                        "path_ends": frozenset((start, end)),
                    }
                    _bump(counters, "sites_considered")
                    if not filter_sites(site, preview):
                        _bump(counters, "sites_skipped")
                        continue
                    combos.append((map1, info1, map2, info2, preview))
                if not combos:
                    continue
                if cache is None:
                    cache = _kekule_cache(mol)
                ends = parents_for_ends(mol, start, end, cache)
                paths: list[tuple[Mol, list[int]]] = []
                for parent in ends.parents:
                    path = alternating_path(
                        _current_bond_map(parent), start, end, neighbors
                    )
                    if path:
                        paths.append((parent, path))
                paths.sort(key=lambda item: len(item[1]))
                if not paths:
                    continue
                for map1, info1, map2, info2, preview in combos:
                    if rings is None and (
                        info1.get("skip_same_rings") or info2.get("skip_same_rings")
                    ):
                        rings = ring_membership(mol)
                    ring_table = rings or {}
                    for parent, path in paths:
                        _bump(counters, "mol_edits")
                        rw = rw_copy(parent)
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
                        products = list(sanitized_fragments(rw, counters).pieces)
                        if not products:
                            continue
                        info = dict(preview)
                        info["path"] = tuple(path)
                        yield ProductsOfReaction(info=info, products=products)
                        break


class Hydroxylation(SmartsReactionRule):
    """Adds a hydroxyl to carbon.

    Both patterns add OH and remove one H. ``[#6h]`` is h=1, 2, or 3;
    ``[#6h2]`` is the subset with at least two hydrogens. The match records
    which of those the atom actually is.
    """

    phase1_sites_on = "atom_hydrogen"
    sites_on = "atom_hydrogen"

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
                ),
                name="h",
            ),
        ),
        (
            "[#6h2:1]>>[*:1]O",
            describe(
                *branches(
                    ({"map": 1, "z": 6, "h": 2}, {"map": 1, "z": 6, "h": 3}),
                    adds="O",
                    removes="H",
                ),
                name="h2",
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

    phase1_sites_on = "atom_hydrogen"
    sites_on = "atom_hydrogen"

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


def _whens(mapno, atomic_nums):
    """One branch constraint per atomic number, for :func:`branches`."""

    return tuple({"map": mapno, "z": int(z)} for z in atomic_nums)


class Dealkylation(SmartsReactionRule):
    """Cleaves a C-N, C-O, C-S, or C-C bond and oxygenates the carbon side.

    The site is both atoms of the broken bond.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="OO", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="OO", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6:1][#6:2]>>(O-[*:1].[*:2])",
            describe(adds="O", cleaves=True, partner="C", site_map=(1, 2)),
        ),
        (
            "[#6h:1][#6:2]>>(O-[*:1].[*:2])",
            describe(
                *branches(
                    (
                        {"map": 1, "z": 6, "h": 1},
                        {"map": 1, "z": 6, "h": 2},
                        {"map": 1, "z": 6, "h": 3},
                    ),
                    adds="O",
                    cleaves=True,
                    partner="C",
                ),
                site_map=(1, 2),
            ),
        ),
        (
            "[#6h:1][#6:2]>>(O=[*:1].[*:2])",
            describe(
                *branches(
                    (
                        {"map": 1, "z": 6, "h": 1},
                        {"map": 1, "z": 6, "h": 2},
                        {"map": 1, "z": 6, "h": 3},
                    ),
                    adds="O",
                    cleaves=True,
                    partner="C",
                ),
                site_map=(1, 2),
            ),
        ),
        (
            "[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])",
            describe(
                *branches(_whens(2, (7, 8, 16)), removes="H", cleaves=True),
                site_map=(1, 2),
            ),
        ),
    )


def _ndealk(smarts: str, leave_count: int | None, **effect) -> tuple[str, PatternInfo]:
    """One N-dealkylation pattern. ``leave_count`` is the named leaving atoms."""

    return (
        smarts,
        describe(
            cleaves=True,
            partner="N",
            leave_count=leave_count,
            site_map=(1, 2),
            **effect,
        ),
    )


class NDealkylation(SmartsReactionRule):
    """Cleaves a carbon-nitrogen bond and oxygenates the carbon side.

    These are the nitrogen rows of dealkylation. A methyl carbon is the whole
    leaving piece (``leave_count`` 1). Any larger alkyl still carries atoms
    the pattern does not name, so ``leave_count`` is None. ``breaks_ring`` is
    filled from the cleaved bond, not stored as a second class.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        _ndealk("[#6H3:1][#7:2]>>([*:2].[*:1](=O)O)", 1, adds="OO"),
        _ndealk("[#6H3:1][#7:2]>>([*:2].[*:1]=O)", 1, adds="O"),
        _ndealk("[#6H3:1][#7:2]>>([*:2].[*:1]-O)", 1, adds="O"),
        _ndealk("[#6H2:1][#7:2]>>([*:2].[*:1](=O)O)", None, adds="OO"),
        _ndealk("[#6H2:1][#7:2]>>([*:2].[*:1]=O)", None, adds="O"),
        _ndealk("[#6H2:1][#7:2]>>([*:2].[*:1]-O)", None, adds="O"),
        _ndealk("[#6H1:1][#7:2]>>([*:2].[*:1]=O)", None, adds="O"),
        _ndealk("[#6H1:1][#7:2]>>([*:2].[*:1]-O)", None, adds="O"),
        _ndealk("[#6H0:1][#7:2]>>([*:2].[*:1]-O)", None, adds="O"),
        _ndealk(
            "[#8H1:3]-[#6:1]-[#7:2]>>([*:3]=[*:1].[*:2])",
            None,
            removes="H",
        ),
    )


class AzoSplitting(SmartsReactionRule):
    """Splits an N=N bond. Both fragments stay.

    The pattern names both nitrogens and no leaving piece, so ``leave_count``
    stays None. ``breaks_ring`` is filled from the cleaved bond. A ring N=N
    and an open azo are the same pattern; a filter reads ``breaks_ring``.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#7:1]=[#7:2]>>[*:1].[*:2]",
            describe(cleaves=True, partner="N", site_map=(1, 2)),
        ),
    )


class BenzodioxoleReduction(SmartsReactionRule):
    """Cleaves both C-O bonds of the methylene in a 1,3-dioxole.

    That carbon is the whole leaving piece, and the pattern names it, so
    ``leave_count`` is 1. ``breaks_ring`` is filled from one cleaved bond.
    Both bonds are in that ring. A filter reads ``leave_count``.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6R:1]-[#8R:2]-[#6H2R:3]-[#8R:4]-[#6R:5]>>([*:1]-[*:2].[*:3].[*:4]-[*:5])",
            describe(
                cleaves=True,
                partner="O",
                leave_count=1,
                site_map=(2, 3),
            ),
        ),
    )


class NitroaromaticReduction(SmartsReactionRule):
    """Cleaves one N-O of a nitro group on a ring carbon, leaving the nitroso.

    That oxygen is the whole leaving piece, and both patterns name it, so
    ``leave_count`` is 1. ``breaks_ring`` is filled from the cleaved bond.
    A filter reads ``leave_count``.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8-1:1]-[#7+1:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
            describe(
                cleaves=True,
                partner="N",
                leave_count=1,
                site_map=(1, 2),
            ),
        ),
        (
            "[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
            describe(
                cleaves=True,
                partner="N",
                leave_count=1,
                site_map=(1, 2),
            ),
        ),
    )


class ThiopheneSulfurOxidation(SmartsReactionRule):
    """Oxidizes the sulfur of a thiophene to the S-oxide.

    The pattern adds oxygen and names no leaving atom, so ``leave_count``
    stays None. It does not cleave. A filter reads ``adds``.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:2]1=[#6:3][#6:4]=[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*&H0&+:1]1[O-]",
            describe(adds="O", symbol="S"),
        ),
    )


_HALIDE = (9, 17, 35, 53, 85)


class Dephosphorylation(SmartsReactionRule):
    """Cleaves an O-P bond of a phosphate. The oxygen stays on the organic fragment."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8:1][#15:2](=[#8:3])([#8:4])[#8:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]",
            describe(*branches(_whens(2, (15,)), cleaves=True)),
        ),
    )


class EpoxideOpening(SmartsReactionRule):
    """Opens an epoxide. One pattern only rearranges bonds; the other also adds OH."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])",
            describe(adds=""),
        ),
        (
            "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)",
            describe(adds="O"),
        ),
    )


class Hydrolysis(SmartsReactionRule):
    """Cleaves the single bond of a carboxylic derivative. One pattern also adds O."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2](O).[*:3])",
            describe(
                *branches(
                    _whens(3, (7, 8, 16)),
                    site_map=2,
                    adds="O",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2].[*:3])",
            describe(
                *branches(
                    _whens(3, (7, 8, 16)),
                    site_map=2,
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
    )


class Dehydration(SmartsReactionRule):
    """Drops an OH, or a carbonyl oxygen, off carbon or nitrogen."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]",
            describe(
                *branches(
                    ({"map": 1, "z": 6}, {"map": 1, "z": 7}),
                    removes="OH",
                    cleaves=True,
                    partner="O",
                )
            ),
        ),
        (
            "[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]",
            describe(removes="OH", cleaves=True, partner="O"),
        ),
        (
            "[#6,#7:1]=[#8:2]>>[*:1].[*:2]",
            describe(
                *branches(
                    ({"map": 1, "z": 6}, {"map": 1, "z": 7}),
                    removes="O",
                    cleaves=True,
                    partner="O",
                )
            ),
        ),
    )


class Hydrogenation(ResonancePairRule):
    """Reduces C#C to C=C, C=C to C-C, and a conjugated pair across the path.

    The double-bond pattern is ``=,:``, so an aromatic bond matches once.
    The pair names ``keep``: nothing changes at the end except the path
    flip, which adds H where a double bond becomes single. A carbon on a
    triple bond is not an end. That reduction is the ``#`` SMARTS.
    Heavy-atom formula is unchanged (``adds`` is ``HH``).
    """

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]#[#6:2]>>[*:1]=[*:2]",
            describe(adds="HH", partner="C"),
        ),
        (
            "[#6:1]=,:[#6:2]>>[*:1]-[*:2]",
            describe(adds="HH", partner="C"),
        ),
    )
    endpoints: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6;!$(*#[#6]):1]",
            describe(adds="H", partner="C", edit="keep"),
        ),
    )


class NitrogenReduction(SmartsReactionRule):
    """Cleaves N-O of nitro, nitroso, and hydroxylamine groups.

    The nitroso pattern uses ``[*:2]``. The old ``[*2]`` string emitted a dummy atom.
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O"),
        ),
        (
            "[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O"),
        ),
        (
            "[#7:1]-[#8:2]>>([*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#7D2:1]=[#8:2]>>([*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O"),
        ),
    )


class OxygenReduction(SmartsReactionRule):
    """Turns C=O / N=O into a single bond, or cleaves a peroxide."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8:1]=[#6,#7:2]>>[*:1]-[*:2]",
            describe(
                *branches(
                    ({"map": 2, "z": 6}, {"map": 2, "z": 7}),
                    adds="HH",
                )
            ),
        ),
        (
            "[#8:1]-[#8:2]>>[*:1].[*:2]",
            describe(cleaves=True, partner="O"),
        ),
    )

# TODO: In chemistry (not current smarts), can ReductiveDehalogenation ever work on
# an aromatic bond? If so, maybe this should be a ResonanceRule.
class ReductiveDehalogenation(SmartsReactionRule):
    """Cleaves a carbon-halogen bond. The second pattern also makes a double bond."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
    )


class SulfurReduction(SmartsReactionRule):
    """Cleaves S=O, S-S, and S-C / S-O single bonds."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#16:1]=[#8:2]>>[*:1].[*:2]",
            describe(removes="O", cleaves=True, partner="O"),
        ),
        (
            "[#16:1]-[#16:2]>>[*:1].[*:2]",
            describe(cleaves=True, partner="S"),
        ),
        (
            "[#16:1]-[#6,#8:2]>>[*:1].[*:2]",
            describe(
                *branches(({"map": 2, "z": 6},), cleaves=True),
                *branches(({"map": 2, "z": 8},), cleaves=True, removes="O"),
            ),
        ),
    )



class Epoxidation(ResonanceRule):
    """Adds an epoxide across a C=C or C=N bond.

    The reactant bond is ``=,:``, so an aromatic bond matches on the parent.
    The reaction runs on the cached kekulé parent for that bond.

    NOTE: Downstream epoxidation model only considers carbone-carbon epoxides.
    Phase 1 model additionally considers carbon-nitrogen.
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]=,:[#6,#7:2]>>[*:1]1-[*:2][O]1",
            describe(
                *branches(
                    ({"map": 2, "z": 6}, {"map": 2, "z": 7}),
                    adds="O",
                )
            ),
        ),
    )


class SulfurOxidation(SmartsReactionRule):
    """Adds oxygen to divalent or tetravalent sulfur (S-oxide, S-OH, or S=O)."""

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#16;v2,v4:1]>>[*&H0&+:1][O-]",
            describe(adds="O", symbol="S"),
        ),
        (
            "[#16;v2,v4:1]>>[*:1][O]",
            describe(adds="O", symbol="S"),
        ),
        (
            "[#16;v2,v4:1]>>[*:1]=O",
            describe(adds="O", symbol="S"),
        ),
    )


class NitrogenOxidation(SmartsReactionRule):
    """N-H to hydroxylamine, primary amine to nitroso, or tertiary N to N-oxide."""

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#7v3h:1]>>[*:1]O",
            describe(
                *branches(
                    ({"map": 1, "z": 7, "h": 1}, {"map": 1, "z": 7, "h": 2}),
                    adds="O",
                )
            ),
        ),
        (
            "[#7v3H2:1]>>[*:1]=O",
            describe(adds="O", h=2, symbol="N"),
        ),
        (
            "[#7v3H0:1]>>[*&H0&+:1][O-]",
            describe(adds="O", h=0, symbol="N"),
        ),
    )


class OxidativeDehalogenation(SmartsReactionRule):
    """Replaces a carbon-bound halogen with OH, carbonyl, or a carboxylic acid."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    adds="O",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    adds="O",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    adds="OO",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    adds="O",
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    adds="OO",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
        (
            "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]",
            describe(
                *branches(
                    _whens(1, _HALIDE),
                    site_map=2,
                    removes_partner=True,
                    adds="OO",
                    cleaves=True,
                ),
                site_map=2,
            ),
        ),
    )


# No adduct molecule for these labels. They stay stars, and only on this rule.
_STAR_ONLY = frozenset({"Protein", "DNA", "Cyanide"})


def _has_star_conjugate(mol: Mol) -> bool:
    return any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())


def _collapse_conjugate_to_star(product: Mol) -> Mol:
    """Replace atoms the reaction added with one ``*`` at each attachment."""

    new_atoms: set[int] = set()
    for atom in product.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.HasProp("react_atom_idx") or atom.HasProp("current_idx"):
            continue
        new_atoms.add(atom.GetIdx())
    if not new_atoms:
        return Mol(product)

    attach: set[int] = set()
    for bond in product.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        begin_new = begin in new_atoms
        end_new = end in new_atoms
        if begin_new == end_new:
            continue
        attach.add(end if begin_new else begin)
    if not attach:
        return Mol(product)

    def _after_remove(idx: int) -> int:
        return idx - sum(1 for removed in new_atoms if removed < idx)

    rw = rw_copy(product)
    for idx in sorted(new_atoms, reverse=True):
        rw.RemoveAtom(idx)
    for parent in sorted(_after_remove(atom) for atom in attach):
        dummy = rw.AddAtom(Atom(0))
        rw.AddBond(parent, dummy, BondType.SINGLE)
    collapsed = rw.GetMol()
    sanitize_catch(collapsed)
    return collapsed


def _apply_star_conjugate(mol: Mol, label: str | None = None) -> Mol:
    collapsed = _collapse_conjugate_to_star(mol)
    if not label:
        return collapsed
    for atom in collapsed.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetProp("atomLabel", label)
    return collapsed


class ConjugationRule(SmartsReactionRule):
    """Attaches an acetyl and, by default, collapses that group to ``*``.

    ``as_star`` and ``star_label`` belong on this class. ``Protein``,
    ``DNA``, and ``Cyanide`` stay stars. The site heteroatom is ``symbol``.
    A filter reads that. This reaction does not cleave.
    """

    as_star: bool = True
    star_label: str | None = None
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]",
            describe(
                *branches(_whens(1, (7, 8, 16)), adds="CCO", removes="H"),
                pin=(1,),
            ),
        ),
    )

    def __init__(
        self,
        as_star: bool | None = None,
        star_label: str | None = None,
        *args,
        **kwargs,
    ):
        if as_star is None:
            as_star = type(self).as_star
        if star_label is None:
            star_label = type(self).star_label
        if star_label in _STAR_ONLY and as_star is False:
            raise ValueError(
                "%r cannot emit a full conjugate structure; use a star adduct"
                % (star_label,)
            )
        if star_label in _STAR_ONLY:
            as_star = True
        self.as_star = as_star
        self.star_label = star_label
        super().__init__(*args, **kwargs)

    def is_terminal_product(self, mol: Mol) -> bool:
        if _has_star_conjugate(mol):
            return True
        return super().is_terminal_product(mol)

    def metabolites(
        self,
        mol: Mol,
        filter_rules=lambda rule, info: True,
        filter_sites=lambda site, info: True,
        context_mol: Mol | None = None,
        **kwargs,
    ):
        for por in super().metabolites(
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            context_mol=context_mol,
            **kwargs,
        ):
            products = por.products
            if self.as_star:
                products = [
                    _apply_star_conjugate(product, self.star_label)
                    for product in products
                ]
            yield ProductsOfReaction(info=por.info, products=products)


class Acetylation(ConjugationRule):
    """Adds an acetyl to OH, NH, and SH.

    Calls :class:`ConjugationRule`. The acetyl SMARTS and the star collapse
    stay there. ``as_star=False`` keeps the acetyl. A filter reads ``symbol``.
    """


class Sulfation(ConjugationRule):
    """Adds a sulfate to an alcohol or a phenol.

    Calls :class:`ConjugationRule`. The sulfate SMARTS live here. The star
    collapse stays there. ``as_star=False`` keeps the sulfate. The oxygen
    has one hydrogen. A filter reads ``partner_h``.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1][#8H1:2]>>[*:1][*:2]S(=O)(=O)O",
            describe(
                *branches(
                    _whens(1, (6,)),
                    site_map=2,
                    adds="SOOO",
                    removes="H",
                ),
                site_map=2,
            ),
        ),
    )


class Glucuronidation(ConjugationRule):
    """Adds a glucuronide to an alcohol, a phenol, or a carboxylic oxygen.

    Calls :class:`ConjugationRule`. The glucuronide SMARTS live here. The
    star collapse stays there. ``as_star=False`` keeps the glucuronide. The
    carbon beside an alcohol oxygen is the partner. A filter reads
    ``partner_h``. The carbonyl pattern matches ``=[#8]`` only.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8H1:1][#6:2]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2])1",
            describe(
                *branches(
                    _whens(2, (6,)),
                    site_map=1,
                    adds="CCCCCCOOOOOO",
                    removes="H",
                ),
                site_map=1,
            ),
        ),
        (
            "[#8H1,#8-:1][#6:2](=[#8:3])[#6:4]>>"
            "O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2](=[#8:3])[*:4])1",
            describe(adds="CCCCCCOOOOOO", site_map=1),
        ),
    )


# Glutathione, cysteine sulfur open at ``{attach}``. Heavy atoms only: C10 N3 O6 S.
_GSH = "C(CC(=O)N[C@@H](CS({attach}))C(=O)NCC(=O)O)[C@@H](C(=O)O)N"
_GSH_ADDS = "CCCCCCCCCCNNNOOOOOOS"


def _gsh(attach: str) -> str:
    return _GSH.format(attach=attach)


class Glutathionation(ConjugationRule):
    """Adds glutathione to a soft electrophile.

    Calls :class:`ConjugationRule`. The glutathione SMARTS live here. The
    star collapse stays there. ``as_star=False`` keeps the peptide. A
    carbon-halogen site names the halogen as ``partner``. A filter reads
    ``partner``. Ring carbons are separate patterns: one query of the three
    ring atoms keeps a single embedding.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6H1:1]1[#8:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6H2:1]1[#8:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6:1]([!#1:4])1[#8:2][#6:3]1>>" + _gsh("[*:1]([*:4])[*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6:1][#9,#17,#35,#53:2]>>" + _gsh("[*:1]"),
            describe(
                *branches(
                    _whens(2, (9, 17, 35, 53)),
                    site_map=1,
                    adds=_GSH_ADDS,
                    removes_partner=True,
                ),
                site_map=1,
            ),
        ),
        (
            "[#16h1:1]>>" + _gsh("[*:1]"),
            describe(adds=_GSH_ADDS, removes="H", site_map=1),
        ),
        (
            "[#6H2:1]=[#6:2]>>" + _gsh("[*:1]-[*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6H1:1]=[#6:2][#6:3]=[#8,#7:4]>>" + _gsh("[*:1][*:2]=[*:3][*:4]"),
            describe(
                *branches(_whens(4, (8, 7)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
            ),
        ),
        (
            "[#6;H1,H2:1]=[#8:2]>>" + _gsh("[*:1]([*:2])"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6H1:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6H2:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6:1]([!#1:4])1[#7:2][#6:3]1>>" + _gsh("[*:1]([*:4])[*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#6:1][#8:2]S(=O)(=O)>>" + _gsh("[*:1]"),
            describe(adds=_GSH_ADDS, site_map=1),
        ),
        (
            "[#7:1]=[#6:2]=[#8,#16:3]>>" + _gsh("[*:2](=[*:3])[*:1]"),
            describe(
                *branches(_whens(3, (8, 16)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
            ),
        ),
    )
