"""Define specific reaction rules."""

from __future__ import annotations

# Standard Library
import copy
import itertools
from collections import defaultdict, deque
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, NamedTuple, TypeAlias, cast

if TYPE_CHECKING:
    # Static re-export so ``from .rules import RuleSet`` types correctly.
    # Runtime still goes through ``__getattr__`` to avoid the cycle.
    from .rulesets import RuleSet as RuleSet

from xenosite.refactor_poc.canonical_plan import (
    CanonicalStep,
    identity_canonical_plan,
    quinone_canonical_plan,
)
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
    Mol,
    RWMol,
    TracingMol,
    _bond_key,
    _current_bond_map,
    aromatic_parent_atoms,
    copy_mol,
    ensure_kekule_parents,
    move_charge_with_bonds,
    parent_for_bond,
    parents_for_ends,
    reaction_from_smarts,
    resonance_bond_maps,
    run_reactants,
    rw_copy,
    sanitize_catch,
    sanitized_fragments,
)
from xenosite.refactor_poc.records import (
    AtomPairOrbitSignature,
    EditCounters,
    Effect,
    EffectField,
    Formula,
    InitializedAtomTrace,
    KekuleParents,
    PairSiteInfo,
    PatternInfo,
    ProductInfo,
    Site,
    SiteInfo,
    SitesOn,
    SmartsSiteInfo,
    Span,
    TraceAddition,
    TraceInfo,
    When,
)

# Dedup key for one SMARTS site across Kekulé forms / equivalent carbons.
# Last field: AtomPairOrbitSignature, or None for a one-atom site.
SiteSignature: TypeAlias = tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int, float], ...],
    int,
    str | None,
    str | None,
    bool,
    bool,
    AtomPairOrbitSignature | None,
]


class _LazyProductInfo(dict):
    """Product ``info`` view. ``csmi`` is ``product.xf.csmi``.

    The hot path does not compute SMILES. Dedup, export, and ``info["csmi"]``
    read ``mol.xf.csmi``, which caches on ``structure["csmi"]``.
    """

    __slots__ = ("_mol",)

    def __init__(self, data: dict, mol: ForestMol):
        super().__init__(data)
        self._mol = mol

    def __getitem__(self, key):
        if key == "csmi":
            return self._mol.xf.csmi
        return super().__getitem__(key)

    def get(self, key, default=None):  # type: ignore[override]
        if key == "csmi":
            return self._mol.xf.csmi
        return super().get(key, default)

    def __contains__(self, key):
        return key == "csmi" or super().__contains__(key)


def _copy_when(raw: When | Mapping[str, int]) -> When:
    copied: When = {}
    mapno = raw.get("map")
    if isinstance(mapno, int):
        copied["map"] = mapno
    atomic = raw.get("z")
    if isinstance(atomic, int):
        copied["z"] = atomic
    hydrogens = raw.get("h")
    if isinstance(hydrogens, int):
        copied["h"] = hydrogens
    return copied


def set_terminal_product(mol: Mol, value: bool = True) -> ForestMol:
    """Mark ``mol`` terminal on the forest. Prefer ``mol.xf._mark_terminal``."""

    held = mol.xf.forestmol
    held.xf._mark_terminal(value)
    return held


class ProductsOfReaction(NamedTuple):
    """One edit from :meth:`ReactionRule.metabolites`, before tracing.

    ``info`` describes the edit. ``products`` are the mols it made.
    :meth:`ReactionRule.metabolize` turns each of those mols into a
    ``(product, info)`` pair and is what callers should use.
    """

    info: SiteInfo
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
    sites_on: SitesOn | None = None

    def _clear_atom_maps(self, mol: Mol) -> Mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    def __init__(
        self,
        name: str | None = None,
        sites_on: SitesOn | None = None,
        longname: str | None = None,
        *args,
        **kwargs,
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

    def canonical_plan(self, mol: Mol, info: SiteInfo) -> tuple[CanonicalStep, ...]:
        """Elementary steps a search should record for this hop.

        The default is identity: this rule is already canonical, so the plan
        is one step at ``info["site"]``. That matches forest
        ``phase1_equivalent`` singletons (including Epoxidation and
        NDealkylation). QuinoneFormation overrides with hydroxylation then
        dehydrogenation. ``Addition.phase1`` is not used; its schema is not
        decided.
        """

        if self.name is None:
            return ()
        return identity_canonical_plan(self.name, mol, info["site"])

    def __call__(
        self, mol: Mol, **kwargs
    ) -> Generator[tuple[TracingMol, ProductInfo], None, None]:
        yield from self.metabolize(mol, **kwargs)

    def __iter__(self) -> Iterator[ReactionRule]:
        return iter([self])

    def metabolize(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        unique_csmi: bool = True,
        **kwargs,
    ) -> Generator[tuple[TracingMol, ProductInfo], None, None]:
        """Apply this rule and yield ``(product, info)`` pairs.

        This is the method callers use. It keeps the invariants below.
        :meth:`metabolites` supplies the chemistry and must not try to.

        The mol that was passed in:

        - Its atoms, bonds, charges, hydrogens, and atom-map numbers are
          unchanged.
        - If it has no ``_forest`` trace, one is installed and left there.
        - If it already has a trace, that trace's depth is not changed.

        Every product:

        - Is a sanitized, connected mol with its own ``_forest``.
        - Has ``atom_trace["depth"]`` one greater than the parent.
        - Has ``atom_trace["formula"]`` equal to the atom counts and formal
          charge of that product, hydrogens included.
        - Records each new transform once, as ``R1``, ``R2``, and so on.
          A new atom's ``added_by`` is that id. The site, the rule
          hierarchy, the resolved effect, the name, ``phase1``, and the
          depth of the site's index frame live under
          ``atom_trace["additions"][id]``. The change in formula lives
          under ``atom_trace["delta_formula"][id]``.
        - Reads ``info["csmi"]`` (or ``product.xf.csmi``)
          as the canonical SMILES of that product; computed on demand and
          cached on the product forest. With ``unique_csmi`` (default), the
          key is ``(rule name, PatternInfo.name | SMARTS, product csmi)`` —
          site topology is not part of it. Prefers ``info["pattern"]["name"]``;
          falls back to the SMARTS string for that pattern. Pair emissions
          without a pattern use ``None`` for the middle field.

        ``filter_rules(mol, rule, pattern_info)`` sees the pattern before a
        match. ``filter_sites(mol, site, info)`` sees the resolved effect
        before an edit. Either may refuse. A refusal edits nothing.
        """
        if mol is None:
            raise ValueError("mol is required")
        # The caller's chemistry is not edited. A mol with no trace gets one,
        # at its current depth, so products can sit one step below it.
        mol = mol.xf.tracing._stamp()
        # Matching, map clearing, and forest stamps happen on a copy.
        mol = _work_copy(mol).xf.tracing._stamp()

        if self.is_terminal_product(mol):
            return

        # Maps must not be present for CanonicalRankAtoms or SMARTS matching.
        self._clear_atom_maps(mol)

        # Product layer: (rule name, PatternInfo.name | SMARTS, csmi).
        # Site topology lives only in unique-edit upstream.
        seen: set[tuple[str, str | None, str]] = set()

        for por in self.metabolites(
            mol,
            filter_rules=filter_rules,
            filter_sites=filter_sites,
            **kwargs,
        ):
            info = por.info
            # react → split → trace. Cleavage pieces are separate mols before
            # tracing; each fragment gets its own atom_trace from the parent
            # (siblings are absent, not deleted from a shared pre-split
            # trace). A dotted leftover is expanded here as a safety net.
            products: list[Mol] = []
            for raw in por.products:
                pieces = list(sanitized_fragments(raw).pieces)
                if not pieces:
                    products = []
                    break
                products.extend(pieces)
            if not products:
                continue

            for p in products:
                p.xf.sanitize()
                # Terminal marking is of_products (reads is_terminal_rule).

            # Stamp + trace + clear structure caches on each fragment.
            finished = mol.xf.of_products(products, info, executed=self)

            # Same (rule, pattern, product csmi) is one outcome. Two sites in
            # one atom class can still be different molecules (ortho / para).
            for n, p in enumerate(finished):
                assert p.xf.tracing.active

                if unique_csmi:
                    key = _unique_csmi_key(info, p.xf.csmi)
                    if key in seen:
                        continue
                    seen.add(key)

                i = cast(
                    ProductInfo,
                    _LazyProductInfo(
                        {
                            **dict(info),
                            "product_index": n,
                            "product_count": len(finished),
                        },
                        p,
                    ),
                )

                yield p, i

    def _top_site(self, site: Site, mol: Mol) -> Site:
        te = mol.xf.topol_equiv
        if isinstance(site, int):
            return te[site]
        if isinstance(site, tuple):
            return frozenset(int(te[s]) for s in site)
        if isinstance(site, frozenset):
            nested: list[frozenset[int]] = []
            flat: list[int] = []
            for item in site:
                if isinstance(item, frozenset):
                    nested.append(frozenset(int(te[index]) for index in item))
                elif isinstance(item, int):
                    flat.append(int(te[item]))
            if nested:
                return frozenset(nested)
            return frozenset(flat)
        return site

    def is_terminal_product(self, mol: Mol) -> bool:
        """True if ``mol`` must not be expanded further in guided path search.

        Reads the forest-level ``is_terminal_product`` flag (survives
        ``clear_structure``). Prefer ``mol.xf.is_terminal`` at call sites that
        already hold a :class:`ForestMol`.
        """

        return mol.xf.is_terminal

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        **kwargs,
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
        - ``products`` is a list of connected mols. Cleavage puts each
          fragment in that list before :meth:`metabolize` runs
          ``forest_trace``. It does not return one mol that is several pieces.

        ``filter_rules(mol, rule, pattern_info)`` is called before a match and
        can see the pattern's ``span``. False means that pattern is skipped.
        ``filter_sites(mol, site, info)`` is called after the effect is resolved
        and before any edit. False means that site is skipped.

        This method must not edit the mol it is given. ``metabolize`` has
        already handed it a copy. It must not attach ``_forest`` to the
        caller's mol, and it must not set product depth. ``metabolize``
        does both after this method returns.
        """
        raise NotImplementedError


FilterRules = Callable[[TracingMol, ReactionRule, PatternInfo], bool]
FilterSites = Callable[[TracingMol, Site, SiteInfo], bool]


def __getattr__(name: str) -> type[ReactionRule]:
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


def _rule_name(rule: ReactionRule | str | None) -> str | None:
    if rule is None:
        return None
    if isinstance(rule, str):
        return rule
    return getattr(rule, "name", type(rule).__name__)


def _work_copy(mol: Mol) -> Mol:
    """A mol the rule may stamp. The caller's object is left alone."""

    return copy_mol(mol)


def stamp_forest_labels(mol: Mol) -> TracingMol:
    """Ensure tracing and stamp labels. Prefer ``mol.xf.tracing._stamp()``.

    Thin shim for tests that still take a bare ``Mol``. ``install_forest`` /
    ``ensure_tracing`` / ``install_forest`` are gone — use ``mol.xf.tracing._stamp``.
    """

    return mol.xf.tracing._stamp()


def install_forest(mol: Mol) -> TracingMol:
    """Deprecated alias of :func:`stamp_forest_labels` (tests)."""

    return stamp_forest_labels(mol)


def reordered_forest_labels(mol: TracingMol) -> None:
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


def _as_site(value: Site | None) -> Site:
    """Narrow a known-index :class:`Site`. Resolve AtomRef leaves first."""

    if isinstance(value, int):
        return value
    if isinstance(value, tuple) and all(isinstance(item, int) for item in value):
        return value
    if isinstance(value, frozenset):
        if all(isinstance(item, int) for item in value):
            return value
        if all(
            isinstance(item, frozenset) and all(isinstance(inner, int) for inner in item)
            for item in value
        ):
            return value
    raise TypeError(value)


def _site_tuple(site: Site) -> Site:
    if isinstance(site, int):
        return (site,)
    if isinstance(site, tuple):
        return tuple(sorted(site))
    sample = next(iter(site), None)
    if isinstance(sample, frozenset):
        return site
    indexes = [item for item in site if isinstance(item, int)]
    return tuple(sorted(indexes))


def _trace_info(info: SiteInfo) -> TraceInfo:
    """Pattern fields worth keeping, without a second copy of the rule object."""

    kept: TraceInfo = {
        "site": info["site"],
        "rule": _rule_name(info["rule"]),
    }
    if "rxn_num" in info:
        kept["rxn_num"] = info["rxn_num"]
    if "ends" in info:
        kept["ends"] = info["ends"]
        kept["end_atoms"] = info["end_atoms"]
        kept["end_maps"] = info["end_maps"]
        kept["path_ends"] = info["path_ends"]
    return kept


def _as_effect(value: Effect | Mapping[str, EffectField] | None) -> Effect:
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
        effect["when"] = _copy_when(when)
    return effect


def forest_trace(
    reactant: Mol,
    product: Mol,
    info: SiteInfo,
    executed: ReactionRule | None = None,
) -> InitializedAtomTrace:
    """Record one transform on the product's atom trace.

    Thin wrapper over ``product.xf.tracing._trace(...)``. Prefer
    ``reactant.xf.of_products(product, info)`` for the finishing path.
    """

    return product.xf.tracing._trace(reactant, info, executed=executed)


def _apply_forest_trace(
    reactant: Mol,
    product: Mol,
    info: SiteInfo,
    executed: ReactionRule | None = None,
) -> InitializedAtomTrace:
    """Record one transform on the product's atom trace.

    A new atom's ``added_by`` is an id such as ``R1``. The site, the rule
    hierarchy, the resolved effect, and the formula change live once, under
    ``atom_trace["additions"][id]``. ``depth`` is the reactant index frame
    the site is written in.
    """

    rule = info["rule"]
    site = _as_site(info["site"])
    parent = reactant.xf.tracing._stamp()
    held = product.xf.forestmol
    trace: InitializedAtomTrace = copy.deepcopy(parent._forest["atom_trace"])
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
    if executed is not None and executed is not rule:
        chain.append(executed)
    if all(rule is not seen for seen in chain):
        chain.append(rule)

    before = trace.get("formula") or reactant.xf.formula
    after = held.xf.formula
    trace["formula"] = after
    trace["delta_formula"][transform_id] = formula_delta(before, after)
    # PatternInfo is stored on info as the same dict the rule holds.
    pattern: PatternInfo | None
    if "pattern" in info:
        pattern = info["pattern"]
    else:
        pattern = None
    addition: TraceAddition = {
        "site": _site_tuple(site),
        "rules": tuple(chain),
        "info": _trace_info(info),
        "effect": _as_effect(info["options"]),
        "name": _rule_name(rule),
        "phase1": None,
        "depth": frame,
        "pattern": pattern,
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


_EFFECT_DEFAULTS: Effect = {
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


def _span(possibilities: Sequence[Effect]) -> Span:
    """Certain values stay bare. Disagreeing values become a tuple."""

    keys: list[str] = []
    for possibility in possibilities:
        for key in possibility:
            if key != "when" and key not in keys:
                keys.append(key)
    raw: dict[str, EffectField | tuple[EffectField, ...]] = {}
    for key in keys:
        values: list[EffectField] = []
        for possibility in possibilities:
            value = cast(EffectField, possibility.get(key, _EFFECT_DEFAULTS.get(key)))
            if value not in values:
                values.append(value)
        raw[key] = values[0] if len(values) == 1 else tuple(values)
    return cast(Span, raw)


def branches(
    whens: Sequence[When],
    site_map: int = 1,
    removes_partner: bool = False,
    **effect: EffectField,
) -> tuple[Effect, ...]:
    """Copy ``effect`` once per ``when``. The site atom and the OR atom differ.

    A ``when`` on ``site_map`` records ``h`` / ``symbol`` (the site was the
    ambiguous atom). A ``when`` on any other map records ``partner`` /
    ``partner_h``.
    """

    out: list[Effect] = []
    for raw in whens:
        when = _copy_when(raw)
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
    **single: EffectField,
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


def may(info: PatternInfo, key: str, value: EffectField = True) -> bool:
    """True if any possibility has this outcome.

    For ``adds`` / ``removes`` / ``needs``, ``value`` may be a substring.
    """

    for possibility in info.get("possibilities") or ():
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


def must(info: PatternInfo, key: str, value: EffectField = True) -> bool:
    """True if every possibility has this outcome."""

    possibilities = info.get("possibilities") or ()
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

    possibilities = info.get("possibilities") or (_as_effect({}),)
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
    when = chosen.get("when")
    branch = when.get("map") if when is not None else None
    if isinstance(branch, int) and branch in mapped and branch != site_map:
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
    rings = mol.xf.rings
    return bool(set(rings.get(left, ())) & set(rings.get(right, ())))


def merge_effects(
    left: Effect, right: Effect, system_aromatic: bool
) -> Effect:
    """One effect for a pair. Per-end detail stays on ``info["ends"]``.

    ``system_aromatic`` is whether this conjugated system contains any
    aromatic atom. A pattern may be able to dearomatize; the resolved
    effect does so only when the system has aromaticity. A system with
    none is the same edit, with ``dearomatizes`` false. That flag is
    what filters read on ``SiteInfo["options"]``.
    """

    needs = (left.get("needs") or "") + (right.get("needs") or "")
    can = left.get("dearomatizes") or right.get("dearomatizes")
    return {
        "adds": left.get("adds", "") + right.get("adds", ""),
        "removes": left.get("removes", "") + right.get("removes", ""),
        "cleaves": bool(left.get("cleaves") or right.get("cleaves")),
        "dearomatizes": bool(can and system_aromatic),
        "methide": bool(left.get("methide")) ^ bool(right.get("methide")),
        "needs": needs,
    }


def _bump(counters: EditCounters | None, name: str, amount: int = 1) -> None:
    if counters is None:
        return
    current = getattr(counters, name)
    if not isinstance(current, int):
        raise TypeError(name)
    setattr(counters, name, current + amount)


def _split_smarts_or(body: str) -> list[str]:
    """Split a SMARTS atom body on top-level ``,``.

    Commas inside parentheses stay put. Those are the OR arms ``&`` would
    otherwise bind to only the last of.
    """

    parts: list[str] = []
    depth = 0
    start = 0
    for index, char in enumerate(body):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == "," and depth == 0:
            parts.append(body[start:index])
            start = index + 1
    parts.append(body[start:])
    return parts


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
        # ``,`` binds looser than ``&``. The isotope has to sit on every
        # alternative, or only the last OR arm is pinned and RunReactants
        # returns a different arm first.
        parts = _split_smarts_or(match.group(1))
        body = ",".join(part + "&" + str(isotope) + "*" for part in parts)
        return "[" + body + ":" + match.group(2) + "]"

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
    counters: EditCounters | None = None,
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
        lifted = rule._lift_forest_labels(mol, prod)
        # Forest clean: one connected mol per piece, never a dotted product.
        pieces = list(sanitized_fragments(lifted, counters).pieces)
        if not pieces:
            return []
        products.extend(pieces)
    return products


class SmartsReactionRule(ReactionRule):
    """Performs reactions specified by SMARTS.

    Each entry is ``(smarts, options)`` with :class:`PatternInfo` options.
    Matches are filtered before ``RunReactants``. One SMARTS is applied once
    per topological site; a later SMARTS with the same formula effect still runs.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = ()

    rxns: list[tuple[str, ChemicalReaction, PatternInfo]]

    def __init__(self, *args, **kwargs) -> None:

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
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        context_mol: Mol | None = None,
        **kwargs,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`ReactionRule.metabolites`.

        ``filter_rules(mol, rule, pattern)`` sees this rule and the pattern,
        including ``span``, before SMARTS runs. ``filter_sites(mol, site, info)``
        sees the resolved effect in ``info["options"]`` before ``RunReactants``.
        ``mol`` is the live tracing parent (``context_mol`` when ``mol`` is a
        kekulé copy). ``context_mol`` is the unsubstituted parent when ``mol``
        is a kekulé copy. Aromatic flags are read from it. ``counters``, when
        passed, records one ``rule_expansions`` per call, one
        ``sites_considered`` per match, ``sites_skipped`` when a filter or a
        topological duplicate refuses the site, and one ``mol_edits`` inside
        :func:`react_at`.
        """

        counters = kwargs.get("counters")
        context = mol if context_mol is None else context_mol
        live = cast(TracingMol, context)
        _bump(counters, "rule_expansions")
        seen: set[SiteSignature] = set()
        ranks = context.xf.topol_equiv

        for work in _kekule_forms(mol):
            for rxn_num, (smarts, _rxn, pattern) in enumerate(self.rxns):
                if not filter_rules(live, self, pattern):
                    continue

                reactant = smarts.split(">>", 1)[0]
                for mapped in work.xf.smarts_matches(reactant):
                    site = _site_indexes(mapped, pattern)
                    if not site:
                        continue
                    effect = resolve_effect(context, mapped, pattern)
                    info: SmartsSiteInfo = {
                        "site": site,
                        "rule": self,
                        "options": effect,
                        "rxn_num": rxn_num,
                        "pattern": pattern,
                    }
                    _bump(counters, "sites_considered")
                    if not filter_sites(live, site, info):
                        _bump(counters, "sites_skipped")
                        continue
                    # Same map roles and the same incident bond orders are one
                    # edit. Equivalent carbons share a rank. Another Kekulé
                    # writing, or swapping which atom is map 1, is not.
                    # A two-atom site also carries its pair orbit, which is
                    # not those ranks.
                    signature = _site_signature(
                        context, work, mapped, ranks, site, rxn_num, effect
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
        **kwargs,
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
    """Shortest alternating path. Either end may hold the double bond."""

    paths = alternating_paths(bond_map, start, end, neighbors)
    if not paths:
        return None
    return min(paths, key=len)


def alternating_paths(
    bond_map: Mapping[tuple[int, int], float],
    start: int,
    end: int,
    neighbors: Mapping[int, Sequence[int]],
) -> list[list[int]]:
    """Every alternating path the phase search can reach.

    A node is expanded once per bond order it still needs, so this is not
    every simple walk. It does keep a longer route that leaves a shared
    atom on the other bond. Either phase may start.
    """

    if start == end or start not in neighbors or end not in neighbors:
        return []
    found: list[list[int]] = []
    seen_paths: set[tuple[int, ...]] = set()
    for first in (2.0, 1.0):
        for path in _alternating_from(bond_map, start, end, neighbors, first):
            key = tuple(path)
            if key in seen_paths:
                continue
            seen_paths.add(key)
            found.append(path)
    return found


def _alternating_from(
    bond_map: Mapping[tuple[int, int], float],
    start: int,
    end: int,
    neighbors: Mapping[int, Sequence[int]],
    first: float,
) -> list[list[int]]:
    queue: deque[tuple[int, float, tuple[int, ...]]] = deque([(start, first, (start,))])
    seen = {(start, int(first))}
    found: list[list[int]] = []
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
                found.append(list(nxt))
                continue
            state = (nbr, int(next_want))
            if state in seen:
                continue
            seen.add(state)
            queue.append((nbr, next_want, nxt))
    return found


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


def _pattern_dedup_token(info: SiteInfo) -> str | None:
    """Prefer ``PatternInfo.name``; else the SMARTS string that owns that pattern.

    Pair emissions have no ``pattern`` → ``None``. Names are the long-term key
    (see TODO); the SMARTS lookup is only the fallback when name is missing.
    """

    if "pattern" not in info:
        return None
    pattern = info["pattern"]
    name = pattern.get("name")
    if name:
        return name
    rule = info["rule"]
    if "rxn_num" in info:
        rxns = getattr(rule, "rxns", None) or ()
        rxn_num = info["rxn_num"]
        if 0 <= rxn_num < len(rxns):
            return rxns[rxn_num][0]
    for entry in getattr(rule, "rxns", ()) or ():
        if len(entry) >= 3 and entry[2] is pattern:
            return entry[0]
    for entry in getattr(rule, "endpoints", ()) or ():
        if len(entry) >= 2 and entry[1] is pattern:
            return entry[0]
    for entry in getattr(rule, "smarts", ()) or ():
        if len(entry) >= 2 and entry[1] is pattern:
            return entry[0]
    return None


def _unique_csmi_key(info: SiteInfo, csmi: str) -> tuple[str, str | None, str]:
    """Product dedup: ``(rule name, PatternInfo.name | SMARTS, csmi)``.

    No site topology. Middle field is ``None`` when the emission has no pattern
    (pair sites).
    """

    rule = info["rule"]
    rule_name = getattr(rule, "name", None) or type(rule).__name__
    return (rule_name, _pattern_dedup_token(info), csmi)


def _site_signature(
    context: Mol,
    work: Mol,
    mapped: Mapping[int, int],
    ranks: dict[int, int],
    site: frozenset[int],
    rxn_num: int,
    effect: Effect,
) -> SiteSignature:
    """Dedup key. Last field is ``((ga, gb), pair_group_id)``, or ``None`` for one atom."""

    return (
        tuple((mapno, ranks[idx]) for mapno, idx in sorted(mapped.items())),
        _incident_orders(work, ranks, mapped),
        rxn_num,
        effect.get("adds"),
        effect.get("removes"),
        bool(effect.get("cleaves")),
        bool(effect.get("dearomatizes")),
        context.xf.atom_pair_orbit_key(site),
    )


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
    aromatic = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
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
    move_charge_with_bonds(rw, before, aromatic)
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
    info: PatternInfo,
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
    info: PatternInfo,
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
    info: PatternInfo,
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
    info: PatternInfo,
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
    if not edit_single_to_double(rw, mapped, {}, rings):
        return False
    rw.GetAtomWithIdx(mapped[2]).SetFormalCharge(1)
    return True


def edit_keep(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: PatternInfo,
    rings: Mapping[int, tuple[tuple[int, ...], ...]],
) -> bool:
    """Leave the endpoint bond alone. The path flip is the reaction."""

    return 1 in mapped


def edit_dealkylate(
    rw: RWMol,
    mapped: Mapping[int, int],
    info: PatternInfo,
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
            PatternInfo,
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
) -> Effect:
    return merge_effects(left, right, dearomatizes)


def _site_atoms(mapped: Mapping[int, int], info: PatternInfo) -> int | None:
    key = info.get("site_map", 1)
    if not isinstance(key, int) or key not in mapped:
        return None
    return mapped[key]


def _kekule_cache(mol: Mol) -> KekuleParents:
    """The dict the resonance rules store. Helpers never touch ``_forest``."""

    forest = mol.xf.forest
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
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        context_mol: Mol | None = None,
        **kwargs,
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
        live = cast(TracingMol, context)
        _bump(counters, "rule_expansions")
        cache = _kekule_cache(mol)
        seen: set[SiteSignature] = set()
        ranks = context.xf.topol_equiv

        for rxn_num, (smarts, _rxn, pattern) in enumerate(self.rxns):
            if not filter_rules(live, self, pattern):
                continue
            reactant = smarts.split(">>", 1)[0]
            for mapped in mol.xf.smarts_matches(reactant):
                site = _site_indexes(mapped, pattern)
                if not site:
                    continue
                effect = resolve_effect(context, mapped, pattern)
                info: SmartsSiteInfo = {
                    "site": site,
                    "rule": self,
                    "options": effect,
                    "rxn_num": rxn_num,
                    "pattern": pattern,
                }
                _bump(counters, "sites_considered")
                if not filter_sites(live, site, info):
                    _bump(counters, "sites_skipped")
                    continue
                work = _reactant_parent(mol, mapped, cache)
                if work is None:
                    _bump(counters, "sites_skipped")
                    continue
                signature = _site_signature(
                    context, work, mapped, ranks, site, rxn_num, effect
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
        filter_rules: FilterRules = lambda mol, rule, info: True,
        filter_sites: FilterSites = lambda mol, site, info: True,
        context_mol: Mol | None = None,
        **kwargs,
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
        counters: EditCounters | None = None,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Endpoint loop.

        ``filter_rules`` drops endpoint patterns before they are matched.
        ``filter_sites`` drops a pair before a cached parent is copied.
        """

        _bump(counters, "rule_expansions")

        live = cast(TracingMol, mol)
        active = [
            (smarts, info)
            for smarts, info in self.endpoints
            if filter_rules(live, self, info)
        ]
        if not active:
            return

        if self.systems == "aromatic":
            systems = mol.xf.aromatic_systems
        else:
            systems = mol.xf.conjugated_systems
        if not systems:
            return

        hits: dict[int, list[tuple[dict[int, int], PatternInfo]]] = defaultdict(list)
        for smarts, info in active:
            for mapped in mol.xf.smarts_matches(smarts):
                hits[mapped[1]].append((mapped, info))
        if len(hits) < 2:
            return

        rings: dict[int, tuple[tuple[int, ...], ...]] | None = None
        cache: KekuleParents | None = None
        resonance_parents: tuple[Mol, ...] | None = None
        for system in systems:
            anchors = [atom for atom in hits if atom in system]
            neighbors = system_neighbors(mol, system)
            # Graph distance is not the alternating path. A shorter even
            # walk must not drop a longer odd alternating path.
            for start, end in itertools.combinations(sorted(anchors), 2):
                combos: list[
                    tuple[
                        dict[int, int],
                        PatternInfo,
                        dict[int, int],
                        PatternInfo,
                        PairSiteInfo,
                    ]
                ] = []
                # Any aromatic atom makes the whole system the dearomatizing
                # case. None means the same edit is not a quinone.
                system_aromatic = any(
                    mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in system
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
                    preview: PairSiteInfo = {
                        "site": site,
                        "rule": self,
                        "options": merge_effects(end1, end2, system_aromatic),
                        "ends": (end1, end2),
                        "end_atoms": (site_a, site_b),
                        "end_maps": (map1, map2),
                        "path_ends": frozenset((start, end)),
                    }
                    _bump(counters, "sites_considered")
                    if not filter_sites(live, site, preview):
                        _bump(counters, "sites_skipped")
                        continue
                    combos.append((map1, info1, map2, info2, preview))
                if not combos:
                    continue
                if self.systems == "aromatic":
                    if cache is None:
                        cache = _kekule_cache(mol)
                    scope = aromatic_parent_atoms(mol, start, end)
                    if scope is None:
                        continue
                    parent_mols = parents_for_ends(
                        mol, start, end, cache, atoms=scope
                    ).parents
                else:
                    # Resonance writings of this molecule. The conjugated
                    # system is kekulized whole, including atoms that are
                    # not themselves aromatic. A carbon-only matching misses
                    # some of those writings.
                    if resonance_parents is None:
                        resonance_parents = tuple(
                            overlay_kekule(mol, bond_map).GetMol()
                            for bond_map in resonance_bond_maps(mol)
                        )
                    parent_mols = resonance_parents
                paths: list[tuple[Mol, list[int]]] = []
                for parent in parent_mols:
                    for path in alternating_paths(
                        _current_bond_map(parent), start, end, neighbors
                    ):
                        paths.append((parent, path))
                paths.sort(key=lambda item: len(item[1]))
                if not paths:
                    continue
                for map1, info1, map2, info2, preview in combos:
                    if rings is None and (
                        info1.get("skip_same_rings") or info2.get("skip_same_rings")
                    ):
                        rings = mol.xf.rings
                    ring_table = rings or {}
                    for parent, path in paths:
                        _bump(counters, "mol_edits")
                        rw = rw_copy(parent)
                        edit1 = EDITS.get(info1.get("edit", ""))
                        edit2 = EDITS.get(info2.get("edit", ""))
                        if edit1 is None or edit2 is None:
                            break
                        if not edit1(rw, map1, info1, ring_table):
                            continue
                        if not edit2(rw, map2, info2, ring_table):
                            continue
                        if not swap_bonds_along_path(rw, path):
                            continue
                        products = list(sanitized_fragments(rw, counters).pieces)
                        if not products:
                            continue
                        yield ProductsOfReaction(info=preview, products=products)


class Hydroxylation(SmartsReactionRule):
    """Adds a hydroxyl to carbon.

    Both patterns add OH and remove one H. ``[#6h]`` is h=1, 2, or 3;
    ``[#6h2,#6h3]`` is the subset with two or three hydrogens. The match
    records which of those the atom actually is.
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
            "[#6h2,#6h3:1]>>[*:1]O",
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
            "[#16v4:1]-[#8H1:2]>>[*:1]=[*:2]",
            describe(removes="HH", name="sulfoxide"),
        ),
        (
            "[#6h:1]-[#8H1:2]>>[*:1]=[*:2]",
            describe(removes="HH", partner="O", name="alcohol"),
        ),
        (
            "[#6h:1]-[#7D1H2,#7D2H1:2]>>[*:1]=[*:2]",
            describe(
                *branches(
                    ({"map": 2, "z": 7, "h": 2}, {"map": 2, "z": 7, "h": 1}),
                    removes="HH",
                ),
                name="amine",
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
                ),
                name="alkyl",
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
                name="phenol_end",
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
                name="amine_end",
            ),
        ),
    )


class QuinoneFormation(ResonancePairRule):
    """Dearomatize an aromatic pair into a quinone, imine, or methide.

    A bare aromatic CH can gain its carbonyl oxygen here (``needs`` ``"O"``).
    The same quinone is hydroxylation of that carbon, then dehydrogenation.
    Use that split when a search should hydroxylate only atoms that still
    need oxygen, and dearomatize only the ring.

    :meth:`canonical_plan` reports that elementary split. Metabolize still
    applies this rule in one hop; the plan is parallel information for search.

    The exocyclic bond is unspecified, so an aromatic ring nitrogen matches
    as well as a single-bonded phenol. ``systems`` is ``conjugated``: every
    conjugated system is walked, aromatic or not. ``options["dearomatizes"]``
    is true only when that system contains an aromatic atom. A non-aromatic
    system is the same edit and is not a quinone.
    """

    systems = "conjugated"

    def canonical_plan(self, mol: Mol, info: SiteInfo) -> tuple[CanonicalStep, ...]:
        """Hydroxylations for missing oxygens, then one dehydrogenation."""

        return quinone_canonical_plan(mol, info)

    endpoints: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6R:1][#8H,#7D1H2,#7D2H1,#6D1H3,#6D2H2,#6D3H1:2]",
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
                name="single_to_double",
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
                name="add_carbonyl_o",
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
                name="replace_halogen",
            ),
        ),
        (
            "[#6H0R:1][#7D3:2]",
            describe(
                partner="N",
                dearomatizes=True,
                edit="iminium",
                site_map=1,
                skip_same_rings=True,
                name="iminium",
            ),
        ),
        (
            "[#6R:1][#7,#8:2][#6:3]",
            describe(
                *branches(
                    ({"map": 2, "z": 7}, {"map": 2, "z": 8}),
                    cleaves=True,
                    dearomatizes=True,
                ),
                edit="dealkylate",
                site_map=1,
                skip_same_rings=True,
                name="dealkylate",
            ),
        ),
    )


def _whens(mapno: int, atomic_nums: Sequence[int]) -> tuple[When, ...]:
    """One branch constraint per atomic number, for :func:`branches`."""

    out: list[When] = []
    for z in atomic_nums:
        when: When = {"map": mapno, "z": int(z)}
        out.append(when)
    return tuple(out)


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
                name="methyl_carboxylic",
            ),
        ),
        (
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methyl_carbonyl",
            ),
        ),
        (
            "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methyl_alcohol",
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="OO", cleaves=True),
                site_map=(1, 2),
                name="methylene_carboxylic",
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methylene_carbonyl",
            ),
        ),
        (
            "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methylene_alcohol",
            ),
        ),
        (
            "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methine_carbonyl",
            ),
        ),
        (
            "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methine_alcohol",
            ),
        ),
        (
            "[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="quaternary_alcohol",
            ),
        ),
        (
            # H0 only: [#6] would also match [#6h] and double-emit the same
            # alcohol under unique_csmi (distinct PatternInfo names, same csmi).
            "[#6H0:1][#6:2]>>(O-[*:1].[*:2])",
            describe(
                adds="O",
                cleaves=True,
                partner="C",
                site_map=(1, 2),
                name="cc_quaternary_alcohol",
            ),
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
                name="cc_alcohol",
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
                name="cc_carbonyl",
            ),
        ),
        (
            "[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])",
            describe(
                *branches(_whens(2, (7, 8, 16)), removes="H", cleaves=True),
                site_map=(1, 2),
                name="hemiaminal",
            ),
        ),
    )


def _ndealk(
    smarts: str, leave_count: int | None, name: str, **effect
) -> tuple[str, PatternInfo]:
    """One N-dealkylation pattern. ``leave_count`` is the named leaving atoms."""

    return (
        smarts,
        describe(
            cleaves=True,
            partner="N",
            leave_count=leave_count,
            site_map=(1, 2),
            name=name,
            **effect,
        ),
    )


class NDealkylation(SmartsReactionRule):
    """Cleaves a carbon-nitrogen bond and oxygenates the carbon side.

    These are the nitrogen rows of dealkylation. A methyl carbon is the whole
    leaving piece (``leave_count`` 1). Any larger alkyl still carries atoms
    the pattern does not name, so ``leave_count`` is None. ``breaks_ring`` is
    filled from the cleaved bond, not stored as a second class.

    Forest ``phase1_steps`` is a degenerate singleton naming this rule; the
    default :meth:`canonical_plan` matches that (not UnstableOxygenation).
    UnstableOxygenation in forest is Dealkylation + OxidativeDehalogenation;
    NDealkylation has its own ruleset.
    """

    smarts: tuple[tuple[str, PatternInfo], ...] = (
        _ndealk(
            "[#6H3:1][#7:2]>>([*:2].[*:1](=O)O)",
            1,
            "methyl_carboxylic",
            adds="OO",
        ),
        _ndealk("[#6H3:1][#7:2]>>([*:2].[*:1]=O)", 1, "methyl_carbonyl", adds="O"),
        _ndealk("[#6H3:1][#7:2]>>([*:2].[*:1]-O)", 1, "methyl_alcohol", adds="O"),
        _ndealk(
            "[#6H2:1][#7:2]>>([*:2].[*:1](=O)O)",
            None,
            "methylene_carboxylic",
            adds="OO",
        ),
        _ndealk(
            "[#6H2:1][#7:2]>>([*:2].[*:1]=O)", None, "methylene_carbonyl", adds="O"
        ),
        _ndealk(
            "[#6H2:1][#7:2]>>([*:2].[*:1]-O)", None, "methylene_alcohol", adds="O"
        ),
        _ndealk("[#6H1:1][#7:2]>>([*:2].[*:1]=O)", None, "methine_carbonyl", adds="O"),
        _ndealk("[#6H1:1][#7:2]>>([*:2].[*:1]-O)", None, "methine_alcohol", adds="O"),
        _ndealk(
            "[#6H0:1][#7:2]>>([*:2].[*:1]-O)", None, "quaternary_alcohol", adds="O"
        ),
        _ndealk(
            "[#8H1:3]-[#6:1]-[#7:2]>>([*:3]=[*:1].[*:2])",
            None,
            "hemiaminal",
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
            describe(cleaves=True, partner="N", site_map=(1, 2), name="azo"),
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
                name="dioxole_methylene",
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
                name="nitro_charged",
            ),
        ),
        (
            "[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
            describe(
                cleaves=True,
                partner="N",
                leave_count=1,
                site_map=(1, 2),
                name="nitro_neutral",
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
            describe(adds="O", symbol="S", name="thiophene_s_oxide"),
        ),
    )


_HALIDE = (9, 17, 35, 53, 85)


class Dephosphorylation(SmartsReactionRule):
    """Cleaves an ester O-P bond of a phosphate. The oxygen stays on the organic fragment.

    Map 1 must be the carbon-bound oxygen. Plain P-OH matches are excluded so
    RDKit uniquify cannot keep a water / methyl-phosphite cleavage instead of
    the ester (see DIVERGENCES.md).
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#8;$([#8][#6]):1][#15:2](=[#8:3])([#8:4])[#8:5]>>"
            "[*:1].[*:2](=[*:3])([*:4])[*:5]",
            describe(
                *branches(_whens(2, (15,)), cleaves=True), name="phosphate_ester"
            ),
        ),
    )


class EpoxideOpening(SmartsReactionRule):
    """Opens an epoxide. One pattern only rearranges bonds; the other also adds OH."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])",
            describe(adds="", name="rearrange"),
        ),
        (
            "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)",
            describe(adds="O", name="hydrate"),
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
                name="add_water",
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
                name="cleave",
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
                ),
                name="alcohol",
            ),
        ),
        (
            "[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]",
            describe(removes="OH", cleaves=True, partner="O", name="beta_elimination"),
        ),
        (
            "[#6,#7:1]=[#8:2]>>[*:1].[*:2]",
            describe(
                *branches(
                    ({"map": 1, "z": 6}, {"map": 1, "z": 7}),
                    removes="O",
                    cleaves=True,
                    partner="O",
                ),
                name="carbonyl",
            ),
        ),
    )


class Hydrogenation(ResonancePairRule):
    """Reduces C#C to C=C, C=C to C-C, and a conjugated pair across the path.

    The double-bond pattern is ``=,:``, so an aromatic bond matches once.
    The pair end is any atom (``[*:1]``) and names ``keep``: nothing changes
    at the end except the path flip, which adds H where a double bond
    becomes single. A carbonyl oxygen is an end, so ``CC=O`` becomes
    ``CCO`` here rather than as an oxygen reduction. There is no
    ``partner``: that field is the methide alkyl, and a carbon partner
    would drop this path.
    Heavy-atom formula is unchanged (``adds`` is ``HH``).
    """

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#6:1]#[#6:2]>>[*:1]=[*:2]",
            describe(adds="HH", name="alkyne"),
        ),
        (
            "[#6:1]=,:[#6:2]>>[*:1]-[*:2]",
            describe(adds="HH", name="alkene"),
        ),
    )
    endpoints: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[*:1]",
            describe(adds="H", edit="keep", name="path_end"),
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
            describe(removes="O", cleaves=True, partner="O", name="nitro_charged"),
        ),
        (
            "[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O", name="nitro_anion"),
        ),
        (
            "[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O", name="nitro_neutral"),
        ),
        (
            "[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O", name="nitro_to_amine"),
        ),
        (
            "[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O", name="nitro_both"),
        ),
        (
            "[#7:1]-[#8:2]>>([*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O", name="hydroxylamine"),
        ),
        (
            "[#7D2:1]=[#8:2]>>([*:1].[*:2])",
            describe(removes="O", cleaves=True, partner="O", name="nitroso"),
        ),
        (
            "[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])",
            describe(removes="OO", cleaves=True, partner="O", name="nitro_both_any"),
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
                ),
                name="carbonyl",
            ),
        ),
        (
            "[#8:1]-[#8:2]>>[*:1].[*:2]",
            describe(cleaves=True, partner="O", name="peroxide"),
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
                name="cleave",
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
                name="alkene",
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
            describe(removes="O", cleaves=True, partner="O", name="sulfoxide"),
        ),
        (
            "[#16:1]-[#16:2]>>[*:1].[*:2]",
            describe(cleaves=True, partner="S", name="disulfide"),
        ),
        (
            "[#16:1]-[#6,#8:2]>>[*:1].[*:2]",
            describe(
                *branches(({"map": 2, "z": 6},), cleaves=True),
                *branches(({"map": 2, "z": 8},), cleaves=True, removes="O"),
                name="thioether",
            ),
        ),
    )



class Epoxidation(ResonanceRule):
    """Adds an epoxide across a C=C or C=N bond.

    The reactant bond is ``=,:``, so an aromatic bond matches on the parent.
    The reaction runs on the cached kekulé parent for that bond.

    Forest ``phase1_steps`` is a degenerate singleton naming this rule; the
    default :meth:`canonical_plan` matches that (not StableOxygenation).
    StableOxygenation is the group that *contains* Epoxidation among peers.

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
                ),
                name="epoxide",
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
            describe(adds="O", symbol="S", name="zwitterion"),
        ),
        (
            "[#16;v2,v4:1]>>[*:1][O]",
            describe(adds="O", symbol="S", name="hydroxy"),
        ),
        (
            "[#16;v2,v4:1]>>[*:1]=O",
            describe(adds="O", symbol="S", name="oxo"),
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
                ),
                name="hydroxylamine",
            ),
        ),
        (
            "[#7v3H2:1]>>[*:1]=O",
            describe(adds="O", h=2, symbol="N", name="nitroso"),
        ),
        (
            "[#7v3H0:1]>>[*&H0&+:1][O-]",
            describe(adds="O", h=0, symbol="N", name="n_oxide"),
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
                name="alcohol",
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
                name="carbonyl",
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
                name="carboxylic",
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
                name="rearrange",
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
                name="gem_carboxylic",
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
                name="gem_hydrate",
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
    A filter reads that. This reaction does not cleave. Products are
    terminal (``is_terminal_rule``): conjugation ends further expansion.
    """

    is_terminal_rule: bool = True
    as_star: bool = True
    star_label: str | None = None
    smarts: tuple[tuple[str, PatternInfo], ...] = (
        (
            "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]",
            describe(
                *branches(_whens(1, (7, 8, 16)), adds="CCO", removes="H"),
                pin=(1,),
                name="acetyl",
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
        filter_rules=lambda mol, rule, info: True,
        filter_sites=lambda mol, site, info: True,
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
                name="alcohol",
            ),
        ),
        (
            # Epoxide on a cyclohexadiene opens to a methyl sulfone.
            "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>"
            "[*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1",
            describe(adds="CSO", removes="O", site_map=4, name="epoxide_methyl_sulfone"),
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
                name="alcohol",
            ),
        ),
        (
            "[#8H1,#8-:1][#6:2](=[#8:3])[#6:4]>>"
            "O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2](=[#8:3])[*:4])1",
            describe(adds="CCCCCCOOOOOO", site_map=1, name="carboxylate"),
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
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_ch"),
        ),
        (
            "[#6H2:1]1[#8:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_ch2"),
        ),
        (
            "[#6:1]([!#1:4])1[#8:2][#6:3]1>>" + _gsh("[*:1]([*:4])[*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_c"),
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
                name="halide",
            ),
        ),
        (
            "[#16h1:1]>>" + _gsh("[*:1]"),
            describe(adds=_GSH_ADDS, removes="H", site_map=1, name="thiol"),
        ),
        (
            "[#6H2:1]=[#6:2]>>" + _gsh("[*:1]-[*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="alkene"),
        ),
        (
            "[#6H1:1]=[#6:2][#6:3]=[#8,#7:4]>>" + _gsh("[*:1][*:2]=[*:3][*:4]"),
            describe(
                *branches(_whens(4, (8, 7)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
                name="michael",
            ),
        ),
        (
            "[#6;H1,H2:1]=[#8:2]>>" + _gsh("[*:1]([*:2])"),
            describe(adds=_GSH_ADDS, site_map=1, name="carbonyl"),
        ),
        (
            "[#6H1:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_ch"),
        ),
        (
            "[#6H2:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_ch2"),
        ),
        (
            "[#6:1]([!#1:4])1[#7:2][#6:3]1>>" + _gsh("[*:1]([*:4])[*:3][*:2]"),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_c"),
        ),
        (
            "[#6:1][#8:2]S(=O)(=O)>>" + _gsh("[*:1]"),
            describe(adds=_GSH_ADDS, site_map=1, name="mesylate"),
        ),
        (
            "[#7:1]=[#6:2]=[#8,#16:3]>>" + _gsh("[*:2](=[*:3])[*:1]"),
            describe(
                *branches(_whens(3, (8, 16)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
                name="isocyanate",
            ),
        ),
    )
