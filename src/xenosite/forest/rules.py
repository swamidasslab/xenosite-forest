"""Define specific reaction rules."""

from __future__ import annotations

# Standard Library
import itertools
import logging
import re
import warnings
from collections import defaultdict, deque
from collections.abc import Callable, Generator, Iterable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, NamedTuple, cast

if TYPE_CHECKING:
    # Static re-export so ``from .rules import RuleSet`` types correctly.
    # Runtime still goes through ``__getattr__`` to avoid the cycle.
    from .rulesets import RuleSet as RuleSet

from xenosite.forest.canonical_plan import (
    CanonicalStep,
    identity_canonical_plan,
    quinone_canonical_plan,
)
from xenosite.forest.forest_copy import copy_mutable
from xenosite.forest.graph_isomorphism import (
    SiteSignature,
    canonical_emitted_sites_requested,
    canonicalize_pair_match,
    canonicalize_smarts_match,
    pair_site_signature,
    site_signature,
)
from xenosite.forest.rdkit_api import (
    ChemicalReaction,
    MolFromSmarts,
    MolToSmarts,
    RenumberAtoms,
)
from xenosite.forest.rdkitutil import (
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
    reaction_from_smirks,
    resonance_bond_maps,
    run_reactants,
    rw_copy,
    sanitize_catch,
    sanitized_fragments,
)
from xenosite.forest.records import (
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
    RuleSiteKind,
    SiteInfo,
    SitesOn,
    Smarts,
    SmirksSiteInfo,
    Smirks,
    Span,
    TraceAddition,
    TraceInfo,
    When,
)


class _LazyProductInfo(dict):
    """Product ``info`` view. ``csmi`` is the emission frozenset of fragment CSMIs.

    The hot path does not compute SMILES until ``info["csmi"]`` / yield dedup
    reads each ``mol.xf.csmi`` (cached on ``structure["csmi"]``).
    """

    __slots__ = ("_mols",)

    def __init__(self, data: dict[str, object], mols: Sequence[ForestMol]):
        super().__init__(data)
        self._mols = list(mols)

    def __getitem__(self, key: str) -> object:
        if key == "csmi":
            return frozenset(m.xf.csmi for m in self._mols)
        return super().__getitem__(key)

    def get(self, key, default=None):  # type: ignore[override]
        if key == "csmi":
            return frozenset(m.xf.csmi for m in self._mols)
        return super().get(key, default)

    def __contains__(self, key: object) -> bool:
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
    :meth:`ReactionRule.metabolize` yields ``(list[TracingMol], info)`` —
    one list per emission (siblings stay together for cleavage).
    """

    info: SiteInfo
    products: list[Mol]


def _accept_all_rules(
    mol: TracingMol, rule: ReactionRule, info: PatternInfo
) -> bool:
    return True


def _accept_all_sites(mol: TracingMol, site: Site, info: SiteInfo) -> bool:
    return True


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
    # Emitted Site shape. Class data — not Generic[SiteT] (heterogeneous
    # RuleSets erase the param; pyright cannot enforce container size).
    # ``"atom"`` singleton frozenset; ``"bond"`` undirected frozenset;
    # ``"directed_bond"``: unique-edit uses map-order tuple; public ``site``
    # is frozenset and ``discovered_site`` is the ordered tuple (same Site
    # union). ``"atom_pair"`` ResonancePair ends only (frozenset).
    site_kind: RuleSiteKind = "atom"
    # Internal SMILES guaranteed to yield metabolites (site_kind meta-test).
    # TODO: expand so examples cover all patterns/whens on this rule.
    _example_substrates: tuple[str, ...] = ()

    def _clear_atom_maps(self, mol: Mol) -> Mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return mol

    def __init__(
        self,
        name: str | None = None,
        sites_on: SitesOn | None = None,
        longname: str | None = None,
        *args: Any,
        **kwargs: Any,
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
        self, mol: Mol, **kwargs: Any
    ) -> Generator[tuple[list[TracingMol], ProductInfo], None, None]:
        yield from self.metabolize(mol, **kwargs)

    def __iter__(self) -> Iterator[ReactionRule]:
        return iter([self])

    def metabolize(
        self,
        mol: Mol,
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        unique_csmi: bool = True,
        **kwargs: Any,
    ) -> Generator[tuple[list[TracingMol], ProductInfo], None, None]:
        """Apply this rule and yield ``(products, info)`` pairs.

        ``products`` is always a list: one mol for non-cleavage, all
        sibling fragments for cleavage. This is the method callers use.
        It keeps the invariants below. :meth:`metabolites` supplies the
        chemistry and must not try to.

        The mol that was passed in:

        - Its atoms, bonds, charges, hydrogens, and atom-map numbers are
          unchanged.
        - If it has no ``_forest`` trace, one is installed and left there.
        - If it already has a trace, that trace's depth is not changed.

        Every product mol in each yielded list:

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
        - Shares emission ``info["csmi"]``: the frozenset of fragment
          canonical SMILES for that yield (``product.xf.csmi`` is still
          each mol's own SMILES). Product CSMI has two roles: **check**
          (always, while site unique-edit is on) warns
          ``SiteDeduplicationWarning`` when a later emission under the
          same ``(rule, pattern)`` repeats an earlier emission's frozenset
          of fragment CSMIs *and* matching site ranks (unique-edit miss);
          **yield** (only when ``unique_csmi``, default on) drops duplicate
          ``(rule, pattern, emission frozenset)`` emissions from the
          stream. Check does not imply drop; ``unique_csmi=False`` still
          checks but yields every emission after site unique-edit. Yield
          keys prefer ``info["pattern"]["name"]``, else the SMARTS string;
          pair emissions without a pattern use ``None`` for the middle
          field. Within one cleavage emission, identical sibling CSMIs
          stay in the list; yield dedup is emission-level only (same
          frozenset as check).

        ``filter_rules(mol, rule, pattern_info)`` sees the pattern before a
        match. ``filter_sites(mol, site, info)`` sees the **discovery** site
        and resolved effect before an edit. Either may refuse. A refusal
        edits nothing.

        ``canonical_emitted_sites`` (opt-in, default off): after a site is
        accepted, isomorphic embeddings are remapped to the lex-smallest
        orbit representative for chemistry and the emitted ``site`` key.
        ``discovered_site`` holds the pre-canonical indexes when they differ.
        Filters must not assume ``info["site"]`` equals discovery — use
        ``discovered_site`` when present. See docs/forest/PAIR_ORBITS.md / HEURISTICS.

        For ``site_kind="directed_bond"``, public ``info["site"]`` is always a
        frozenset (API stays unordered-site shaped). Orientation lives on
        ``info["discovered_site"]`` as the ordered map-order tuple (same
        ``Site`` union — no extra field). With ``canonical_emitted_sites``,
        ``site`` is the frozenset of the lex representative and
        ``discovered_site`` remains the directed discovery tuple. Unique-edit
        keys directed ``MapRankKey`` before presentation.
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

        # Yield layer (unique_csmi): (rule, PatternInfo.name | SMARTS, S).
        # Site topology lives only in unique-edit upstream. S is the emission
        # frozenset of fragment CSMIs (same as check).
        seen_yield: set[tuple[str, str | None, frozenset[str]]] = set()
        # Check layer: keepers are (emission CSMI frozenset, site ranks) under
        # (rule, pattern). Same S + equal ranks → SiteDeduplicationWarning (miss).
        # Check always runs (site unique-edit is on); does not imply drop.
        seen_emissions: dict[
            tuple[str, str | None], list[tuple[frozenset[str], tuple[int, ...]]]
        ] = {}

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
            if "discovered_site" in info:
                from xenosite.forest.rdkitutil import (
                    restamp_product_forest_last_layer,
                )

                for p in finished:
                    # of_products cleared product cache — lex reps live on parent.
                    restamp_product_forest_last_layer(p, parent=mol)

            # Cached xf.csmi after of_products — one read per fragment.
            fragment_csmis = [p.xf.csmi for p in finished]
            emission_csmi = frozenset(fragment_csmis)
            site_ranks = _site_ranks_for_csmi_warn(info, mol)
            rule_pat = (
                _rule_dedup_name(info["rule"]),
                _pattern_dedup_token(info),
            )
            prior = seen_emissions.setdefault(rule_pat, [])
            if any(
                emission_csmi == kept_s and site_ranks == kept_r
                for kept_s, kept_r in prior
            ):
                # Unique-edit miss: same emission set + ranks as a keeper.
                _report_csmi_dedup_drop(
                    mol, info, next(iter(emission_csmi))
                )
            else:
                # Quiet unequal-rank iso / leaving-group / cleavage siblings.
                prior.append((emission_csmi, site_ranks))

            for p in finished:
                assert p.xf.tracing.active

            if unique_csmi:
                key = _unique_csmi_key(info, emission_csmi)
                if key in seen_yield:
                    continue
                seen_yield.add(key)

            i = cast(
                ProductInfo,
                _LazyProductInfo(dict(info), finished),
            )
            yield finished, i

    def _top_site(self, site: Site, mol: Mol) -> Site:
        te = mol.xf.topol_equiv
        if isinstance(site, int):
            return te[site]
        if isinstance(site, tuple):
            # Preserve map order for directed_bond emission.
            return tuple(int(te[s]) for s in site)
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
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        **kwargs: Any,
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
    """Normalize a site for undirected trace bookkeeping (sorted ends).

    Do **not** use for ``discovered_site`` on ``directed_bond`` — that field
    keeps map order. Callers that need order preserved should copy the tuple.
    """

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
    if "discovered_site" in info:
        # Preserve tuple order (directed_bond orientation).
        kept["discovered_site"] = info["discovered_site"]
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
    trace: InitializedAtomTrace = copy_mutable(parent._forest["atom_trace"])
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
    if "discovered_site" in info:
        # Preserve map order for directed_bond (do not sort via _site_tuple).
        disc = info["discovered_site"]
        addition["discovered_site"] = (
            disc if isinstance(disc, tuple) else _site_tuple(disc)
        )
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
    **effect: Any,
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
    swap_group: str | None = None,
    **single: EffectField,
) -> PatternInfo:
    """Build a :class:`PatternInfo`.

    One outcome: ``describe(adds="O", removes="H")``.
    Several: ``describe(*branches(...), edit="dealkylate")``.
    ``name`` distinguishes this pattern from the others on the same rule.
    ``swap_group`` is optional; omit when it equals ``name`` (the default).
    Set it only when interchangeable ends group differently from ``name``.
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
    if swap_group is not None:
        info["swap_group"] = swap_group
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

    ``dearomatizes`` on PatternInfo is capability (same as ResonancePair ends).
    The resolved ``options["dearomatizes"]`` is true only when a site-map atom
    is aromatic — mirroring ``merge_effects(..., system_aromatic)``.
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
    if effect.get("dearomatizes"):
        effect["dearomatizes"] = _site_map_aromatic(mol, mapped, site_map)
    return effect


def _site_map_aromatic(
    mol: Mol, mapped: Mapping[int, int], site_map: int | tuple[int, ...]
) -> bool:
    """True when any PatternInfo ``site_map`` atom is aromatic on ``mol``."""

    maps: tuple[int, ...]
    if isinstance(site_map, int):
        maps = (site_map,)
    else:
        maps = tuple(site_map)
    for mapno in maps:
        idx = mapped.get(mapno)
        if isinstance(idx, int) and mol.GetAtomWithIdx(idx).GetIsAromatic():
            return True
    return False


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


def _isotope_smirks(smirks: Smirks, pin: tuple[int, ...]) -> Smirks:
    """Restrict a SMIRKS reaction to the pinned maps.

    The matched atoms are stamped with isotope ``8000 + map number``. The
    isotope is applied to the whole atom query, so an OR list does not keep
    the label on only its first alternative. Those atoms are written first:
    ``SubstructMatch`` starts at query atom 0 and does not reorder.
    """

    import re

    wanted = set(pin)

    def repl(match: re.Match[str]) -> str:
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

    rewritten = re.sub(r"\[([^\[\]]*):(\d+)\]", repl, smirks)
    return _isotope_atoms_first(Smirks(rewritten), wanted)


def _isotope_atoms_first(smirks: Smirks, pin: set[int]) -> Smirks:
    """Put each isotope-bearing reactant atom before the others."""

    if ">>" not in smirks or not pin:
        return smirks
    reactant, product = smirks.split(">>", 1)
    query = MolFromSmarts(reactant)
    if query is None:
        return smirks
    pinned = [
        idx
        for _mapno, idx in sorted(
            (atom.GetAtomMapNum(), atom.GetIdx())
            for atom in query.GetAtoms()
            if atom.GetAtomMapNum() in pin
        )
    ]
    if not pinned:
        return smirks
    rest = [idx for idx in range(query.GetNumAtoms()) if idx not in set(pinned)]
    order = pinned + rest
    if order == list(range(query.GetNumAtoms())):
        return smirks
    rewritten = MolToSmarts(RenumberAtoms(query, order))
    if MolFromSmarts(rewritten) is None:
        return smirks
    return Smirks(rewritten + ">>" + product)


def _reactant_smarts(smirks: Smirks) -> Smarts:
    """Match-only reactant side of a :class:`Smirks`."""

    return Smarts(smirks.split(">>", 1)[0])


def _site_indexes(
    mapped: Mapping[int, int],
    pattern: PatternInfo,
    *,
    site_kind: RuleSiteKind = "atom",
) -> Site:
    """Atom indexes the pattern calls the site. Defaults to map 1.

    ``directed_bond`` keeps ``site_map`` order as a tuple (map 1 first when
    that is the chemically distinct end) for unique-edit. Public yield
    coerces via :func:`_present_site_info` (frozenset ``site`` + ordered
    ``discovered_site``). Other kinds are already a frozenset (undirected
    ``bond`` / singleton ``atom``).
    """

    key = pattern.get("site_map", 1)
    if isinstance(key, (list, tuple)):
        idxs = tuple(mapped[k] for k in key if k in mapped)
    elif key in mapped:
        idxs = (mapped[key],)
    else:
        idxs = tuple(mapped.values())
    if not idxs:
        return () if site_kind == "directed_bond" else frozenset()
    if site_kind == "directed_bond":
        return idxs
    return frozenset(idxs)


def _present_site_info(
    info: SmirksSiteInfo,
    emit_site: Site,
    discovery: Site,
    *,
    site_kind: RuleSiteKind,
) -> SmirksSiteInfo:
    """Yield-only presentation of ``site`` / ``discovered_site``.

    Unique-edit ``seen`` must already have keyed on ``discovery`` (directed
    tuple + ``MapRankKey`` for ``directed_bond``).

    For ``directed_bond``:
    - ``site`` = ``frozenset(emit_site)`` (lex-canonical when
      ``canonical_emitted_sites`` remapped ``emit_site``)
    - ``discovered_site`` = ordered discovery tuple (orientation; always set)

    For other kinds, leave canonical ``discovered_site`` semantics unchanged
    (only when remap already set it on ``info``); ensure ``site`` is
    ``emit_site``.
    """

    if site_kind == "directed_bond" and isinstance(discovery, tuple):
        public_emit = (
            frozenset(emit_site) if isinstance(emit_site, tuple) else emit_site
        )
        return {
            **info,
            "site": public_emit,
            "discovered_site": discovery,
        }
    if info.get("site") is emit_site:
        return info
    return {**info, "site": emit_site}


def react_at(
    rule: SmirksReactionRule,
    smirks: Smirks,
    mol: Mol,
    mapped: Mapping[int, int],
    counters: EditCounters | None = None,
    pin: tuple[int, ...] | None = None,
) -> list[Mol]:
    """Run ``smirks`` on one match. One call is one ``mol_edits``.

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
        product_sets = run_reactants(_isotope_smirks(smirks, chosen), stamped)
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


class SmirksReactionRule(ReactionRule):
    """Performs reactions specified by SMIRKS (reaction SMARTS).

    Each entry is ``(smirks, options)`` with :class:`PatternInfo` options.
    Values are :class:`~xenosite.forest.records.Smirks` (``reactant>>product``).
    Matches are filtered before ``RunReactants``. One SMIRKS is applied once
    per topological site; a later SMIRKS with the same formula effect still runs.
    """

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = ()

    rxns: list[tuple[Smirks, ChemicalReaction, PatternInfo]]

    def __init__(self, *args: Any, **kwargs: Any) -> None:

        super().__init__(*args, **kwargs)

        self.rxns = [
            (smirks, self._smirks2rxns(smirks, **kwargs), opt)
            for smirks, opt in self.smirks
        ]
        patterns = [opt for _smirks, _rxn, opt in self.rxns]
        endpoints = getattr(self, "endpoints", ()) or ()
        patterns.extend(opt for _smarts, opt in endpoints)
        _assign_pattern_names(patterns)

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        context_mol: Mol | None = None,
        **kwargs: Any,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`ReactionRule.metabolites`.

        ``filter_rules(mol, rule, pattern)`` sees this rule and the pattern,
        including ``span``, before SMARTS runs. ``filter_sites(mol, site, info)``
        sees the resolved effect in ``info["options"]`` before ``RunReactants``.
        ``mol`` is the live tracing parent. ``context_mol``, when set, is an
        alternate unsubstituted parent for aromatic flags / ranks. ``counters``,
        when passed, records one ``rule_expansions`` per call, one
        ``sites_considered`` per match, ``sites_skipped`` when a filter or a
        topological duplicate refuses the site, and one ``mol_edits`` inside
        :func:`react_at`.
        """

        counters = kwargs.get("counters")
        want_canonical = canonical_emitted_sites_requested(kwargs)
        context = mol if context_mol is None else context_mol
        live = cast(TracingMol, context)
        _bump(counters, "rule_expansions")
        seen: set[SiteSignature] = set()
        ranks = context.xf.topol_equiv

        for rxn_num, (smirks, _rxn, pattern) in enumerate(self.rxns):
            if not filter_rules(live, self, pattern):
                continue

            reactant = _reactant_smarts(smirks)
            for mapped in mol.xf.smarts_matches(reactant):
                site = _site_indexes(mapped, pattern, site_kind=self.site_kind)
                if not site:
                    continue
                effect = resolve_effect(context, mapped, pattern)
                info: SmirksSiteInfo = {
                    "site": site,
                    "rule": self,
                    "options": effect,
                    "rxn_num": rxn_num,
                    "pattern": pattern,
                }
                _bump(counters, "sites_considered")
                # Filters see the discovery site (not the lex representative).
                if not filter_sites(live, site, info):
                    _bump(counters, "sites_skipped")
                    continue
                signature = site_signature(
                    context,
                    mol,
                    mapped,
                    ranks,
                    site,
                    rxn_num,
                    effect,
                    site_kind=self.site_kind,
                )
                if signature in seen:
                    _bump(counters, "sites_skipped")
                    continue
                emit_mapped: dict[int, int] = dict(mapped)
                emit_site = site
                if want_canonical:
                    remapped = canonicalize_smarts_match(
                        context,
                        mapped,
                        site,
                        parent=context,
                    )
                    if remapped is None:
                        _bump(counters, "sites_skipped")
                        continue
                    emit_mapped, emit_site = remapped
                    if emit_site != site:
                        info = {
                            **info,
                            "site": emit_site,
                            "discovered_site": site,
                        }
                    else:
                        info = {**info, "site": emit_site}
                products = react_at(
                    self,
                    smirks,
                    mol,
                    emit_mapped,
                    counters,
                    pattern.get("pin"),
                )
                if not products:
                    continue
                # Unique-edit ``seen`` already keyed on directed signature.
                # Presentation: frozenset site; directed_bond orientation on
                # discovered_site (ordered tuple).
                seen.add(signature)
                info = _present_site_info(
                    info,
                    emit_site,
                    site,
                    site_kind=self.site_kind,
                )
                yield ProductsOfReaction(info=info, products=products)

    def _smirks2rxns(
        self,
        smirks: Smirks,
        use_implicit_properties: bool = False,
        **kwargs: Any,
    ) -> ChemicalReaction:
        """Converts SMIRKS reactions to RDKit reactions."""
        return reaction_from_smirks(smirks)

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
    """Shortest double-first alternating path. Either end may hold the opening double."""

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
    """Every double-first alternating path between ``start`` and ``end``.

    The first bond must be double from one endpoint (then alternate). Either
    endpoint may hold that opening double — try ``start→end`` and
    ``end→start``, and orient results ``start→end``. Paths start on a double
    only (not a single). Keep longer odd routes that leave a shared atom on
    the other bond (node expanded once per needed bond order).
    """

    if start == end or start not in neighbors or end not in neighbors:
        return []
    found: list[list[int]] = []
    seen_paths: set[tuple[int, ...]] = set()
    for origin, target, reverse in (
        (start, end, False),
        (end, start, True),
    ):
        for path in _alternating_from(bond_map, origin, target, neighbors, 2.0):
            oriented = list(reversed(path)) if reverse else path
            key = tuple(oriented)
            if key in seen_paths:
                continue
            seen_paths.add(key)
            found.append(oriented)
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
    for entry in getattr(rule, "smirks", ()) or ():
        if len(entry) >= 2 and entry[1] is pattern:
            return entry[0]
    return None


def _unique_csmi_key(
    info: SiteInfo, emission_csmi: frozenset[str]
) -> tuple[str, str | None, frozenset[str]]:
    """Emission dedup: ``(rule name, PatternInfo.name | SMARTS, CSMI frozenset)``.

    No site topology. Middle field is ``None`` when the emission has no pattern
    (pair sites). The frozenset matches the check-layer emission set.
    """

    rule = info["rule"]
    rule_name = getattr(rule, "name", None) or type(rule).__name__
    return (rule_name, _pattern_dedup_token(info), emission_csmi)


def _site_ranks_for_csmi_warn(info: SiteInfo, mol: Mol) -> tuple[int, ...]:
    """Topological ranks of the emission site (for SiteDeduplicationWarning gating).

    Directed discovery tuples keep map order; frozenset sites are sorted.
    Matching ranks on a full-duplicate CSMI drop means unique-edit should have
    merged the embeddings; unequal ranks are product isomorphism only.
    """

    ranks = mol.xf.topol_equiv
    site = info.get("discovered_site", info.get("site"))
    if isinstance(site, int):
        return (int(ranks[site]),)
    if isinstance(site, tuple):
        return tuple(int(ranks[i]) for i in site)
    if isinstance(site, frozenset):
        return tuple(sorted(int(ranks[i]) for i in site if isinstance(i, int)))
    return ()


class SiteDeduplicationWarning(UserWarning):
    """Unique-edit miss: later emission repeats a keeper's CSMI set + site ranks."""


def _rule_dedup_name(rule: object) -> str:
    return getattr(rule, "name", None) or type(rule).__name__


def _csmi_dedup_warning_message(rule_name: str) -> str:
    """Generic text keyed by rule so stdlib once-per-message emits once per rule."""

    return (
        f"CSMI dedup triggered for rule {rule_name}: canonization / unique-edit "
        "filtering is off or incomplete (unique-edit/orbit missed a duplicate "
        "product). Raise log level to INFO on this logger for per-drop detail."
    )


_logger = logging.getLogger(__name__)


def _report_csmi_dedup_drop(substrate: Mol, info: SiteInfo, product_csmi: str) -> None:
    """Warn once-per-rule (generic) and log INFO on a unique-edit miss check hit."""

    rule_name = _rule_dedup_name(info["rule"])
    warnings.warn(
        _csmi_dedup_warning_message(rule_name),
        SiteDeduplicationWarning,
        stacklevel=2,
    )
    _logger.info(
        "CSMI dedup drop: substrate=%s rule=%s site=%s pattern=%s product=%s",
        substrate.xf.csmi,
        rule_name,
        info.get("site"),
        _pattern_dedup_token(info),
        product_csmi,
    )


def _report_redundant_rules_drop(
    substrate: Mol,
    *,
    kept_rule: str,
    dropped_info: SiteInfo,
    product_csmi: str,
) -> None:
    """INFO when a later child rule repeats an earlier child's product CSMI.

    Overlapping coverage is expected on some substrates; not a unique-edit /
    ``SiteDeduplicationWarning`` miss. No ``warnings.warn`` on this path.
    """

    dropped_rule = _rule_dedup_name(dropped_info["rule"])
    _logger.info(
        "Redundant rules drop: substrate=%s kept_rule=%s dropped_rule=%s "
        "site=%s pattern=%s product=%s",
        substrate.xf.csmi,
        kept_rule,
        dropped_rule,
        dropped_info.get("site"),
        _pattern_dedup_token(dropped_info),
        product_csmi,
    )


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
    if "cache" not in forest:
        raise KeyError("cache")
    structure = forest["cache"]
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


def _smirks_mapped_bond_order(
    smirks: Smirks, left_map: int = 1, right_map: int = 2
) -> float | None:
    """Bond order the reactant SMARTS asks for between two atom maps.

    Reads the reactant template only. ``=`` / ``=,:`` → 2.0; ``-`` / ``-,:`` /
    unspecified (single-or-aromatic) → 1.0. ``None`` when the maps or bond are
    missing.
    """

    reactant = _reactant_smarts(smirks)
    query = MolFromSmarts(reactant)
    if query is None:
        return None
    idxs = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in query.GetAtoms()
        if atom.GetAtomMapNum()
    }
    if left_map not in idxs or right_map not in idxs:
        return None
    bond = query.GetBondBetweenAtoms(idxs[left_map], idxs[right_map])
    if bond is None:
        return None
    return float(bond.GetBondTypeAsDouble())


def _reactant_parent(
    mol: Mol,
    mapped: dict[int, int],
    cache: KekuleParents,
    smirks: Smirks | None = None,
) -> Mol | None:
    """Kekulé parent whose bond orders match an aromatic hit. Else ``mol``.

    On an aromatic bond, the order is the one the reactant SMARTS implies
    between maps 1 and 2 (single or double), not a hard-coded double. That
    lets ring-open / S-oxide patterns pick the single-bond Kekulé parent while
    epoxidation still picks the double.
    """

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
        if smirks is not None:
            implied = _smirks_mapped_bond_order(smirks)
            if implied is not None:
                order = implied
    elif begin.GetIsAromatic() or end.GetIsAromatic():
        order = bond.GetBondTypeAsDouble()
    else:
        return mol
    ensure_kekule_parents(mol, left, right, cache)
    parent = parent_for_bond(cache, left, right, order)
    # No assignment with that order (charged rings, awkward systems): keep the
    # aromatic parent so the site is not dropped.
    return mol if parent is None else parent


class ResonanceRule(SmirksReactionRule):
    """Match once on the aromatic parent, then react on a cached kekulé parent.

    The reactant SMARTS should match aromatic bonds (``=,:``, ``-,:``, or
    unspecified). One parent is cached per assignment of the conjugated system
    that contains the match; ``_reactant_parent`` selects the assignment where
    maps 1–2 have the bond order the SMARTS implies. Other systems stay
    aromatic. The dict is stored on the molecule's forest. The helpers that
    fill it do not read ``_forest``.
    """

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        context_mol: Mol | None = None,
        **kwargs: Any,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Same contract as :meth:`SmirksReactionRule.metabolites`.

        SMARTS runs once on this mol (patterns should match aromatic bonds).
        A hit picks the cached Kekulé parent where maps 1–2 have the bond order
        the reactant SMARTS implies; if none exists, the aromatic mol is kept.
        Unique-edit signatures use aromatic ``incident_orders``. ``context_mol``
        keeps the original aromatic flags for :func:`resolve_effect`.
        """

        if not self.rxns:
            return
        counters = kwargs.get("counters")
        want_canonical = canonical_emitted_sites_requested(kwargs)
        context = mol if context_mol is None else context_mol
        live = cast(TracingMol, context)
        _bump(counters, "rule_expansions")
        cache = _kekule_cache(mol)
        seen: set[SiteSignature] = set()
        ranks = context.xf.topol_equiv

        for rxn_num, (smirks, _rxn, pattern) in enumerate(self.rxns):
            if not filter_rules(live, self, pattern):
                continue
            reactant = _reactant_smarts(smirks)
            for mapped in mol.xf.smarts_matches(reactant):
                site = _site_indexes(mapped, pattern, site_kind=self.site_kind)
                if not site:
                    continue
                effect = resolve_effect(context, mapped, pattern)
                info: SmirksSiteInfo = {
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
                work = _reactant_parent(mol, mapped, cache, smirks)
                if work is None:
                    _bump(counters, "sites_skipped")
                    continue
                # Unique-edit ranks / incident_orders stay on the aromatic
                # parent so Kekulé bond-order flips do not split equivalent sites.
                signature = site_signature(
                    context,
                    mol,
                    mapped,
                    ranks,
                    site,
                    rxn_num,
                    effect,
                    site_kind=self.site_kind,
                )
                if signature in seen:
                    _bump(counters, "sites_skipped")
                    continue
                emit_mapped: dict[int, int] = dict(mapped)
                emit_site = site
                if want_canonical:
                    remapped = canonicalize_smarts_match(
                        context,
                        mapped,
                        site,
                        parent=context,
                    )
                    if remapped is None:
                        _bump(counters, "sites_skipped")
                        continue
                    emit_mapped, emit_site = remapped
                    if emit_site != site:
                        info = {
                            **info,
                            "site": emit_site,
                            "discovered_site": site,
                        }
                    else:
                        info = {**info, "site": emit_site}
                    work = _reactant_parent(mol, emit_mapped, cache, smirks)
                    if work is None:
                        _bump(counters, "sites_skipped")
                        continue
                products = react_at(
                    self,
                    smirks,
                    work,
                    emit_mapped,
                    counters,
                    pattern.get("pin"),
                )
                if not products:
                    continue
                # Unique-edit ``seen`` already keyed on directed signature.
                # Presentation: frozenset site; directed_bond orientation on
                # discovered_site (ordered tuple).
                seen.add(signature)
                info = _present_site_info(
                    info,
                    emit_site,
                    site,
                    site_kind=self.site_kind,
                )
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

    endpoints: tuple[tuple[Smarts, PatternInfo], ...] = ()
    systems: str = "conjugated"

    def metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        context_mol: Mol | None = None,
        **kwargs: Any,
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
            mol,
            filter_rules,
            filter_sites,
            counters=kwargs.get("counters"),
            canonical_emitted_sites=canonical_emitted_sites_requested(kwargs),
        )

    def pair_metabolites(
        self,
        mol: Mol,
        filter_rules: FilterRules,
        filter_sites: FilterSites,
        counters: EditCounters | None = None,
        *,
        canonical_emitted_sites: bool = False,
    ) -> Generator[ProductsOfReaction, None, None]:
        """Endpoint loop.

        ``filter_rules`` drops endpoint patterns before they are matched.
        ``filter_sites`` drops a pair before a cached parent is copied
        (discovery site). With ``canonical_emitted_sites``, chemistry and the
        emitted ``site`` key use the lex representative; ``discovered_site``
        records the pre-canonical indexes when they differ.
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

        ranks = mol.xf.topol_equiv
        seen: set[tuple] = set()
        rings: dict[int, tuple[tuple[int, ...], ...]] | None = None
        cache: KekuleParents | None = None
        resonance_parents: tuple[Mol, ...] | None = None
        path_cache: dict[tuple[int, int], list[tuple[Mol, list[int]]]] = {}
        for system in systems:
            anchors = [atom for atom in hits if atom in system]
            neighbors = system_neighbors(mol, system)
            # Enumerate every anchor pair (not BFS-odd graph distance alone):
            # a shorter even walk must not hide a longer odd alternating path
            # (rings). Apply only odd bond-count paths below — even-length
            # walks are not valid ResonancePair flips.
            for start, end in itertools.combinations(sorted(anchors), 2):
                combos: list[
                    tuple[
                        dict[int, int],
                        PatternInfo,
                        dict[int, int],
                        PatternInfo,
                        PairSiteInfo,
                        tuple,
                        int,
                        int,
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
                    # At most one methide end is data: two methide ends do not
                    # resolve (docs/forest/DROPPED.md / data-not-branches). Not a search
                    # filter — the pair is never built.
                    if end1.get("methide") and end2.get("methide"):
                        continue
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
                    # Filters see discovery indexes.
                    if not filter_sites(live, site, preview):
                        _bump(counters, "sites_skipped")
                        continue
                    signature = pair_site_signature(
                        mol,
                        ranks,
                        map1,
                        map2,
                        site_a,
                        site_b,
                        info1,
                        info2,
                        preview,
                    )
                    if signature in seen:
                        _bump(counters, "sites_skipped")
                        continue
                    emit_map1: dict[int, int] = dict(map1)
                    emit_map2: dict[int, int] = dict(map2)
                    emit_a, emit_b = site_a, site_b
                    emit_start, emit_end = start, end
                    if canonical_emitted_sites:
                        remapped = canonicalize_pair_match(
                            mol,
                            map1,
                            map2,
                            site_a,
                            site_b,
                            info1,
                            info2,
                            effect1=end1,
                            effect2=end2,
                            parent=mol,
                        )
                        if remapped is None:
                            _bump(counters, "sites_skipped")
                            continue
                        emit_map1, emit_map2, emit_a, emit_b = remapped
                        emit_start = emit_map1[1]
                        emit_end = emit_map2[1]
                        emit_site = frozenset((emit_a, emit_b))
                        preview = {
                            **preview,
                            "site": emit_site,
                            "end_atoms": (emit_a, emit_b),
                            "end_maps": (emit_map1, emit_map2),
                            "path_ends": frozenset((emit_start, emit_end)),
                        }
                        if emit_site != site:
                            preview = {**preview, "discovered_site": site}
                    # Reserve the unique-edit slot once per swappable/ordered site.
                    seen.add(signature)
                    combos.append(
                        (
                            emit_map1,
                            info1,
                            emit_map2,
                            info2,
                            preview,
                            signature,
                            emit_start,
                            emit_end,
                        )
                    )
                if not combos:
                    continue
                for (
                    map1,
                    info1,
                    map2,
                    info2,
                    preview,
                    signature,
                    path_start,
                    path_end,
                ) in combos:
                    path_key = (path_start, path_end)
                    if path_key not in path_cache:
                        if self.systems == "aromatic":
                            if cache is None:
                                cache = _kekule_cache(mol)
                            scope = aromatic_parent_atoms(
                                mol, path_start, path_end
                            )
                            if scope is None:
                                path_cache[path_key] = []
                            else:
                                parent_mols = parents_for_ends(
                                    mol,
                                    path_start,
                                    path_end,
                                    cache,
                                    atoms=scope,
                                ).parents
                                paths: list[tuple[Mol, list[int]]] = []
                                for parent in parent_mols:
                                    for path in alternating_paths(
                                        _current_bond_map(parent),
                                        path_start,
                                        path_end,
                                        neighbors,
                                    ):
                                        paths.append((parent, path))
                                paths.sort(key=lambda item: len(item[1]))
                                path_cache[path_key] = paths
                        else:
                            if resonance_parents is None:
                                resonance_parents = tuple(
                                    overlay_kekule(mol, bond_map).GetMol()
                                    for bond_map in resonance_bond_maps(mol)
                                )
                            paths = []
                            for parent in resonance_parents:
                                for path in alternating_paths(
                                    _current_bond_map(parent),
                                    path_start,
                                    path_end,
                                    neighbors,
                                ):
                                    paths.append((parent, path))
                            paths.sort(key=lambda item: len(item[1]))
                            path_cache[path_key] = paths
                    paths = [
                        item
                        for item in path_cache[path_key]
                        if (len(item[1]) - 1) % 2 == 1
                    ]
                    if not paths:
                        continue
                    if rings is None and (
                        info1.get("skip_same_rings") or info2.get("skip_same_rings")
                    ):
                        rings = mol.xf.rings
                    ring_table = rings or {}
                    emitted_csmi: set[str] = set()
                    _bump(counters, "mol_edits")
                    for parent, path in paths:
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
                        fresh = []
                        for product in products:
                            csmi = product.xf.csmi
                            if csmi in emitted_csmi:
                                continue
                            emitted_csmi.add(csmi)
                            fresh.append(product)
                        if fresh:
                            yield ProductsOfReaction(info=preview, products=fresh)


class Hydroxylation(SmirksReactionRule):
    """Adds a hydroxyl to carbon.

    Patterns partition by H count so the same product is not emitted twice
    under ``unique_csmi`` (distinct ``PatternInfo.name``, same csmi). ``h`` is
    exactly one hydrogen; ``h2`` is the widened ``[#6h2,#6h3]`` OR for two or
    three. ``when`` records which branch the atom actually is.
    """

    phase1_sites_on = "atom_hydrogen"
    sites_on = "atom_hydrogen"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO', 'c1ccccc1')

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6h1:1]>>[*:1]O"),
            describe(
                *branches(({"map": 1, "z": 6, "h": 1},), adds="O", removes="H"),
                name="h",
            ),
        ),
        (
            Smirks("[#6h2,#6h3:1]>>[*:1]O"),
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

    Sites are unordered atom pairs (``site_kind="atom_pair"``): one-bond SMARTS emit
    both bond endpoints; path emissions are the two end atoms. Unique-edit
    uses atom–atom pair orbits. Same PatternInfo role → unordered; different
    roles (phenol vs amine) → ordered so ``(a,b)`` ≠ ``(b,a)``.
    """

    phase1_sites_on = "atom_hydrogen"
    sites_on = "atom_pairs"
    site_kind: RuleSiteKind = "atom_pair"
    _example_substrates: tuple[str, ...] = ('CCO', 'Oc1ccc(O)cc1')

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#16v4:1]-[#8H1:2]>>[*:1]=[*:2]"),
            describe(removes="HH", name="sulfoxide", site_map=(1, 2)),
        ),
        (
            Smirks("[#6h:1]-[#8H1:2]>>[*:1]=[*:2]"),
            describe(removes="HH", partner="O", name="alcohol", site_map=(1, 2)),
        ),
        (
            Smirks("[#6h:1]-[#7D1H2,#7D2H1:2]>>[*:1]=[*:2]"),
            describe(
                *branches(
                    ({"map": 2, "z": 7, "h": 2}, {"map": 2, "z": 7, "h": 1}),
                    removes="HH",
                ),
                name="amine",
                site_map=(1, 2),
            ),
        ),
        (
            Smirks("[#6h:1]-[#6D1H3,#6D2H2,#6D3H1:2]>>[*:1]=[*:2]"),
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
                site_map=(1, 2),
            ),
        ),
    )

    endpoints: tuple[tuple[Smarts, PatternInfo], ...] = (
        (
            Smarts("[#6:1]-[#8H:2]"),
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
            Smarts("[#6:1]-[#7D1H2,#7D2H1:2]"),
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
        (
            Smarts("[#6:1]-[#6D1H3,#6D2H2,#6D3H1:2]"),
            describe(
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
                name="methide_end",
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
    sites_on = "atom_pairs"
    site_kind: RuleSiteKind = "atom_pair"
    _example_substrates: tuple[str, ...] = ('c1ccccc1', 'Oc1ccccc1')

    def canonical_plan(self, mol: Mol, info: SiteInfo) -> tuple[CanonicalStep, ...]:
        """Hydroxylations for missing oxygens, then one dehydrogenation."""

        return quinone_canonical_plan(mol, info)

    endpoints: tuple[tuple[Smarts, PatternInfo], ...] = (
        (
            Smarts("[#6R:1][#8H,#7D1H2,#7D2H1,#6D1H3,#6D2H2,#6D3H1:2]"),
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
            Smarts("[#6D2H1;R:1]"),
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
            Smarts("[#6H0R:1]-[F,Cl,Br,I:2]"),
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
            Smarts("[#6H0R:1][#7D3:2]"),
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
            Smarts("[#6R:1][#7,#8:2][#6:3]"),
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


class Dealkylation(ResonanceRule):
    """Cleaves a C-N, C-O, C-S, or C-C bond and oxygenates the carbon side.

    Internally the site is both atoms of the broken bond as an ordered tuple
    (map 1 = oxygenated carbon, map 2 = heteroatom partner) for unique-edit.
    Public ``info["site"]`` is a frozenset; orientation is on
    ``info["discovered_site"]`` (ordered tuple). Aromatic hits react on the
    Kekulé parent where that bond is single (SMARTS-implied order).
    ``site_kind="directed_bond"``: unique-edit keeps directed MapRankKey
    because map 1 is chemically distinct.
    """
    sites_on = "bonds"
    site_kind: RuleSiteKind = "directed_bond"
    _example_substrates: tuple[str, ...] = ('CCO', 'COc1ccccc1')


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="OO", cleaves=True),
                site_map=(1, 2),
                name="methyl_carboxylic",
            ),
        ),
        (
            Smirks("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methyl_carbonyl",
            ),
        ),
        (
            Smirks("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methyl_alcohol",
            ),
        ),
        (
            Smirks("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="OO", cleaves=True),
                site_map=(1, 2),
                name="methylene_carboxylic",
            ),
        ),
        (
            Smirks("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methylene_carbonyl",
            ),
        ),
        (
            Smirks("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methylene_alcohol",
            ),
        ),
        (
            Smirks("[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methine_carbonyl",
            ),
        ),
        (
            Smirks("[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="methine_alcohol",
            ),
        ),
        (
            Smirks("[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)"),
            describe(
                *branches(_whens(2, (7, 8, 16)), adds="O", cleaves=True),
                site_map=(1, 2),
                name="quaternary_alcohol",
            ),
        ),
        (
            # H0 only: [#6] would also match [#6h] and double-emit the same
            # alcohol under unique_csmi (distinct PatternInfo names, same csmi).
            Smirks("[#6H0:1][#6:2]>>(O-[*:1].[*:2])"),
            describe(
                adds="O",
                cleaves=True,
                partner="C",
                site_map=(1, 2),
                name="cc_quaternary_alcohol",
            ),
        ),
        (
            Smirks("[#6h:1][#6:2]>>(O-[*:1].[*:2])"),
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
            Smirks("[#6h:1][#6:2]>>(O=[*:1].[*:2])"),
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
            Smirks("[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])"),
            describe(
                *branches(_whens(2, (7, 8, 16)), removes="H", cleaves=True),
                site_map=(1, 2),
                name="hemiaminal",
            ),
        ),
    )


def _ndealk(
    smirks: Smirks, leave_count: int | None, name: str, **effect: Any
) -> tuple[Smirks, PatternInfo]:
    """One N-dealkylation pattern. ``leave_count`` is the named leaving atoms."""

    return (
        smirks,
        describe(
            cleaves=True,
            partner="N",
            leave_count=leave_count,
            site_map=(1, 2),
            name=name,
            **effect,
        ),
    )


class NDealkylation(ResonanceRule):
    """Cleaves a carbon-nitrogen bond and oxygenates the carbon side.

    These are the nitrogen rows of dealkylation. A methyl carbon is the whole
    leaving piece (``leave_count`` 1). Any larger alkyl still carries atoms
    the pattern does not name, so ``leave_count`` is None. ``breaks_ring`` is
    filled from the cleaved bond, not stored as a second class.

    Aromatic C–N hits (e.g. pyridine ring-open) match on the aromatic parent
    and react on the Kekulé parent where that bond is single — same parenting
    as :class:`Dealkylation`. ``site_kind="directed_bond"``: unique-edit uses
    ordered ``(carbon, nitrogen)`` (map 1 = oxygenated carbon); public ``site``
    is a frozenset and ``discovered_site`` holds that ordered tuple.

    Forest ``phase1_steps`` is a degenerate singleton naming this rule; the
    default :meth:`canonical_plan` matches that (not UnstableOxygenation).
    UnstableOxygenation in forest is Dealkylation + OxidativeDehalogenation;
    NDealkylation has its own ruleset.
    """
    sites_on = "bonds"
    site_kind: RuleSiteKind = "directed_bond"
    _example_substrates: tuple[str, ...] = ("CCN", "CN(C)C", "c1ccncc1")


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        _ndealk(
            Smirks("[#6H3:1][#7:2]>>([*:2].[*:1](=O)O)"),
            1,
            "methyl_carboxylic",
            adds="OO",
        ),
        _ndealk(Smirks("[#6H3:1][#7:2]>>([*:2].[*:1]=O)"), 1, "methyl_carbonyl", adds="O"),
        _ndealk(Smirks("[#6H3:1][#7:2]>>([*:2].[*:1]-O)"), 1, "methyl_alcohol", adds="O"),
        _ndealk(
            Smirks("[#6H2:1][#7:2]>>([*:2].[*:1](=O)O)"),
            None,
            "methylene_carboxylic",
            adds="OO",
        ),
        _ndealk(
            Smirks("[#6H2:1][#7:2]>>([*:2].[*:1]=O)"), None, "methylene_carbonyl", adds="O"
        ),
        _ndealk(
            Smirks("[#6H2:1][#7:2]>>([*:2].[*:1]-O)"), None, "methylene_alcohol", adds="O"
        ),
        _ndealk(Smirks("[#6H1:1][#7:2]>>([*:2].[*:1]=O)"), None, "methine_carbonyl", adds="O"),
        _ndealk(Smirks("[#6H1:1][#7:2]>>([*:2].[*:1]-O)"), None, "methine_alcohol", adds="O"),
        _ndealk(
            Smirks("[#6H0:1][#7:2]>>([*:2].[*:1]-O)"), None, "quaternary_alcohol", adds="O"
        ),
        _ndealk(
            Smirks("[#8H1:3]-[#6:1]-[#7:2]>>([*:3]=[*:1].[*:2])"),
            None,
            "hemiaminal",
            removes="H",
        ),
    )


class AzoSplitting(ResonanceRule):
    """Splits an N=N bond. Both fragments stay.

    The pattern names both nitrogens and no leaving piece, so ``leave_count``
    stays None. ``breaks_ring`` is filled from the cleaved bond. A ring N=N
    and an open azo are the same pattern; a filter reads ``breaks_ring``.
    ``=,:`` matches aromatic ring N=N; the Kekulé parent keeps that bond double.
    """
    sites_on = "bonds"
    site_kind: RuleSiteKind = "bond"
    _example_substrates: tuple[str, ...] = ('N=Nc1ccccc1',)


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#7:1]=,:[#7:2]>>[*:1].[*:2]"),
            describe(cleaves=True, partner="N", site_map=(1, 2), name="azo"),
        ),
    )


class BenzodioxoleReduction(SmirksReactionRule):
    """Cleaves both C-O bonds of the methylene in a 1,3-dioxole.

    That carbon is the whole leaving piece, and the pattern names it, so
    ``leave_count`` is 1. ``breaks_ring`` is filled from one cleaved bond.
    Both bonds are in that ring. A filter reads ``leave_count``.
    """
    sites_on = "bonds"
    site_kind: RuleSiteKind = "directed_bond"
    _example_substrates: tuple[str, ...] = ('c1ccc2c(c1)OCO2',)


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6R:1]-[#8R:2]-[#6H2R:3]-[#8R:4]-[#6R:5]>>([*:1]-[*:2].[*:3].[*:4]-[*:5])"),
            describe(
                cleaves=True,
                partner="O",
                leave_count=1,
                site_map=(2, 3),
                name="dioxole_methylene",
            ),
        ),
    )


class NitroaromaticReduction(SmirksReactionRule):
    """Cleaves one N-O of a nitro group on a ring carbon, leaving the nitroso.

    That oxygen is the whole leaving piece, and both patterns name it, so
    ``leave_count`` is 1. ``breaks_ring`` is filled from the cleaved bond.
    A filter reads ``leave_count``.
    """
    sites_on = "bonds"
    site_kind: RuleSiteKind = "directed_bond"
    _example_substrates: tuple[str, ...] = ('[O-][N+](=O)c1ccccc1',)


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#8-1:1]-[#7+1:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]"),
            describe(
                cleaves=True,
                partner="N",
                leave_count=1,
                site_map=(1, 2),
                name="nitro_charged",
            ),
        ),
        (
            Smirks("[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]"),
            describe(
                cleaves=True,
                partner="N",
                leave_count=1,
                site_map=(1, 2),
                name="nitro_neutral",
            ),
        ),
    )


class ThiopheneSulfurOxidation(ResonanceRule):
    """Oxidizes the sulfur of a thiophene to the S-oxide.

    The pattern adds oxygen and names no leaving atom, so ``leave_count``
    stays None. It does not cleave. A filter reads ``adds``. ``=,:`` matches
    the aromatic ring; maps 1–2 are the S–C single in the Kekulé parent.
    """
    sites_on = "atoms"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('c1ccsc1',)


    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6:2]1=,:[#6:3][#6:4]=,:[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*&H0&+:1]1[O-]"),
            describe(adds="O", symbol="S", name="thiophene_s_oxide"),
        ),
    )


_HALIDE = (9, 17, 35, 53, 85)


class Dephosphorylation(SmirksReactionRule):
    """Cleaves an ester O-P bond of a phosphate. The oxygen stays on the organic fragment.

    Map 1 must be the carbon-bound oxygen. Plain P-OH matches are excluded so
    RDKit uniquify cannot keep a water / methyl-phosphite cleavage instead of
    the ester (see docs/forest/DIVERGENCES.md).
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('COP(=O)(O)O',)
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks(
                "[#8;$([#8][#6]):1][#15:2](=[#8:3])([#8:4])[#8:5]>>"
                "[*:1].[*:2](=[*:3])([*:4])[*:5]"
            ),
            describe(
                *branches(_whens(2, (15,)), cleaves=True), name="phosphate_ester"
            ),
        ),
    )


class EpoxideOpening(SmirksReactionRule):
    """Opens an epoxide. One pattern only rearranges bonds; the other also adds OH."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('C1OC1', 'c1ccccc1C1CO1')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])"),
            describe(adds="", name="rearrange"),
        ),
        (
            Smirks("[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)"),
            describe(adds="O", name="hydrate"),
        ),
    )


class Hydrolysis(SmirksReactionRule):
    """Cleaves the single bond of a carboxylic derivative. One pattern also adds O."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CC(=O)OC',)
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2](O).[*:3])"),
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
            Smirks("[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2].[*:3])"),
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


class Dehydration(SmirksReactionRule):
    """Drops an OH, or a carbonyl oxygen, off carbon or nitrogen."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO',)
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]"),
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
            Smirks("[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]"),
            describe(removes="OH", cleaves=True, partner="O", name="beta_elimination"),
        ),
        (
            Smirks("[#6,#7:1]=[#8:2]>>[*:1].[*:2]"),
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

    ``path_end`` declares ``dearomatizes=True`` (capability): the path flip on
    an aromatic system is a dearomatizing reduction. ``merge_effects`` resolves
    that against ``system_aromatic``, same as Dehydrogenation / QF ends.
    ``filter_sites`` then reads ``options["adds"]`` / ``dearomatizes`` against
    ``atom_diff`` — not a Hydrogenation-named branch.
    """

    phase1_sites_on = "atoms"
    sites_on = "atom_pairs"
    site_kind: RuleSiteKind = "atom_pair"
    _example_substrates: tuple[str, ...] = ('C=C', 'C#C')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6:1]#[#6:2]>>[*:1]=[*:2]"),
            describe(adds="HH", name="alkyne", site_map=(1, 2)),
        ),
        (
            Smirks("[#6:1]=,:[#6:2]>>[*:1]-[*:2]"),
            # Span ``dearomatizes`` is capability for filter_rules. Leaving it
            # false here: aliphatic C=C→CC must not be refused when the target
            # has no loses_aromaticity. Aromatic one-bond hits still add HH;
            # ``filter_sites`` reads adds vs h_delta. Path dearomatizing
            # reductions are the ``path_end`` capability below.
            describe(adds="HH", name="alkene", site_map=(1, 2)),
        ),
    )
    endpoints: tuple[tuple[Smarts, PatternInfo], ...] = (
        (
            Smarts("[*:1]"),
            # Adds H (reduction). Capability to dearomatize on aromatic
            # systems; merge_effects resolves against system_aromatic.
            describe(adds="H", dearomatizes=True, edit="keep", name="path_end"),
        ),
    )


class TautomerRule(ResonancePairRule):
    """Stub for long-range tautomerization (forest ``Tautomerization``).

    Raises ``NotImplementedError`` from :meth:`metabolites`. Captures the
    design so it is not lost in chat history.

    **What forest did.** ``Tautomerization`` was a ``ResonanceRule`` that
    walked ``resonate_with_pair_paths``, then extended each alternating path
    by one H-bearing neighbor (``[#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#8H]``)
    and flipped bonds along that full path. Net heavy-atom formula and H
    count were unchanged (enol↔keto style). It did not use SMARTS
    ``RunReactants``; the edit was path swap only. Example pairs in forest
    tests: ``O=C1CCCCC1`` ↔ ``OC1=CCCCC1``, long-range ``ClCC=CC=CC=CC=CO``
    ↔ ``ClCCC=CC=CC=CC=O``.

    **What this is not.** Separate deferred work (TODO / LOG): tautomer
    *SMARTS matching* — match a tautomer of a reactant pattern via
    RDKit ``TautomerQuery`` (or equivalent), then infer the tautomerization
    on hit. Preferred direction only; Status: not decided; wait for mature
    tests. That matching helper is orthogonal to a first-class tautomer
    *rule* that emits tautomer metabolites.

    **Stub.** Not in PhaseOne. Patternless until chemistry and
    unique-edit (ordered vs unordered ends) are decided. When implemented,
    prefer data on ``PatternInfo`` / endpoints over a silent search branch
    (see docs/forest/HEURISTICS.md).
    """
    sites_on = "atom_pairs"
    site_kind: RuleSiteKind = "atom_pair"
    _example_substrates: tuple[str, ...] = ()


    name = "TautomerRule"
    longname = "Tautomerization"
    endpoints = ()
    smirks = ()

    def metabolites(self, mol: Mol, *args: Any, **kwargs: Any):  # type: ignore[override]
        raise NotImplementedError(
            "TautomerRule is a design stub; see class docstring and "
            "docs/forest/HEURISTICS.md / TODO.md (tautomer SMARTS vs tautomer rule)"
        )


class NitrogenReduction(ResonanceRule):
    """Cleaves N-O of nitro, nitroso, and hydroxylamine groups.

    The nitroso pattern uses ``[*:2]``. The old ``[*2]`` string emitted a dummy atom.
    Hydroxylamine uses ``-,:`` so aromatic N–O (isoxazole / benzisoxazole) matches;
    the Kekulé parent keeps that bond single.
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCNO', '[O-][N+](=O)c1ccccc1')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])"),
            describe(removes="O", cleaves=True, partner="O", name="nitro_charged"),
        ),
        (
            Smirks("[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])"),
            describe(removes="O", cleaves=True, partner="O", name="nitro_anion"),
        ),
        (
            Smirks("[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])"),
            describe(removes="O", cleaves=True, partner="O", name="nitro_neutral"),
        ),
        (
            Smirks("[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])"),
            describe(removes="OO", cleaves=True, partner="O", name="nitro_to_amine"),
        ),
        (
            Smirks("[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])"),
            describe(removes="OO", cleaves=True, partner="O", name="nitro_both"),
        ),
        (
            Smirks("[#7:1]-,:[#8:2]>>([*:1].[*:2])"),
            describe(removes="O", cleaves=True, partner="O", name="hydroxylamine"),
        ),
        (
            Smirks("[#7D2:1]=[#8:2]>>([*:1].[*:2])"),
            describe(removes="O", cleaves=True, partner="O", name="nitroso"),
        ),
        (
            Smirks("[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])"),
            describe(removes="OO", cleaves=True, partner="O", name="nitro_both_any"),
        ),
    )


class OxygenReduction(SmirksReactionRule):
    """Turns C=O / N=O into a single bond, or cleaves a peroxide."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CC(=O)OC', 'CC=O')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#8:1]=[#6,#7:2]>>[*:1]-[*:2]"),
            describe(
                *branches(
                    ({"map": 2, "z": 6}, {"map": 2, "z": 7}),
                    adds="HH",
                ),
                name="carbonyl",
            ),
        ),
        (
            Smirks("[#8:1]-[#8:2]>>[*:1].[*:2]"),
            describe(cleaves=True, partner="O", name="peroxide"),
        ),
    )

# TODO: In chemistry (not current smarts), can ReductiveDehalogenation ever work on
# an aromatic bond? If so, maybe this should be a ResonanceRule.
class ReductiveDehalogenation(SmirksReactionRule):
    """Cleaves a carbon-halogen bond. The second pattern also makes a double bond."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCCl', 'Clc1ccccc1')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]"),
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


class SulfurReduction(SmirksReactionRule):
    """Cleaves S=O, S-S, and S-C / S-O single bonds."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CS(=O)C', 'CCSO')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#16:1]=[#8:2]>>[*:1].[*:2]"),
            describe(removes="O", cleaves=True, partner="O", name="sulfoxide"),
        ),
        (
            Smirks("[#16:1]-[#16:2]>>[*:1].[*:2]"),
            describe(cleaves=True, partner="S", name="disulfide"),
        ),
        (
            Smirks("[#16:1]-[#6,#8:2]>>[*:1].[*:2]"),
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
    The reaction runs on the cached kekulé parent for that bond. The site is
    the undirected bond (``site_kind="bond"``); unique-edit keys unordered
    bond-end ranks so symmetry-related embeddings collapse.

    Forest ``phase1_steps`` is a degenerate singleton naming this rule; the
    default :meth:`canonical_plan` matches that (not StableOxygenation).
    StableOxygenation is the group that *contains* Epoxidation among peers.

    NOTE: Downstream epoxidation model only considers carbone-carbon epoxides.
    Phase 1 model additionally considers carbon-nitrogen.
    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "bond"
    _example_substrates: tuple[str, ...] = ('C=C', 'c1ccccc1')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6:1]=,:[#6,#7:2]>>[*:1]1-[*:2][O]1"),
            describe(
                *branches(
                    ({"map": 2, "z": 6}, {"map": 2, "z": 7}),
                    adds="O",
                ),
                name="epoxide",
                site_map=(1, 2),
            ),
        ),
    )


class SulfurOxidation(SmirksReactionRule):
    """Adds oxygen to divalent or tetravalent sulfur (S-oxide, S-OH, or S=O)."""

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCS', 'CSC')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#16;v2,v4:1]>>[*&H0&+:1][O-]"),
            describe(adds="O", symbol="S", name="zwitterion"),
        ),
        (
            Smirks("[#16;v2,v4:1]>>[*:1][O]"),
            describe(adds="O", symbol="S", name="hydroxy"),
        ),
        (
            Smirks("[#16;v2,v4:1]>>[*:1]=O"),
            describe(adds="O", symbol="S", name="oxo"),
        ),
    )


class NitrogenOxidation(SmirksReactionRule):
    """N-H to hydroxylamine, primary amine to nitroso, or tertiary N to N-oxide."""

    phase1_sites_on = "atoms"
    sites_on = "atoms"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCN', 'CN(C)C')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#7v3h:1]>>[*:1]O"),
            describe(
                *branches(
                    ({"map": 1, "z": 7, "h": 1}, {"map": 1, "z": 7, "h": 2}),
                    adds="O",
                ),
                name="hydroxylamine",
            ),
        ),
        (
            Smirks("[#7v3H2:1]>>[*:1]=O"),
            describe(adds="O", h=2, symbol="N", name="nitroso"),
        ),
        (
            Smirks("[#7v3H0:1]>>[*&H0&+:1][O-]"),
            describe(adds="O", h=0, symbol="N", name="n_oxide"),
        ),
    )


class OxidativeDehalogenation(SmirksReactionRule):
    """Replaces a carbon-bound halogen with OH, carbonyl, or a carboxylic acid."""

    phase1_sites_on = "bonds"
    sites_on = "bonds"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCCl', 'Clc1ccccc1')
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]"),
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
            Smirks("[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]"),
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


class ConjugationRule(SmirksReactionRule):
    """Attaches an acetyl and, by default, collapses that group to ``*``.

    ``as_star`` and ``star_label`` belong on this class. ``Protein``,
    ``DNA``, and ``Cyanide`` stay stars. The site heteroatom is ``symbol``.
    A filter reads that. This reaction does not cleave. Products are
    terminal (``is_terminal_rule``): conjugation ends further expansion.
    """
    sites_on = "atoms"
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO',)


    is_terminal_rule: bool = True
    as_star: bool = True
    star_label: str | None = None
    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]"),
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
        *args: Any,
        **kwargs: Any,
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
        filter_rules: FilterRules = _accept_all_rules,
        filter_sites: FilterSites = _accept_all_sites,
        context_mol: Mol | None = None,
        **kwargs: Any,
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
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO', 'Nc1ccccc1')


class Sulfation(ConjugationRule):
    """Adds a sulfate to an alcohol or a phenol.

    Calls :class:`ConjugationRule`. The sulfate SMARTS live here. The star
    collapse stays there. ``as_star=False`` keeps the sulfate. The oxygen
    has one hydrogen. A filter reads ``partner_h``.
    """
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO', 'Oc1ccccc1')

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6:1][#8H1:2]>>[*:1][*:2]S(=O)(=O)O"),
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
            Smirks(
                "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>"
                "[*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1"
            ),
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
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('CCO', 'Oc1ccccc1')

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#8H1:1][#6:2]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2])1"),
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
            Smirks(
                "[#8H1,#8-:1][#6:2](=[#8:3])[#6:4]>>"
                "O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2](=[#8:3])[*:4])1"
            ),
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
    site_kind: RuleSiteKind = "atom"
    _example_substrates: tuple[str, ...] = ('C=C', 'C1OC1')

    smirks: tuple[tuple[Smirks, PatternInfo], ...] = (
        (
            Smirks("[#6H1:1]1[#8:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]")),
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_ch"),
        ),
        (
            Smirks("[#6H2:1]1[#8:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]")),
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_ch2"),
        ),
        (
            # H0 only: bare [#6]([!#1]) also matches [#6H1] with a substituent
            # (e.g. styrene oxide) and double-emits the same GSH adduct.
            Smirks(
                "[#6H0:1]([!#1:4])1[#8:2][#6:3]1>>"
                + _gsh("[*:1]([*:4])[*:3][*:2]")
            ),
            describe(adds=_GSH_ADDS, site_map=1, name="epoxide_c"),
        ),
        (
            Smirks("[#6:1][#9,#17,#35,#53:2]>>" + _gsh("[*:1]")),
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
            Smirks("[#16h1:1]>>" + _gsh("[*:1]")),
            describe(adds=_GSH_ADDS, removes="H", site_map=1, name="thiol"),
        ),
        (
            Smirks("[#6H2:1]=[#6:2]>>" + _gsh("[*:1]-[*:2]")),
            describe(adds=_GSH_ADDS, site_map=1, name="alkene"),
        ),
        (
            Smirks(
                "[#6H1:1]=[#6:2][#6:3]=[#8,#7:4]>>"
                + _gsh("[*:1][*:2]=[*:3][*:4]")
            ),
            describe(
                *branches(_whens(4, (8, 7)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
                name="michael",
            ),
        ),
        (
            Smirks("[#6;H1,H2:1]=[#8:2]>>" + _gsh("[*:1]([*:2])")),
            describe(adds=_GSH_ADDS, site_map=1, name="carbonyl"),
        ),
        (
            Smirks("[#6H1:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]")),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_ch"),
        ),
        (
            Smirks("[#6H2:1]1[#7:2][#6:3]1>>" + _gsh("[*:1][*:3][*:2]")),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_ch2"),
        ),
        (
            # H0 only: same partition as epoxide_c (see above).
            Smirks(
                "[#6H0:1]([!#1:4])1[#7:2][#6:3]1>>"
                + _gsh("[*:1]([*:4])[*:3][*:2]")
            ),
            describe(adds=_GSH_ADDS, site_map=1, name="aziridine_c"),
        ),
        (
            Smirks("[#6:1][#8:2]S(=O)(=O)>>" + _gsh("[*:1]")),
            describe(adds=_GSH_ADDS, site_map=1, name="mesylate"),
        ),
        (
            Smirks(
                "[#7:1]=[#6:2]=[#8,#16:3]>>" + _gsh("[*:2](=[*:3])[*:1]")
            ),
            describe(
                *branches(_whens(3, (8, 16)), site_map=1, adds=_GSH_ADDS),
                site_map=1,
                name="isocyanate",
            ),
        ),
    )
