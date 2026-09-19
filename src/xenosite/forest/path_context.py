"""MCS-backed context for guided reactant→product path search.

Uses RDKit FindMCS with bond-order tolerance so aromatic parents still match
Kekulé / dearomatized / quinoid products.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional, Sequence, Union

from rdkit import Chem
from rdkit.Chem import rdFMCS


def _formula(mol) -> Counter:
    return Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def oxygen_deficit(mol, target) -> int:
    """How many more oxygen heavy-atoms ``target`` has than ``mol``."""
    return int(_formula(target).get("O", 0) - _formula(mol).get("O", 0))


def _oxygen_neighbor_count(mol, idx: int) -> int:
    atom = mol.GetAtomWithIdx(int(idx))
    return sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8)


def _mcs_oxygen_images(ctx) -> dict:
    """Reactant atom → target atoms it maps to, across every MCS match.

    Cached on ``ctx``. An atom missing from the map is in no embedding.
    """
    cached = getattr(ctx, "_mcs_oxygen_images", None)
    if cached is not None:
        return cached
    images: dict = {}
    reactant = getattr(ctx, "reactant", None)
    target = getattr(ctx, "target", None)
    mcs = getattr(ctx, "mcs_mol", None)
    if reactant is not None and target is not None and mcs is not None:
        try:
            n_query = int(mcs.GetNumAtoms())
        except Exception:
            n_query = 0
        if n_query > 0:
            try:
                r_matches = reactant.GetSubstructMatches(mcs)
                t_matches = target.GetSubstructMatches(mcs)
            except Exception:
                r_matches, t_matches = (), ()
            for rm in r_matches:
                for tm in t_matches:
                    if len(rm) != len(tm):
                        continue
                    for pos, ri in enumerate(rm):
                        images.setdefault(int(ri), set()).add(int(tm[pos]))
    if not images:
        for ri, ti in (getattr(ctx, "r_to_t", None) or {}).items():
            images.setdefault(int(ri), set()).add(int(ti))
    try:
        ctx._mcs_oxygen_images = images
    except Exception:
        pass
    return images


def atom_has_enough_oxygens(mol, atom_idx, ctx) -> bool:
    """True if another oxygen on ``atom_idx`` cannot close a mapped gap.

    Every MCS image of a mapped atom already has at least as many oxygen
    neighbors. An unmapped atom that already carries an oxygen is also
    enough: the target did not keep that atom as a place that still needs one.
    """
    if mol is None or ctx is None:
        return False
    try:
        idx = int(atom_idx)
        n_here = _oxygen_neighbor_count(mol, idx)
    except Exception:
        return False
    images = _mcs_oxygen_images(ctx).get(idx)
    target = getattr(ctx, "target", None)
    if not images or target is None:
        return n_here >= 1
    try:
        return all(_oxygen_neighbor_count(target, ti) <= n_here for ti in images)
    except Exception:
        return False


def heavy_formula_equal(mol, target) -> bool:
    return _formula(mol) == _formula(target)


def _aromatic_bond_count(mol) -> int:
    return sum(1 for b in mol.GetBonds() if b.GetIsAromatic())


@dataclass(frozen=True)
class FormulaHint:
    """One concrete heavy-atom formula effect (product − reactant).

    - ``delta``: element count changes for a single non-fragmenting product.
    - ``cleave``: may split / shrink the kept piece (target must be smaller).
    - Empty ``delta`` and ``cleave=False`` ⇒ formula-neutral (bond / H only).
    """

    delta: Mapping[str, int] = field(default_factory=dict)
    cleave: bool = False

    def compatible(self, mol, target) -> bool:
        """True if this effect could move ``mol`` toward ``target``."""
        try:
            fm, ft = _formula(mol), _formula(target)
        except Exception:
            return True
        if self.cleave:
            try:
                return target.GetNumHeavyAtoms() < mol.GetNumHeavyAtoms()
            except Exception:
                return True
        if not self.delta:
            return fm == ft
        for el, change in self.delta.items():
            need = int(ft.get(el, 0) - fm.get(el, 0))
            if change > 0 and need <= 0:
                return False
            if change < 0 and need >= 0:
                return False
        return True


@dataclass(frozen=True)
class FormulaAny:
    """Match may realize **any** of these effects (static polymorphism).

    Pre-match SMARTS skip is sound iff **no** choice is compatible.
    Typical for cleavage+oxidation SMARTS (e.g. dealkylation).
    """

    choices: tuple


@dataclass(frozen=True)
class FormulaMatch:
    """Match-dependent effect with an explicit pre-match possibility set.

    ``resolve(mol, match)`` returns a :class:`FormulaHint` or sequence of them
    for a concrete hit. When ``match`` is omitted (before ``RunReactants`` /
    without a query hit), ``possible`` is used so skips stay sound.
    """

    possible: tuple
    resolve: Callable  # (mol, match) -> FormulaHint | Sequence[FormulaHint]


# Common Phase I effects (share instances across rules).
ADD_O = FormulaHint(delta={"O": 1})
REMOVE_O = FormulaHint(delta={"O": -1})
NEUTRAL = FormulaHint()
CLEAVE = FormulaHint(cleave=True)

# Shorthand for cleavage SMARTS that may also oxygenate a fragment.
CLEAVE_OR_ADD_O = FormulaAny(choices=(CLEAVE, ADD_O))

FormulaSpec = Union[FormulaHint, FormulaAny, FormulaMatch, Sequence, Callable, None]


def expand_formula_effects(
    spec: FormulaSpec, mol=None, match=None
) -> tuple:
    """Resolve a ``formula_hint`` option to concrete :class:`FormulaHint`\\ s.

    ``match`` is optional context (query hit, mapids, …). When omitted, match
    resolvers must return their **possibility set** (or be permissive).

    Accepted ``spec`` forms (also nested under :class:`FormulaAny`):

    - :class:`FormulaHint` — static
    - :class:`FormulaAny` / ``(hint, hint, …)`` — any-of polymorphism
    - :class:`FormulaMatch` — ``possible`` pre-match; ``resolve(mol, match)`` after
    - ``callable(mol, match)`` — same contract as ``FormulaMatch.resolve``;
      if ``match is None`` and the callable returns ``None``, treat as unknown
      (caller should not skip)
    """
    if spec is None:
        return ()
    if isinstance(spec, FormulaHint):
        return (spec,)
    if isinstance(spec, FormulaAny):
        out = []
        for c in spec.choices:
            out.extend(expand_formula_effects(c, mol, match))
        return tuple(out)
    if isinstance(spec, FormulaMatch):
        if match is None:
            out = []
            for c in spec.possible:
                out.extend(expand_formula_effects(c, mol, None))
            return tuple(out)
        return expand_formula_effects(spec.resolve(mol, match), mol, match)
    if callable(spec):
        try:
            resolved = spec(mol, match)
        except TypeError:
            # Older resolvers that only accept (mol,)
            try:
                resolved = spec(mol)
            except Exception:
                return ()
        except Exception:
            return ()
        if resolved is None:
            return ()
        return expand_formula_effects(resolved, mol, match)
    if isinstance(spec, (list, tuple)) and not isinstance(spec, FormulaHint):
        # Sequence of specs → any-of
        out = []
        for c in spec:
            out.extend(expand_formula_effects(c, mol, match))
        return tuple(out)
    raise TypeError("unknown formula_hint spec %r" % (type(spec),))


def formula_compatible(spec: FormulaSpec, mol, target, match=None) -> bool:
    """True if ``spec`` might advance ``mol`` toward ``target``.

    Sound for pre-match skips: False only when every expanded hint is
    incompatible. Unknown / empty expansion ⇒ True (do not skip).
    """
    hints = expand_formula_effects(spec, mol, match)
    if not hints:
        return True
    return any(h.compatible(mol, target) for h in hints if isinstance(h, FormulaHint))


def hint_includes_cleave(spec: FormulaSpec) -> bool:
    """True if ``spec`` can realize a cleaving / fragmenting effect."""
    for h in expand_formula_effects(spec, None, None):
        if isinstance(h, FormulaHint) and h.cleave:
            return True
    return False


def any_hint_compatible(
    hints: Optional[Sequence], mol, target, match=None
) -> bool:
    """True if any declared hint/spec is compatible (or hints are unset)."""
    if not hints:
        return True
    for h in hints:
        if h is None:
            return True
        if formula_compatible(h, mol, target, match=match):
            return True
    return False


def _site_atom_indices(site) -> frozenset:
    if isinstance(site, tuple) and len(site) >= 2 and isinstance(site[0], str):
        site = site[1]
    out = []
    for i in site:
        try:
            out.append(int(i))
        except (TypeError, ValueError):
            # AtomRef / deferred refs: skip non-int indices for graph tests.
            idx = getattr(i, "idx", None)
            if idx is None:
                idx = getattr(i, "origin", None)
            if idx is not None:
                try:
                    out.append(int(idx))
                except (TypeError, ValueError):
                    pass
    return frozenset(out)


def bond_bridge_partitions(mol, a: int, b: int):
    """Atom sets on each side of bond ``a–b`` if it is a bridge; else ``None``.

    Uses a graph flood that never crosses the bond — no ``RunReactants``, no
    fragment mols. When the bond lies in a cycle, each flood reaches the other
    endpoint via another path and this returns ``None`` (ring-open / non-bridge).
    """
    a, b = int(a), int(b)
    if mol.GetBondBetweenAtoms(a, b) is None:
        return None

    def flood(start, blocked):
        seen = {start}
        stack = [start]
        while stack:
            i = stack.pop()
            atom = mol.GetAtomWithIdx(i)
            for nbr in atom.GetNeighbors():
                j = nbr.GetIdx()
                if j == blocked or j in seen:
                    continue
                seen.add(j)
                stack.append(j)
        return frozenset(seen)

    sa, sb = flood(a, b), flood(b, a)
    # Ring / non-bridge: the other endpoint is reachable without crossing.
    if a in sb or b in sa or (sa & sb):
        return None
    return sa, sb


def cleavage_site_bond(mol, site):
    """Return ``(a, b)`` for the bonded pair in a 2-atom cleavage site, else ``None``."""
    atoms = sorted(_site_atom_indices(site))
    if len(atoms) < 2:
        return None
    for i, a in enumerate(atoms):
        for b in atoms[i + 1 :]:
            if mol.GetBondBetweenAtoms(a, b) is not None:
                return a, b
    return None


def cleavage_site_may_reach(mol, site, target, ctx, *, oxygen_slack: int = 1) -> bool:
    """True if cutting ``site`` could leave an MCS embedding of T on one piece.

    Cheap structural test (SMILES graph + :class:`PathContext` MCS) — does **not**
    materialize cleavage products.

    Uses **every** reactant MCS embedding on ``ctx`` (full-size and smaller
    secondary placements). A site is kept if any embedding's conserved core
    (minus cleaved site atoms) lies on one bridge side with enough heavy atoms.

    Ring-opening N–C bonds are never pruned (N-dealk ring open / intermediate).
    Other ring-opens stay only when the ring intersects some embedding.
    """
    if ctx is None:
        return True
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        cons = getattr(ctx, "conserved_r_atoms", None)
        if cons:
            embeddings = [frozenset(cons)]
    if not embeddings:
        return True
    try:
        t_n = int(target.GetNumHeavyAtoms())
    except Exception:
        return True

    site_atoms = _site_atom_indices(site)
    bond = cleavage_site_bond(mol, site)
    if bond is None:
        return True
    a, b = bond
    parts = bond_bridge_partitions(mol, a, b)
    if parts is None:
        # N-dealk ring open: never prune.
        if any(mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 7 for i in site_atoms):
            return True
        # Other ring-opens: only full-size MCS embeddings (remainder / smaller
        # placements often latch onto unrelated rings and re-open fanout).
        full = max((len(e) for e in embeddings), default=0)
        try:
            ri = mol.GetRingInfo()
            rings = [set(ring) for ring in ri.AtomRings()]
        except Exception:
            return True
        for cons in embeddings:
            if len(cons) < full:
                continue
            for ring in rings:
                if a in ring and b in ring and (ring & set(cons)):
                    return True
        return False

    for cons in embeddings:
        cons = set(cons)
        # Site atoms are the reaction center; MCS may map a leaving heteroatom
        # onto T (e.g. ester O ↔ acid OH). Require the rest on one fragment.
        cons_core = cons - set(site_atoms)
        for side in parts:
            if cons_core and not (cons_core <= side):
                continue
            if not cons_core and not (cons & side):
                continue
            heavy = sum(
                1 for i in side if mol.GetAtomWithIdx(int(i)).GetAtomicNum() > 1
            )
            if heavy + int(oxygen_slack) >= t_n:
                return True
    return False


def mcs_embedding_disagreement(ctx) -> frozenset:
    """Atoms in some MCS embedding but not all (asymmetric / unequal placements)."""
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if len(embeddings) < 2:
        return frozenset()
    full = max(len(e) for e in embeddings)
    full_embs = [set(e) for e in embeddings if len(e) == full]
    if len(full_embs) < 2:
        ordered = sorted(embeddings, key=len, reverse=True)
        if len(ordered) < 2:
            return frozenset()
        a, b = set(ordered[0]), set(ordered[1])
        return frozenset((a | b) - (a & b))
    inter = set.intersection(*full_embs)
    return frozenset(set.union(*full_embs) - inter)


def _bonds_chemically_equal(rb, tb) -> bool:
    """True if reactant/target bonds agree ignoring FindMCS CompareAny tolerance."""
    if tb is None:
        return False
    if bool(rb.GetIsAromatic()) != bool(tb.GetIsAromatic()):
        return False
    if rb.GetIsAromatic():
        return True
    return rb.GetBondType() == tb.GetBondType()


def mcs_chem_disagree_bonds(ctx) -> frozenset:
    """Reactant bonds whose MCS-mapped endpoints disagree chemically with T.

    FindMCS uses ``BondCompare.CompareAny``, so aromatic↔Kekulé / missing ring
    bonds still place atoms in the conserved core. Those bonds are *not* true
    interior — ring-opens and bond-order edits live here.
    """
    r = getattr(ctx, "reactant", None)
    t = getattr(ctx, "target", None)
    r_to_t = getattr(ctx, "r_to_t", None) or {}
    if r is None or t is None or not r_to_t:
        return frozenset()
    out = []
    for b in r.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a not in r_to_t or c not in r_to_t:
            continue
        tb = t.GetBondBetweenAtoms(int(r_to_t[a]), int(r_to_t[c]))
        if not _bonds_chemically_equal(b, tb):
            out.append(frozenset((int(a), int(c))))
    return frozenset(out)


def mcs_chem_disagree_atoms(ctx) -> frozenset:
    """Endpoints of :func:`mcs_chem_disagree_bonds` (match-but-different region)."""
    atoms = set()
    for bond in mcs_chem_disagree_bonds(ctx):
        atoms |= set(bond)
    return frozenset(atoms)


def _embedding_heavy(mol, embedding) -> int:
    return sum(
        1 for i in embedding if mol.GetAtomWithIdx(int(i)).GetAtomicNum() > 1
    )


def match_covers_large_part(ctx) -> bool:
    """True when the largest MCS embedding is a substantial piece of both mols.

    Same size bar as cleavage interior trust: at least half the smaller mol,
    and at least 3 heavy atoms.
    """
    mol = getattr(ctx, "reactant", None)
    target = getattr(ctx, "target", None)
    if mol is None or target is None:
        return False
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        cons = getattr(ctx, "conserved_r_atoms", None)
        embeddings = [frozenset(cons)] if cons else []
    if not embeddings:
        return False
    try:
        n = min(int(mol.GetNumHeavyAtoms()), int(target.GetNumHeavyAtoms()))
    except Exception:
        return False
    best = max(_embedding_heavy(mol, emb) for emb in embeddings)
    return best >= max(3, n // 2)


def attachment_boundary_atoms(ctx) -> frozenset:
    """Mapped reactant atoms whose target image is bonded to an unmapped atom.

    Those are the sites where the target still has atoms the match did not
    place — the boundary new atoms have to be added on.
    """
    target = getattr(ctx, "target", None)
    r_to_t = getattr(ctx, "r_to_t", None) or {}
    t_only = set(getattr(ctx, "t_only_atoms", ()) or ())
    if target is None or not r_to_t or not t_only:
        return frozenset()
    out = []
    for ri, ti in r_to_t.items():
        try:
            atom = target.GetAtomWithIdx(int(ti))
        except Exception:
            continue
        if any(nbr.GetIdx() in t_only for nbr in atom.GetNeighbors()):
            out.append(int(ri))
    return frozenset(out)


def aromatic_mismatch_atoms(mol, ctx) -> frozenset:
    """Mapped atoms that are aromatic on ``mol`` and not on the target."""
    target = getattr(ctx, "target", None)
    r_to_t = getattr(ctx, "r_to_t", None) or {}
    if mol is None or target is None or not r_to_t:
        return frozenset()
    out = []
    for ri, ti in r_to_t.items():
        try:
            src = mol.GetAtomWithIdx(int(ri))
            dst = target.GetAtomWithIdx(int(ti))
        except Exception:
            continue
        if src.GetIsAromatic() and not dst.GetIsAromatic():
            out.append(int(ri))
    return frozenset(out)


def dearomatization_systems(mol, ctx):
    """Aromatic systems that the match says must lose aromaticity.

    Empty unless the MCS covers a large part of both mols and at least half
    of a system (and at least 3 atoms) is aromatic here and not on the target.
    """
    if not match_covers_large_part(ctx):
        return ()
    mismatch = aromatic_mismatch_atoms(mol, ctx)
    if len(mismatch) < 3:
        return ()
    from .base import AromaticSystems

    kept = []
    for system in AromaticSystems().systems(mol):
        system = frozenset(int(i) for i in system)
        hit = system & mismatch
        if len(hit) >= max(3, len(system) // 2):
            kept.append(system)
    return tuple(kept)


def site_match_boundary_rank(ctx, atoms) -> tuple:
    """Sort key for expansion sites. Lower is tried first.

    Prefer atoms on the attachment boundary (unmapped target neighbors), then
    atoms whose matched bonds disagree with the target. A site that also edits
    atoms outside that boundary sorts later — it is not the match-boundary edit.
    """
    idxs = []
    for atom in atoms or ():
        try:
            idxs.append(int(atom))
        except (TypeError, ValueError):
            continue
    if not idxs:
        return (1, 1, 1)
    atoms_fs = frozenset(idxs)
    disagree = mcs_chem_disagree_atoms(ctx)
    attach = attachment_boundary_atoms(ctx)
    boundary = disagree | attach
    on_attach = 0 if atoms_fs & attach else 1
    on_disagree = 0 if atoms_fs & disagree else 1
    off = 1 if (atoms_fs - boundary) else 0
    return (on_attach, on_disagree, off)


def mcs_shared_core(ctx) -> frozenset:
    """Intersection of full-size embeddings (agreed conserved core)."""
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        cons = getattr(ctx, "conserved_r_atoms", None)
        return frozenset(cons) if cons else frozenset()
    full = max(len(e) for e in embeddings)
    full_embs = [set(e) for e in embeddings if len(e) == full]
    if not full_embs:
        return frozenset(embeddings[0])
    if len(full_embs) == 1:
        return frozenset(full_embs[0])
    return frozenset(set.intersection(*full_embs))


def cleavage_site_on_mcs_frontier(mol, site, ctx) -> bool:
    """True if ``site`` sits on an MCS embedding boundary or disagreement edge.

    "Disagreement" includes multi-embedding placement diffs **and** bonds that
    MCS matched under CompareAny but that differ in order / aromaticity / are
    missing in T (ring-open loci). Used for **priority** (try these first).
    Hard-dropping non-frontier sites from Required is *not* always sound —
    multi-hop peels and partial MCS can need non-frontier cuts. See
    :func:`cleavage_site_safe_to_drop_required`.
    """
    if ctx is None:
        return True
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        cons = getattr(ctx, "conserved_r_atoms", None)
        embeddings = [frozenset(cons)] if cons else []
    if not embeddings:
        return True

    site_atoms = _site_atom_indices(site)
    bond = cleavage_site_bond(mol, site)
    if bond is not None:
        parts = bond_bridge_partitions(mol, bond[0], bond[1])
        if parts is None and any(
            mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 7 for i in site_atoms
        ):
            return True
        # Match-but-different under CompareAny (incl. ring bond missing in T).
        if frozenset(bond) in mcs_chem_disagree_bonds(ctx):
            return True
        chem_atoms = mcs_chem_disagree_atoms(ctx)
        if bond[0] in chem_atoms or bond[1] in chem_atoms:
            return True

    disagree = mcs_embedding_disagreement(ctx)
    if site_atoms & disagree:
        return True
    if bond is not None and (bond[0] in disagree or bond[1] in disagree):
        return True

    for emb in embeddings:
        emb = set(emb)
        if bond is not None:
            a, b = bond
            if (a in emb) != (b in emb):
                return True
        if site_atoms & emb and site_atoms - emb:
            return True
        if bond is not None:
            for endpoint in bond:
                if endpoint not in emb:
                    continue
                atom = mol.GetAtomWithIdx(int(endpoint))
                for nbr in atom.GetNeighbors():
                    j = nbr.GetIdx()
                    if j not in emb and j in site_atoms:
                        return True
    return False


def cleavage_site_deep_interior(mol, site, ctx) -> bool:
    """True if the cleaved bond lies strictly inside the agreed MCS core.

    Both endpoints are in :func:`mcs_shared_core` and have no neighbor outside
    that core. The bond itself must chemically agree with T (not in
    :func:`mcs_chem_disagree_bonds`) — CompareAny can park ring-open / Kekulé
    diffs inside the atom core. Cutting a true interior bond can only damage
    the conserved scaffold — safe to drop from Required ``sites_toward``.
    """
    if ctx is None:
        return False
    core = set(mcs_shared_core(ctx))
    # Threshold: core must be substantial vs target (avoid tiny/noisy MCS).
    try:
        t_n = int(ctx.target.GetNumHeavyAtoms()) if ctx.target is not None else 0
    except Exception:
        t_n = 0
    if t_n and len(core) + 1 < max(3, t_n // 2):
        return False  # MCS too weak to trust interior/exterior geometry

    bond = cleavage_site_bond(mol, site)
    if bond is None:
        return False
    a, b = bond
    if a not in core or b not in core:
        return False
    # Not interior if the mapped bond differs / is absent in T.
    if frozenset((a, b)) in mcs_chem_disagree_bonds(ctx):
        return False
    chem_atoms = mcs_chem_disagree_atoms(ctx)
    if a in chem_atoms or b in chem_atoms:
        return False
    for endpoint in (a, b):
        atom = mol.GetAtomWithIdx(int(endpoint))
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core:
                return False
    return True


# When |largest embedding| + slack >= |T|, a single frontier cleave can yield T;
# deep-exterior sites are then droppable from Required (still not for multi-hop
# ring-open paths — those sites are not deep-exterior).
_SINGLE_CLEAVE_OXYGEN_SLACK = 2


def cleavage_site_deep_exterior(mol, site, ctx) -> bool:
    """True if the cleaved bond is outside every MCS embedding and not adjacent."""
    if ctx is None:
        return False
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        return False
    covered = set().union(*embeddings)
    bond = cleavage_site_bond(mol, site)
    if bond is None:
        return False
    a, b = bond
    if a in covered or b in covered:
        return False
    for endpoint in (a, b):
        atom = mol.GetAtomWithIdx(int(endpoint))
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in covered:
                return False
    return True


def mcs_suggests_single_cleave(ctx, target, *, oxygen_slack: int = _SINGLE_CLEAVE_OXYGEN_SLACK) -> bool:
    """True when the largest MCS embedding is already about as big as T."""
    if ctx is None or target is None:
        return False
    embeddings = list(getattr(ctx, "conserved_r_embeddings", ()) or ())
    if not embeddings:
        return False
    try:
        t_n = int(target.GetNumHeavyAtoms())
        mol = ctx.reactant
        best = 0
        for emb in embeddings:
            heavy = sum(
                1 for i in emb if mol.GetAtomWithIdx(int(i)).GetAtomicNum() > 1
            )
            if heavy > best:
                best = heavy
        return best + int(oxygen_slack) >= t_n
    except Exception:
        return False


def cleavage_site_safe_to_drop_required(mol, site, ctx, target=None) -> bool:
    """When True, omit ``site`` from Required ``sites_toward`` (still may be Maybe).

    Safe thresholds (sound under current MCS use):

    1. **Deep interior** of the shared full-size MCS core — cutting cannot
       separate T from discarded atoms; only damages the scaffold.
    2. **Deep exterior** *and* :func:`mcs_suggests_single_cleave` — T already
       matches one embedding's size, so Required cleavage must touch that
       frontier; far exterior peels are not Required.

    Everything else is only **reordered** via :func:`cleavage_site_priority`,
    not dropped. Ring-opens and disagreement edges are never dropped here.
    """
    if ctx is None:
        return False
    if cleavage_site_on_mcs_frontier(mol, site, ctx):
        return False
    if cleavage_site_deep_interior(mol, site, ctx):
        return True
    tgt = target if target is not None else getattr(ctx, "target", None)
    if mcs_suggests_single_cleave(ctx, tgt) and cleavage_site_deep_exterior(
        mol, site, ctx
    ):
        return True
    return False


def cleavage_site_priority(mol, site, ctx) -> tuple:
    """Sort key for Required sites (lower = try first).

    Prefer multi-embedding disagreement edges, then MCS geometric boundaries.
    CompareAny chem-disagree bonds stay **frontier** (never hard-dropped) but
    are not priority-boosted: on a large parent vs small T they are often
    ring-open red herrings while Required peels are exterior chain cuts.
    """
    if ctx is None:
        return (1, 1)
    disagree = mcs_embedding_disagreement(ctx)
    site_atoms = _site_atom_indices(site)
    on_disagree = 0 if (site_atoms & disagree) else 1
    on_frontier = 0 if cleavage_site_on_mcs_frontier(mol, site, ctx) else 1
    return (on_disagree, on_frontier)


def heavy_element_symbols(mol) -> set:
    return {a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1}


# Common Forest organic / halide set — anything else is "exotic" (metals, …).
_ORGANIC_ATOMIC_NUMS = frozenset({5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53})


def exotic_element_symbols(mol) -> set:
    """Heavy atoms outside the usual organic/halide palette (e.g. metals)."""
    out = set()
    for a in mol.GetAtoms():
        z = a.GetAtomicNum()
        if z > 1 and z not in _ORGANIC_ATOMIC_NUMS:
            out.add(a.GetSymbol())
    return out


def elements_may_add_from_hints(hints) -> Optional[frozenset]:
    """Element symbols formula hints may introduce; ``None`` if unknown/undeclared."""
    if not hints:
        return None
    out = set()
    saw = False
    for spec in hints:
        if spec is None:
            return None
        for h in expand_formula_effects(spec, None, None):
            if not isinstance(h, FormulaHint):
                continue
            saw = True
            for el, delta in (h.delta or {}).items():
                if int(delta) > 0:
                    out.add(str(el))
    return frozenset(out) if saw else None


def unreachable_new_elements(rules, reactant, target) -> frozenset:
    """New heavy elements in ``target`` that no rule can introduce.

    Exotic atoms (metals, …) always fail fast unless a rule explicitly lists them
    in its formula-hint adds. For ordinary C/N/O/…, fail only when every active
    rule declares hints and none can add the missing element.
    """
    try:
        need = heavy_element_symbols(target) - heavy_element_symbols(reactant)
    except Exception:
        return frozenset()
    if not need:
        return frozenset()

    known_add = set()
    any_unknown = False
    for rule in rules or ():
        try:
            adds = rule.elements_may_add()
        except Exception:
            any_unknown = True
            continue
        if adds is None:
            any_unknown = True
            continue
        known_add |= set(adds)

    exotic_need = set()
    try:
        exotic_need = exotic_element_symbols(target) - exotic_element_symbols(reactant)
    except Exception:
        exotic_need = set()
    # Metals / exotics: unreachable unless explicitly addable.
    blocked = frozenset(el for el in exotic_need if el not in known_add)

    if any_unknown:
        return blocked
    return frozenset(need - known_add) | blocked


def _topo_rank_signature(mol, atom_idxs) -> frozenset:
    """Canonical-rank multiset for an atom set (symmetry fingerprint)."""
    try:
        ranks = Chem.CanonicalRankAtoms(mol)
        return frozenset(int(ranks[int(i)]) for i in atom_idxs)
    except Exception:
        return frozenset(int(i) for i in atom_idxs)


def _dedupe_symmetric_embeddings(mol, embeddings):
    """Drop topo-equivalent embeddings; keep asymmetric placements distinct."""
    out = []
    seen = set()
    for emb in embeddings:
        sig = _topo_rank_signature(mol, emb)
        if sig in seen:
            continue
        seen.add(sig)
        out.append(frozenset(int(i) for i in emb))
    return out


def _mol_without_atoms(mol, drop):
    """Copy of ``mol`` with ``drop`` atom indices removed; map new→old index."""
    drop = {int(i) for i in drop}
    if not drop:
        return Chem.Mol(mol), {i: i for i in range(mol.GetNumAtoms())}
    rw = Chem.RWMol(Chem.Mol(mol))
    for idx in sorted(drop, reverse=True):
        if 0 <= idx < rw.GetNumAtoms():
            rw.RemoveAtom(idx)
    rem = rw.GetMol()
    try:
        Chem.SanitizeMol(rem, catchErrors=True)
    except Exception:
        pass
    # Rebuild old-index map: after reverse deletions, surviving atoms keep order
    old_idxs = [i for i in range(mol.GetNumAtoms()) if i not in drop]
    new_to_old = {new: old for new, old in enumerate(old_idxs)}
    return rem, new_to_old


def _find_mcs_smarts(mol_a, mol_b) -> str:
    params = rdFMCS.MCSParameters()
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareAny
    params.Timeout = 5
    try:
        params.BondCompareParameters.CompleteRingsOnly = False
        params.BondCompareParameters.MatchFusedRings = True
    except Exception:
        pass
    result = rdFMCS.FindMCS([mol_a, mol_b], params)
    return result.smartsString if result and result.smartsString else ""


def _collect_conserved_r_embeddings(r, t, mcs_mol, mcs_smarts, *, max_rounds: int = 4):
    """All placements of T onto R: full MCS matches, then smaller remainders.

    Secondary rounds delete atoms already covered by a prior embedding and
    re-run FindMCS so a second (possibly smaller) match is not truncated away.
    Topo-symmetric duplicates of the same size class are collapsed.
    """
    raw = []
    seen = set()

    def add(idxs):
        key = frozenset(int(i) for i in idxs)
        if len(key) < 2 or key in seen:
            return False
        seen.add(key)
        raw.append(key)
        return True

    if mcs_mol is not None and mcs_mol.GetNumAtoms() > 0:
        for m in r.GetSubstructMatches(mcs_mol):
            add(m)

    # Smaller / alternate placements on the uncovered remainder.
    covered = set().union(*raw) if raw else set()
    for _ in range(max_rounds - 1):
        if len(covered) >= r.GetNumAtoms() - 1:
            break
        rem, new_to_old = _mol_without_atoms(r, covered)
        if rem.GetNumHeavyAtoms() < 2:
            break
        smarts = _find_mcs_smarts(rem, t)
        if not smarts:
            break
        q = Chem.MolFromSmarts(smarts)
        if q is None or q.GetNumAtoms() < 2:
            break
        matches = rem.GetSubstructMatches(q)
        if not matches:
            break
        grew = False
        for m in matches:
            orig = frozenset(new_to_old[int(i)] for i in m)
            if add(orig):
                grew = True
                covered |= orig
        if not grew:
            # Block the first remainder match's atoms and retry once more shape.
            covered |= frozenset(new_to_old[int(i)] for i in matches[0])
    return tuple(_dedupe_symmetric_embeddings(r, raw))


@dataclass
class PathContext:
    """Maps between reactant ``R`` and target ``T`` via aromatic-tolerant MCS."""

    reactant: object
    target: object
    mcs_smarts: str = ""
    mcs_mol: object = None
    r_to_t: dict = field(default_factory=dict)
    t_to_r: dict = field(default_factory=dict)
    conserved_r_atoms: frozenset = field(default_factory=frozenset)
    conserved_t_atoms: frozenset = field(default_factory=frozenset)
    conserved_r_embeddings: tuple = field(default_factory=tuple)
    r_only_atoms: frozenset = field(default_factory=frozenset)
    t_only_atoms: frozenset = field(default_factory=frozenset)
    formula_delta: dict = field(default_factory=dict)
    dearomatization_delta: int = 0

    @classmethod
    def from_mols(cls, reactant, target) -> PathContext:
        r = Chem.Mol(reactant)
        t = Chem.Mol(target)
        Chem.SanitizeMol(r, catchErrors=True)
        Chem.SanitizeMol(t, catchErrors=True)

        smarts = _find_mcs_smarts(r, t)
        mcs_mol = Chem.MolFromSmarts(smarts) if smarts else None

        embeddings = _collect_conserved_r_embeddings(r, t, mcs_mol, smarts)

        r_to_t: dict = {}
        t_to_r: dict = {}
        # Primary mapping: largest embedding paired with first T match of that MCS.
        if mcs_mol is not None and mcs_mol.GetNumAtoms() > 0 and embeddings:
            primary = max(embeddings, key=len)
            t_matches = t.GetSubstructMatches(mcs_mol)
            r_match = None
            for m in r.GetSubstructMatches(mcs_mol):
                if frozenset(int(i) for i in m) == primary:
                    r_match = m
                    break
            if r_match is None and r.GetSubstructMatch(mcs_mol):
                r_match = r.GetSubstructMatch(mcs_mol)
            t_match = t_matches[0] if t_matches else t.GetSubstructMatch(mcs_mol)
            if r_match and t_match and len(r_match) == len(t_match):
                for ri, ti in zip(r_match, t_match):
                    r_to_t[int(ri)] = int(ti)
                    t_to_r[int(ti)] = int(ri)

        # Union of all embeddings (asymmetric placements included).
        conserved_r = frozenset().union(*embeddings) if embeddings else frozenset(r_to_t)
        conserved_t = frozenset(t_to_r)
        r_only = frozenset(range(r.GetNumAtoms())) - conserved_r
        t_only = frozenset(range(t.GetNumAtoms())) - conserved_t

        fr, ft = _formula(r), _formula(t)
        delta = {el: ft[el] - fr[el] for el in set(fr) | set(ft) if ft[el] != fr[el]}

        return cls(
            reactant=r,
            target=t,
            mcs_smarts=smarts,
            mcs_mol=mcs_mol,
            r_to_t=r_to_t,
            t_to_r=t_to_r,
            conserved_r_atoms=conserved_r,
            conserved_t_atoms=conserved_t,
            conserved_r_embeddings=embeddings,
            r_only_atoms=r_only,
            t_only_atoms=t_only,
            formula_delta=delta,
            dearomatization_delta=_aromatic_bond_count(t) - _aromatic_bond_count(r),
        )

    def distance_proxy(self, mol) -> int:
        """Smaller is closer: unmatched heavy atoms vs target formula."""
        fm, ft = _formula(mol), _formula(self.target)
        return sum(abs(ft[el] - fm[el]) for el in set(fm) | set(ft))
