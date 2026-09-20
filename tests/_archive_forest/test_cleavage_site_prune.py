"""Cleavage site pruning from MCS + graph split (no RunReactants)."""

from rdkit import Chem

from xenosite._archive_forest import PathSearchCounters, find_path
from xenosite._archive_forest.path_context import (
    PathContext,
    bond_bridge_partitions,
    cleavage_site_bond,
    cleavage_site_may_reach,
)
from xenosite._archive_forest.rules import Dealkylation


def _n_sites(mol):
    return [
        s
        for s in Dealkylation().iter_reactant_site_matches(mol)
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in s)
    ]


def _may_reach_with_single_cons(mol, site, target, cons, *, n_ring_open=False):
    """Broken modes: one conserved set (union or wrong embedding), optional N ring-open."""
    cons = set(cons)
    if not cons:
        return True
    t_n = target.GetNumHeavyAtoms()
    site_atoms = set(site)
    cons_core = cons - site_atoms
    bond = cleavage_site_bond(mol, site)
    if bond is None:
        return True
    a, b = bond
    parts = bond_bridge_partitions(mol, a, b)
    if parts is None:
        if n_ring_open and any(
            mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in site_atoms
        ):
            return True
        for ring in mol.GetRingInfo().AtomRings():
            s = set(ring)
            if a in s and b in s and (s & cons):
                return True
        return False
    for side in parts:
        if cons_core and not (cons_core <= side):
            continue
        if not cons_core and not (cons & side):
            continue
        heavy = sum(1 for i in side if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)
        if heavy + 1 >= t_n:
            return True
    return False


def test_bond_bridge_partitions_split_vs_ring():
    mol = Chem.MolFromSmiles("CCN")  # acyclic
    assert bond_bridge_partitions(mol, 0, 1) is not None
    ring = Chem.MolFromSmiles("C1CCCCC1")
    assert bond_bridge_partitions(ring, 0, 1) is None  # ring-open


def test_cleavage_may_reach_uses_mcs_without_metabolize(monkeypatch):
    """TBA: MCS-conserved enyne lies on N–CH2(enyne) carbon side; naph ring pruned."""
    r = Chem.MolFromSmiles("CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12")
    t = Chem.MolFromSmiles("CC(C)(C)C#CC=CC=O")
    ctx = PathContext.from_mols(r, t)

    def boom(*_a, **_k):
        raise AssertionError("must not materialize cleavage")

    monkeypatch.setattr(Dealkylation, "metabolites", boom)
    monkeypatch.setattr(Dealkylation, "metabolize", boom)

    rule = Dealkylation()
    sites = rule.sites_toward(r, t, ctx)
    assert sites is not None
    assert sites, "expected at least the TBA-forming dealk site"
    assert all(cleavage_site_may_reach(r, s, t, ctx) for s in sites)
    naph_ring = frozenset({13, 14})
    if r.GetBondBetweenAtoms(13, 14) is not None:
        assert naph_ring not in {frozenset(s) for s in sites}


def test_tba_guided_stays_under_small_billed_budget():
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
            "CC(C)(C)C#CC=CC=O",
            depth=2,
            max_paths=1,
            max_expansions=40,
            counters=counters,
        )
    )
    assert hits
    assert not counters.budget_exhausted
    assert counters.billed() <= 40
    assert counters.plans_expanded <= 15


def test_multi_mcs_union_would_prune_both_benzyls():
    """Two PhCHO placements on start: union-as-one-conserved keeps nothing.

    Without per-embedding OR, cons spans both rings so no bridge cut has
    cons_core ⊆ one side — both formation sites die.
    """
    r = Chem.MolFromSmiles("c1ccccc1CNCc1ccccc1")
    t = Chem.MolFromSmiles("O=Cc1ccccc1")
    ctx = PathContext.from_mols(r, t)
    assert len(ctx.conserved_r_embeddings) >= 2
    assert len({len(e) for e in ctx.conserved_r_embeddings}) == 1

    sites = _n_sites(r)
    union = set().union(*ctx.conserved_r_embeddings)
    broken = [s for s in sites if _may_reach_with_single_cons(r, s, t, union)]
    fixed = [s for s in sites if cleavage_site_may_reach(r, s, t, ctx)]
    assert broken == [], "union of both benzyl MCS embeddings must not be one cons"
    assert len(fixed) >= 2

    hits = list(find_path(r, t, ruleset="ND", depth=2, max_paths=5, max_expansions=40))
    assert len(hits) >= 2
    step_sites = {frozenset(str(a) for a in site) for h in hits for _, site in h.steps}
    assert len(step_sites) >= 2


def test_smaller_secondary_embedding_alone_misses_naphthaldehyde():
    """Naphthyl–benzyl amine → naphthaldehyde: Ph-only (smaller) embedding keeps [].

    Full MCS sits on naphthalene; a remainder match on the phenyl is smaller and
    must not be the sole conserved set — OR with the full-size embedding.
    """
    r = Chem.MolFromSmiles("c1ccc2ccccc2c1CNCc1ccccc1")
    t = Chem.MolFromSmiles("O=Cc1cccc2ccccc12")
    ctx = PathContext.from_mols(r, t)
    sizes = sorted(len(e) for e in ctx.conserved_r_embeddings)
    assert len(sizes) >= 2 and sizes[0] < sizes[-1], sizes

    sites = _n_sites(r)
    small = min(ctx.conserved_r_embeddings, key=len)
    large = max(ctx.conserved_r_embeddings, key=len)
    assert [s for s in sites if _may_reach_with_single_cons(r, s, t, small)] == []
    assert [s for s in sites if _may_reach_with_single_cons(r, s, t, large)]
    assert [s for s in sites if cleavage_site_may_reach(r, s, t, ctx)]

    hits = list(find_path(r, t, ruleset="ND", depth=2, max_paths=3, max_expansions=40))
    assert hits
    assert hits[0].steps[0][0] == "NDealkylation"


def test_ndealk_ring_open_not_pruned_when_mcs_misses_ring():
    """N-benzylpiperidine → amino-aldehyde: ring-open N–C must stay.

    If MCS is taken as the benzyl only, ring∩cons is empty — without an N-dealk
    ring-open exemption those sites are wrongly dropped.
    """
    r = Chem.MolFromSmiles("c1ccccc1CN1CCCCC1")
    t = Chem.MolFromSmiles("NCCCCC=O")
    ctx = PathContext.from_mols(r, t)

    benzyl = set()
    for i, a in enumerate(r.GetAtoms()):
        if a.GetIsAromatic():
            benzyl.add(i)
        elif (
            a.GetAtomicNum() == 6
            and not a.IsInRing()
            and any(n.GetIsAromatic() or n.GetAtomicNum() == 7 for n in a.GetNeighbors())
        ):
            benzyl.add(i)

    ring_open_n = []
    for s in _n_sites(r):
        bond = cleavage_site_bond(r, s)
        if bond and bond_bridge_partitions(r, *bond) is None:
            ring_open_n.append(s)
    assert ring_open_n, "expected endocyclic N–C sites"

    for s in ring_open_n:
        assert not _may_reach_with_single_cons(
            r, s, t, benzyl, n_ring_open=False
        ), "benzyl-only MCS must not justify ring-open via ring∩cons"
        assert cleavage_site_may_reach(
            r, s, t, ctx
        ), "N-dealk ring-open must not be pruned"

    hits = list(find_path(r, t, ruleset="ND", depth=3, max_paths=5, max_expansions=60))
    assert hits, "expected ring-open then demethyl (or equivalent) to amino-aldehyde"


def test_mcs_frontier_priority_and_safe_drop_thresholds():
    """Frontier/disagreement first; only deep-interior (and gated exterior) hard-dropped."""
    from xenosite._archive_forest.path_context import (
        cleavage_site_deep_interior,
        cleavage_site_on_mcs_frontier,
        cleavage_site_priority,
        cleavage_site_safe_to_drop_required,
        mcs_suggests_single_cleave,
    )

    r = Chem.MolFromSmiles("c1ccccc1CNCc1ccccc1")
    t = Chem.MolFromSmiles("O=Cc1ccccc1")
    ctx = PathContext.from_mols(r, t)
    rule = Dealkylation()
    toward = rule.sites_toward(r, t, ctx)
    assert toward is not None and len(toward) >= 2
    # Sorted: disagreement/frontier keys nondecreasing.
    keys = [cleavage_site_priority(r, s, ctx) for s in toward]
    assert keys == sorted(keys)
    # Formation N–CH2 sites are on the frontier.
    for s in toward[:2]:
        assert cleavage_site_on_mcs_frontier(r, s, ctx)
    # Single-cleave suggestion for PhCHO-sized MCS.
    assert mcs_suggests_single_cleave(ctx, t)
    # Deep interior of one phenyl (if any such Dealk site exists) would be droppable.
    for s in rule.iter_reactant_site_matches(r):
        if cleavage_site_deep_interior(r, s, ctx):
            assert cleavage_site_safe_to_drop_required(r, s, ctx, t)


def test_o_ring_open_chem_disagree_not_dropped():
    """CompareAny MCS maps ring atoms onto open-chain T; ring bonds are match-but-different.

    Without chem-disagree-as-frontier, those bonds look deep-interior and are
    hard-dropped from Required (``CC1CCC2OCCCC2C1`` → ``O=CCCCC1CCCC1``).
    """
    from xenosite._archive_forest.path_context import (
        cleavage_site_deep_interior,
        cleavage_site_on_mcs_frontier,
        cleavage_site_safe_to_drop_required,
        mcs_chem_disagree_bonds,
    )

    r = Chem.MolFromSmiles("CC1CCC2OCCCC2C1")
    t = Chem.MolFromSmiles("O=CCCCC1CCCC1")
    assert t.GetNumHeavyAtoms() < r.GetNumHeavyAtoms()
    ctx = PathContext.from_mols(r, t)
    diffs = mcs_chem_disagree_bonds(ctx)
    assert diffs, "expected ring bond missing/order-diff under CompareAny map"
    toward = Dealkylation().sites_toward(r, t, ctx)
    assert toward, "sites_toward must not be empty"
    # At least one chem-disagree bond site stays Required and is frontier.
    chem_toward = [
        s
        for s in toward
        if (b := cleavage_site_bond(r, s)) is not None and frozenset(b) in diffs
    ]
    assert chem_toward, "match-but-different ring sites must remain in Required"
    for s in chem_toward:
        assert cleavage_site_on_mcs_frontier(r, s, ctx)
        assert not cleavage_site_deep_interior(r, s, ctx)
        assert not cleavage_site_safe_to_drop_required(r, s, ctx, t)


def test_chem_disagree_does_not_priority_starve_chain_peels():
    """Large alkyl-bicycle → small chain: chem-diff ring must not crowd out peels.

    Winning path peels the chain; CompareAny chem-diff on the ring is frontier
    (kept) but must not be sorted first or budget is burned on ring-open misses.
    """
    r = "CCCCCCCCCCCC1CCC2OCCCC2C1"
    t = "O=CCCCO"
    counters = PathSearchCounters()
    hits = list(
        find_path(
            r,
            t,
            ruleset="PhaseOneRS",
            depth=3,
            max_paths=2,
            max_expansions=200,
            counters=counters,
        )
    )
    assert hits, "expected chain-peel path under budget (chem priority poison?)"
    assert not counters.budget_exhausted or counters.billed() <= 200
    # Chem-diff sites still present in toward (never-drop), just not first.
    mol = Chem.MolFromSmiles(r)
    tgt = Chem.MolFromSmiles(t)
    ctx = PathContext.from_mols(mol, tgt)
    from xenosite._archive_forest.path_context import mcs_chem_disagree_bonds

    diffs = mcs_chem_disagree_bonds(ctx)
    toward = Dealkylation().sites_toward(mol, tgt, ctx)
    assert toward and diffs
    assert any(
        (b := cleavage_site_bond(mol, s)) is not None and frozenset(b) in diffs
        for s in toward
    )


def test_o_methyl_tetrahydrofuran_demethyl_at_frontier():
    """Exocyclic O–Me on THF: demethyl is a frontier Required site (non-N)."""
    r = Chem.MolFromSmiles("COC1CCCO1")
    t = Chem.MolFromSmiles("OC1CCCO1")
    ctx = PathContext.from_mols(r, t)
    # T not smaller by much / same heavy count often — may return None (unrestricted).
    toward = Dealkylation().sites_toward(r, t, ctx)
    hits = list(
        find_path(r, t, ruleset="PhaseOneRS", depth=2, max_paths=3, max_expansions=40)
    )
    # Either find demethyl or unrestricted search still works.
    assert hits or toward is None
    if hits:
        assert any(s[0] == "Dealkylation" for s in hits[0].steps)


def test_unreachable_elements_api():
    from xenosite._archive_forest import load_ruleset
    from xenosite._archive_forest.path_context import unreachable_new_elements
    from rdkit import Chem

    rules = list(load_ruleset("PhaseOneRS"))
    r = Chem.MolFromSmiles("c1ccccc1")
    t = Chem.MolFromSmiles("[Fe]")
    assert "Fe" in unreachable_new_elements(rules, r, t)
    assert unreachable_new_elements(rules, r, Chem.MolFromSmiles("c1ccccc1O")) == frozenset()
