"""PathContext MCS helpers and guided find_path gold / counter tests."""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest import PathSearchCounters, find_path, load_ruleset
from xenosite.forest.path_context import PathContext
from xenosite.forest.rules import (
    Acetylation,
    Dealkylation,
    Hydroxylation,
    NDealkylation,
)
from xenosite.forest.utils import canon_smi, unmapped_smiles


def _smi(s):
    return Chem.MolFromSmiles(s)


def test_path_context_apap_napqi():
    ctx = PathContext.from_mols(
        _smi("CC(=O)Nc1ccc(O)cc1"),
        _smi("CC(=O)N=C1C=CC(=O)C=C1"),
    )
    assert ctx.mcs_smarts
    assert len(ctx.conserved_r_atoms) >= 8
    assert ctx.dearomatization_delta != 0 or ctx.formula_delta.get("H", 0) == 0


def test_path_context_benzene_benzoquinone_aromatic_tolerant():
    ctx = PathContext.from_mols(_smi("c1ccccc1"), _smi("O=C1C=CC(=O)C=C1"))
    assert len(ctx.conserved_r_atoms) == 6
    assert ctx.formula_delta.get("O", 0) == 2
    assert ctx.dearomatization_delta < 0


def test_path_context_etoh_acetaldehyde():
    ctx = PathContext.from_mols(_smi("CCO"), _smi("CC=O"))
    assert len(ctx.conserved_r_atoms) == 3
    assert not ctx.t_only_atoms


def test_path_context_aromatic_vs_kekule_parent():
    arom = _smi("c1ccccc1O")
    kekule = _smi("C1=CC=C(O)C=C1")
    # Product dearomatized quinone-like
    prod = _smi("O=C1C=CC(=O)C=C1")
    ctx_a = PathContext.from_mols(arom, prod)
    ctx_k = PathContext.from_mols(kekule, prod)
    assert len(ctx_a.conserved_r_atoms) >= 5
    assert len(ctx_k.conserved_r_atoms) >= 5


def test_include_exclude_sites_metabolize():
    mol = _smi("CC")
    rule = Hydroxylation()
    all_sites = [
        frozenset(s[1]) for s, _ in rule.metabolize(mol, tag_atoms=False)
    ]
    assert all_sites
    only = all_sites[0]
    filtered = list(
        rule.metabolize(mol, tag_atoms=False, include_sites=[only])
    )
    assert filtered
    assert all(frozenset(s[1]) == only for s, _ in filtered)
    excluded = list(
        rule.metabolize(mol, tag_atoms=False, exclude_sites=[only])
    )
    assert all(frozenset(s[1]) != only for s, _ in excluded)


def test_conjugation_is_terminal_and_redundant_peers():
    a = Acetylation()
    b = Acetylation()
    peers = [a, b]
    assert not a.is_redundant(peers)
    assert b.is_redundant(peers)
    # star product is terminal
    mol = _smi("CCO")
    site, products = next(a.metabolize(mol, tag_atoms=False))
    assert a.is_terminal_product(products[0])


def test_ndealkylation_redundant_when_dealkylation_peer():
    nd = NDealkylation()
    dealk = Dealkylation()
    assert nd.is_redundant([dealk, nd])
    assert not nd.is_redundant([nd])


def test_find_path_etoh_crosscheck_bfs_dfs():
    """Guided gold path must also be reachable by classic BFS and DFS.

    Classic BFS cross-checks stay at depth ≤ 2 (depth ≥ 3 is too expensive).
    """
    r, t = _smi("CCO"), _smi("CC=O")
    classic = load_ruleset("PhaseOneRS")
    bfs_hits = list(classic.find_path(r, t, depth=2, search="bfs"))
    dfs_hits = list(classic.find_path(r, t, depth=2, search="dfs"))
    guided = list(
        find_path(
            r, t, ruleset="PhaseOneRS", depth=2, maybe_prefixes=False, max_paths=3
        )
    )
    assert bfs_hits and dfs_hits and guided


def test_find_path_apap_napqi():
    """APAP→NAPQI via guided (depth 3 OK); classic cross-check at depth 1 only."""
    counters = PathSearchCounters()
    hits = list(
        find_path(
            "CC(=O)Nc1ccc(O)cc1",
            "CC(=O)N=C1C=CC(=O)C=C1",
            ruleset="PhaseOneRS",
            depth=3,
            maybe_prefixes=False,
            expand_phase1_plans=True,
            counters=counters,
            max_paths=2,
        )
    )
    assert hits
    assert canon_smi(hits[0][0][-1]) == canon_smi("CC(=O)N=C1C=CC(=O)C=C1")
    # Direct DH path is depth 1 — never use classic BFS at depth 3 here.
    classic = list(
        load_ruleset("PhaseOneRS").find_path(
            _smi("CC(=O)Nc1ccc(O)cc1"),
            _smi("CC(=O)N=C1C=CC(=O)C=C1"),
            depth=1,
            search="bfs",
        )
    )
    assert classic


def test_enumerate_for_path_prunes_topo_orbits():
    """Benzene hydroxylation: one topo orbit (metabolize distinct sites), not six."""
    from xenosite.forest.path_context import PathContext
    from xenosite.forest.rules import Hydroxylation

    mol = _smi("c1ccccc1")
    ctx = PathContext.from_mols(mol, _smi("O=C1C=CC(=O)C=C1"))
    items = list(
        Hydroxylation().enumerate_for_path(mol, ctx, expand_phase1_plans=True)
    )
    assert len(items) == 1


def test_find_path_benzene_bq_qf_vs_phaseone():
    """With QF: one expansion. PhaseOneRS (OH→OH→DH) stays a small handful."""
    qf_c = PathSearchCounters()
    qf_hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset="QF",
            depth=3,
            maybe_prefixes=False,
            counters=qf_c,
            max_paths=1,
        )
    )
    assert qf_hits
    assert qf_c.rule_expansions == 1

    po_c = PathSearchCounters()
    po_hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset="PhaseOneRS",
            depth=3,
            maybe_prefixes=False,
            counters=po_c,
            max_paths=1,
        )
    )
    assert po_hits
    assert po_c.rule_expansions > qf_c.rule_expansions
    assert po_c.rule_expansions <= 30  # formula/MCS could_help; not hundreds


# Terbinafine → TBF-A (6,6-dimethylhept-2-en-4-ynal). Stereo not required —
# guided matching strips double-bond stereo.
TERBINAFINE = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
TBF_A = "CC(C)(C)C#CC=CC=O"


def test_terbinafine_path_context_prunes_naphthalene_side():
    ctx = PathContext.from_mols(_smi(TERBINAFINE), _smi(TBF_A))
    assert ctx.r_only_atoms  # naphthalene retained on R only
    assert ctx.conserved_r_atoms  # enyne chain conserved


def test_find_path_default_phaseone_qf_no_ruleset_pick():
    """Default PhaseOneQF: quinones stay cheap; TBA still found via Dealkylation."""
    bq = PathSearchCounters()
    bq_hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            depth=3,
            max_paths=1,
            max_expansions=20,
            counters=bq,
        )
    )
    assert bq_hits
    assert bq.rule_expansions == 1  # QF first, not Phase I OH orbit search

    tba = PathSearchCounters()
    tba_hits = list(
        find_path(
            TERBINAFINE,
            TBF_A,
            depth=2,
            max_paths=1,
            max_expansions=40,
            counters=tba,
        )
    )
    assert tba_hits
    outcome = tba_hits[0]
    assert [s[0] for s in outcome.steps] == ["Dealkylation"]
    assert outcome.maybe  # CleavageSide bag for naphthyl-amine leaving group
    # Preceding demethyl (same N as TBA-forming dealk) passes Maybe.
    assert outcome.allows("NDealkylation", frozenset({0, 1}))
    assert tba.rule_expansions <= 10
    assert tba.billed() <= 40
    assert not tba.budget_exhausted
