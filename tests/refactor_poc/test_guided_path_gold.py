"""Guided find_path gold cases only (no PathContext MCS / forest internals).

Port of the metabolite golds in ``tests/test_guided_path.py``. Skip
PathContext helpers, dearomatization site-rank filters, and
``enumerate_for_path``. Peer ``is_redundant`` is forest-only (not in poc);
assert terminal conjugations instead.
"""

from __future__ import annotations

from xenosite.refactor_poc.find_path import PathCounters, find_path
from xenosite.refactor_poc.rdkit_api import MolFromSmiles
from xenosite.refactor_poc.rules import (
    Acetylation,
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
)
from xenosite.refactor_poc.rulesets import PhaseOne, RuleSet

from .helpers import canon

TERBINAFINE = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
TBF_A = "CC(C)(C)C#CC=CC=O"


def test_find_path_apap_napqi():
    counters = PathCounters()
    hits = list(
        find_path(
            "CC(=O)Nc1ccc(O)cc1",
            "CC(=O)N=C1C=CC(=O)C=C1",
            ruleset=PhaseOne,
            max_paths=2,
            max_nodes=80,
            counters=counters,
        )
    )
    assert hits
    assert hits[0].smiles == canon("CC(=O)N=C1C=CC(=O)C=C1")
    assert counters.nodes <= 80


def test_find_path_ring_end_sites_are_expanded():
    """Hydroquinone and APAP are found at the OH / NH ends, inside budget."""

    qf = RuleSet(
        (QuinoneFormation, Hydroxylation, Dehydrogenation, Dealkylation),
        name="ringEndQf",
    )
    counters = PathCounters()
    hits = list(
        find_path(
            "Oc1ccc(O)cc1",
            "O=C1C=CC(=O)C=C1",
            ruleset=qf,
            max_paths=1,
            max_nodes=120,
            counters=counters,
        )
    )
    assert hits, (
        f"billed={counters.billed} edits={counters.mol_edits} nodes={counters.nodes}"
    )
    assert counters.sites_considered >= 1
    assert counters.nodes <= 120
    assert hits[0].smiles == canon("O=C1C=CC(=O)C=C1")

    phase1 = RuleSet((Hydroxylation, Dehydrogenation), name="ringEndP1")
    counters = PathCounters()
    hits = list(
        find_path(
            "CC(=O)Nc1ccc(O)cc1",
            "CC(=O)N=C1C=CC(=O)C=C1",
            ruleset=phase1,
            max_paths=1,
            max_nodes=40,
            counters=counters,
        )
    )
    assert hits, (
        f"billed={counters.billed} edits={counters.mol_edits} nodes={counters.nodes}"
    )
    assert counters.sites_considered >= 1
    assert hits[0].smiles == canon("CC(=O)N=C1C=CC(=O)C=C1")


def test_find_path_benzene_bq_qf_vs_phaseone():
    """QF is one hop. PhaseOne includes QF, so both stay cheap (not hundreds)."""

    qf_c = PathCounters()
    qf_hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset=RuleSet((QuinoneFormation,), name="QF"),
            max_paths=1,
            max_nodes=40,
            counters=qf_c,
        )
    )
    assert qf_hits
    assert qf_c.mol_edits == 1

    po_c = PathCounters()
    po_hits = list(
        find_path(
            "c1ccccc1",
            "O=C1C=CC(=O)C=C1",
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=80,
            counters=po_c,
        )
    )
    assert po_hits
    # Forest PhaseOneRS lacked QF and spent more OH orbits; poc PhaseOne
    # includes QuinoneFormation, so billed stays a small handful.
    assert po_c.mol_edits >= qf_c.mol_edits
    assert po_c.billed <= 30
    assert po_c.nodes <= 80


def test_conjugation_is_terminal():
    """Forest peer ``is_redundant`` is not in poc; terminal conjugates still stop."""

    a = Acetylation()
    mol = MolFromSmiles("CCO")
    product, _info = next(a.metabolize(mol))
    assert a.is_terminal_product(product)
    assert product.xf.is_terminal
    assert list(a.metabolize(product)) == []


def test_find_path_default_phaseone_tba():
    """PhaseOne: TBA via Dealkylation with CleavageSide bag; demethyl passes Maybe."""

    tba = PathCounters()
    tba_hits = list(
        find_path(
            TERBINAFINE,
            TBF_A,
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=40,
            counters=tba,
        )
    )
    assert tba_hits
    outcome = tba_hits[0]
    assert [s.rule for s in outcome.plan.children] == ["Dealkylation"]
    assert outcome.maybe
    assert outcome.allows("NDealkylation", frozenset({0, 1}))
    assert tba.nodes <= 40
    assert tba.billed <= 40
