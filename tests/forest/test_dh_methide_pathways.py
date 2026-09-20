"""Methide is always in rule data; refuse via filter_sites (no pathways= opt-in).

Adapted from ``tests/test_dh_methide_pathways.py``. Forest gated methide with
``pathways=("methide",)`` — approved drop (docs/forest/DROPPED.md). Forest offers methide on
Dehydrogenation / QuinoneFormation PatternInfo. At most one methide end is
resolved (both-methide pairs are not built). Callers who want none pass
``filter_sites`` that reads ``options["methide"]``.
"""

from __future__ import annotations

from xenosite.forest.find_path import PathCounters, find_path
from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import Dehydrogenation, Hydroxylation, QuinoneFormation
from xenosite.forest.rulesets import PhaseOne, RuleSet

from .helpers import canon, product_smiles

_O_QM = "C=C1C=CC=CC1=O"


def _no_methide(mol, site, info) -> bool:
    return not (info.get("options") or {}).get("methide")


def test_methide_is_always_in_rule_data():
    """No pathways= flag: PatternInfo possibilities declare methide."""

    dh_names = {
        p.get("name")
        for _smarts, p in Dehydrogenation.endpoints
        if any(poss.get("methide") for poss in (p.get("possibilities") or ()))
    }
    assert "methide_end" in dh_names

    qf_methide = any(
        poss.get("methide")
        for _smarts, p in QuinoneFormation.endpoints
        for poss in (p.get("possibilities") or ())
    )
    assert qf_methide


def test_dh_methide_emits_o_quinone_methide_from_o_cresol():
    assert canon(_O_QM) in product_smiles(Dehydrogenation(), "Oc1ccccc1C")


def test_filter_sites_refuses_methide():
    """Caller who wants no methides filters the effect — not a pathways= branch."""

    assert canon(_O_QM) not in {
        p.xf.csmi
        for products, _info in Dehydrogenation().metabolize(
            MolFromSmiles("Oc1ccccc1C"), filter_sites=_no_methide
        )
        for p in products
    }
    assert canon(_O_QM) not in {
        p.xf.csmi
        for products, _info in QuinoneFormation().metabolize(
            MolFromSmiles("Cc1ccccc1"), filter_sites=_no_methide
        )
        for p in products
    }


def test_at_most_one_methide_end():
    """Xylene: no bis-methide pair; every methide product has exactly one end."""

    assert product_smiles(Dehydrogenation(), "Cc1ccccc1C") == set()

    for products, info in QuinoneFormation().metabolize(MolFromSmiles("Cc1ccccc1C")):
        opts = info.get("options") or {}
        ends = info.get("ends") or ()
        n_methide = sum(1 for end in ends if end.get("methide"))
        assert n_methide <= 1, info["csmi"]
        if opts.get("methide"):
            assert n_methide == 1, info["csmi"]
        for product in products:
            assert product.xf.csmi != canon("C=c1ccccc1=C")


def test_find_path_o_cresol_to_o_qm_via_dehydrogenation():
    counters = PathCounters()
    hits = list(
        find_path(
            "Oc1ccccc1C",
            _O_QM,
            ruleset=RuleSet((Dehydrogenation,), name="DH"),
            max_paths=1,
            max_nodes=40,
            counters=counters,
        )
    )
    assert hits
    assert hits[0].smiles == canon(_O_QM)
    assert counters.nodes <= 40


def test_find_path_toluene_to_o_qm():
    """Toluene → o-QM (PhaseOne / QF). Forest used OH+DH with methide opt-in.

    Atom-diff oxygen filters block the OH→DH walk on toluene when filters are
    on; PhaseOne still finds via QuinoneFormation. Direct OH+DH works with
    ``use_filters=False`` (chemistry present; heuristic gap is separate).
    """

    counters = PathCounters()
    hits = list(
        find_path(
            "Cc1ccccc1",
            _O_QM,
            ruleset=PhaseOne,
            max_paths=1,
            max_nodes=80,
            counters=counters,
        )
    )
    assert hits
    assert hits[0].smiles == canon(_O_QM)
    assert counters.nodes <= 80

    qf = list(
        find_path(
            "Cc1ccccc1",
            _O_QM,
            ruleset=RuleSet((QuinoneFormation,), name="QF"),
            max_paths=1,
            max_nodes=40,
        )
    )
    assert qf and qf[0].smiles == canon(_O_QM)

    # Chemistry of the forest gold walk without atom-diff filters.
    oh_dh = list(
        find_path(
            "Cc1ccccc1",
            _O_QM,
            ruleset=RuleSet((Hydroxylation, Dehydrogenation), name="t"),
            max_paths=1,
            max_nodes=40,
            use_filters=False,
        )
    )
    assert oh_dh
    assert [s.rule for s in oh_dh[0].plan.children] == [
        "Hydroxylation",
        "Dehydrogenation",
    ]
