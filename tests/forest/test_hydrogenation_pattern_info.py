"""Hydrogenation PatternInfo: adds H, dearomatize capability, filter semantics.

Hydrogenation **adds** H (reduction / saturation). Dehydrogenation **removes** H.
Do not conflate them.

Locks data-not-branches: capability on PatternInfo; ``merge_effects`` /
``filter_sites`` read ``adds`` / ``dearomatizes`` against ``atom_diff``.
"""

from __future__ import annotations

from xenosite.forest.find_path import PathCounters, _filters, atom_diff, find_path
from xenosite.forest.rdkitutil import as_mol, copy_mol, wipe_forest
from xenosite.forest.rules import Hydrogenation
from xenosite.forest.rulesets import PhaseOne


def test_hydrogenation_path_end_adds_h_and_dearomatize_capability():
    end = Hydrogenation.endpoints[0][1]
    span = end.get("span") or {}
    assert "H" in (span.get("adds") or "")
    assert not (span.get("removes") or "")
    assert span.get("dearomatizes") is True


def test_hydrogenation_alkene_adds_hh_without_span_dearomatize():
    """Alkene SMARTS must keep span.dearomatizes false so aliphatic C=C→CC
    is not refused by filter_rules when loses_aromaticity is empty."""

    _smirks, info = next(
        item for item in Hydrogenation.smirks if item[1].get("name") == "alkene"
    )
    span = info.get("span") or {}
    assert span.get("adds") == "HH"
    assert span.get("dearomatizes") is False


def test_aromatic_path_hydrogenation_resolves_dearomatizes_true():
    """Path reduction on benzene: adds H and resolved dearomatizes=True."""

    mol = as_mol("c1ccccc1")
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in Hydrogenation().metabolize(mol)
        if "ends" in info
    }
    adds = {
        info["options"].get("adds")
        for _products, info in Hydrogenation().metabolize(mol)
        if "ends" in info
    }
    assert True in flags
    assert any(isinstance(a, str) and "H" in a for a in adds)


def test_aliphatic_path_hydrogenation_resolves_dearomatizes_false():
    mol = as_mol("C=CC=C")
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in Hydrogenation().metabolize(mol)
        if "ends" in info
    }
    assert flags == {False} or flags == set()


def test_adds_h_filter_refuses_when_no_atom_gains_hydrogen():
    """Oxidative MeOPhOH→hydroxyQ: no atom gains H → Hydrogenation skipped.

    Dearomatize-alone would still allow H (target is non-aromatic); reading
    ``adds`` vs on-demand H delta is what blocks reductive saturation.
    """

    from xenosite.forest.find_path import any_h_gain

    reactant = wipe_forest(copy_mol(as_mol("COc1ccc(O)cc1"))).xf.tracing._stamp()
    target = as_mol("O=C1C=C(O)C(=O)C(O)=C1")
    diff = atom_diff(reactant, target)
    assert diff.loses_aromaticity  # dearomatize filter alone would allow
    assert not any_h_gain(diff)
    fr, filter_sites = _filters(diff, True)
    kept = list(
        Hydrogenation().metabolites(
            reactant, filter_rules=fr, filter_sites=filter_sites
        )
    )
    assert kept == []


def test_adds_h_filter_allows_alkene_when_target_gains_hydrogen():
    from xenosite.forest.find_path import any_h_gain

    reactant = wipe_forest(copy_mol(as_mol("C=C"))).xf.tracing._stamp()
    target = as_mol("CC")
    diff = atom_diff(reactant, target)
    assert any_h_gain(diff)
    fr, fs = _filters(diff, True)
    hits = list(Hydrogenation().metabolites(reactant, filter_rules=fr, filter_sites=fs))
    assert hits, "alkene→alkane must keep Hydrogenation under adds-H filter"


def test_meoph_oh_hydroxyq_bill_after_adds_h_filter():
    """Regression: bill ≪ ~900 once reductive H is filtered on oxidative target."""

    counters = PathCounters()
    hits = list(
        find_path(
            "COc1ccc(O)cc1",
            "O=C1C=C(O)C(=O)C(O)=C1",
            ruleset=PhaseOne,
            counters=counters,
            max_paths=1,
            max_nodes=800,
        )
    )
    assert hits
    assert counters.billed < 250, (
        f"expected bill≪900 after adds-H filter; got {counters.billed} "
        f"(nd={counters.nodes} ed={counters.mol_edits})"
    )
