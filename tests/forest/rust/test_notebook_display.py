"""Notebook display helpers for forest pyclasses."""

from __future__ import annotations

from xenosite.forest import ForestMol, find_path, hydroxylation, phase_one, product_graph_bfs


def test_ruleset_str_is_indented_hierarchy():
    text = str(phase_one())
    assert "\n" in text
    assert "Hydroxylation" in text or "members" in text


def test_forest_mol_str_uses_trace_not_steps():
    mol = ForestMol("CCO")
    text = str(mol)
    assert "ForestMol" in text
    assert "trace" in text
    assert "survivors=" in text
    assert "born=" in text
    assert "tags" in text
    stamp_end, survivors, born, untagged = mol.trace_counts()
    assert stamp_end >= 1
    assert survivors + born + untagged == len(mol.atom_tags())
    assert len(mol.stamp_origins()) == 3


def test_network_str_summary():
    net = product_graph_bfs("CC", target="CCO", max_nodes=32, max_depth=3)
    text = str(net)
    assert "MetabolicNetwork" in text
    assert "root:" in text


def test_path_outcome_html_prefers_trace_and_plan():
    hits, _ = find_path("CC", "CCO", max_paths=1, max_nodes=64)
    assert hits
    hit = hits[0]
    html = hit._repr_html_()
    assert "PathOutcome" in html
    assert "trace" in html
    assert "StepPlan" in html
    # No per-hop intermediate dump.
    assert "hop 0:" not in html
    mol = hit.mol
    assert "survivors=" in str(mol)


def test_repr_html_ruleset_and_emission_uses_product_trace():
    rs = hydroxylation()
    assert "pre" in rs._repr_html_()
    mol = ForestMol("c1ccccc1")
    emissions = rs.metabolize(mol)
    assert emissions
    ehtml = emissions[0]._repr_html_()
    assert "Emission" in ehtml
    assert "trace" in ehtml or "product" in ehtml


def test_forest_mol_wrapper_repr_html():
    mol = ForestMol("CCO")
    html = mol._repr_html_()
    assert "ForestMol" in html
    assert "trace" in html
