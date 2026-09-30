"""Notebook display helpers for forest pyclasses."""

from __future__ import annotations

from xenosite.forest import (
    ForestMol,
    find_path,
    find_path_partial,
    hydroxylation,
    phase_one,
    product_graph_bfs,
)


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
    hits, counters, net = find_path("CC", "CCO", max_paths=1, max_nodes=64)
    assert hits
    assert net.n_edges() >= 1
    hit = hits[0]
    html = hit._repr_html_()
    assert "PathOutcome" in html
    assert "hop" in html or "start" in html
    assert "SOM" in html or "@" in html
    hops = hit.hops()
    assert hops
    assert "reactant" in hops[0]
    assert hops[0]["reactant"]
    # Site indexes the reactant, not only the final product.
    assert hops[0]["site"] is not None
    mol = hit.mol
    assert "survivors=" in str(mol)
    assert "PathCounters" in str(counters)
    assert "billed=" in str(counters)
    assert "PathCounters" in counters._repr_html_()


def test_repr_html_ruleset_and_emission_uses_product_trace():
    rs = hydroxylation()
    assert "pre" in rs._repr_html_()
    mol = ForestMol("c1ccccc1")
    emissions = rs.metabolize(mol)
    assert emissions
    e = emissions[0]
    ehtml = e._repr_html_()
    assert "Emission" in ehtml
    assert "reactant" in ehtml.lower() or "SOM" in ehtml
    # Site is on the substrate, not mis-applied only to products.
    assert e.reactant.csmi == mol.csmi
    assert e.site == int(e.site)


def test_forest_mol_wrapper_repr_html():
    mol = ForestMol("CCO")
    html = mol._repr_html_()
    assert "ForestMol" in html
    assert "trace" in html


def test_graph_node_and_edge_html_use_mol_trace():
    net = product_graph_bfs("CC", target="CCO", max_nodes=32, max_depth=3)
    assert net.n_nodes() >= 1
    node = net[0]
    nhtml = node._repr_html_()
    assert "GraphNode" in nhtml
    assert "trace" in nhtml
    assert "GraphNode" in str(node)
    # Find a child with an inbound edge.
    for i in range(net.n_nodes()):
        n = net[i]
        if n.n_inbound() > 0:
            edge = n.inbound_edge(0)
            ehtml = edge._repr_html_()
            assert "GraphEdge" in ehtml
            assert "SOM" in ehtml or "reactant" in ehtml or "start" in ehtml
            # Site marks parent, not only kept child.
            assert edge.parent_mol.csmi
            assert edge.site == int(edge.site)
            break
    else:
        raise AssertionError("expected at least one inbound edge in BFS graph")


def test_partial_outcome_html_traces_mol():
    exact, partials, _counters, net = find_path_partial(
        "c1ccccc1", "CCO", max_paths=1, max_nodes=48
    )
    assert net.n_nodes() >= 1
    # Prefer a partial when present; otherwise a hit still exercises plan HTML.
    if partials:
        part = partials[0]
        html = part._repr_html_()
        assert "PartialOutcome" in html
        assert "trace" in html or "hop" in html or "start" in html
        assert "residual_cost" in html
        assert part.residual_cost >= 0
    else:
        assert exact
        assert "PathOutcome" in exact[0]._repr_html_()


def test_step_plan_linearizations_api():
    hits, _, _net = find_path("CC", "CCO", max_paths=1, max_nodes=64)
    plan = hits[0].plan
    n = plan.n_linearizations()
    lins = plan.linearizations()
    assert len(lins) == min(n, 64)
    assert n >= 1
    text = str(plan)
    assert "StepPlan" in text
    assert "linearization" in text
