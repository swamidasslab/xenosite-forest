"""GraphNode / GraphEdge views on MetabolicNetwork."""

from __future__ import annotations

from xenosite.forest import product_graph_bfs


def test_network_indexing_and_inbound_edges():
    net = product_graph_bfs("CC", target="CCO", max_nodes=32, max_depth=3)
    assert len(net) >= 2
    root = net[0]
    assert root.index == 0
    assert root.csmi
    assert root.mol.csmi == root.csmi
    assert "expanded" in root or root.expanded in (True, False)
    root["note"] = "root"
    assert root["note"] == "root"
    assert root.get_attr("note") == "root"
    d = root.to_dict()
    assert d["index"] == 0
    assert d["attrs"]["note"] == "root"
    ethanol_i = net.index_of(
        next(c for c in net.children(root.csmi) if c != root.csmi)
    )
    assert ethanol_i is not None
    child = net[ethanol_i]
    assert child.n_inbound() >= 1
    edge = child.inbound_edge(0)
    assert edge.rule == "Hydroxylation"
    assert edge.parent_index == 0
    assert edge.kept_mol.csmi == child.csmi
    assert "rule" in edge
    edge["custom"] = 1
    assert edge["custom"] == 1
    assert edge.to_dict()["attrs"]["custom"] == 1
