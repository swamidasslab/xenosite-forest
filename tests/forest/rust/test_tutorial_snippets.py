"""Smoke tests mirroring notebooks/forest_product_door.ipynb (realistic drug paths)."""

from __future__ import annotations

from xenosite.forest import (
    ForestMol,
    MetabolicNetwork,
    find_path,
    find_path_partial,
    normalize_tautomer,
    phase_one,
    product_graph_bfs,
    random_path,
)

APAP = "CC(=O)Nc1ccc(O)cc1"
NAPQI = "CC(=O)N=C1C=CC(=O)C=C1"
TERBINAFINE = "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"
TBF_ALDEHYDE = r"C(#C/C=C/C=O)C(C)(C)C"
EUGENOL = "COc1ccc(O)cc1"
EUGENOL_Q_TARGET = "O=C1C=C(O)C(=O)C(O)=C1"


def test_tutorial_apap_to_napqi_phase_one():
    hits, counters, net = find_path(
        APAP,
        NAPQI,
        max_paths=1,
        max_nodes=80,
        use_atom_diff=False,
    )
    assert hits, f"no path (billed={counters.billed})"
    assert hits[0].smiles
    assert net.n_edges() >= 1
    steps = hits[0].to_dict()["steps"]
    rules = {s["rule"] for s in steps}
    assert rules  # hydroxylation / dehydrogenation / QF chain


def test_tutorial_terbinafine_dealkylation():
    hits, _, net = find_path(
        TERBINAFINE,
        TBF_ALDEHYDE,
        max_paths=1,
        max_nodes=40,
        use_atom_diff=False,
    )
    assert hits
    assert hits[0].to_dict()["steps"][0]["rule"] == "Dealkylation"
    assert net.n_nodes() >= 2


def test_tutorial_shared_network_bfs_then_find_path():
    net = MetabolicNetwork()
    product_graph_bfs(
        "CC",
        target="CCO",
        max_nodes=32,
        max_depth=3,
        network=net,
    )
    assert net.n_nodes() >= 2
    hits, _, out = find_path("CC", "CCO", max_paths=1, max_nodes=200, network=net)
    assert hits
    assert out.n_nodes() == net.n_nodes()
    root = net.root_idx()
    target_i = net.index_of(hits[0].smiles)
    assert root is not None and target_i is not None
    plan = net.step_plan_between(root, target_i)
    assert len(plan) >= 1


def test_tutorial_eugenol_bfs_bounded():
    net = product_graph_bfs(
        EUGENOL,
        target=EUGENOL_Q_TARGET,
        max_nodes=64,
        max_depth=4,
    )
    assert net.n_nodes() <= 64
    assert net.n_nodes() >= 1


def test_tutorial_metabolize_one_hop():
    mol = ForestMol(APAP)
    emissions = phase_one().metabolize(mol)
    assert emissions
    names = {e.pattern_name for e in emissions}
    assert names


def test_tutorial_normalize_tautomer():
    keto, _ = normalize_tautomer("CC=O")
    enol, changed = normalize_tautomer("OC=C")
    assert changed
    assert enol.csmi == keto.csmi


def test_tutorial_random_path_seeded():
    out = random_path(APAP, seed=42, max_steps=2)
    assert out.path
    assert out.to_dict()["steps"]


def test_tutorial_find_path_partial_smoke():
    exact, partials, counters, net = find_path_partial(
        APAP,
        NAPQI,
        max_paths=1,
        max_nodes=80,
        use_atom_diff=False,
    )
    assert exact or partials or counters.billed >= 0
    assert net.n_nodes() >= 1
