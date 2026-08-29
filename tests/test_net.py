import pytest

nx = pytest.importorskip("networkx")

from xenosite.forest.net import MetaboliteNetwork, metabolites
from xenosite.forest.phaseone import PhaseOneRS


def test_expand_and_paths():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)
    assert "CCO" in net
    paths = list(net.paths("CCO", "C=CO", extra_depth=1))
    assert paths
    assert paths[0][0] == "CCO"
    assert paths[0][-1] == "C=CO"


def test_rxn_map_pairs_smiles_output_order():
    rxns = list(metabolites(PhaseOneRS, "CCO"))
    assert rxns
    for rxn in rxns:
        frm, to = rxn.map
        assert len(frm) == len(to)
        assert len(frm) >= 1
        assert list(frm) == sorted(frm)
