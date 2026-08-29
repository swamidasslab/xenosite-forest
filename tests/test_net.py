import pytest

nx = pytest.importorskip("networkx")

from xenosite.metabolite.net import MetaboliteNetwork
from xenosite.metabolite.phaseone import PhaseOneRS


def test_expand_and_paths():
    net = MetaboliteNetwork("CCO")
    net.expand(PhaseOneRS)
    assert "CCO" in net
    paths = list(net.paths("CCO", "C=CO", extra_depth=1))
    assert paths
    assert paths[0][0] == "CCO"
    assert paths[0][-1] == "C=CO"
