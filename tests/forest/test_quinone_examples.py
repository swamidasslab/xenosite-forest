"""Port of ``tests/test_quinone.py``: each quinone example must appear as a product.

SMILES stay in ``test_quinone.examples``. The fused oxazine is a shared miss
with forest, asserted as absence rather than an xfail.
"""

from __future__ import annotations

import pytest
from test_quinone import examples

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import QuinoneFormation

from .helpers import emits_product

_CASES = [(name, reactant, product, site) for name, reactant, product, site in examples]

# Forest also does not open this fused N-alkyl oxazine.
_FOREST_MISS = (
    "CN1CCOc2c1cc(cc2)Nc1ncc2n(n1)c(cc2)c1cccc(c1)S(=O)(=O)NC(C)(C)C"
)


@pytest.mark.parametrize(
    "name, reactant, product, site",
    _CASES,
    ids=[name for name, *_rest in _CASES],
)
def test_quinone_emits_historical_product(name, reactant, product, site):
    if reactant == _FOREST_MISS:
        # Forest also misses this fused N-alkyl oxazine.
        assert not emits_product(QuinoneFormation(), reactant, product)
        return

    if emits_product(QuinoneFormation(), reactant, product):
        return

    assert False, f"Failed to find {product} in {reactant} ({name})"


def test_nonaromatic_conjugate_is_not_dearomatizing():
    """A ring enol pair is the same edit. No aromatic atom, so not a quinone."""

    mol = MolFromSmiles("OC1=CCCC=C1O")
    assert mol is not None
    flags = {
        info["options"].get("dearomatizes")
        for _product, info in QuinoneFormation().metabolize(mol)
    }
    assert flags == {False}


def test_aromatic_conjugate_is_dearomatizing():
    mol = MolFromSmiles("Oc1ccc(O)cc1")
    assert mol is not None
    flags = {
        info["options"].get("dearomatizes")
        for _product, info in QuinoneFormation().metabolize(mol)
    }
    assert True in flags
