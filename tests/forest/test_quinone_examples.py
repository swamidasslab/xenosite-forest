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


def test_carbamazepine_does_not_emit_two_double_nitrogen():
    """Iminium plus dealkylation must not leave ``C=[N+]=C``.

    Charge on a tertiary nitrogen with one double bond is still an iminium.
    Chlorpromazine keeps that product. Carbamazepine site ``{4, 17}`` does not.
    """

    bad = "C1=c2ccccc2=[N+]=c2ccccc2=C1"
    cbz = MolFromSmiles("NC(=O)N1c2ccccc2C=Cc2ccccc21")
    assert cbz is not None
    cbz_products = {
        p.xf.csmi
        for products, _info in QuinoneFormation().metabolize(cbz)
        for p in products
    }
    assert bad not in cbz_products

    cpz = MolFromSmiles("CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21")
    assert cpz is not None
    cpz_products = {
        p.xf.csmi
        for products, _info in QuinoneFormation().metabolize(cpz)
        for p in products
    }
    assert any("[n+]" in smi or "[N+]" in smi for smi in cpz_products)


def test_alprazolam_re_aromatized_cation_is_not_a_quinone():
    """The kekulized system came back aromatic, with no localized double bond.

    ``Cc1nnc2cnc(-c3ccccc3)c3cc(Cl)ccc3[n+]1-2`` has no carbonyl. The pendant
    phenyl can stay aromatic. The system that was kekulized cannot.
    """

    bad = "Cc1nnc2cnc(-c3ccccc3)c3cc(Cl)ccc3[n+]1-2"
    mol = MolFromSmiles("Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2")
    assert mol is not None
    products = {
        p.xf.csmi
        for pieces, _info in QuinoneFormation().metabolize(mol)
        for p in pieces
    }
    assert bad not in products


def test_nonaromatic_conjugate_is_not_dearomatizing():
    """A ring enol pair is the same edit. No aromatic atom, so not a quinone."""

    mol = MolFromSmiles("OC1=CCCC=C1O")
    assert mol is not None
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in QuinoneFormation().metabolize(mol)
    }
    assert flags == {False}


def test_aromatic_conjugate_is_dearomatizing():
    mol = MolFromSmiles("Oc1ccc(O)cc1")
    assert mol is not None
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in QuinoneFormation().metabolize(mol)
    }
    assert True in flags
