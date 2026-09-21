"""Quinone goldens the fragment cleanup used to drop.

The charged fused product is not one of them. Two saturated goldens from the
downstream report are not products of this rule in the archive either, so
they are not locked here.
"""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.rules import QuinoneFormation

_PARENT = "CS(=O)(=O)c1ccc(-c2cn3ccccc3n2)cc1"
_GOLDENS = (
    "CS(=O)(=O)c1ccc(C2=CN=CC(=O)C=CC=N2)cc1",
    "CS(=O)(=O)c1ccc(C2=CN=CC=CC(=O)C=N2)cc1",
    "CS(=O)(=O)c1ccc(C2=NC=CC=CC(=O)N=C2)cc1",
)
_CHARGED = "CS(=O)(=O)c1ccc(C2=Nc3cccc[n+]3C2=O)cc1"


def _products(smiles: str) -> set[str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    found: set[str] = set()
    for batch in QuinoneFormation().metabolites(mol):
        for product in batch.products:
            fragments = Chem.GetMolFrags(product, asMols=True)
            largest = max(fragments, key=lambda item: item.GetNumAtoms())
            found.add(Chem.MolToSmiles(largest, canonical=True, isomericSmiles=False))
    return found


def test_imidazole_pyridine_ring_opened_quinones() -> None:
    got = _products(_PARENT)
    for golden in _GOLDENS:
        assert golden in got
    assert _CHARGED not in got
