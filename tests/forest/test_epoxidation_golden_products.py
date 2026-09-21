"""Neutral epoxides that the cached kekulé parent dropped or charged.

These are structure failures against the archive, not atom-index shifts.
Progressions and index-only rows are not locked here.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.rules import Epoxidation


def _products(smiles: str) -> set[str]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    found: set[str] = set()
    for batch in Epoxidation().metabolites(mol):
        for product in batch.products:
            fragments = Chem.GetMolFrags(product, asMols=True)
            largest = max(fragments, key=lambda item: item.GetNumAtoms())
            found.add(Chem.MolToSmiles(largest, canonical=True, isomericSmiles=False))
    return found


@pytest.mark.parametrize(
    "smiles,golden",
    [
        ("CSc1ccc2sc(N)nc2c1", "CSC1=CC2OC23SC(N)=NC3=C1"),
        ("O=S(=O)(O)Oc1ccc(O)cc1", "O=S(=O)(O)OC1=CC=C(O)C2OC12"),
        ("CC(C)(C)c1ccc(S(N)(=O)=O)cc1", "CC(C)(C)C1=CC=C(S(N)(=O)=O)C2OC12"),
        (
            "NC(Cc1ccc([N+](=O)[O-])cc1)C(=O)O",
            "NC(CC1=CC=C([N+](=O)[O-])C2OC12)C(=O)O",
        ),
        ("Cc1ncc([N+](=O)[O-])n1CCO", "CC1=NC2OC2([N+](=O)[O-])N1CCO"),
        ("Cc1ncc([N+](=O)[O-])n1CCO", "CC12ON1C=C([N+](=O)[O-])N2CCO"),
        ("CC(=O)c1ccc2oc(=O)ccc2c1", "CC(=O)c1ccc2c(c1)C1OC1C(=O)O2"),
    ],
)
def test_epoxidation_keeps_neutral_golden(smiles: str, golden: str) -> None:
    assert golden in _products(smiles)
