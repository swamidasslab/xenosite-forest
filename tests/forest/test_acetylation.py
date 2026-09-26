"""Acetylation organic-subset product covers aliphatic and aromatic heteroatoms.

Chematic expand of product ``[#6](=[#8])[#6]`` misses (parity #5). Workaround
``C(=O)C`` is aliphatic-only on the *product*; reactant SMARTS keep ``#`` so
both aliphatic (``O``/``N``/``S``) and aromatic (``n``) specialize branches
must still acetylate.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from xenosite.forest.rules import Acetylation


def _canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


@pytest.mark.parametrize(
    "substrate,want",
    [
        # aliphatic heteroatom
        ("CCO", "CC(=O)OCC"),
        ("CCN", "CCNC(C)=O"),
        ("CS", "CSC(C)=O"),
        # aryl-attached aliphatic heteroatom
        ("Oc1ccccc1", "CC(=O)Oc1ccccc1"),
        ("Nc1ccccc1", "CC(=O)Nc1ccccc1"),
        ("Sc1ccccc1", "CC(=O)Sc1ccccc1"),
        # aromatic heteroatom (pyrrole NH → specialize ``n``)
        ("[nH]1cccc1", "CC(=O)n1cccc1"),
    ],
)
def test_acetylation_covers_aliphatic_and_aromatic_heteroatom_branches(
    substrate: str, want: str
) -> None:
    mol = Chem.MolFromSmiles(substrate)
    assert mol is not None
    products = {
        Chem.MolToSmiles(p)
        for products, _info in Acetylation(as_star=False).metabolize(mol)
        for p in products
    }
    assert _canon(want) in products, products


def test_acetylation_example_substrates_include_both_attachment_classes() -> None:
    examples = Acetylation._example_substrates
    assert "CCO" in examples  # aliphatic OH
    assert any("c" in s and ("N" in s or "n" in s or "O" in s) for s in examples)
