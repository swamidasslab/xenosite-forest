"""Acetylation: product `#` apply + aliphatic/aromatic heteroatom branches.

Rust expands product ``[#6](=[#8])[#6]`` to organic aliphatic then aromatic
spellings (chematic bracket expand misses). Catalog keeps ``#``. Reactant
``#`` must hit aliphatic and aromatic heteroatoms (pyrrole ``n``).
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
        ("CCO", "CC(=O)OCC"),
        ("CCN", "CCNC(C)=O"),
        ("CS", "CSC(C)=O"),
        ("Oc1ccccc1", "CC(=O)Oc1ccccc1"),
        ("Nc1ccccc1", "CC(=O)Nc1ccccc1"),
        ("Sc1ccccc1", "CC(=O)Sc1ccccc1"),
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


def test_acetylation_catalog_keeps_atomic_product_smarts() -> None:
    """Product side stays ``#``; apply-layer expand covers chematic gap."""

    text, _info = Acetylation.smirks[0]
    assert "[#6](=[#8])[#6]" in str(text)
    assert "C(=O)C" not in str(text)


def test_acetylation_example_substrates_include_both_attachment_classes() -> None:
    examples = Acetylation._example_substrates
    assert "CCO" in examples
    assert any("c" in s and ("N" in s or "n" in s or "O" in s) for s in examples)
