"""Hydrogenation SMARTS unique-edit uses undirected bond ranks.

``site_kind="atom_pair"`` is for ResonancePair path ends. One-bond alkene /
alkyne SMARTS still go through ``site_signature``; directed MapRankKey used to
split symmetry-related aromatic bonds (phenol ortho/meta) that yield the same
dihydro product, and CSMI then warned. Undirected ``bond_rank_key`` for
``atom_pair`` SMARTS collapses those embeddings before CSMI.
"""

from __future__ import annotations

import warnings

import pytest
from rdkit import Chem

from xenosite.forest.rules import Hydrogenation, SiteDeduplicationWarning


@pytest.mark.parametrize(
    "smiles",
    [
        "Oc1ccccc1",
        "COc1ccccc1",
        "Clc1ccccc1",
        "Nc1ccccc1",
        "Oc1ccc(O)cc1",
        "C=CC=C",
        "CCc1ccccc1",
        "CC(=O)c1ccccc1",
    ],
)
def test_hydrogenation_unique_edit_no_csmi(smiles: str) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", SiteDeduplicationWarning)
        list(Hydrogenation().metabolize(Chem.MolFromSmiles(smiles)))
    csmi = [w for w in caught if issubclass(w.category, SiteDeduplicationWarning)]
    assert not csmi, f"{smiles}: {[str(w.message) for w in csmi]}"


def test_hydrogenation_declares_atom_pair() -> None:
    assert Hydrogenation().site_kind == "atom_pair"
