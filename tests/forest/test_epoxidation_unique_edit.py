"""Epoxidation unique-edit collapses undirected bond embeddings.

Phenol (and peers) used to emit symmetry-related aromatic bonds as distinct
sites because directed MapRankKey ``((mapno, rank), …)`` distinguished map
1↔2 placement while the epoxide product is the same. CSMI then dropped the
duplicates. ``site_kind="bond"`` keys unique-edit by sorted site ranks
(undirected bond ends) so those embeddings merge before CSMI.
"""

from __future__ import annotations

import warnings

import pytest
from rdkit import Chem

from xenosite.forest.rules import CsmiDedupWarning, Epoxidation


@pytest.mark.parametrize(
    "smiles,n_keep",
    [
        ("Oc1ccccc1", 3),  # ipso-ortho, ortho-meta, meta-para
        ("C=Cc1ccccc1", 4),  # vinyl + three aromatic classes
        ("c1ccccc1", 1),
        ("Oc1ccc(O)cc1", 2),
    ],
)
def test_epoxidation_unique_edit_no_csmi(smiles: str, n_keep: int) -> None:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", CsmiDedupWarning)
        products = list(Epoxidation().metabolize(Chem.MolFromSmiles(smiles)))
    csmi = [w for w in caught if issubclass(w.category, CsmiDedupWarning)]
    assert not csmi, f"{smiles}: unexpected CSMI drops {[str(w.message) for w in csmi]}"
    assert len(products) == n_keep
    for _product, info in products:
        assert len(info["site"]) == 2, info["site"]


def test_epoxidation_declares_bond_site() -> None:
    rule = Epoxidation()
    assert rule.site_kind == "bond"
    _smarts, pattern = rule.smarts[0]
    assert pattern.get("site_map") == (1, 2)
