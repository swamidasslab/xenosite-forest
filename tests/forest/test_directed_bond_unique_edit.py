"""Directed-bond unique-edit keeps chemically distinct orientations.

``site_kind="directed_bond"`` keys unique-edit by ordered MapRankKey
(map 1 = oxygenated carbon). Undirected ``bond_rank_key`` would merge
anisole ring-open regioisomers that yield different products.
"""

from __future__ import annotations

import warnings
from collections import defaultdict

from rdkit import Chem

from xenosite.forest.graph_isomorphism import bond_rank_key, map_rank_key
from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.rules import (
    CsmiDedupWarning,
    Dealkylation,
    NDealkylation,
    ResonanceRule,
    _site_indexes,
)


def test_dealkylation_declares_directed_bond() -> None:
    assert Dealkylation.site_kind == "directed_bond"
    assert NDealkylation.site_kind == "directed_bond"
    assert issubclass(NDealkylation, ResonanceRule)


def test_anisole_directed_bond_keeps_ring_open_regioisomers() -> None:
    """Opposite map orientations on the same undirected bond → distinct products."""

    mol = MolFromSmiles("COc1ccccc1")
    ranks = mol.xf.topol_equiv
    rule = Dealkylation()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", CsmiDedupWarning)
        products = list(rule.metabolize(mol))
    assert not [
        w for w in caught if issubclass(w.category, CsmiDedupWarning)
    ], "anisole Dealkylation should not CSMI-drop under directed_bond"

    alcohols = {
        Chem.MolToSmiles(p)
        for p, info in products
        if (info.get("pattern") or {}).get("name") == "cc_alcohol"
    }
    # HEURISTICS exemplar: undirected ranks would collapse these regioisomers.
    assert alcohols >= {
        "COC=CC=CC=CO",
        "C=CC=CC(=CO)OC",
        "C=CC=C(C=CO)OC",
        "C=CC(=CC=CO)OC",
        "C=C(C=CC=CO)OC",
    }

    # Same undirected bond ends, two directed MapRankKeys (map 1↔2).
    by_und: dict[tuple[int, ...], set[tuple]] = defaultdict(set)
    for smarts, _rxn, pattern in rule.rxns:
        if pattern.get("name") != "cc_alcohol":
            continue
        reactant = smarts.split(">>", 1)[0]
        for mapped in mol.xf.smarts_matches(reactant):
            site = _site_indexes(mapped, pattern)
            by_und[bond_rank_key(ranks, site)].add(map_rank_key(ranks, mapped))
    assert any(len(dirs) >= 2 for dirs in by_und.values()), by_und


def test_ndealkylation_pyridine_ring_open_and_no_csmi() -> None:
    """Aromatic C–N matches need ResonanceRule parenting (same as Dealkylation)."""

    mol = MolFromSmiles("c1ccncc1")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", CsmiDedupWarning)
        products = list(NDealkylation().metabolize(mol))
    assert not [
        w for w in caught if issubclass(w.category, CsmiDedupWarning)
    ]
    csmi = {Chem.MolToSmiles(p) for p, _info in products}
    assert "N=CC=CC=CO" in csmi
    assert "N=CC=CC=C=O" in csmi
    for _p, info in products:
        assert len(info["site"]) == 2
        assert any(
            mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in info["site"]
        )
