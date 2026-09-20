"""Directed-bond unique-edit keeps chemically distinct orientations.

``site_kind="directed_bond"`` keys unique-edit by ordered MapRankKey
(map 1 = oxygenated carbon). Undirected ``bond_rank_key`` would merge
anisole ring-open regioisomers that yield different products.

Public ``info["site"]`` is a frozenset; orientation is on
``info["discovered_site"]`` (ordered tuple). Unique-edit is unaffected
(yield-only presentation).
"""

from __future__ import annotations

import warnings
from collections import defaultdict

from rdkit import Chem

from xenosite.forest.graph_isomorphism import bond_rank_key, map_rank_key
from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.rules import (
    Dealkylation,
    NDealkylation,
    ResonanceRule,
    SiteDeduplicationWarning,
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
        warnings.simplefilter("always")
        products = list(rule.metabolize(mol))
    assert not [
        w for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    ], "anisole Dealkylation should not CSMI-drop under directed_bond"
    assert not caught, f"unexpected warnings on default anisole Dealk: {caught}"

    alcohols = {
        Chem.MolToSmiles(p)
        for pl, info in products
        for p in pl
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

    # Same undirected bond ends, two directed MapRankKeys (map 1↔2) at match.
    by_und: dict[tuple[int, ...], set[tuple]] = defaultdict(set)
    for smarts, _rxn, pattern in rule.rxns:
        if pattern.get("name") != "cc_alcohol":
            continue
        reactant = smarts.split(">>", 1)[0]
        for mapped in mol.xf.smarts_matches(reactant):
            site = _site_indexes(mapped, pattern, site_kind=rule.site_kind)
            by_und[bond_rank_key(ranks, site)].add(map_rank_key(ranks, mapped))
    assert any(len(dirs) >= 2 for dirs in by_und.values()), by_und

    for _pl, info in products:
        assert isinstance(info["site"], frozenset), info["site"]
        disc = info["discovered_site"]
        assert isinstance(disc, tuple), disc
        assert frozenset(disc) == info["site"]


def test_anisole_discovered_site_orientation_vs_frozenset_site() -> None:
    """Products keep frozenset site; discovered_site carries map order."""

    mol = MolFromSmiles("COc1ccccc1")
    products = list(Dealkylation().metabolize(mol))
    alcohol_rows = [
        (Chem.MolToSmiles(p), info["site"], info["discovered_site"])
        for pl, info in products
        for p in pl
        if (info.get("pattern") or {}).get("name") == "cc_alcohol"
    ]
    assert alcohol_rows
    assert all(isinstance(site, frozenset) for _s, site, _d in alcohol_rows)
    assert all(isinstance(disc, tuple) for _s, _site, disc in alcohol_rows)

    # Distinct products (regioisomers) — directed unique-edit preserved.
    smiles = {s for s, _site, _d in alcohol_rows}
    assert len(smiles) >= 5

    # discovered_site tuples are ordered; at least one is not sorted
    # (map 1 ≠ lower index) or we still see map1=carbon chemistry via atoms.
    for _smi, site, disc in alcohol_rows:
        assert frozenset(disc) == site
        a0, a1 = disc
        # map 1 = oxygenated carbon for cc_alcohol
        assert mol.GetAtomWithIdx(a0).GetAtomicNum() == 6


def test_ndealkylation_pyridine_frozenset_site_directed_discovered() -> None:
    """Aromatic C–N: frozenset site; discovered_site (C, N) ordered."""

    mol = MolFromSmiles("c1ccncc1")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        products = list(NDealkylation().metabolize(mol))
    assert not [
        w for w in caught if issubclass(w.category, SiteDeduplicationWarning)
    ]
    assert not caught, f"unexpected warnings on pyridine NDealk: {caught}"
    csmi = {Chem.MolToSmiles(p) for pl, _info in products for p in pl}
    assert "N=CC=CC=CO" in csmi
    assert "N=CC=CC=C=O" in csmi
    for _pl, info in products:
        site = info["site"]
        disc = info["discovered_site"]
        assert isinstance(site, frozenset), site
        assert isinstance(disc, tuple), disc
        assert len(site) == 2 and len(disc) == 2
        assert frozenset(disc) == site
        assert any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in site)
        # map 1 = oxygenated carbon, map 2 = nitrogen
        assert mol.GetAtomWithIdx(disc[0]).GetAtomicNum() == 6
        assert mol.GetAtomWithIdx(disc[1]).GetAtomicNum() == 7
