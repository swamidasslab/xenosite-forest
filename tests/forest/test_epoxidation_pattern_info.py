"""Epoxidation PatternInfo: dearomatizes capability (parity with Rust catalog).

Catalog ``dearomatizes`` is capability; ``resolve_effect`` narrows it to site-map
aromaticity (same shape as Hydrogenation path_end / ResonancePair ends).
"""

from __future__ import annotations

from xenosite.forest.rdkitutil import as_mol
from xenosite.forest.rules import Epoxidation


def test_epoxidation_declares_dearomatizes_capability():
    _smirks, info = Epoxidation.smirks[0]
    span = info.get("span") or {}
    assert span.get("adds") == "O"
    assert span.get("dearomatizes") is True
    for arm in info.get("possibilities") or ():
        assert arm.get("dearomatizes") is True, arm


def test_aromatic_epoxidation_resolves_dearomatizes_true():
    mol = as_mol("c1ccccc1")
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in Epoxidation().metabolize(mol)
    }
    assert flags == {True}


def test_aliphatic_epoxidation_resolves_dearomatizes_false():
    mol = as_mol("C=C")
    flags = {
        info["options"].get("dearomatizes")
        for _products, info in Epoxidation().metabolize(mol)
    }
    assert flags == {False}
