"""Rust-door conjugation products: star + CX atomLabel pass-through."""

from __future__ import annotations

import pytest

from xenosite.forest import ForestMol, leaf_rule, reactivity

pytestmark = pytest.mark.skipif(
    __import__("xenosite.forest", fromlist=["available"]).available() is False,
    reason="Rust extension not built",
)

# leaf name → hardcoded CX conjugate label + probe substrate
_CASES = [
    ("Acetylation", "Ac", "CCO"),
    ("Sulfation", "SO3", "CCO"),
    ("Glucuronidation", "GlcA", "c1ccccc1O"),
    ("Glutathionation", "GSH", "C1OC1c1ccccc1"),
    ("GSH", "GSH", "C1OC1c1ccccc1"),
    ("Protein", "Protein", "C1OC1c1ccccc1"),
    ("DNA", "DNA", "C1OC1c1ccccc1"),
    ("Cyanide", "Cyanide", "C1OC1c1ccccc1"),
]


@pytest.mark.parametrize("leaf,label,smi", _CASES, ids=[c[0] for c in _CASES])
def test_conjugation_product_csmis_are_cx_labeled_stars(leaf, label, smi):
    rule = leaf_rule(leaf)
    mol = ForestMol(smi)
    rows = list(rule.metabolize(mol))
    assert rows, f"{leaf}: no products on {smi}"
    for emission in rows:
        products = emission.product_csmis()
        assert products, f"{leaf}: empty product_csmis"
        for cx in products:
            bare = cx.split()[0]
            assert "*" in bare, f"{leaf}: missing * in {cx!r}"
            assert label in cx, f"{leaf}: missing label {label!r} in {cx!r}"
            assert "|" in cx, f"{leaf}: expected CX block in {cx!r}"


def test_reactivity_catalog_products_carry_head_labels():
    mol = ForestMol("C1OC1c1ccccc1")
    rows = list(reactivity().metabolize(mol))
    assert rows
    labels_seen = set()
    for emission in rows:
        for cx in emission.product_csmis():
            assert "*" in cx.split()[0]
            for lab in ("GSH", "Protein", "DNA", "Cyanide"):
                if lab in cx:
                    labels_seen.add(lab)
    assert labels_seen & {"GSH", "Protein", "DNA", "Cyanide"}
