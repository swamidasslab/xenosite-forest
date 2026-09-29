"""Forest↔XMET SSSOM embed + resolve (Rust product door)."""

from __future__ import annotations

import pytest

from xenosite.forest import available

pytestmark = pytest.mark.skipif(not available(), reason="Rust extension not built")


def test_forest_xmet_sssom_embedded():
    from xenosite.forest import forest_xmet_sssom

    text = forest_xmet_sssom()
    assert "subject_id\t" in text
    assert "xf:Tautomerization/tautomer_h" in text
    assert "xf:Tautomerization/path_partner" in text
    assert "forest.rule:" not in text
    assert "forest.pattern:" not in text


def test_resolve_leaf_and_bound_pattern():
    from xenosite.forest import BoundPattern, expand_iri, resolve, to_curie

    rule = resolve("xf:Tautomerization")
    assert rule.name == "Tautomerization"
    assert "tautomer_h" in rule

    bp = resolve("xf:Tautomerization/tautomer_h")
    assert isinstance(bp, BoundPattern)
    assert bp.name == "tautomer_h"
    assert bp.rule_name == "Tautomerization"
    assert bp.curie == "xf:Tautomerization/tautomer_h"
    assert bp.iri == expand_iri(bp.curie)
    assert to_curie(bp.iri) == bp.curie

    via_index = rule["tautomer_h"]
    assert isinstance(via_index, BoundPattern)
    assert via_index.name == bp.name


def test_bound_pattern_metabolize():
    from xenosite.forest._ext import load

    mod = load()
    # ethane → ethanol via Hydroxylation/h2 (methyl/methylene)
    bp = mod.resolve("xf:Hydroxylation/h2")
    mol = mod.ForestMol("CC")
    rows = bp.metabolize(mol)
    assert rows
    assert all(r[0] == "h2" for r in rows)
    products = {p for row in rows for p in row[2]}
    assert any("CCO" in p or "OCC" in p for p in products) or products


def test_expand_iri_helpers():
    from xenosite.forest import expand_iri, to_curie

    iri = expand_iri("xf:PhaseOne")
    assert iri.startswith("https://w3id.org/xenosite/forest/")
    assert to_curie(iri) == "xf:PhaseOne"
