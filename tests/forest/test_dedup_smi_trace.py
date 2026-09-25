"""atom_trace caches dedup_smi per depth frame (survives clear_structure)."""

from __future__ import annotations

from rdkit import Chem

from xenosite.forest.rules import Hydroxylation


def test_stamp_caches_depth0_dedup_smi():
    mol = Chem.MolFromSmiles("CC")
    stamped = mol.xf.tracing._stamp()
    frames = stamped._forest["atom_trace"]["dedup_smi"]
    assert frames == ["CC"]
    assert stamped.xf.tracing.dedup_smi == "CC"
    assert stamped.xf.tracing.dedup_smi_at(0) == "CC"
    assert stamped.xf.tracing.dedup_smi_at(1) is None


def test_metabolize_appends_product_frame_and_survives_clear():
    mol = Chem.MolFromSmiles("CC")
    products = list(Hydroxylation().metabolize(mol))
    assert products
    finished, _info = products[0]
    product = finished[0]
    trace = product._forest["atom_trace"]
    assert trace["depth"] == 1
    assert len(trace["dedup_smi"]) == 2
    assert trace["dedup_smi"][0] == "CC"
    assert trace["dedup_smi"][1] == product.xf.tracing.dedup_smi
    assert product.xf.tracing.dedup_smi == "CCO"
    # of_products cleared structure caches; display csmi is not re-seeded
    # until read — dedup key stays on the trace.
    assert "csmi" not in (product._forest.get("cache") or {})
    assert product.xf.tracing.dedup_smi == "CCO"


def test_dedup_smi_at_parent_frame_after_hop():
    mol = Chem.MolFromSmiles("CC")
    finished, _ = next(Hydroxylation().metabolize(mol))
    product = finished[0]
    assert product.xf.tracing.dedup_smi_at(0) == "CC"
    assert product.xf.tracing.dedup_smi_at(1) == "CCO"
