"""Topological emission dedup should match the XenoSite UI identity key.

UI key: pathway + product SMILES + sorted topological ranks of site atoms
(``CanonicalRankAtoms(..., breakTies=False)`` / predict ``atoms.cipRank``).

Collapse only when that full key matches — never drop a chemically distinct
product SMILES just because it sits on a symmetry-equivalent site.
"""

from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite._archive_forest import rules
from xenosite._archive_forest.base import AtomTracker, can_smi_set


def _product_smiles(mets) -> frozenset:
    mets = [m for m in mets if m]
    if not mets:
        return frozenset()
    return frozenset(can_smi_set(mets))


def _ui_identity_keys(rule, mol):
    """UI-style identity keys over raw ``metabolites()`` (no topo filter)."""
    ranks = AtomTracker.topol_equiv(mol)
    keys = set()
    for (raw_name, atoms), mets in rule.metabolites(mol):
        prods = _product_smiles(mets)
        if not prods:
            continue
        pathway = raw_name.split("_", 1)[0]
        rank_key = tuple(sorted(ranks[i] for i in atoms))
        keys.add((pathway, rank_key, prods))
    return keys


def _all_product_smiles(iterable) -> set[str]:
    out: set[str] = set()
    for _site, mets in iterable:
        out.update(_product_smiles(mets))
    return out


def _metabolize_rows(rule, mol):
    return list(rule.metabolize(mol, tag_atoms=False))


def test_phenol_epoxidation_metabolize_matches_ui_topo_dedup():
    """Ortho/meta/para: equivalent ring sites with the same epoxide collapse."""
    mol = MolFromSmiles("Oc1ccccc1")
    rule = rules.Epoxidation()
    keys = _ui_identity_keys(rule, mol)
    rows = _metabolize_rows(rule, mol)
    assert len(keys) == 3
    assert len(rows) == len(keys)


def test_naphthalene_quinone_metabolize_matches_ui_topo_dedup():
    mol = MolFromSmiles("c1ccc2ccccc2c1")
    rule = rules.QuinoneFormation()
    keys = _ui_identity_keys(rule, mol)
    rows = _metabolize_rows(rule, mol)
    assert len(keys) < 22  # today's leaky metabolize count
    assert len(rows) == len(keys)


def test_distinct_product_smiles_are_never_dropped():
    """Topo collapse must not discard a unique product SMILES."""
    for smi, rule in [
        ("Oc1ccccc1", rules.Epoxidation()),
        ("c1ccc2ccccc2c1", rules.QuinoneFormation()),
        ("c1ccccc1", rules.Epoxidation()),
        ("CN(C)Cc1ccccc1", rules.Dealkylation()),
    ]:
        mol = MolFromSmiles(smi)
        raw = _all_product_smiles(rule.metabolites(mol))
        emitted = _all_product_smiles(rule.metabolize(mol, tag_atoms=False))
        assert raw
        assert emitted == raw, f"{smi} {rule.__class__.__name__}: dropped {raw - emitted}"


def test_distinct_products_at_same_topo_class_are_kept():
    """UI keeps different product SMILES even when site ranks match."""
    mol = MolFromSmiles("c1ccc2ccccc2c1")
    rule = rules.QuinoneFormation()
    keys = _ui_identity_keys(rule, mol)
    by_ranks: dict = {}
    for _pathway, ranks, prods in keys:
        by_ranks.setdefault(ranks, set()).add(prods)
    assert any(len(prods) > 1 for prods in by_ranks.values())
    assert len(_metabolize_rows(rule, mol)) == len(keys)


def test_benzene_epoxidation_single_topo_class():
    mol = MolFromSmiles("c1ccccc1")
    rule = rules.Epoxidation()
    keys = _ui_identity_keys(rule, mol)
    assert len(keys) == 1
    assert len(_metabolize_rows(rule, mol)) == 1
