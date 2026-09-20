"""Unique-edit / unique_csmi never drop distinct product SMILES.

Adapted from ``tests/test_topo_emission.py``. Forest keyed
``(rule, topo ranks, product SMILES)``. Live forest splits that into unique-edit
(site ranks + pair orbit) and ``unique_csmi`` ``(rule, pattern, csmi)`` —
see docs/forest/DIVERGENCES.md (symmetry collapse approved; product csmi key approved).
"""

from __future__ import annotations

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import Dealkylation, Epoxidation, QuinoneFormation


def test_phenol_epoxidation_symmetry_collapse():
    """Ortho/meta/para: three distinct epoxides under approved symmetry collapse."""

    mol = MolFromSmiles("Oc1ccccc1")
    rows = list(Epoxidation().metabolize(mol))
    assert len(rows) == 3
    assert len({p.xf.csmi for _pl, _ in rows for p in _pl}) == 3


def test_naphthalene_quinone_symmetry_collapse():
    """Naphthalene QF: unique-edit + unique_csmi keep a small orbit set."""

    mol = MolFromSmiles("c1ccc2ccccc2c1")
    rows = list(QuinoneFormation().metabolize(mol))
    assert len(rows) < 22
    assert len(rows) == len({p.xf.csmi for _pl, _ in rows for p in _pl})


def test_unique_csmi_never_drops_distinct_smiles():
    """``unique_csmi`` may collapse duplicate CSMIs, never a distinct structure."""

    for smi, rule in [
        ("Oc1ccccc1", Epoxidation()),
        ("c1ccc2ccccc2c1", QuinoneFormation()),
        ("c1ccccc1", Epoxidation()),
        ("CN(C)Cc1ccccc1", Dealkylation()),
    ]:
        mol = MolFromSmiles(smi)
        with_dedup = {
            p.xf.csmi for pl, _ in rule.metabolize(mol, unique_csmi=True) for p in pl
        }
        without = {
            p.xf.csmi for pl, _ in rule.metabolize(mol, unique_csmi=False) for p in pl
        }
        assert without
        assert with_dedup == without, (
            f"{smi} {rule.__class__.__name__}: dropped {without - with_dedup}"
        )


def test_benzene_epoxidation_single_topo_class():
    mol = MolFromSmiles("c1ccccc1")
    rows = list(Epoxidation().metabolize(mol))
    assert len(rows) == 1
