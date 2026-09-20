"""Port of ``tests/test_ndealkylation.py`` without forest ``load_ruleset``.

Diphenhydramine SMILES is imported from the forest suite. ND ruleset load
and UO.Dealkylation load are skipped (unfinished public surface).
"""

from __future__ import annotations

from test_ndealkylation import DIPHENHYDRAMINE

from xenosite.forest.rdkit_api import Mol, MolFromSmiles
from xenosite.forest.rules import Dealkylation, NDealkylation

from .helpers import canon


def _site_has_nitrogen(mol: Mol, site) -> bool:
    atoms = (site,) if isinstance(site, int) else tuple(site)
    return any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in atoms)


def _sites_and_products(rule, mol: Mol):
    out = []
    for products, info in rule.metabolize(mol):
        site = info["site"]
        site_fs = frozenset((site,) if isinstance(site, int) else site)
        out.append((site_fs, frozenset(canon(p) for p in products)))
    return out


def test_dealkylation_still_emits_non_nitrogen_sites_on_diphenhydramine():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    sites = {s for s, _ in _sites_and_products(Dealkylation(), mol)}
    assert sites
    non_n = [s for s in sites if not _site_has_nitrogen(mol, s)]
    assert non_n, "Dealkylation should still emit at least one non-N site"


def test_ndealkylation_only_emits_nitrogen_containing_sites():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    pairs = _sites_and_products(NDealkylation(), mol)
    assert pairs, "expected at least one N-dealkylation product"
    for site, _ in pairs:
        assert _site_has_nitrogen(mol, site), f"non-N site leaked: {sorted(site)}"


def test_nd_products_are_subset_of_dealkylation():
    """Every ND product SMILES must also appear under Dealkylation.

    Site indexes need not match across rules; product chemistry is the check.
    """

    mol = MolFromSmiles(DIPHENHYDRAMINE)
    dealk = {p for _site, products in _sites_and_products(Dealkylation(), mol) for p in products}
    nd = {p for _site, products in _sites_and_products(NDealkylation(), mol) for p in products}
    assert nd
    assert nd <= dealk
    assert len(nd) < len(dealk)


def test_ndealkylation_short_circuits_when_molecule_has_no_nitrogen():
    mol = MolFromSmiles("CCOC")
    assert not any(a.GetAtomicNum() == 7 for a in mol.GetAtoms())
    assert _sites_and_products(Dealkylation(), mol), "Dealkylation should still fire on ethers"
    assert _sites_and_products(NDealkylation(), mol) == []


# Skipped: test_nd_ruleset_loadable_and_nitrogen_only — needs forest load_ruleset("ND").
