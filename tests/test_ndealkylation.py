"""NDealkylation / ND: Dealkylation post-filtered to nitrogen-containing sites."""

from __future__ import annotations

from rdkit.Chem.rdmolfiles import MolFromSmiles

from xenosite.forest import load_ruleset, rules
from xenosite.forest.base import can_smi

# Amine + ether (diphenhydramine): N–C and O–C / C–C dealkylation sites.
DIPHENHYDRAMINE = "CN(C)CCOC(c1ccccc1)c2ccccc2"


def _site_has_nitrogen(mol, site) -> bool:
    return any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in site)


def _sites_and_products(rule_or_rs, mol, *, via_ruleset: bool = False):
    """Collect (site frozenset, frozenset of canonical product SMILES)."""
    out = []
    if via_ruleset:
        iterable = rule_or_rs.metabolites(mol, unique=True)
        for (rxnname, site), metabolites in iterable:
            smiles = frozenset(can_smi(rdmol=m)[0] for m in metabolites if m)
            out.append((frozenset(site), smiles))
    else:
        for (_rxn, site), metabolites in rule_or_rs.metabolites(mol):
            smiles = frozenset(can_smi(rdmol=m)[0] for m in metabolites if m)
            out.append((frozenset(site), smiles))
    return out


def test_dealkylation_still_emits_non_nitrogen_sites_on_diphenhydramine():
    """Regression: broad Dealkylation / UO must keep O–C and C–C cleavage."""
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    sites = {s for s, _ in _sites_and_products(rules.Dealkylation(), mol)}
    assert sites
    non_n = [s for s in sites if not _site_has_nitrogen(mol, s)]
    assert non_n, "Dealkylation should still emit at least one non-N site"

    uo = load_ruleset("UO.Dealkylation")
    uo_sites = {s for s, _ in _sites_and_products(uo, mol, via_ruleset=True)}
    assert any(not _site_has_nitrogen(mol, s) for s in uo_sites)


def test_ndealkylation_only_emits_nitrogen_containing_sites():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    rule = rules.NDealkylation()
    pairs = _sites_and_products(rule, mol)
    assert pairs, "expected at least one N-dealkylation product"
    for site, _ in pairs:
        assert _site_has_nitrogen(mol, site), f"non-N site leaked: {sorted(site)}"


def test_nd_ruleset_loadable_and_nitrogen_only():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    rs = load_ruleset("ND")
    pairs = _sites_and_products(rs, mol, via_ruleset=True)
    assert pairs
    for site, _ in pairs:
        assert _site_has_nitrogen(mol, site), f"ND leaked non-N site: {sorted(site)}"
    assert any(r.name == "NDealkylation" for r in rs)


def test_nd_products_are_subset_of_dealkylation():
    mol = MolFromSmiles(DIPHENHYDRAMINE)
    dealk = {
        (site, products)
        for site, products in _sites_and_products(rules.Dealkylation(), mol)
    }
    nd = {
        (site, products)
        for site, products in _sites_and_products(rules.NDealkylation(), mol)
    }
    assert nd
    assert nd <= dealk
    assert len(nd) < len(dealk)


def test_ndealkylation_short_circuits_when_molecule_has_no_nitrogen():
    """No N → no N-dealkylation work (Dealkylation may still cleave O/C)."""
    mol = MolFromSmiles("CCOC")  # ether only
    assert not any(a.GetAtomicNum() == 7 for a in mol.GetAtoms())
    dealk_pairs = _sites_and_products(rules.Dealkylation(), mol)
    assert dealk_pairs, "Dealkylation should still fire on ethers"
    assert _sites_and_products(rules.NDealkylation(), mol) == []
    assert _sites_and_products(load_ruleset("ND"), mol, via_ruleset=True) == []
