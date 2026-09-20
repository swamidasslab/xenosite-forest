"""Dehydrogenation unique-edit: bond_atom orbits + ordered vs unordered ends."""

from __future__ import annotations

from xenosite.forest.graph_isomorphism import (
    bond_atom_orbit_key,
    ends_swappable,
    pair_orbit,
    pair_site_signature,
    resolved_swap_group,
)
from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.records import BondAtomOrbitSignature, BondAtomPairOrbitSignature
from xenosite.forest.rules import Dehydrogenation, QuinoneFormation


def _mol(smi: str):
    mol = MolFromSmiles(smi)
    assert mol is not None
    return mol


def test_dehydrogenation_declares_bond_atom_unique_orbit():
    assert Dehydrogenation.unique_orbit == "bond_atom"
    assert QuinoneFormation.unique_orbit == "atom_atom"


def test_hydroquinone_pair_emits_once_not_per_resonance_parent():
    """Symmetric phenol+phenol: resonance parents must not double-emit."""

    mol = _mol("Oc1ccc(O)cc1")
    rows = list(Dehydrogenation().metabolites(mol))
    pair_rows = [r for r in rows if "ends" in r.info]
    assert len(pair_rows) == 1
    assert [p.xf.csmi for p in pair_rows[0].products] == ["O=C1C=CC(=O)C=C1"]


def test_aminophenol_asymmetric_ends_emit_once():
    """Phenol+amine is an ordered role pair; still one unique-edit."""

    mol = _mol("Nc1ccc(O)cc1")
    pair_rows = [r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info]
    assert len(pair_rows) == 1
    assert pair_rows[0].products[0].xf.csmi == "N=C1C=CC(=O)C=C1"


def test_end_roles_swappable_vs_ordered():
    phenol = Dehydrogenation.endpoints[0][1]
    amine = Dehydrogenation.endpoints[1][1]
    assert resolved_swap_group(phenol) == "phenol_end"
    assert resolved_swap_group(amine) == "amine_end"
    assert ends_swappable(phenol, phenol)
    assert not ends_swappable(phenol, amine)
    # QF: distinct names are never swappable; same add_carbonyl_o is.
    add_o = next(i for _, i in QuinoneFormation.endpoints if i.get("name") == "add_carbonyl_o")
    std = next(i for _, i in QuinoneFormation.endpoints if i.get("name") == "single_to_double")
    assert ends_swappable(add_o, add_o)
    assert not ends_swappable(add_o, std)


def test_asymmetric_pair_orbit_is_name_ordered_not_arg_ordered():
    """Different PatternInfo roles → ordered by canonical name; arg-swap stable."""

    mol = _mol("Nc1ccc(O)cc1")
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    site_a, site_b = por.info["end_atoms"]
    phenol = Dehydrogenation.endpoints[0][1]
    amine = Dehydrogenation.endpoints[1][1]

    def info_for_site(site_atom: int) -> object:
        atom = mol.GetAtomWithIdx(site_atom)
        return amine if atom.GetAtomicNum() == 7 else phenol

    i1, i2 = info_for_site(site_a), info_for_site(site_b)
    assert {i1.get("name"), i2.get("name")} == {"phenol_end", "amine_end"}
    assert not ends_swappable(i1, i2)
    orbit = pair_orbit(
        mol,
        map1,
        map2,
        site_a,
        site_b,
        frozenset({site_a, site_b}),
        i1,
        i2,
        unique_orbit="bond_atom",
    )
    assert isinstance(orbit, BondAtomPairOrbitSignature)
    assert orbit.ordered is True
    assert all(isinstance(x, BondAtomOrbitSignature) for x in orbit.ends)
    reversed_orbit = pair_orbit(
        mol,
        map2,
        map1,
        site_b,
        site_a,
        frozenset({site_a, site_b}),
        i2,
        i1,
        unique_orbit="bond_atom",
    )
    # Same chemical pairing → identical signature under argument swap.
    assert orbit == reversed_orbit


def test_symmetric_pair_orbit_is_unordered():
    mol = _mol("Oc1ccc(O)cc1")
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    site_a, site_b = por.info["end_atoms"]
    phenol = Dehydrogenation.endpoints[0][1]
    left = pair_orbit(
        mol,
        map1,
        map2,
        site_a,
        site_b,
        frozenset({site_a, site_b}),
        phenol,
        phenol,
        unique_orbit="bond_atom",
    )
    right = pair_orbit(
        mol,
        map2,
        map1,
        site_b,
        site_a,
        frozenset({site_a, site_b}),
        phenol,
        phenol,
        unique_orbit="bond_atom",
    )
    assert left == right
    assert isinstance(left, BondAtomPairOrbitSignature)
    assert left.ordered is False


def test_one_bond_alcohol_uses_bond_atom_orbit_key():
    """SMARTS alcohol DH: site is one carbon; unique-edit carries bond_atom."""

    mol = _mol("CCO")
    dh = Dehydrogenation()
    alcohol = [r for r in dh.metabolites(mol) if r.info.get("pattern", {}).get("name") == "alcohol"]
    assert alcohol
    # Bond between C(map1) and O(map2); site atom is map1 carbon.
    mapped = mol.xf.smarts_matches("[#6h:1]-[#8H1:2]")[0]
    bond = mol.GetBondBetweenAtoms(mapped[1], mapped[2])
    assert bond is not None
    sig = bond_atom_orbit_key(mol, bond.GetIdx(), mapped[1])
    assert isinstance(sig, BondAtomOrbitSignature)


def test_pair_site_signature_order_invariant_for_asymmetric():
    """Argument swap of the same amine/phenol pairing keeps one signature."""

    mol = _mol("Nc1ccc(O)cc1")
    ranks = mol.xf.topol_equiv
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    site_a, site_b = por.info["end_atoms"]
    phenol = Dehydrogenation.endpoints[0][1]
    amine = Dehydrogenation.endpoints[1][1]

    def info_for(site: int):
        return amine if mol.GetAtomWithIdx(site).GetAtomicNum() == 7 else phenol

    i_a, i_b = info_for(site_a), info_for(site_b)
    sig_ab = pair_site_signature(
        mol,
        ranks,
        map1,
        map2,
        site_a,
        site_b,
        i_a,
        i_b,
        por.info,
        unique_orbit="bond_atom",
    )
    sig_ba = pair_site_signature(
        mol,
        ranks,
        map2,
        map1,
        site_b,
        site_a,
        i_b,
        i_a,
        por.info,
        unique_orbit="bond_atom",
    )
    assert sig_ab == sig_ba
    assert ends_swappable(phenol, phenol)
    assert not ends_swappable(phenol, amine)
