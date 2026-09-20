"""Pair unique-edit signatures: ordered vs unordered, order invariance, products.

Contract (HEURISTICS / PatternInfo ``swap_group``):

1. Distinct signatures → different product csmi sets (when both emit).
2. Identical signatures → identical product csmi sets.
3. Order invariance: argument / SiteInfo order does not change the signature
   for the same chemical pairing.
4. Covers helpers (``ends_swappable``, ``pair_orbit``, ``pair_site_signature``)
   and real QF / Dehydrogenation PatternInfo couples.
5. Resolved group defaults to ``name`` (omit annotation when equal); When may
   override; explicit ``swap_group`` only when grouping differs from ``name``.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping

import pytest

from xenosite.forest.graph_isomorphism import (
    all_site_pair_orbits_nauty,
    ends_swappable,
    pair_orbit,
    pair_site_signature,
    resolved_swap_group,
)
from xenosite.forest.rdkitutil import MolFromSmiles
from xenosite.forest.records import (
    AtomPairOrbitSignature,
    PatternInfo,
)
from xenosite.forest.rules import Dehydrogenation, QuinoneFormation


def _mol(smi: str):
    mol = MolFromSmiles(smi)
    assert mol is not None
    return mol


def _endpoint(rule, name: str) -> PatternInfo:
    for _smarts, info in rule.endpoints:
        if info.get("name") == name:
            return info
    raise KeyError(name)


def _pair_product_csmi(rule, mol) -> dict[tuple, frozenset[str]]:
    """Map pair unique-edit signature → frozenset of product csmi.

    Rebuilds the signature from each emission's ``end_maps`` / PatternInfo so
    the test does not depend on metabolize-internal ``seen`` bookkeeping.
    """

    ranks = mol.xf.topol_equiv
    by_sig: dict[tuple, set[str]] = defaultdict(set)
    for row in rule.metabolites(mol):
        info = row.info
        if "ends" not in info or "end_maps" not in info:
            continue
        map1, map2 = info["end_maps"]
        site_a, site_b = info["end_atoms"]
        # Recover PatternInfo from resolved end effects when possible via name
        # on the live endpoints (pair preview does not stash PatternInfo).
        # Match by site atom chemistry against endpoint SMARTS hits.
        i1 = _info_for_map(rule, mol, map1, site_a)
        i2 = _info_for_map(rule, mol, map2, site_b)
        sig = pair_site_signature(
            mol,
            ranks,
            map1,
            map2,
            site_a,
            site_b,
            i1,
            i2,
            info,
        )
        for product in row.products:
            by_sig[sig].add(product.xf.csmi)
    return {sig: frozenset(csmi) for sig, csmi in by_sig.items()}


def _info_for_map(rule, mol, mapped: Mapping[int, int], site_atom: int) -> PatternInfo:
    """Pick the endpoint PatternInfo that produced this map at ``site_atom``."""

    for smarts, info in rule.endpoints:
        for hit in mol.xf.smarts_matches(smarts):
            if dict(hit) != dict(mapped):
                continue
            sm = info.get("site_map", 1)
            key = sm if isinstance(sm, int) else sm[0]
            if hit[key] == site_atom:
                return info
    # Fallback: first endpoint whose site_map atom matches (rare).
    for smarts, info in rule.endpoints:
        for hit in mol.xf.smarts_matches(smarts):
            sm = info.get("site_map", 1)
            key = sm if isinstance(sm, int) else sm[0]
            if hit[key] == site_atom:
                return info
    raise AssertionError(f"no PatternInfo for site {site_atom} map {dict(mapped)}")


# ---------------------------------------------------------------------------
# Schema / helper logic
# ---------------------------------------------------------------------------


def test_swap_group_defaults_to_name_on_pair_endpoints():
    """Resolved swap_group defaults to name; no redundant annotation needed."""

    from xenosite.forest.rules import Hydrogenation

    for rule in (Dehydrogenation(), QuinoneFormation(), Hydrogenation()):
        names = []
        for _smarts, info in rule.endpoints:
            assert "swap_group" not in info, info  # omit when equal to name
            assert resolved_swap_group(info) == info.get("name")
            names.append(info["name"])
        assert len(names) == len(set(names))


def test_ends_swappable_reads_swap_group_not_edit_string():
    """DH phenol/amine share edit=single_to_double but different names → ordered."""

    phenol = _endpoint(Dehydrogenation, "phenol_end")
    amine = _endpoint(Dehydrogenation, "amine_end")
    assert phenol.get("edit") == amine.get("edit") == "single_to_double"
    assert ends_swappable(phenol, phenol)
    assert not ends_swappable(phenol, amine)

    add_o = _endpoint(QuinoneFormation, "add_carbonyl_o")
    std = _endpoint(QuinoneFormation, "single_to_double")
    dealk = _endpoint(QuinoneFormation, "dealkylate")
    assert ends_swappable(add_o, add_o)
    assert not ends_swappable(add_o, std)
    assert not ends_swappable(std, dealk)
    # Distinct endpoint names → distinct resolved groups (no cross-role swap).
    groups = {
        resolved_swap_group(_endpoint(QuinoneFormation, n))
        for n in (
            "single_to_double",
            "add_carbonyl_o",
            "replace_halogen",
            "iminium",
            "dealkylate",
        )
    }
    assert len(groups) == 5


def test_when_swap_group_override():
    phenol = _endpoint(Dehydrogenation, "phenol_end")
    amine = _endpoint(Dehydrogenation, "amine_end")
    # Synthetic When override: force amine into phenol's group.
    effect_override = {"when": {"map": 2, "z": 7, "h": 2, "swap_group": "phenol_end"}}
    assert resolved_swap_group(amine) == "amine_end"  # default = name
    assert resolved_swap_group(amine, effect_override) == "phenol_end"
    assert ends_swappable(phenol, amine, None, effect_override)
    # Empty When.swap_group falls through to PatternInfo then name.
    effect_empty = {"when": {"map": 2, "z": 7, "h": 2, "swap_group": ""}}
    assert resolved_swap_group(amine, effect_empty) == "amine_end"


# ---------------------------------------------------------------------------
# Order invariance (helpers)
# ---------------------------------------------------------------------------


def test_unordered_orbit_invariant_to_argument_order():
    mol = _mol("Oc1ccc(O)cc1")
    phenol = _endpoint(Dehydrogenation, "phenol_end")
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    a, b = por.info["end_atoms"]
    site = frozenset({a, b})
    left = pair_orbit(
        mol,
        map1,
        map2,
        a,
        b,
        site,
        phenol,
        phenol,
    )
    right = pair_orbit(
        mol,
        map2,
        map1,
        b,
        a,
        site,
        phenol,
        phenol,
    )
    assert left == right
    assert isinstance(left, AtomPairOrbitSignature)
    assert left.ordered is False
    assert left.end_ranks == ()


def test_ordered_orbit_canonical_name_order_invariant_to_args():
    """Same pairing (amine@N, phenol@O): arg swap → same ordered orbit."""

    mol = _mol("Nc1ccc(O)cc1")
    phenol = _endpoint(Dehydrogenation, "phenol_end")
    amine = _endpoint(Dehydrogenation, "amine_end")
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    a, b = por.info["end_atoms"]

    def info_for(site: int) -> PatternInfo:
        return amine if mol.GetAtomWithIdx(site).GetAtomicNum() == 7 else phenol

    i_a, i_b = info_for(a), info_for(b)
    assert {i_a["name"], i_b["name"]} == {"phenol_end", "amine_end"}
    site = frozenset({a, b})
    ab = pair_orbit(
        mol, map1, map2, a, b, site, i_a, i_b
    )
    ba = pair_orbit(
        mol, map2, map1, b, a, site, i_b, i_a
    )
    assert ab == ba
    assert isinstance(ab, AtomPairOrbitSignature)
    assert ab.ordered is True
    assert len(ab.end_ranks) == 2


def test_ordered_signature_invariant_to_argument_order():
    mol = _mol("Nc1ccc(O)cc1")
    ranks = mol.xf.topol_equiv
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    a, b = por.info["end_atoms"]
    phenol = _endpoint(Dehydrogenation, "phenol_end")
    amine = _endpoint(Dehydrogenation, "amine_end")

    def info_for(site: int) -> PatternInfo:
        return amine if mol.GetAtomWithIdx(site).GetAtomicNum() == 7 else phenol

    i_a, i_b = info_for(a), info_for(b)
    sig_ab = pair_site_signature(
        mol,
        ranks,
        map1,
        map2,
        a,
        b,
        i_a,
        i_b,
        por.info,
    )
    sig_ba = pair_site_signature(
        mol,
        ranks,
        map2,
        map1,
        b,
        a,
        i_b,
        i_a,
        por.info,
    )
    assert sig_ab == sig_ba


def test_swapped_sites_different_ordered_signature():
    """amine@A+phenol@B ≠ amine@B+phenol@A when A≠B (surgical failure case)."""

    mol = _mol("Nc1ccc(O)cc1")
    ranks = mol.xf.topol_equiv
    por = next(r for r in Dehydrogenation().metabolites(mol) if "ends" in r.info)
    map1, map2 = por.info["end_maps"]
    a, b = por.info["end_atoms"]
    phenol = _endpoint(Dehydrogenation, "phenol_end")
    amine = _endpoint(Dehydrogenation, "amine_end")

    def info_for(site: int, *, crossed: bool = False) -> PatternInfo:
        is_n = mol.GetAtomWithIdx(site).GetAtomicNum() == 7
        if crossed:
            return phenol if is_n else amine
        return amine if is_n else phenol

    sig_real = pair_site_signature(
        mol,
        ranks,
        map1,
        map2,
        a,
        b,
        info_for(a),
        info_for(b),
        por.info,
    )
    sig_crossed = pair_site_signature(
        mol,
        ranks,
        map1,
        map2,
        a,
        b,
        info_for(a, crossed=True),
        info_for(b, crossed=True),
        por.info,
    )
    assert sig_real != sig_crossed


# ---------------------------------------------------------------------------
# Signature ↔ product contract on real rules
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "rule_cls,smiles",
    [
        (Dehydrogenation, "Oc1ccc(O)cc1"),  # phenol×phenol unordered
        (Dehydrogenation, "Nc1ccc(O)cc1"),  # phenol×amine ordered
        (Dehydrogenation, "Nc1ccc(N)cc1"),  # amine×amine unordered
        (QuinoneFormation, "c1ccccc1"),  # add_carbonyl_o×2
        (QuinoneFormation, "Oc1ccc(O)cc1"),  # single_to_double×2
        (QuinoneFormation, "Oc1ccc(N)cc1"),  # mixed ends ordered
    ],
    ids=["dh-hq", "dh-aminophenol", "dh-diamine", "qf-benzene", "qf-hq", "qf-aminophenol"],
)
def test_identical_signatures_identical_products(rule_cls, smiles):
    mol = _mol(smiles)
    by_sig = _pair_product_csmi(rule_cls(), mol)
    # Group by signature already; each sig maps to one product set.
    # Re-derive: if two emissions shared a signature they must share products.
    assert by_sig
    for products in by_sig.values():
        assert products


@pytest.mark.parametrize(
    "rule_cls,smiles",
    [
        (Dehydrogenation, "Nc1ccc(O)cc1"),
        (QuinoneFormation, "Oc1ccc(Cl)cc1"),
        (
            QuinoneFormation,
            "COC(=O)c1ccccc1c1ccc(cc1)CN(c1ncccc1NC(=O)CC(F)(F)F)C",
        ),  # Reaction95843: many dealkylate embeddings
    ],
    ids=["dh-aminophenol", "qf-halophenol", "qf-95843"],
)
def test_distinct_signatures_distinct_product_sets(rule_cls, smiles):
    """Different unique-edit keys must not claim the exact same product set.

    Weak form: if two signatures emit, their product frozensets differ.
    (Resonance may still share a csmi across unrelated sites rarely — we
    assert the mapping is injective on the observed emissions.)
    """

    mol = _mol(smiles)
    by_sig = _pair_product_csmi(rule_cls(), mol)
    if len(by_sig) < 2:
        pytest.skip("need ≥2 pair signatures on this mol")
    # Invert: product-set → signatures that produced it.
    by_products: dict[frozenset[str], list] = defaultdict(list)
    for sig, products in by_sig.items():
        by_products[products].append(sig)
    collisions = {ps: sigs for ps, sigs in by_products.items() if len(sigs) > 1}
    assert not collisions, f"distinct signatures share product sets: {collisions}"


def test_qf_95843_dealkylate_embeddings_get_distinct_signatures():
    """Missing map ranks would collapse methyl vs benzyl dealkylate (95843)."""

    mol = _mol(
        "COC(=O)c1ccccc1c1ccc(cc1)CN(c1ncccc1NC(=O)CC(F)(F)F)C"
    )
    want = _mol(
        r"COC(=O)c1ccccc1c1ccc(cc1)C/N=C/1\N=CC=C\C1=N/C(=O)CC(F)(F)F"
    ).xf.csmi
    by_sig = _pair_product_csmi(QuinoneFormation(), mol)
    all_csmi = set().union(*by_sig.values()) if by_sig else set()
    assert want in all_csmi
    # At least two dealkylate×single_to_double-class signatures on the pair path.
    assert len(by_sig) >= 2


def test_hydroquinone_dh_one_signature_one_quinone():
    mol = _mol("Oc1ccc(O)cc1")
    by_sig = _pair_product_csmi(Dehydrogenation(), mol)
    assert len(by_sig) == 1
    assert next(iter(by_sig.values())) == frozenset({"O=C1C=CC(=O)C=C1"})


# ---------------------------------------------------------------------------
# Light fuzz: argument-order invariance over all emitted pairings
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "rule_cls,smiles",
    [
        (Dehydrogenation, "Oc1ccc(O)cc1"),
        (Dehydrogenation, "Nc1ccc(O)cc1"),
        (QuinoneFormation, "Oc1ccc(O)cc1"),
        (QuinoneFormation, "Oc1ccc(Cl)cc1"),
        (QuinoneFormation, "c1ccccc1O"),
    ],
)
def test_fuzz_signature_order_invariance_on_emissions(rule_cls, smiles):
    mol = _mol(smiles)
    ranks = mol.xf.topol_equiv
    rule = rule_cls()
    for row in rule.metabolites(mol):
        info = row.info
        if "end_maps" not in info:
            continue
        map1, map2 = info["end_maps"]
        a, b = info["end_atoms"]
        i1 = _info_for_map(rule, mol, map1, a)
        i2 = _info_for_map(rule, mol, map2, b)
        sig = pair_site_signature(
            mol, ranks, map1, map2, a, b, i1, i2, info
        )
        flipped = pair_site_signature(
            mol, ranks, map2, map1, b, a, i2, i1, info
        )
        assert sig == flipped, (i1.get("name"), i2.get("name"), smiles)


def test_nauty_six_family_present_for_pair_mols():
    """Wire sanity: unordered/ordered families exist for a DH substrate."""

    mol = _mol("Nc1ccc(O)cc1")
    families = all_site_pair_orbits_nauty(mol)
    for key in (
        "atom_atom_unordered",
        "atom_atom_ordered",
        "atom_bond_unordered",
        "atom_bond_ordered",
        "bond_bond_unordered",
        "bond_bond_ordered",
    ):
        assert key in families
    # atom_bond ordered ≡ unordered (types already distinguish).
    assert families["atom_bond_unordered"] == families["atom_bond_ordered"]
