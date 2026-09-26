"""``Effect.exclusive_partner``: bridging N/O must not be shared across pair ends.

Chemistry sets the bit (C14); the couple loop reads it. Not a global
shared-map refuse — methide alkyl (partner C) stays allowed to share.
"""

from __future__ import annotations

from xenosite.forest.rdkit_api import MolFromSmiles
from xenosite.forest.rules import (
    QuinoneFormation,
    _shared_exclusive_partner,
    _site_atoms,
    resolve_effect,
)


def _pair_emissions_sharing_exclusive_partner(smiles: str) -> list[str]:
    """Return human-readable couples that would share an exclusive partner."""

    mol = MolFromSmiles(smiles)
    assert mol is not None
    rule = QuinoneFormation()
    bad: list[str] = []
    for _products, info in rule.metabolize(mol):
        if "ends" not in info or "end_maps" not in info:
            continue
        map1, map2 = info["end_maps"]
        end1, end2 = info["ends"]
        # Rebuild PatternInfo stubs from ends is unnecessary: gate already ran.
        # Assert no survivor still shares exclusive partner atoms.
        info1 = info2 = None
        for smarts, pattern in rule.endpoints:
            for mapped in mol.xf.smarts_matches(smarts):
                if dict(mapped) == dict(map1):
                    info1 = pattern
                if dict(mapped) == dict(map2):
                    info2 = pattern
        if info1 is None or info2 is None:
            continue
        if _shared_exclusive_partner(map1, info1, end1, map2, info2, end2):
            bad.append(
                f"{info1.get('name')}×{info2.get('name')} "
                f"ends={info.get('end_atoms')}"
            )
    return bad


def test_bridging_n_emissions_do_not_share_exclusive_partner():
    """N-methyl diphenylamine: both aryl carbons match; shared N refused."""

    bad = _pair_emissions_sharing_exclusive_partner("c1ccc(N(C)c2ccccc2)cc1")
    assert not bad, "emitted couples still share exclusive partner:\n" + "\n".join(
        bad
    )


def test_dihydroacridine_bridging_n_refused():
    bad = _pair_emissions_sharing_exclusive_partner("c1ccc2c(c1)Nc1ccccc1C2")
    assert not bad


def test_phenoxazine_bridging_n_refused():
    bad = _pair_emissions_sharing_exclusive_partner("c1ccc2c(c1)Nc1ccccc1O2")
    assert not bad


def test_catechol_identical_o_partners_still_emit():
    """Distinct O atoms (identical role) are not exclusive-partner shares."""

    mol = MolFromSmiles("Oc1ccccc1O")
    rule = QuinoneFormation()
    pairs = [info for _, info in rule.metabolize(mol) if "ends" in info]
    assert pairs, "expected ortho catechol pair emissions"
    # At least one couple with partner O on both ends (distinct atoms).
    saw = False
    for info in pairs:
        e0, e1 = info["ends"]
        if e0.get("partner") == "O" and e1.get("partner") == "O":
            saw = True
            assert not (
                e0.get("exclusive_partner")
                and e1.get("exclusive_partner")
                and _site_atoms(info["end_maps"][0], {"site_map": 1})
                == _site_atoms(info["end_maps"][1], {"site_map": 1})
            )
    assert saw


def test_methide_shared_alkyl_not_gated_by_exclusive_partner():
    """Shared CH2-bridge methide (partner C) is outside C14 — bit stays off."""

    mol = MolFromSmiles("c1ccc2c(c1)Cc1ccccc1O2")  # xanthene
    rule = QuinoneFormation()
    for smarts, info in rule.endpoints:
        if info.get("name") != "single_to_double":
            continue
        for mapped in mol.xf.smarts_matches(smarts):
            effect = resolve_effect(mol, mapped, info)
            if effect.get("partner") == "C":
                assert not effect.get("exclusive_partner")


def test_qf_hetero_branches_declare_exclusive_partner():
    rule = QuinoneFormation()
    by_name = {info.get("name"): info for _, info in rule.endpoints}
    for name in ("iminium", "dealkylate", "single_to_double", "replace_halogen"):
        info = by_name[name]
        arms = [
            p
            for p in (info.get("possibilities") or ())
            if p.get("exclusive_partner")
            or (name == "iminium" and info.get("span", {}).get("exclusive_partner"))
        ]
        if name == "iminium":
            # Single-possibility pattern: bit on the only arm / span.
            poss = info.get("possibilities") or ()
            assert any(p.get("exclusive_partner") for p in poss) or info.get(
                "span", {}
            ).get("exclusive_partner")
        elif name == "single_to_double":
            # O/N arms yes; methide C arms no.
            assert any(p.get("exclusive_partner") for p in info["possibilities"])
            assert any(
                p.get("methide") and not p.get("exclusive_partner")
                for p in info["possibilities"]
            )
        else:
            assert any(p.get("exclusive_partner") for p in info["possibilities"]), name
