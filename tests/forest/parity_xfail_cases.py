"""Known leaf product-parity gaps under C18 (chemistry-first plan).

Parametric :func:`test_leaf_rule_product_parity` marks these ``xfail`` until
Phase 1–3 close them against **approved chemistry**, not literal Python match.

Product / charge / tautomer **standardization** belongs in Rust as explicit,
auditable transforms — not RDKit ``SanitizeMol``, Chematic black-box
reparsing, or silent unique_csmi (C11 / C18).
"""

from __future__ import annotations

# (rule_name, substrate_smiles) → short reason (C18 tag).
# Keep SMILES exactly as in the parametric corpus / PARITY_REMAINING.md.
PARITY_PRODUCT_XFAIL: dict[tuple[str, str], str] = {
    # --- Phase 1 (implement next) ---
    (
        "Dealkylation",
        "[O-][N+](=O)c1ccccc1",
    ): "C18 P1: nitro charge form — normalize to [N+](=O)[O-] in Rust",
    (
        "Dealkylation",
        "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",
    ): "C18 P1: nitro charge form (chloramphenicol)",
    (
        "Dealkylation",
        "O=Nc1ccccc1",
    ): "C18 P1: nitrosobenzene ring-open extras (Rust-ahead)",
    (
        "Dealkylation",
        "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
    ): "C18 P1/P3: olsalazine Dealk ring-open extras",
    (
        "Dealkylation",
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    ): "C18 P1: SMx C–N cleavage Rust-ahead",
    (
        "Dealkylation",
        "c1ccc2nnccc2c1",
    ): "C18 P1: cinnoline ring-open keep N=N (not NN)",
    (
        "NDealkylation",
        "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",
    ): "C18 P1: nitro leave charge form",
    (
        "NDealkylation",
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    ): "C18 P1: SMx C–N cleavage Rust-ahead",
    (
        "NDealkylation",
        "c1ccc2nnccc2c1",
    ): "C18 P1: cinnoline ring-open N=N",
    (
        "Dehydration",
        "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
    ): "C18 P1: Py ketene-like junk; Rust refuse",
    (
        "Dehydration",
        "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",
    ): "C18 P1: Py quinoid junk; Rust refuse",
    (
        "Dehydrogenation",
        "c1ccc2c(c1)Nc1ccccc1C2",
    ): "C18 P1: dihydroacridine→acridine under DH (Rust-ahead)",
    (
        "Hydrogenation",
        "O=C=Nc1ccccc1",
    ): "C18 P1: cumulated×ring H refuse (Rust-ahead)",
    (
        "Hydrogenation",
        "S=C=Nc1ccccc1",
    ): "C18 P1: cumulated×ring H refuse (Rust-ahead)",
    (
        "Hydrogenation",
        "N=C=Nc1ccccc1",
    ): "C18 P1: cumulated×ring H refuse (Rust-ahead)",
    (
        "Hydrogenation",
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    ): "C18 P1: refuse generic S=O hydrogenation",
    (
        "QuinoneFormation",
        "c1ccc(N(C)c2ccccc2)cc1",
    ): "C18 P1: Ph2NMe enumerate Me and Ph carbon leaves",
    # --- Phase 2 (policy / schema) ---
    (
        "Dealkylation",
        "CC(=O)Oc1ccccc1C(=O)O",
    ): "C18 P2/C13: aspirin equivalent-embedding (scissile bond)",
    (
        "Dehydrogenation",
        "C1=CC2OC2C=C1",
    ): "C18 P2/C13: benzene-oxide DH event identity",
    (
        "Hydrogenation",
        "[O-][N+](=O)c1ccccc1",
    ): "C18 P2: nitro out of generic H → NitrogenReduction",
    (
        "Hydrogenation",
        "O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",
    ): "C18 P2: nitro out of generic H (chloramphenicol)",
    (
        "NitrogenReduction",
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    ): "C18 P2: isoxazole tautomer canon (Rust, not sanitize)",
    (
        "QuinoneFormation",
        "c1ccc2c(c1)Nc1ccccc1C2",
    ): "C18 P2: acridine is DH-only, not QF",
    # --- Phase 3 (deferred multi-ring) ---
    (
        "Dehydrogenation",
        "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
    ): "C18 P3: olsalazine linked-π DH — defer",
    (
        "Hydrogenation",
        "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
    ): "C18 P3: olsalazine no generic path across azo — defer",
    (
        "Hydrogenation",
        "c1ccc2nnccc2c1",
    ): "C18 P3: cinnoline multi-ring H — defer (P1 owns ring-open N=N)",
    (
        "QuinoneFormation",
        "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O",
    ): "C18 P3: olsalazine cross-azo QF — defer",
    (
        "QuinoneFormation",
        "C1=CC2OC2C=C1",
    ): "C18 P3: benzene-oxide QF until quinonoid contract — defer",
    (
        "QuinoneFormation",
        "c1ccc2nnccc2c1",
    ): "C18 P3: cinnoline multi-ring QF — defer",
}


def parity_xfail_reason(rule_name: str, smiles: str) -> str | None:
    """Return xfail reason if ``(rule, smiles)`` is a known C18 gap."""

    return PARITY_PRODUCT_XFAIL.get((rule_name, smiles))
