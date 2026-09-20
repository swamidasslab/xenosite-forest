# Cross-rule redundant pairs

Inventory of PhaseOne rule pairs that emit the same product CSMI on some
substrates. Undirected. First rule in PhaseOne order keeps the yield;
the later rule is logged at INFO under `RuleSet.metabolize` with
`unique_csmi` on (no Warning emission). This is overlapping coverage, not a
unique-edit / `SiteDeduplicationWarning` bug.

No pair is a full-rule perfect subset. Closest near-subset: OxygenReduction
carbonyl applicability ⊆ Hydrogenation near that region (peroxide chemistry
stays exclusive to OxygenReduction).

Status of treating these as permanent expected overlap: **not decided** —
record only. Do not invent a “less likely” schema field to silence them
(see `.cursor/rules/data-not-branches.mdc`). Overlaps are INFO only, so they
do not fail the forest WAE gate; `SiteDeduplicationWarning` stays an error.

| # | Pair | Example substrate | Notes |
|---|------|-------------------|-------|
| 1 | Dehydration ↔ NitrogenReduction | nitrobenzene | Often shared leaving `O` |
| 2 | Dehydrogenation ↔ QuinoneFormation | hydroquinone | Same aromatization product |
| 3 | Hydrogenation ↔ OxygenReduction | acetophenone | OR carbonyl ⊆ H near-subset; peroxide exclusive to OR |
| 4 | Hydroxylation ↔ OxidativeDehalogenation | chlorobenzene | Phenol from either path |
| 5 | OxidativeDehalogenation ↔ ReductiveDehalogenation | chlorobenzene | Often leaving-group `Cl` only — chemically spurious overlap |
| 6 | Dehydration ↔ Hydrolysis | histidine | Seen under forest WAE |

Panel census snapshot (`artifacts/census_csmi_vs_redundant.out`): pairs 1–5
appear in the metabolize panel; pair 6 shows up on broader WAE runs.

Related: two dedup layers in `HEURISTICS.md` (site unique-edit vs product
CSMI); `DIVERGENCES.md` product-csmi vs archive.
