# TODO

All items below target the **Rust** door (`crates/xenosite-forest` +
`xenosite.forest` wrapper). Do not implement them in `xenosite.forest.native`
— see `docs/forest/NATIVE.md`. No CLI work for now.

## Now

- Re-enable Keep-H / Effect-formula materialize filters + parity harness
  formula gates (`allow_formula_delta_mismatch`, SiteDeduplicationWarning)
- **Do not remove** `inventory_fully_covered_in_sssom` xfail until the **next
  version bump** (map new Effect SMARTS-split PatternInfo names into
  `mappings/xmet-forest.sssom.tsv` as part of that bump)

## After parity + catalog

- Pair materialize: stop minting path-end iminium/sulfinic junk; collapse
  `accept_pair_product` → `accept_product` (HEURISTICS: soft failure today)
- Re-enable parity-chase soft filters: Rust Keep-H / Effect-formula drop in
  `materialize_pair_mols` + `DeferredSite::materialize_mols`; drop
  `allow_formula_delta_mismatch` and `SiteDeduplicationWarning` ignore on
  `test_rule_parity_fuzz.py` (conftest formula collector gate back on)
- Tautomerization: cut find_path bill on tacrine further; normalize_tautomer
  stays opt-in (default off)
- Validate walk-history residual bags on MetX hard misses
- find_path: mapping argument between target and reactant
- StepSequence / StepPlan / PathOutcome sampling doors
- XMET SSSOM / chemist-home lint; Conjugation catalog when Phase II settles

## Later

- Optimize `alternating_paths`
- PatternInfo coverage gaps (product-side SMARTS; `pin`; `breaks_ring`)
- `find_path_to_MS1` / MS2; rank hits via `xenosite-predict`
- Deferred: tautomer SMARTS / TautomerQuery (not decided)
- SMARTS least-common atoms first
