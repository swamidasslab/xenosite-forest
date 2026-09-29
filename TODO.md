# TODO

All items below target the **Rust** door (`crates/xenosite-forest` +
`xenosite.forest` wrapper). Do not implement them in `xenosite.forest.native`
— see `docs/forest/NATIVE.md`. No CLI work for now (no Rust CLI; do not grow
the temporary native console script).

## Now (find_path + sampling)

- find_path_partial: like find_path when rules cannot reach target; report closest reach and what was missed (may replace find_path); Expand is now one DeferredSite stream — next, peel child-enqueue / heap update out of the PathSearch loop so partial can share it
- find_path: mapping argument between target and reactant
- find_path: return graph alongside paths; emit graph even when path not found (crate has `product_graph`; public API still `(hits, counters)`)
- Collect hard metabolism-DB cases where find_path / find_path_partial do not fail fast; harden fail-fast (includes MeOPhOH seen-after-`mol_edits` residual)
- StepSequence: ordered container of steps that **acts like a Step** (composite). Same `apply` (intermediates + final, `PathwayOptions`, ensure-tags-never-overwrite) for length 1..n or `random_path`
- StepPlan (`Deps` in Rust): compact set of linearizations under constraints; `linearizations()` yields `StepSequence`s then `apply`
- PathOutcome: randomly sample a StepSequence from a plan; `apply` it on a molecule (same as above)

## Next (conjugation / reactivity)

- Preserve text labels on stars through the pipeline; wrap chematics where labels drop and reapply
- Normalize conjugated/adduct mols (collapse glucuronide, GSH, etc. to `*` with normalized labels; standardize common/hinted text ids); wire into target-seeking
- Validate conjugation rules with examples + review; ship a prebuilt Reactivity ruleset (protein, DNA, cyanide, GSH)

## Later

- Optimize `alternating_paths` (path clone per BFS step; odd-ring / 2-colorable decomposition)
- PatternInfo coverage gaps still worth asserting (product-side SMARTS; `pin`; `skip_same_rings` / `edit`; `breaks_ring` + `partner`/`partner_h` from `resolve_effect`)
- `find_path_to_MS1` / then `find_path_to_MS2` (isotope-aware m/z within tolerance; MS fragmentation rules)
- Rank PathOutcome hits by likelihood via `xenosite-predict` (cleavage-side bags; no side-reaction enum)
- Deferred: tautomer SMARTS / TautomerQuery (Status: not decided) — orthogonal to `TautomerRule` stub
- SMARTS least-common atoms first
