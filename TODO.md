# TODO

## Now (find_path + sampling)

- find_path_partial: like find_path when rules cannot reach target; report closest reach and what was missed (may replace find_path); DRY/clarity refactor of shared search so top-level find_path stays easy to follow
- find_path: mapping argument between target and reactant
- find_path: return graph alongside paths; emit graph even when path not found (crate has `product_graph`; public API still `(hits, counters)`)
- Collect hard metabolism-DB cases where find_path / find_path_partial do not fail fast; harden fail-fast (includes MeOPhOH seen-after-`mol_edits` residual)
- random_path: reactant + seed + ruleset (or PathOutcome) → randomly apply rules; return applied steps, final mol, PatternInfos (tests + sampling)
- Linearization.apply: richer apply-to-molecule output in the random_path shape (apply already returns mols)
- PathOutcome: randomly sample a linearization; apply a random linearization to a molecule (same shape as random_path)

## Next (conjugation / reactivity)

- Preserve text labels on stars through the pipeline; wrap chematics where labels drop and reapply
- Normalize conjugated/adduct mols (collapse glucuronide, GSH, etc. to `*` with normalized labels; standardize common/hinted text ids); wire into target-seeking
- Validate conjugation rules with examples + review; ship a prebuilt Reactivity ruleset (protein, DNA, cyanide, GSH)

## Later

- Ship a Rust CLI entry; retire temporary `xenosite.forest.native.cli:main`
- Optimize `alternating_paths` (path clone per BFS step; odd-ring / 2-colorable decomposition)
- PatternInfo coverage gaps still worth asserting (product-side SMARTS; `pin`; `skip_same_rings` / `edit`; `breaks_ring` + `partner`/`partner_h` from `resolve_effect`)
- `find_path_to_MS1` / then `find_path_to_MS2` (isotope-aware m/z within tolerance; MS fragmentation rules)
- Rank PathOutcome hits by likelihood via `xenosite-predict` (cleavage-side bags; no side-reaction enum)
- Deferred: tautomer SMARTS / TautomerQuery (Status: not decided) — orthogonal to `TautomerRule` stub
- SMARTS least-common atoms first
