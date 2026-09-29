# TODO

All items below target the **Rust** door (`crates/xenosite-forest` +
`xenosite.forest` wrapper). Do not implement them in `xenosite.forest.native`
— see `docs/forest/NATIVE.md`. No CLI work for now (no Rust CLI; do not grow
the temporary native console script).

## Now (find_path + sampling)

- Tautomerization: cut find_path bill on tacrine further (fanout / unique-edit); imine↔amine coverage fixed; normalize_tautomer door shipped (chematic pick + Forest adopt)
- Validate walk-history residual bags on MetX hard misses (chromenone OH thrash, S-ox, arene-epoxide) + +GSH before coding invariant-leftover cuts
- Boost / `stop_after_sealed_basins` stay parked until bag history is measured (early-stop hurt closest)
- find_path: mapping argument between target and reactant
- StepSequence: ordered container of steps that **acts like a Step** (composite). Same `apply` (intermediates + final, `PathwayOptions`, ensure-tags-never-overwrite) for length 1..n or `random_path`
- StepPlan (`Deps` in Rust): compact set of linearizations under constraints; `linearizations()` yields `StepSequence`s then `apply`
- PathOutcome: randomly sample a StepSequence from a plan; `apply` it on a molecule (same as above)

## Next (conjugation / reactivity)

- Preserve text labels on stars through the pipeline; wrap chematics where labels drop and reapply
- Normalize conjugated/adduct mols (collapse glucuronide, GSH, etc. to `*` with normalized labels; standardize common/hinted text ids); wire into target-seeking
- Validate conjugation rules with examples + review; ship a prebuilt Reactivity ruleset (protein, DNA, cyanide, GSH)

## Mapping / XMET (upstream)

- Strip Forest/Rainbow identity language from 29 allowlisted chemist homes in tagger `xmet.yaml` (see `test_xmet_definition_lint.py`); drop IDs from allowlist as they clear
- Consider distinct chemist homes / always_with for `Tautomerization/tautomer_h` vs `path_partner` (both map to `xmet:4000186` today)
- Add a correct `Conjugation` catalog + SSSOM rows when Phase II composition is settled (CJ rows removed)
- On forest release: snapshot SSSOM into xenosite-xmet (`docs/release.md`)

## Later

- Optimize `alternating_paths` (path clone per BFS step; odd-ring / 2-colorable decomposition)
- PatternInfo coverage gaps still worth asserting (product-side SMARTS; `pin`; `skip_same_rings` / `edit`; `breaks_ring` + `partner`/`partner_h` from `resolve_effect`)
- `find_path_to_MS1` / then `find_path_to_MS2` (isotope-aware m/z within tolerance; MS fragmentation rules)
- Rank PathOutcome hits by likelihood via `xenosite-predict` (cleavage-side bags; no side-reaction enum)
- Deferred: tautomer SMARTS / TautomerQuery (Status: not decided) — orthogonal to Default `Tautomerization`
- SMARTS least-common atoms first
