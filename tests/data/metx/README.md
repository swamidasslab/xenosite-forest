# MetXBioDB hard-case bins

Provenance: **MetXBioDB Phase I** reactant→product pairs
(`artifacts/metx_phase1_pairs.tsv`), scored with nostereo call-site
`find_path_partial` (`metx_hard_cases` example).

Each TSV starts with `# commit=` / `# scan=` so you know which tree produced
the labels. Regenerate:

```bash
cargo run -p xenosite-forest --example metx_curate_bins --release -- \
  artifacts/metx_hard_misses_nostereo.tsv \
  artifacts/metx_phase1_pairs.tsv \
  artifacts/metx_hard_misses.tsv \
  tests/data/metx
```

| File | Bin |
|------|-----|
| [`metx_thrash_rules_cover.tsv`](metx_thrash_rules_cover.tsv) | Thrash/framing but rules should cover |
| [`metx_thrash_rules_gap.tsv`](metx_thrash_rules_gap.tsv) | Thrash; DB saturation; rules cannot cover |
| [`metx_quiet_rules_gap.tsv`](metx_quiet_rules_gap.tsv) | No thrash; rule/chemistry gap |
| [`metx_near_miss_progress.tsv`](metx_near_miss_progress.tsv) | Missed goal; large residual drop + `path_to_closest` |
| [`metx_hard_case_bins.tsv`](metx_hard_case_bins.tsv) | Combined |
