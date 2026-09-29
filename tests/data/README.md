# Test data libraries (`tests/data/`)

Shared probe / case lists for forest tests. Prefer a file here over burying
SMILES tuples inside a single test module. Rust may `include_str!` these paths;
Python loaders live next to the tests that need them (thin wrappers only).

## Files

| File | Why it matters |
|------|----------------|
| [`coverage_substrates.txt`](coverage_substrates.txt) | **PatternInfo / When / Effect coverage + native↔Rust product parity.** Sections `[library]` (phase-I / formula / conjugation probes) and `[pattern]` (rare OR/When fillers). Expand `[pattern]` when a catalog arm stays mute. Do **not** grow per-leaf `LEAF_EXAMPLE_SUBSTRATES` / `_example_substrates` for this. Loaders: `crates/…/substrate_library.rs`, `tests/forest/native/substrate_library.py`. |
| [`quick_substrates.txt`](quick_substrates.txt) | **Small parametrize slice** of the coverage library for cheap meta-tests (`test_canonical_plan`, etc.). Every line must also appear in `coverage_substrates.txt` `[library]`. |
| [`phase1_smarts_probes.txt`](phase1_smarts_probes.txt) | **Collective SMARTS-family hit list** for QuinoneFormation / Dehydrogenation (and related Phase I) query coverage. Named rows (`id\tsmiles`). Used historically by legacy `test_phase1_steps` / formula-hint suites. |
| [`conjugation_probes.txt`](conjugation_probes.txt) | **Conjugation SMARTS site probes** — one SMILES per reactive class (phenol, acid, epoxide, Michael acceptor, …). Test wiring still binds rule factory + rxn index; this file owns the molecule list so it is not duplicated across legacy/native conjugate tests. |
| [`bfs_fuzz_corpus.txt`](bfs_fuzz_corpus.txt) | **Drug-like RDKit crashers + suite anchors** for bounded BFS/DFS fuzz (no crash, no dotted product). Keep when a new field crash appears. |
| [`find_path_bench_cases.txt`](find_path_bench_cases.txt) | **find_path wall / bill benches** — `label\treactant\ttarget` rows, sections `[mid]` and `[larger]`. Shared by native + rust `bench_find_path_h2h` / profile scripts so doors do not drift apart. |

## Not moved here (on purpose)

| Location | Why it stays |
|----------|----------------|
| `LEAF_EXAMPLE_SUBSTRATES` in `rules.rs` / native `_example_substrates` | **Rule metadata** for site_kind / emit smoke (`seal_leaf`), not a coverage pool. Keep short. |
| Per-test goldens (`test_*_golden_products.py`, `test_guided_path_gold.py`, `gold_csmi.rs`) | Locked **product expectations** for one chemistry story, not a reusable probe library. |
| Conjugation probe **rows** `(id, rule, smiles, rxn_index)` | Rule/rxn binding is test logic; only the SMILES set is shared data. |
| Hypothesis strategies / inline 2–4 SMILES in a unit test | Ephemeral fixtures; no second SoT. |

## Rules

1. One SoT per concern — do not fork a second copy in a test module.
2. Header comment on every file: purpose, who loads it, when to expand.
3. After moving a list, delete the buried constant and load from here.
4. Coverage mute → expand `coverage_substrates.txt` `[pattern]`, not leaf examples.
