# Rust derisk crate

`crates/xenosite-forest` is a **WASM-clean** chemistry door for a future port of live `xenosite.forest`. It is not the Python package.

## Why this crate exists

The Python engine is RDKit + pynauty. This crate checks the seams that would block a full port:

- chematic SMILES / RDKit-parity aromaticity / SMARTS maps / one-match SMIRKS
- unique-edit (topological ranks)
- canonaut pair orbits (benzene meta ≠ para)
- kekulé parents for ResonanceRule
- valence gate (two-double nitrogen)
- hydroquinone pair-path edit
- `ForestMol` owns the chematic mol plus caches (no `_forest` / `xf` facade)
- kekulé parents: one conjugated system, on demand, shared `Rc` across relatives; an edit that changes a system’s shape starts a new bag
- PyO3 `#[pyclass]` wrap of `ForestMol` (native) and wasm-bindgen JS class (WASM)
- `RuleSet` of `PatternInfo` records, with `FilterRules` / `FilterSites` as `impl Fn` or `Box<dyn Fn>`
- `wasm32-unknown-unknown` (no C FFI: no `inchi` native, no canonaut `c-nauty-bench`)

Aromaticity is chematic’s RDKit-parity engine (`apply_aromaticity_rdkit_parity_experimental`). Canonical SMILES is chematic’s, not RDKit’s (`C(C)O` vs `CCO`). Tests compare `canon_of` identities, not a spelling. Kekulé parents stamp integer orders onto **one** conjugated system; other systems stay aromatic. Chematic SMARTS treats bare `h` / `;h` as inert — digitize like Hydroxylation (`[#6h1,#6h2,#6h3]` not `[#6h]`).

## `ForestMol`

Python hangs `_forest` on a foreign RDKit `Mol` and mints `xf` on every read. Rust does not need that. [`ForestMol`](../../crates/xenosite-forest/src/forest_mol.rs) **is** the object:

- Owns a chematic `Molecule`.
- Structure answers (`csmi`, formula, ranks, SMARTS) live on the object as `Rc<RefCell<Structure>>`, filled on first read.
- Kekulé assignments live on the same object as `Rc<RefCell<KekuleCache>>`, keyed by system atom set plus a fingerprint of aromatic/bond shape.
- Tags (`Tag`) are a sidecar parallel to atom index, not `atom_map` and not on chematic `Atom`. `from_apply` remaps them when a correspondence is supplied. Chematic's public apply/write do not return that map (see [Atom identity](#atom-identity-chematic-has-no-public-correspondence)).
- `copy_mol` shares both caches. `edit_copy` / `product` (after an edit) start a **new** structure bag and **keep** the kekulé `Rc`. The first relative to fill a system shares it with every relative whose key still matches. An edit that changes kekulization of a system is a new key and an empty bag.
- `clear_structure` drops `csmi`/formula/SMARTS, not the kekulé `Rc`.

There is no `has_forest`, no `wipe_forest`, and no `xf()` facade.

## Wrapping the class (PyO3 and wasm-bindgen)

The Rust `ForestMol` is the payload. Bindings do **not** reimplement the cache.

**PyO3** (`--features python`, native only — CPython C-API, not WASM):

```rust
#[pyclass(name = "ForestMol", unsendable)]
pub struct PyForestMol {
    inner: ForestMol,              // the payload
    csmi: Option<Py<PyString>>,    // interned so `mol.csmi is mol.csmi`
}
```

- `#[pyclass]` makes a Python type whose instance **is** that struct.
- `#[new]` is `__init__`. `#[getter]` is a Python property. `#[pymethods]` are methods.
- `unsendable`: `RefCell` cache + chematic mol are not `Send`.
- Nested wrap: `Formula` is its own `#[pyclass]`.
- Returning a fresh `String` would allocate a new Python `str` each time and break `is`. The wrapper intern/caches `Py<PyString>`.

```bash
cargo test -p xenosite-forest --features python
```

`python3-dev` (libpython) is required to link the test binary. `#[pyclass]` is CPython; it is not WASM.

**RuleSet / candidates:** Primary walk is [`RuleSet::candidates`](../../crates/xenosite-forest/src/ruleset.rs)
— site–pattern–[`ParentRef`](../../crates/xenosite-forest/src/candidate.rs) triples
without applying edits. A search reads `PatternInfo` / `Effect` to filter, then
[`Candidate::materialize`](../../crates/xenosite-forest/src/candidate.rs) only for
survivors. Filter closures on `metabolize` remain for Python parity; they are
not required. `metabolites(mol)` materializes everything. ResonancePair uses
[`PairCandidate`](../../crates/xenosite-forest/src/pair_edit.rs) the same way
(merged `Effect` for filtering; path flip on materialize).

**RuleSet / closures (optional):** Python `FilterRules` / `FilterSites` are
`Callable`. Rust can still pass `impl Fn` on `RuleSet::metabolize`, or
`BoxedFilters` (`Box<dyn Fn>`). Prefer filtering candidates by reading data.

Compose the set in Python once (`RuleSet([PatternInfo(...), ...])` or `RuleSet.compose([hydroxylation, dealkylation])`). That copies pattern data into the Rust payload. Later `metabolize(mol)` / `candidates(mol)` passes handles only.

**RuleSet namespaces:** nested sets stay nested (`compose` does not flatten). Each emission carries a leaf-first `rule_path` (emitting rule, then each containing set), matching Python `info["rule"]`. Children run with `unique_csmi=false` so alternate rules bubble; the caller's `unique_csmi` is the cross-child CSMI layer. Filters see the leaf set, not the outer compose container.

**Rule catalog:** every concrete Python reaction rule is a leaf `RuleSet` in
[`rules.rs`](../../crates/xenosite-forest/src/rules.rs) (`phase_one`,
`default_ruleset`, `all_rules`). Pair-endpoint patterns metabolize through
the ResonancePair door (`Edit::PairEndpoint` + path flip). Atom/bond SMIRKS
use Kekulé `reactant_parent` when the matched bond is aromatic.

```python
rs = RuleSet([
    PatternInfo("h", "[#6h1:1]", "hydroxyl", "O", "H"),
    PatternInfo("h2", "[#6h2,#6h3:1]", "hydroxyl", "O", "H"),
], "Hydroxylation")
rs.metabolize(ForestMol("c1ccccc1"))  # no per-search marshal of patterns
```

`python` and `wasm` are mutually exclusive (both want the `cdylib`).

**wasm-bindgen** (`--features wasm`): the same payload, JS classes named `ForestMol` and `RuleSet`. JS does not get `filter_rules` callbacks on this door; compose patterns in Rust/JS data, then `metabolize`.


## Build

Native tests:

```bash
cargo test -p xenosite-forest
cargo test -p xenosite-forest --features python
```

Browser WASM (the whole rlib + cdylib, including hydroxylation and the JS `ForestMol` class):

```bash
rustup target add wasm32-unknown-unknown
cargo build -p xenosite-forest --target wasm32-unknown-unknown --features wasm --release
```

`.cargo/config.toml` sets `getrandom_backend="wasm_js"` so canonaut’s `rand` compiles on `wasm32-unknown-unknown`. The crate does not call C.

## WASM size

`[profile.wasm-size]` in the workspace `Cargo.toml` is opt-level `z`, fat LTO, one codegen unit, stripped symbols.

Measured `xenosite_forest.wasm` (`--features wasm`):

| Build | Raw | gzip -9 |
| ----- | --- | ------- |
| debug | 39 MiB | — |
| `--release` | 1.07 MiB (1 126 646 B) | 300 KiB |
| `--profile wasm-size` | 782 KiB (801 241 B) | 217 KiB |
| then `wasm-opt -Oz` | 645 KiB (660 208 B) | 213 KiB |

```bash
cargo build -p xenosite-forest --target wasm32-unknown-unknown --features wasm --profile wasm-size
wasm-opt -Oz --enable-bulk-memory --enable-sign-ext --enable-mutable-globals \
  --enable-nontrapping-float-to-int --enable-multivalue --enable-reference-types \
  target/wasm32-unknown-unknown/wasm-size/xenosite_forest.wasm -o xenosite_forest.wasm
```

`wasm-opt` is optional (binaryen). gzip of the cargo-only wasm-size artifact is already close to gzip of the wasm-opt output.

## Wheel size

The PyPI package is still pure Python (hatchling + RDKit). The native door wheel is a separate maturin artifact from `crates/xenosite-forest` (`xenosite-forest-native`, import `xenosite_forest`). Same size profile as WASM.

Measured CPython 3.12 linux x86_64 (`--features python,extension-module`):

| Build | Wheel (zip) | `.so` uncompressed |
| ----- | ----------- | ------------------ |
| `--release` | 576 KiB (589 373 B) | 1.37 MiB |
| `--profile wasm-size` | 347 KiB (355 072 B) | 741 KiB |

```bash
cd crates/xenosite-forest
maturin build --profile wasm-size
```

The size-profile wheel is already zip-compressed, so compare it to gzip WASM (~217 KiB), not to the raw `.wasm`. Native x86_64 + PyO3 is a bit larger than the WASM module; it is still well under a megabyte. This is not a manylinux2014 auditwheel rebuild — tags follow the builder (`manylinux_2_34` here).

## Speed (measured on this machine)

RDKit is already C++. Chematic is not a free 10× on SMARTS/SMIRKS. The gain is dropping Python around unique-edit, nauty glue, forest copy, and `find_path`.

Same five molecules, release Rust vs live Python/RDKit/pynauty (this VM):

| Step | Rust | Python | Ratio |
| ---- | ---- | ------ | ----- |
| Parse benzene | 1.7 µs | 15 µs | ~9× |
| Uncached `csmi` | 8–17 µs | 9–21 µs | ~1× |
| SMARTS `[#6h1:1]` | 0.9–1.4 µs | 0.34 µs | RDKit **~3× faster** |
| Hydroxylation (all unique sites) | 13–154 µs | 24–100 µs `RunReactants` | ~1×; see below |
| Full `Hydroxylation.metabolize` | (not in crate) | 0.5–4.3 ms | Python+trace+copy |
| Nauty generators / unordered pair orbits | 7–30 µs | 53–107 µs gens; 0.2–1.0 ms six-family | ~4–6× gens |
| Anisole dealkylation apply | 8.3 µs | 6.9 µs | RDKit slightly faster |

The 13–154 µs hydroxylation range is **how many unique carbons**, not a slower edit. Benzene and hydroquinone have 1 site (~12 µs). Phenol 3, anisole 4, ibuprofen 7. Each product pays `canon_smiles` (write, reparse, aromatize, write again). Per product is ~11–14 µs on the small aromatics and ~22 µs on ibuprofen. Unique-edit itself stays 3–7 µs.

`find_path` hits today are **8–200 ms** ([PERFORMANCE.md](PERFORMANCE.md)). Cheap one-step paths are already RDKit-bound; a chematic door alone would not make those 10×. MeOPhOH-style bills spend that time in **many** Python `metabolize` / unique-edit / nauty calls — that is where a full Rust port would show.

Estimate, not a promise:

- **Door only (PyO3, Python still runs `find_path`)**: ~1–2×. FFI still pays per call; SMARTS may regress.
- **Full native `find_path` + unique-edit**: typically **~3–5×**, more like **5–10×** when pair-orbit unique-edit dominates, little when a single RDKit reaction is the whole path.
- **BFS-style enum** (Python loop over hundreds of `metabolize`): **~10×** is plausible; live search already avoided that factorial blow-up.
- **WASM in the browser**: native-rust numbers times ~1.5–3×, plus JIT warmup. Still faster than Python loops; not faster than RDKit VF2.

Reproduce the Rust column: `cargo run -p xenosite-forest --example door_bench --release`.

## Atom identity (vendored chematic patch)

Chematic on crates.io still has no atom userdata and no public SMILES visit
order. This repo vendors `chematic` @ `v1.0.21` as a **sparse submodule**
(`vendor/chematic`) and applies
[`patches/chematic-v1.0.21-atom-tag-visit-order.patch`](../../patches/chematic-v1.0.21-atom-tag-visit-order.patch)
via `./scripts/vendor-chematic.sh`.

The patch adds:

- `Atom.tag: Option<u32>` — non-chemical; copied by clone / apply / fragments;
  ignored by SMILES write and canon; not emitted as `:n`. `atom_map` is still
  cleared on apply.
- `write_with_order` / `canonical_smiles_with_order` — string plus DFS visit
  order (`order[k]` = mol index of the k-th atom in the string).

`ForestMol` still keeps a tag sidecar today (derisk continuity). Chematic
`Atom.tag` is the path for production atom-trace: stamp before apply, read
after. Visit order covers write/parse without isotope probes. POC:
[`AtomTracker`](../../crates/xenosite-forest/src/atom_tracker.rs)
(`stamp` / `src_to_new` / `adopt_born` / `write_parse`).

Apply does not yet return `src_to_new`; with tags on the atom that map is
optional. `fragments` still drops its private `old_to_new`; tags survive the
clone into each fragment.

Do **not** call private chematic writers or `build_product` maps, and do not
commit a patched submodule tree — only the pin SHA and the patch file.

## Not in this crate

Full Python `atom_diff` parity (leave_count / bond_raises site filters on every
edge case), archive `Deps` apply. A provisional
[`atom_diff`](../../crates/xenosite-forest/src/atom_diff.rs) door gates
candidates via MCS + effect fields (`find_path_diff` / `use_atom_diff`).
[`find_path`](../../crates/xenosite-forest/src/find_path.rs) walks carry tagged
[`ForestMol`](../../crates/xenosite-forest/src/forest_mol.rs) (structure/`csmi`
cache; `Atom.tag` synced from the sidecar through `adopt_product`). Outcomes
emit elementary [`CanonicalStep`](../../crates/xenosite-forest/src/canonical_plan.rs)
plans; quinone-shaped leaves set [`PlanKind::HydroxylationThenDehydrogenation`]
on the `RuleSet` (data, not a name branch). Default lazy (and eager) closer
[`try_atom_diff_for_child`](../../crates/xenosite-forest/src/atom_diff.rs) /
[`atom_diff_for_child`](../../crates/xenosite-forest/src/atom_diff.rs) via
tag-lift (+ greedy OH extend); full MCS only when lift fails. Nested `RuleSet`
namespaces are in the crate as a door; they are not yet the live Python search.

## find_path H2H (Rust vs Python live)

Same PhaseOne cases / `max_nodes=800` / `max_paths=1` / best-of-5 as
[`bench_find_path_h2h.py`](../../tests/forest/bench_find_path_h2h.py).

```bash
cargo run -p xenosite-forest --example find_path_bench --release
cargo run -p xenosite-forest --example find_path_bench --release -- --larger
cargo run -p xenosite-forest --example find_path_bench --release -- --hard
# unfiltered (slow) is opt-in:
cargo run -p xenosite-forest --example find_path_bench --release -- --hard --nofilter
cargo run -p xenosite-forest --example find_path_profile --release
uv run python tests/forest/bench_find_path_rust_h2h.py
uv run python tests/forest/bench_find_path_rust_h2h.py --larger
uv run python tests/forest/bench_find_path_rust_h2h.py --hard
```

Default bench is **filter-only** (atom_diff on). `--budget-secs N` stops remaining
rows once the suite wall exceeds N (default 30s filter / 60s with `--nofilter`).
`find_path_profile` instruments the tagged `ForestMol` walk (not a CSMI-string
re-parse door).

Raw: `artifacts/bench_find_path_{rust,python}_{mid,larger,hard}.out`.

Filter parity (`BondCompare::Any` + `match_bonds=false` MCS, multi-placement
merge, `leave_count` on `Effect`, site/pattern/pair-end could-help, QF
`methide_end` + `partner: C`, cleavage-first order, atom_diff `closer`)
closes the mid-size miss gap.

### Mid-size / multi-edit

| Case | Python live | Rust atom_diff | Rust no filter |
| --- | --- | --- | --- |
| eugenol→allyl-Q | **ok** 0.030s · bill=20 | **ok** 0.008s · bill=16 | MISS · bill≈7k |
| dimethoxy-PEA→catechol | **ok** 0.011s · bill=10 | **ok** 0.003s · bill=12 | **ok** 0.039s · bill=391 |
| MeOPhOH→hydroxyQ | **ok** 0.31s · bill=226 | **ok** 0.080s · bill=315 | MISS · bill≈29k |
| TBA→aldehyde | **ok** 0.018s · bill=5 | **ok** 0.014s · bill=5 | **ok** 0.19s · bill=168 |
| 2-MeO-naph→1,2-NQ | **ok** 0.010s · bill=7 | **ok** 0.008s · bill=7 | MISS · bill≈48k |
| **TOTAL wall** | **0.38s (5/5)** | **0.11s (5/5)** | **8s (2/5)** |

### Larger (HA≈17–26)

| Case | Python live | Rust atom_diff | Rust no filter |
| --- | --- | --- | --- |
| tBu-bis-ND→dialdehyde | **ok** 0.086s · bill=36 | **ok** 0.027s · bill=33 | **ok** 0.20s · bill=294 |
| macrocycle-ND→aminoK | **ok** 0.16s · bill=26 | **ok** 0.079s · bill=37 | **ok** 1.16s · bill=1109 |
| tribenzyl→PhCHO | **ok** 0.006s · bill=4 | **ok** 0.026s · bill=5 | **ok** 0.17s · bill=84 |
| triPh-butyl→OH | **ok** 0.025s · bill=5 | **ok** 0.031s · bill=5 | **ok** 0.24s · bill=122 |
| MeO-diphenyl→catechol | **ok** 0.042s · bill=7 | **ok** 0.018s · bill=9 | **ok** 0.15s · bill=150 |
| **TOTAL wall** | **0.32s (5/5)** | **0.18s (5/5)** | **1.92s (5/5)** |

### Hard (HA≈14–20, ≥3–8 hops)

Multi-dealkylation + hydroxylation + DH/QF on bigger rings. Without filters,
Rust misses most rows under `max_nodes=800`.

| Case | Python live | Rust atom_diff | Rust no filter |
| --- | --- | --- | --- |
| trimethoxy-PEA→catechol | **ok** 0.12s · 5 steps · bill=60 | **ok** 0.035s · 5 steps · bill=57 | MISS · bill≈43k |
| eugenol-MeO→allylQ | **ok** 0.067s · 5 steps · bill=39 | **ok** 0.020s · 4 steps · bill=40 | **ok** 1.6s · bill≈17k |
| bisMeO-naph→1,2NQ | **ok** 0.062s · 5 steps · bill=27 | **ok** 0.054s · 4 steps · bill=34 | MISS · bill≈47k |
| tetraMeO-biphenyl→tetraOH | **ok** 0.13s · 4 steps · bill=52 | **ok** 0.075s · 4 steps · bill=54 | **ok** 16s · bill≈37k |
| veratrole-allyl→allylQ | **ok** 1.16s · 7 steps · bill=873 | **ok** 0.045s · 6 steps · bill=105 | MISS · bill≈41k |
| tetraMeO-naph→polyOH-NQ | **ok** 5.33s · 8 steps · bill=1311 | **ok** 1.37s · 7 steps · bill=767 | MISS · bill≈54k |
| **TOTAL wall** | **6.86s (6/6)** | **1.60s (6/6)** | **55s (2/6)** |

**Read:** with filters on, Rust hits every H2H row and is **faster than Python** on
mid (0.11s vs 0.38s), larger, and hard (1.60s vs 6.86s). Walks carry tagged
`ForestMol` + `CanonicalStep` plans; default lazy closer tag-lifts MCS at
enqueue (full MCS only on lift miss / root).

### Profile: where hard wall goes (not copies)

```bash
cargo run -p xenosite-forest --example find_path_profile --release
```

Eager closer (MCS on every emit) spent **~62–66%** of hard wall on per-product
child `atom_diff`. Unit: `Molecule.clone` ≈0.3 µs; `atom_diff` ≈21 ms —
clone elision is noise.

### Lazy closer (landed)

`FindPathConfig::lazy_closer` (default on with `use_atom_diff`): enqueue kept
fragments with tag-lifted child diff when possible; on pop verify
`cost() < parent_cost` before expand (full MCS only if lift missed). Same gate
as eager, deferred so siblings never popped skip work.

Hard H2H (release, best-of-5), atom_diff on:

| | eager closer | lazy closer |
| --- | ---: | ---: |
| TOTAL wall | 1.97s | **1.60s** |
| tetraMeO-naph | 1.72s · bill=806 | **1.37s · bill=767** |
| veratrole-allyl | 0.051s · bill=115 | **0.045s · bill=105** |

Mid atom_diff TOTAL **0.11s** (was ~0.23s). Do not HA-gate at enqueue —
oxidation can raise HA distance while lowering cost.

### Tagged walks + step plans (landed)

`find_path` walks carry [`ForestMol`](../../crates/xenosite-forest/src/forest_mol.rs)
(not bare SMILES). Edits use `materialize_mols` → `adopt_product` so tags and
the structure/`csmi` cache stay on the walk. `PathOutcome.plan` is a
`Vec<CanonicalStep>`; quinone leaves expand via `RuleSet::plan_kind`. Default
lazy closer (and eager) uses `try_atom_diff_for_child` / `atom_diff_for_child`
(tag-lift + OH extend) instead of a fresh MCS when tags allow; full MCS only
on lift miss.
