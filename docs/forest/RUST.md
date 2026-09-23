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

Aromaticity is chematic’s RDKit-parity engine (`apply_aromaticity_rdkit_parity_experimental`). Canonical SMILES is chematic’s, not RDKit’s (`C(C)O` vs `CCO`). Tests compare `canon_of` identities, not a spelling. Kekulé parents stamp integer orders onto **one** conjugated system; other systems stay aromatic.

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

**RuleSet / closures:** Python `FilterRules` / `FilterSites` are `Callable`. Rust uses `impl Fn` on `RuleSet::metabolize`, or `BoxedFilters` (`Box<dyn Fn>`) when a search stores them. Capture-by-move; the boxed form is `'static`. Built-in filters should be Rust functions that read `PatternInfo` / `Effect` (methide is an effect field). A Python lambda still crosses the GIL and marshals a `PatternInfo` per call.

Compose the set in Python once (`RuleSet([PatternInfo(...), ...])` or `RuleSet.compose([hydroxylation, dealkylation])`). That copies pattern data into the Rust payload. Later `metabolize(mol)` passes handles only.

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
after. Visit order covers write/parse without isotope probes.

Apply does not yet return `src_to_new`; with tags on the atom that map is
optional. `fragments` still drops its private `old_to_new`; tags survive the
clone into each fragment.

Do **not** call private chematic writers or `build_product` maps, and do not
commit a patched submodule tree — only the pin SHA and the patch file.

## Not in this crate

Full `find_path`, every Phase I rule, production atom-trace on `Atom.tag`.
A `RuleSet` of `PatternInfo` plus closures is in the crate as a door; it is
not the live Python `RuleSet` / `find_path` filters.
