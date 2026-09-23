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
- `_forest` cache / `xf` facade (same-object answers; `copy_mol` shares cache; `rw_copy` does not)
- PyO3 `#[pyclass]` wrap of `ForestMol` (native) and wasm-bindgen JS class (WASM)
- `wasm32-unknown-unknown` (no C FFI: no `inchi` native, no canonaut `c-nauty-bench`)

Aromaticity is chematic’s RDKit-parity engine (`apply_aromaticity_rdkit_parity_experimental`). Canonical SMILES is chematic’s, not RDKit’s (`C(C)O` vs `CCO`). Tests compare `canon_of` identities, not a spelling. Kekulé parents are integer-order graphs with aromatic atom flags cleared so aliphatic SMIRKS can match.

## Cache (`ForestMol` / `xf`)

Python stores answers on `mol._forest["cache"]` and mints a new `Xf` on every `mol.xf` read ([`docs/forest/XF.md`](XF.md)). The Rust door copies that:

- `ForestMol::parse` does **not** install a forest (`has_forest` is false).
- `mol.xf()` is a facade. Answers live on `Forest` (`cache` is `Rc<RefCell<Structure>>`).
- First `csmi` / `formula` / `topol_equiv` / `smarts_matches` fills the structure bag; later reads return the same `Rc`.
- `copy_mol` keeps `cache` by identity. `rw_copy` / `from_molecule` carry no forest (the copy is about to be edited).
- `clear_structure` drops the bag; labels would stay (trace not in this crate yet). `wipe_forest` removes the forest.

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

**wasm-bindgen** (`--features wasm`): the same payload, a JS class named `ForestMol`. `python` and `wasm` are mutually exclusive (both want the `cdylib`).


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

## Not in this crate

Full `find_path`, every Phase I rule, atom-trace, RuleSet. Those wait on these tests staying green.
