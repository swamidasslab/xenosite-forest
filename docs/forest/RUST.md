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
- `wasm32-unknown-unknown` (no C FFI: no `inchi` native, no canonaut `c-nauty-bench`)

Aromaticity is chematic’s RDKit-parity engine (`apply_aromaticity_rdkit_parity_experimental`). Canonical SMILES is chematic’s, not RDKit’s (`C(C)O` vs `CCO`). Tests compare `canon_of` identities, not a spelling. Kekulé parents are integer-order graphs with aromatic atom flags cleared so aliphatic SMIRKS can match.

## Build

Native tests:

```bash
cargo test -p xenosite-forest
```

Browser WASM (the whole rlib + cdylib, including hydroxylation):

```bash
rustup target add wasm32-unknown-unknown
cargo build -p xenosite-forest --target wasm32-unknown-unknown --features wasm --release
```

`.cargo/config.toml` sets `getrandom_backend="wasm_js"` so canonaut’s `rand` compiles on `wasm32-unknown-unknown`. The crate does not call C.

## Not in this crate

Full `find_path`, every Phase I rule, atom-trace, RuleSet, Python bindings. Those wait on these tests staying green.
