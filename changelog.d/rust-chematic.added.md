A WASM-clean Rust derisk crate (`crates/xenosite-forest`) on chematic and canonaut: SMARTS/SMIRKS, unique-edit, pair orbits, per-system on-demand kekulé cache shared across relatives, `ForestMol` owning caches (no `_forest`/`xf` hack), atom tags as a sidecar remapped through apply and SMILES write (not `atom_map`, not `canonical_atom_order`), PyO3 class wrap, `RuleSet` / `PatternInfo` with `FilterRules` closures, `wasm32-unknown-unknown`, and a size-tuned native wheel.

Chematic is vendored as a sparse submodule at `v1.0.21` with a tiny local patch (`Atom.tag`, `write_with_order` / `canonical_smiles_with_order`). See `patches/README.md` and `./scripts/vendor-chematic.sh`.

`AtomTracker` POC follows atoms via chematic `Atom.tag` through SMIRKS apply and SMILES write/parse (no isotope probe, no sidecar).
