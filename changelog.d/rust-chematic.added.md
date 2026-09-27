A WASM-clean Rust derisk crate (`crates/xenosite-forest`) on chematic and canonaut: SMARTS/SMIRKS, unique-edit, pair orbits, per-system on-demand kekulé cache shared across relatives, `ForestMol` owning caches (no `_forest`/`xf` hack), atom tags as a sidecar remapped through apply and SMILES write (not `atom_map`, not `canonical_atom_order`), PyO3 class wrap, `RuleSet` / `PatternInfo` with `FilterRules` closures, `wasm32-unknown-unknown`, and a size-tuned native wheel.

Depends on crates.io `chematic` **1.0.27+** for molecule tags (`set_tag` /
`atom_tag`) and SMILES visit-order helpers (`write_with_atom_order` /
`canonical_smiles_with_atom_order`). Forest `Tag` stays 0-based; sync uses
`tag.0 + 1`. See `docs/forest/RUST.md` § Atom identity.

`AtomTracker` follows atoms via chematic molecule tags through SMIRKS apply
and SMILES write/parse (no isotope probe required).

**Pair door:** ResonancePair endpoints (`Edit::PairEndpoint`) metabolize via
alternating-path flip on Kekulé forms. Hydroquinone → quinone and benzene
QF `add_carbonyl_o`×2 are covered. SMIRKS apply uses a Kekulé
`reactant_parent` when maps 1–2 are aromatic.
