//! Metabolic Forest chemistry door on chematic + canonaut.
//!
//! Derisk crate for a full Rust port. The algorithm lives in the Python
//! package; this crate proves the chemistry seams compile and behave on
//! chematic, including `wasm32-unknown-unknown`.

#[cfg(all(feature = "python", feature = "wasm"))]
compile_error!(
    "features `python` and `wasm` are mutually exclusive (CPython vs wasm-bindgen cdylib)"
);

pub mod forest;
pub mod hydroxylation;
pub mod kekule;
pub mod mol;
pub mod orbits;
pub mod pair_edit;
pub mod smarts;
pub mod smirks;
pub mod unique_edit;
pub mod valence;
pub mod xf;

#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python_api;

#[cfg(feature = "wasm")]
mod wasm_api;

pub use forest::{Forest, Formula, Structure, empty_forest, forest_copy};
pub use hydroxylation::hydroxylate;
pub use mol::{ForestError, Molecule, canon_of, canon_smiles, parse_mol, ranks};
pub use orbits::{atom_pair_orbit_id, unordered_atom_pair_orbit_sizes};
pub use pair_edit::dehydrogenate_hydroquinone;
pub use smarts::smarts_matches;
pub use smirks::apply_smirks_at;
pub use unique_edit::unique_atom_sites;
pub use valence::{accept_product, nitrogen_two_doubles};
pub use xf::{ForestMol, Xf};
