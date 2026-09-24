//! Metabolic Forest chemistry door on chematic + canonaut.
//!
//! Derisk crate for a full Rust port. The algorithm lives in the Python
//! package; this crate proves the chemistry seams compile and behave on
//! chematic, including `wasm32-unknown-unknown`.

#[cfg(all(feature = "python", feature = "wasm"))]
compile_error!(
    "features `python` and `wasm` are mutually exclusive (CPython vs wasm-bindgen cdylib)"
);

pub mod atom_diff;
pub mod atom_tracker;
pub mod candidate;
pub mod chematic_vendor;
pub mod find_path;
pub mod forest;
pub mod forest_mol;
pub mod hydroxylation;
pub mod kekule;
pub mod labels;
pub mod mol;
pub mod orbits;
pub mod pair_edit;
pub mod pattern;
pub mod rules;
pub mod ruleset;
pub mod smarts;
pub mod smirks;
pub mod unique_edit;
pub mod valence;

#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python_api;

#[cfg(feature = "wasm")]
mod wasm_api;

pub use atom_diff::{
    AtomDiff, atom_diff, candidate_could_help, candidate_could_help_on, candidate_order_key,
    keep_against_diff, pair_could_help, pattern_could_help, pattern_could_help_mol,
};
pub use atom_tracker::{AtomTracker, tags_agree_elements};
pub use candidate::{Candidate, ParentRef};
pub use find_path::{
    FindPathConfig, PathCounters, PathOutcome, PathStep, find_path, find_path_default,
    find_path_diff, find_path_with, find_path_with_filters,
};
pub use forest::{Formula, Structure, molecule_formula};
pub use forest_mol::ForestMol;
pub use hydroxylation::{hydroxylate, hydroxylation};
pub use labels::Tag;
pub use mol::{ForestError, Molecule, canon_of, canon_smiles, parse_mol, ranks};
pub use orbits::{atom_pair_orbit_id, unordered_atom_pair_orbit_sizes};
pub use pair_edit::{PairCandidate, dehydrogenate_hydroquinone};
pub use pattern::{Edit, Effect, Emission, PatternInfo, SiteInfo, SiteKind};
pub use rules::{all_rules, catalog_names, default_ruleset, phase_one};
pub use ruleset::{
    BoxedFilters, FilterRules, FilterSites, RuleMember, RuleSet, accept_all_rules,
    accept_all_sites, o_dealkylation,
};
pub use smarts::smarts_matches;
pub use smirks::apply_smirks_at;
pub use unique_edit::unique_atom_sites;
pub use valence::{accept_product, nitrogen_two_doubles};
