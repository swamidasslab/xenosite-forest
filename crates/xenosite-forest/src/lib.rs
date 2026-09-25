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
pub mod canonical_plan;
pub mod chematic_vendor;
pub mod cleavage_graph;
pub mod find_path;
pub mod forest;
pub mod forest_mol;
pub mod formula_check;
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
pub mod stream;
pub mod unique_edit;
pub mod valence;

#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python_api;

#[cfg(feature = "wasm")]
mod wasm_api;

pub use atom_diff::{
    AtomDiff, added_heavy_atoms, atom_diff, atom_diff_after_cleavage, atom_diff_for_child,
    atom_diff_from_mappings, candidate_could_help, candidate_could_help_on, candidate_order_key,
    candidate_order_key_on, extend_mapping_for_added, keep_against_diff, lift_mappings,
    pair_could_help, pattern_could_help, pattern_could_help_mol, site_h_progress,
    try_atom_diff_for_child, try_lift_cleaved_child,
};
pub use atom_tracker::{AtomTracker, tags_agree_elements};
pub use candidate::{Candidate, ParentRef};
pub use canonical_plan::{
    CanonicalPlanFn, CanonicalStep, CleavageSide, Deps, Linearization, Maybe, PlanAtom, Step,
    align_deps_indices, as_deps, bind_deps, canonical_dependency_edges, identity_canonical_plan,
    identity_plan, identity_plan_with_orbit, plan_for_leaf, quinone_canonical_plan, steps_for_leaf,
    transitive_closure_masks,
};
pub use cleavage_graph::{
    CleavageArm, CleavageGraph, CleavageGraphConfig, CleavageGraphStats, CleavageLayer,
    CleavageNode, CleavageOr, CleavageSeed, CleavageSeedHop, cleavage_first_seeds,
    cleavage_first_seeds_smiles, cleavage_graph_stats, cleavage_layer, cleavage_product_graph,
    fold_cleavage_arms,
};
pub use find_path::{
    FindPath, FindPathConfig, FindPathFilters, OpenFindPath, PathCounters, PathOutcome, PathStep,
    find_path, find_path_default, find_path_diff, find_path_with, find_path_with_filters,
};
pub use forest::{Formula, Structure, formula_delta, molecule_formula};
pub use formula_check::check_effect_delta_formula;
pub use forest_mol::ForestMol;
pub use hydroxylation::{hydroxylate, hydroxylation};
pub use labels::Tag;
pub use mol::{
    ForestError, Molecule, canon_of, canon_smiles, parse_mol, ranks, stable_csmi_key,
    stable_csmi_key_of,
};
pub use orbits::{
    AtomBondGenerator, atom_bond_generators, atom_orbit, atom_orbit_with_gens, atom_pair_orbit_id,
    atom_pair_orbit_id_with_gens, atoms_orbit_with_gens, unordered_atom_pair_groups_with_gens,
    unordered_atom_pair_orbit_sizes,
};
pub use pair_edit::{PairCandidate, dehydrogenate_hydroquinone};
pub use pattern::{
    CleaveFoldKey, CleaveSideSig, Edit, Effect, Emission, PatternInfo, SiteInfo, SiteKind, When,
    bag_counts, bag_delta_formula, compose_delta_formula, leave_ch2, leave_me, leave_o,
    merge_delta_formula, named_leave_formula,
};
pub use rules::{all_rules, catalog_names, default_ruleset, leaf_rule, phase_one};
pub use ruleset::{
    BoxedFilters, FilterRules, FilterSites, RuleMember, RuleSet, accept_all_rules,
    accept_all_sites, o_dealkylation,
};
pub use smarts::smarts_matches;
pub use smirks::apply_smirks_at;
pub use stream::{Candidates, Metabolize, PairCandidates};
pub use unique_edit::{
    UniqueSite, same_site_orbit, unique_atom_sites, unique_atom_sites_with_orbits,
};
pub use valence::{accept_product, nitrogen_two_doubles};
