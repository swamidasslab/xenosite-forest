//! Verify site shell deltas against the actual edit when applied.
//!
//! For each live-gate survivor: materialize → tag-align reactant↔product →
//! compare (1) forecast `at_sites` vs target with (2) actual reactant→product
//! shell L1 at the same tagged atoms.
//!
//! ```text
//! cargo run -p xenosite-forest --example verify_site_shells --release
//! ```

use std::collections::BTreeMap;

use xenosite_forest::{
    AtomNeighborhood, ForestMol, aligned_shells, atom_diff, atom_neighborhood,
    candidate_could_help_on, format_shell, molecule_shells, pair_could_help, parse_mol, phase_one,
    shell_l1,
};

const CASES: &[(&str, &str, &str)] = &[
    ("CC→CCO", "CC", "CCO"),
    ("CC→CC=O", "CC", "CC=O"),
    ("anisole→phenol", "COc1ccccc1", "Oc1ccccc1"),
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    (
        "MeOPhOH→hydroxyQ",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
    (
        "eugenol→allylQ",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
    ("hydroquinone→quinone", "Oc1ccc(O)cc1", "O=C1C=CC(=O)C=C1"),
];

fn site_atoms_cand(c: &xenosite_forest::Candidate) -> Vec<usize> {
    let mut atoms: Vec<usize> = c
        .pattern
        .site_map
        .iter()
        .filter_map(|m| c.mapped.get(m).copied())
        .collect();
    if atoms.is_empty() {
        atoms.push(c.site);
    }
    for &i in &c.orbit {
        atoms.push(i);
    }
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn env_line(e: &AtomNeighborhood) -> String {
    format!(
        "n0{{{}}} n1{{{}}} n2{{{}}} Σ={}",
        format_shell(&e.n0),
        format_shell(&e.n1),
        format_shell(&e.n2),
        e.abs_delta()
    )
}

/// Actual reactant→product neighborhood delta at a reactant atom, via tags.
fn actual_delta(parent: &ForestMol, child: &ForestMol, r_idx: usize) -> Option<AtomNeighborhood> {
    let tag = parent.tag_of(r_idx)?;
    let c_idx = child.index_of(tag)?;
    let from = atom_neighborhood(parent.mol(), r_idx);
    let to = atom_neighborhood(child.mol(), c_idx);
    Some(neighborhood_delta(&to, &from))
}

fn shell_delta_map(
    to: &BTreeMap<String, i32>,
    from: &BTreeMap<String, i32>,
) -> BTreeMap<String, i32> {
    let mut out = BTreeMap::new();
    for k in to.keys().chain(from.keys()) {
        let d = to.get(k).copied().unwrap_or(0) - from.get(k).copied().unwrap_or(0);
        if d != 0 {
            out.insert(k.clone(), d);
        }
    }
    out
}

fn neighborhood_delta(to: &AtomNeighborhood, from: &AtomNeighborhood) -> AtomNeighborhood {
    AtomNeighborhood {
        aromatic: to.aromatic - from.aromatic,
        n0: shell_delta_map(&to.n0, &from.n0),
        n1: shell_delta_map(&to.n1, &from.n1),
        n2: shell_delta_map(&to.n2, &from.n2),
    }
}

fn main() {
    println!(
        "verify_site_shells — forecast (reactant→target at_sites) vs actual (reactant→product)\n"
    );

    let mut n_applied = 0usize;
    let mut n_site_atoms = 0usize;
    let mut exact_match = 0usize;
    let mut forecast_only = 0usize; // forecast nonzero, actual zero/missing
    let mut actual_only = 0usize; // actual nonzero, forecast empty/missing
    let mut both_nonzero_diff = 0usize;
    let mut cost_drop = 0usize;
    let mut cost_flat = 0usize;
    let mut cost_up = 0usize;
    let mut hop0_forecast_zero = 0usize;

    for (name, reactant, target) in CASES {
        let parent = ForestMol::parse(reactant).unwrap();
        let tgt = parse_mol(target).unwrap();
        let set = phase_one();
        let forecast_align = aligned_shells(parent.mol(), &tgt);
        let parent_cost = forecast_align.without_unchanged().cost();
        let ad = atom_diff(parent.mol(), &tgt);

        println!("=== {name}  parent_shell_cost={parent_cost} ===");

        // Single-site candidates
        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let mut checked = 0usize;
        for c in &cands {
            if !candidate_could_help_on(c, &ad, Some(parent.mol()), Some(&tgt)) {
                continue;
            }
            let Ok(pieces) = c.materialize_mols(parent.mol()) else {
                continue;
            };
            if pieces.is_empty() {
                continue;
            }
            // Prefer product that drops shell cost to target, else first piece.
            let mut best: Option<(ForestMol, usize)> = None;
            for piece in pieces {
                let child = parent.adopt_product(piece);
                let child_cost = aligned_shells(child.mol(), &tgt).without_unchanged().cost();
                best = Some(match best {
                    None => (child, child_cost),
                    Some((_b, bc)) if child_cost < bc => (child, child_cost),
                    Some(prev) => prev,
                });
            }
            let Some((child, child_cost)) = best else {
                continue;
            };
            n_applied += 1;
            checked += 1;
            match child_cost.cmp(&parent_cost) {
                std::cmp::Ordering::Less => cost_drop += 1,
                std::cmp::Ordering::Equal => cost_flat += 1,
                std::cmp::Ordering::Greater => cost_up += 1,
            }

            let atoms = site_atoms_cand(c);
            let forecast = forecast_align.at_sites(&atoms);
            let fcost = forecast.cost();
            if fcost == 0 {
                hop0_forecast_zero += 1;
            }

            println!(
                "  {} site={:?}  forecast_Σ={fcost}  shell_cost {}→{} ({})",
                c.pattern.name,
                atoms,
                parent_cost,
                child_cost,
                if child_cost < parent_cost {
                    "↓"
                } else if child_cost > parent_cost {
                    "↑"
                } else {
                    "="
                }
            );

            for &r in &atoms {
                n_site_atoms += 1;
                let f_env = forecast.atoms.get(&r);
                let a_env = actual_delta(&parent, &child, r);
                match (f_env, a_env) {
                    (None, None) => {
                        exact_match += 1; // both absent / unchanged
                    }
                    (None, Some(a)) if a.is_unchanged() => exact_match += 1,
                    (None, Some(a)) => {
                        actual_only += 1;
                        println!("    r{r} ACTUAL_ONLY  actual={}", env_line(&a));
                    }
                    (Some(f), None) => {
                        // Atom left (cleavage) — forecast had a delta, product lost the atom.
                        if f.is_unchanged() {
                            exact_match += 1;
                        } else {
                            forecast_only += 1;
                            println!(
                                "    r{r} FORECAST_ONLY (atom gone on product)  forecast={}",
                                env_line(f)
                            );
                        }
                    }
                    (Some(f), Some(a)) => {
                        if f == &a {
                            exact_match += 1;
                        } else if f.is_unchanged() && a.is_unchanged() {
                            exact_match += 1;
                        } else if f.is_unchanged() {
                            actual_only += 1;
                            println!("    r{r} ACTUAL_ONLY  actual={}", env_line(&a));
                        } else if a.is_unchanged() {
                            forecast_only += 1;
                            println!("    r{r} FORECAST_ONLY  forecast={}", env_line(f));
                        } else {
                            both_nonzero_diff += 1;
                            // How much of the forecast appears in the actual edit?
                            let overlap = shell_l1(&f.n1, &a.n1)
                                + shell_l1(&f.n2, &a.n2)
                                + shell_l1(&f.n0, &a.n0);
                            println!(
                                "    r{r} DIFF  forecast={}  actual={}  |f−a|={overlap}",
                                env_line(f),
                                env_line(&a)
                            );
                        }
                    }
                }
            }
            if checked >= 8 {
                break; // cap per case
            }
        }

        // Pairs (DH / QF) — sample a few gate survivors
        let pairs = set
            .pair_candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let mut pair_checked = 0usize;
        for p in &pairs {
            if !pair_could_help(p, &ad, parent.mol(), &tgt) {
                continue;
            }
            let Ok(pieces) = p.materialize_mols(parent.mol()) else {
                continue;
            };
            if pieces.is_empty() {
                continue;
            }
            let child = parent.adopt_product(pieces[0].clone());
            let child_cost = aligned_shells(child.mol(), &tgt).without_unchanged().cost();
            n_applied += 1;
            pair_checked += 1;
            match child_cost.cmp(&parent_cost) {
                std::cmp::Ordering::Less => cost_drop += 1,
                std::cmp::Ordering::Equal => cost_flat += 1,
                std::cmp::Ordering::Greater => cost_up += 1,
            }
            let atoms = p.plan_site_atoms();
            let forecast = forecast_align.at_sites(&atoms);
            let fcost = forecast.cost();
            if fcost == 0 {
                hop0_forecast_zero += 1;
            }
            println!(
                "  PAIR {} atoms={:?} forecast_Σ={fcost}  shell_cost {}→{} ({})",
                p.pattern_name,
                atoms,
                parent_cost,
                child_cost,
                if child_cost < parent_cost {
                    "↓"
                } else if child_cost > parent_cost {
                    "↑"
                } else {
                    "="
                }
            );
            for &r in &atoms {
                n_site_atoms += 1;
                let f_env = forecast.atoms.get(&r);
                let a_env = actual_delta(&parent, &child, r);
                match (f_env, a_env.as_ref()) {
                    (Some(f), Some(a)) if f == a => exact_match += 1,
                    (Some(f), Some(a)) if !f.is_unchanged() && !a.is_unchanged() => {
                        both_nonzero_diff += 1;
                        println!(
                            "    r{r} DIFF  forecast={}  actual={}",
                            env_line(f),
                            env_line(a)
                        );
                    }
                    (Some(f), _) if !f.is_unchanged() => {
                        forecast_only += 1;
                        println!("    r{r} FORECAST_ONLY  {}", env_line(f));
                    }
                    (_, Some(a)) if !a.is_unchanged() => {
                        actual_only += 1;
                        println!("    r{r} ACTUAL_ONLY  {}", env_line(a));
                    }
                    _ => exact_match += 1,
                }
            }
            if pair_checked >= 4 {
                break;
            }
        }
        println!();
    }

    println!("## Summary");
    println!("  applied edits: {n_applied}");
    println!("  site-atom checks: {n_site_atoms}");
    println!("  exact match (forecast==actual or both empty): {exact_match}");
    println!("  both nonzero but differ: {both_nonzero_diff}");
    println!("  forecast only (target residual ≠ applied edit / atom left): {forecast_only}");
    println!("  actual only (edit changed site; forecast silent): {actual_only}");
    println!("  gate survivors with forecast_Σ=0: {hop0_forecast_zero}/{n_applied}");
    println!("  full shell_cost to target after apply: ↓{cost_drop} ={cost_flat} ↑{cost_up}");

    let _ = (molecule_shells, neighborhood_delta);
}
