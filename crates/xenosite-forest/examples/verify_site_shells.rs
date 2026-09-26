//! Verify site shell bags: forecast (reactant→target at site) vs applied edit.
//!
//! Cost at a site is Σ |current + δ − target| (`site_shell_cost`). Bags must
//! match exactly; mismatches are reported (Error mode aborts the example).
//! Close pairs use the joint site atom list.
//!
//! ```text
//! cargo run -p xenosite-forest --example verify_site_shells --release
//! ```

use xenosite_forest::{
    ForestMol, SiteShellCheck, aligned_shells, atom_diff, candidate_could_help_on,
    check_site_shell_bags, edit_shells, edit_site_bag, forecast_site_bag, molecule_shells,
    pair_could_help, parse_mol, phase_one, site_shell_cost,
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

fn main() {
    println!("verify_site_shells — bag match + site_shell_cost = Σ|current+δ−target|\n");

    let mut n_applied = 0usize;
    let mut n_match = 0usize;
    let mut n_mismatch = 0usize;
    let mut n_inconsistent = 0usize;
    let mut cost_now_pos = 0usize;
    let mut cost_edit_zero = 0usize;
    let mut cost_edit_pos = 0usize;
    let mut residual_drop = 0usize;
    let mut residual_flat = 0usize;
    let mut residual_up = 0usize;

    for (name, reactant, target) in CASES {
        let parent = ForestMol::parse(reactant).unwrap();
        let tgt = parse_mol(target).unwrap();
        let set = phase_one();
        let forecast_align = aligned_shells(parent.mol(), &tgt);
        let map = atom_diff(parent.mol(), &tgt).mapping;
        let cur = molecule_shells(parent.mol());
        let tgt_shells = molecule_shells(&tgt);
        let parent_residual = site_shell_cost(
            &cur,
            None,
            &tgt_shells,
            &map,
            &cur.atoms.keys().copied().collect::<Vec<_>>(),
        );

        println!("=== {name}  parent_residual={parent_residual:.4} ===");

        let cands = set
            .candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let mut checked = 0usize;
        for c in &cands {
            if !candidate_could_help_on(
                c,
                &atom_diff(parent.mol(), &tgt),
                Some(parent.mol()),
                Some(&tgt),
            ) {
                continue;
            }
            let Ok(pieces) = c.materialize_mols(parent.mol()) else {
                continue;
            };
            if pieces.is_empty() {
                continue;
            }
            let mut best: Option<(xenosite_forest::ForestMol, f64)> = None;
            for piece in pieces {
                let child = parent.adopt_product(piece);
                let child_shells = molecule_shells(child.mol());
                let child_map = atom_diff(child.mol(), &tgt).mapping;
                let child_residual = site_shell_cost(
                    &child_shells,
                    None,
                    &tgt_shells,
                    &child_map,
                    &child_shells.atoms.keys().copied().collect::<Vec<_>>(),
                );
                best = Some(match best {
                    None => (child, child_residual),
                    Some((_b, bc)) if child_residual < bc => (child, child_residual),
                    Some(prev) => prev,
                });
            }
            let Some((child, child_residual)) = best else {
                continue;
            };
            n_applied += 1;
            checked += 1;
            match child_residual
                .partial_cmp(&parent_residual)
                .unwrap_or(std::cmp::Ordering::Equal)
            {
                std::cmp::Ordering::Less => residual_drop += 1,
                std::cmp::Ordering::Equal => residual_flat += 1,
                std::cmp::Ordering::Greater => residual_up += 1,
            }

            let atoms = site_atoms_cand(c);
            let forecast = forecast_site_bag(&forecast_align, &atoms);
            let actual = edit_site_bag(&parent, &child, &atoms);
            let now = site_shell_cost(&cur, None, &tgt_shells, &map, &atoms);
            let edit = edit_shells(&parent, &child);
            let after = site_shell_cost(&cur, Some(&edit), &tgt_shells, &map, &atoms);
            if now > 0.0 {
                cost_now_pos += 1;
            }
            if after < 1e-9 {
                cost_edit_zero += 1;
            } else {
                cost_edit_pos += 1;
            }

            let label = format!("{} / {}", name, c.pattern.name);
            let bags_eq = forecast.matches(&actual);
            let done = after < 1e-9;
            if bags_eq {
                n_match += 1;
                println!(
                    "  OK  {} site={:?}  now={now:.4} after_δ={after:.4}  residual {parent_residual:.4}→{child_residual:.4}",
                    c.pattern.name, atoms
                );
            } else {
                n_mismatch += 1;
                let _ = check_site_shell_bags(&label, &forecast, &actual, SiteShellCheck::Warn);
                println!(
                    "  DIFF  {} site={:?}  bag_l1={}  now={now:.4} after_δ={after:.4}  cleaved {}→{} kept {}→{}",
                    c.pattern.name,
                    atoms,
                    forecast.l1(&actual),
                    forecast.cleaved,
                    actual.cleaved,
                    forecast.kept.len(),
                    actual.kept.len()
                );
            }
            // Residual bag ≡ edit bag iff the edit completes the site.
            if bags_eq != done {
                n_inconsistent += 1;
                println!("    INCONSISTENT bags_eq={bags_eq} after_δ≈0={done}");
            }
            if checked >= 8 {
                break;
            }
        }

        let pairs = set
            .pair_candidates(parent.mol())
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let ad = atom_diff(parent.mol(), &tgt);
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
            let child_shells = molecule_shells(child.mol());
            let child_map = atom_diff(child.mol(), &tgt).mapping;
            let child_residual = site_shell_cost(
                &child_shells,
                None,
                &tgt_shells,
                &child_map,
                &child_shells.atoms.keys().copied().collect::<Vec<_>>(),
            );
            n_applied += 1;
            pair_checked += 1;
            match child_residual
                .partial_cmp(&parent_residual)
                .unwrap_or(std::cmp::Ordering::Equal)
            {
                std::cmp::Ordering::Less => residual_drop += 1,
                std::cmp::Ordering::Equal => residual_flat += 1,
                std::cmp::Ordering::Greater => residual_up += 1,
            }
            // Joint ends — close pairs share shells.
            let atoms = p.plan_site_atoms();
            let forecast = forecast_site_bag(&forecast_align, &atoms);
            let actual = edit_site_bag(&parent, &child, &atoms);
            let now = site_shell_cost(&cur, None, &tgt_shells, &map, &atoms);
            let edit = edit_shells(&parent, &child);
            let after = site_shell_cost(&cur, Some(&edit), &tgt_shells, &map, &atoms);
            if now > 0.0 {
                cost_now_pos += 1;
            }
            if after < 1e-9 {
                cost_edit_zero += 1;
            } else {
                cost_edit_pos += 1;
            }
            let label = format!("{} / PAIR {}", name, p.pattern_name);
            let bags_eq = forecast.matches(&actual);
            let done = after < 1e-9;
            if bags_eq {
                n_match += 1;
                println!(
                    "  OK  PAIR {} atoms={:?}  now={now:.4} after_δ={after:.4}  residual {parent_residual:.4}→{child_residual:.4}",
                    p.pattern_name, atoms
                );
            } else {
                n_mismatch += 1;
                let _ = check_site_shell_bags(&label, &forecast, &actual, SiteShellCheck::Warn);
                println!(
                    "  DIFF  PAIR {} atoms={:?}  bag_l1={}  now={now:.4} after_δ={after:.4}  cleaved {}→{} kept {}→{}",
                    p.pattern_name,
                    atoms,
                    forecast.l1(&actual),
                    forecast.cleaved,
                    actual.cleaved,
                    forecast.kept.len(),
                    actual.kept.len()
                );
            }
            if bags_eq != done {
                n_inconsistent += 1;
                println!("    INCONSISTENT bags_eq={bags_eq} after_δ≈0={done}");
            }
            if pair_checked >= 4 {
                break;
            }
        }
        println!();
    }

    println!("## Summary");
    println!("  applied edits: {n_applied}");
    println!("  bag match (residual≡edit at site): {n_match}");
    println!("  bag diff (partial / off-target edit): {n_mismatch}");
    println!("  inconsistent (bags_eq XOR after_δ≈0): {n_inconsistent}");
    println!("  site now>0: {cost_now_pos}/{n_applied}");
    println!("  site after_δ==0 (edit completes site): {cost_edit_zero}");
    println!("  site after_δ>0 (partial / off-target): {cost_edit_pos}");
    println!("  full residual after apply: ↓{residual_drop} ={residual_flat} ↑{residual_up}");
    if n_inconsistent > 0 {
        eprintln!("verify_site_shells: {n_inconsistent} bag/cost inconsistenc(ies)");
        std::process::exit(1);
    }
}
