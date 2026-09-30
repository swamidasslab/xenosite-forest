//! Unreachable conjugate on a *combinatorially improving* PhaseOne path.
//!
//! Not OH-specific: any target that admits many cost-dropping PhaseOne edits
//! (dealk, OH, DH, …) but ends in GSH (no PhaseOne rule) burns max_nodes.
//! Also runs [`find_path_partial`] to measure closest residual flush.
//!
//! ```text
//! cargo run -p xenosite-forest --example unreachable_conjugate_probe --release -- 800
//! ```

use std::time::Instant;

use xenosite_forest::{
    FindPathConfig, MetabolicNetwork, PathCounters, canon_of, find_path_partial, find_path_with,
    phase_one,
};

fn run_exact(name: &str, reactant: &str, target: &str, config: FindPathConfig) {
    if let Err(e) = canon_of(target) {
        eprintln!("{name}: target canon failed: {e}");
        return;
    }
    let set = phase_one();
    let mut counters = PathCounters::default();
    let t0 = Instant::now();
    let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
        .unwrap()
        .collect_all()
        .unwrap();
    let secs = t0.elapsed().as_secs_f64();
    println!(
        "{name:<48} hit={:<5} secs={secs:7.3} nodes={:<5} edits={:<7} bill={:<7} timed_out={}",
        !hits.is_empty(),
        counters.nodes,
        counters.mol_edits,
        counters.billed(),
        counters.timed_out
    );
}

fn run_partial(name: &str, reactant: &str, target: &str, config: FindPathConfig) {
    if let Err(e) = canon_of(target) {
        eprintln!("{name}: target canon failed: {e}");
        return;
    }
    let set = phase_one();
    let mut counters = PathCounters::default();
    let mut net = MetabolicNetwork::new();
    let t0 = Instant::now();
    let out = find_path_partial(
        reactant,
        target,
        &set,
        &mut counters,
        config,
        Some(&mut net),
        |_| true,
    )
    .unwrap();
    let secs = t0.elapsed().as_secs_f64();
    let best = out
        .partials
        .first()
        .map(|p| {
            format!(
                "cost={} cats={:?} smiles={}",
                p.residual.cost, p.residual.categories, p.smiles
            )
        })
        .unwrap_or_else(|| "none".into());
    println!(
        "{name:<48} exact={} partials={} secs={secs:7.3} nodes={:<5} edits={:<7} bill={:<7} net_nodes={} best={best}",
        out.exact.len(),
        out.partials.len(),
        counters.nodes,
        counters.mol_edits,
        counters.billed(),
        net.n_nodes(),
    );
}

fn main() {
    let max_nodes = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(800usize);
    let config = FindPathConfig {
        max_paths: 2,
        max_nodes,
        use_atom_diff: true,
        lazy_closer: true,
        ..FindPathConfig::default()
    };
    println!("many improving PhaseOne edits + unreachable GSH  max_nodes={max_nodes}\n");

    let bp_r = "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC";
    let bp_oh = "Oc1ccc(-c2ccc(O)c(O)c2)cc1O";
    let bp_gsh =
        "N[C@@H](CCC(=O)N[C@@H](CSc1c(O)ccc(-c2ccc(O)c(O)c2)c1O)C(=O)NCC(=O)O)C(=O)O";

    println!("=== find_path (exact-only) ===");
    run_exact("CTRL tetraMeO-BP→tetraOH", bp_r, bp_oh, config);
    run_exact("tetraMeO-BP→tetraOH+GSH", bp_r, bp_gsh, config);
    run_exact(
        "CTRL eugenol→allylQ",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
        config,
    );
    run_exact(
        "eugenol→allylQ+GSH",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O",
        config,
    );
    run_exact(
        "CTRL MeOPhOH→hydroxyQ",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
        config,
    );
    run_exact(
        "MeOPhOH→hydroxyQ+GSH",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O",
        config,
    );

    println!("\n=== find_path_partial (exact + closest flush) ===");
    run_partial("CTRL tetraMeO-BP→tetraOH", bp_r, bp_oh, config);
    run_partial("tetraMeO-BP→tetraOH+GSH", bp_r, bp_gsh, config);
    run_partial(
        "eugenol→allylQ+GSH",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O",
        config,
    );
    run_partial(
        "MeOPhOH→hydroxyQ+GSH",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O",
        config,
    );
}
