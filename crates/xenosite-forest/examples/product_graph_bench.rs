//! Measure general product-graph size on mid / hard PhaseOne scaffolds.
//!
//! ```text
//! cargo run -p xenosite-forest --example product_graph_bench --release
//! cargo run -p xenosite-forest --example product_graph_bench --release -- --hard
//! ```

use std::env;
use std::time::Instant;

use xenosite_forest::{phase_one, product_graph_stats};

const MID: &[(&str, &str, &str)] = &[
    (
        "eugenol→allyl-quinone",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    (
        "MeOPhOH→hydroxyquinone",
        "COc1ccc(O)cc1",
        "O=C1C=C(O)C(=O)C(O)=C1",
    ),
    ("anisole→phenol", "COc1ccccc1", "Oc1ccccc1"),
];

const HARD: &[(&str, &str, &str)] = &[
    (
        "trimethoxy-PEA→catechol",
        "COc1cc(OC)c(OC)c(CCN)c1",
        "NCCc1cc(O)c(O)c(O)c1",
    ),
    (
        "eugenol-MeO→allylQ",
        "COc1cc(CC=C)cc(OC)c1O",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "tetraMeO-biphenyl→tetraOH",
        "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC",
        "Oc1ccc(-c2ccc(O)c(O)c2)cc1O",
    ),
    (
        "veratrole-allyl→allylQ",
        "COc1ccc(CC=C)c(OC)c1OC",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
];

fn main() {
    let hard = env::args().any(|a| a == "--hard");
    let cases = if hard { HARD } else { MID };
    let set = phase_one();
    let max_nodes = 256;
    let max_depth = 6;
    println!(
        "product graph (all edits; MCS expand gate)  set={}  max_nodes={max_nodes}  max_depth={max_depth}\n",
        if hard { "HARD" } else { "MID" }
    );
    println!(
        "{:<28} {:>5} {:>5} {:>5} {:>6} {:>8}",
        "case", "nodes", "edges", "rules", "reach", "seconds"
    );
    for &(name, start, target) in cases {
        let t0 = Instant::now();
        let (stats, _) =
            product_graph_stats(start, Some(target), &set, max_nodes, max_depth).unwrap();
        let secs = t0.elapsed().as_secs_f64();
        println!(
            "{:<28} {:>5} {:>5} {:>5} {:>6} {:>8.3}",
            name,
            stats.n_nodes,
            stats.n_edges,
            stats.n_rule_patterns,
            if stats.reaches_target { "yes" } else { "no" },
            secs
        );
    }
}
