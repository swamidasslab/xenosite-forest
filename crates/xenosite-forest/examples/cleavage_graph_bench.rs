//! Measure cleavage product graph size on LARGER / HARD scaffolds.
//!
//! Both sides of each cleavage are nodes; arms fold into Or by fragment multiset.
//!
//! ```text
//! cargo run -p xenosite-forest --example cleavage_graph_bench --release
//! cargo run -p xenosite-forest --example cleavage_graph_bench --release -- --hard
//! ```

use std::env;
use std::time::Instant;

use xenosite_forest::{cleavage_graph_stats, phase_one};

const LARGER: &[(&str, &str, &str)] = &[
    (
        "tBu-bis-ND→dialdehyde",
        "CN(C)Cc1ccc(CN(C)Cc2ccc(C(C)(C)C)cc2)cc1",
        "O=Cc1ccc(C=O)cc1",
    ),
    (
        "tribenzyl→PhCHO",
        "N(Cc1ccccc1)(Cc1ccccc1)Cc1ccccc1",
        "O=Cc1ccccc1",
    ),
    (
        "MeO-diphenyl→catechol",
        "COc1ccc(Cc2ccc(OC)cc2)cc1",
        "Oc1ccc(Cc2ccc(O)cc2)cc1",
    ),
    (
        "dimethoxy-PEA→catechol",
        "COc1ccc(CCN)cc1OC",
        "NCCc1ccc(O)c(O)c1",
    ),
    (
        "eugenol→allyl-quinone",
        "COc1ccc(CC=C)cc1O",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
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
    (
        "tetraMeO-naph→polyOH-NQ",
        "COc1cc(OC)c2c(OC)cc(OC)cc2c1",
        "O=C1C=C(O)C(=O)c2c(O)cc(O)cc12",
    ),
];

fn main() {
    let hard = env::args().any(|a| a == "--hard");
    let cases = if hard { HARD } else { LARGER };
    let set = phase_one();
    println!(
        "cleavage product graph (both sides kept, Or by fragment multiset)  set={}\n",
        if hard { "HARD" } else { "LARGER" }
    );
    println!(
        "{:<28} {:>6} {:>5} {:>5} {:>6} {:>8} {:>8} {:>8} {:>7}",
        "case", "nodes", "ors", "arms", "fanin", "1h_prod", "1h_arms", "1h_frag", "ms"
    );
    for (name, start, target) in cases {
        let t0 = Instant::now();
        let (stats, _) = cleavage_graph_stats(start, Some(target), &set).unwrap();
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        println!(
            "{:<28} {:>6} {:>5} {:>5} {:>6} {:>8} {:>8} {:>8} {:>6.0}",
            name,
            stats.nodes,
            stats.ors,
            stats.arms,
            stats.max_or_fanin,
            stats.one_hop_products,
            stats.one_hop_arms,
            stats.one_hop_fragments,
            ms
        );
    }
}
