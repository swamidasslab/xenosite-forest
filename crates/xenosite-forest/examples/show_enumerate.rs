//! Demo: stream BFS/DFS metabolites with PathInfo.
//!
//! cargo run -p xenosite-forest --example show_enumerate --release

use xenosite_forest::{EnumConfig, PathInfo, bfs, dfs, enumerate_metabolites, phase_one};

fn show_path(path: &PathInfo) -> String {
    path.hops
        .iter()
        .map(|h| format!("{}:{}", h.rule, h.pattern_name))
        .collect::<Vec<_>>()
        .join(" → ")
}

fn main() {
    let set = phase_one();

    println!("=== BFS ethane depth=2 (hydroxylation via PhaseOne) ===");
    let mut n = 0usize;
    for hit in bfs("CC", &set, 2).unwrap() {
        let hit = hit.unwrap();
        n += 1;
        if n <= 12 {
            println!(
                "  d{} {}  path=[{}]",
                hit.depth(),
                hit.smiles(),
                show_path(&hit.path)
            );
        }
    }
    println!("  total={n}");

    println!("\n=== DFS vs BFS order on propane depth=2 (first 8) ===");
    print!("BFS depths: ");
    for (i, hit) in bfs("CCC", &set, 2).unwrap().enumerate() {
        if i >= 8 {
            break;
        }
        print!("{} ", hit.unwrap().depth());
    }
    println!();
    print!("DFS depths: ");
    for (i, hit) in dfs("CCC", &set, 2).unwrap().enumerate() {
        if i >= 8 {
            break;
        }
        print!("{} ", hit.unwrap().depth());
    }
    println!();

    println!("\n=== Anisole PhaseOne depth=2 (capped) ===");
    let hits: Vec<_> = enumerate_metabolites(
        "COc1ccccc1",
        &set,
        EnumConfig::bfs(2).with_max_nodes(40),
    )
    .unwrap()
    .map(|h| h.unwrap())
    .collect();
    println!("  yielded={} (max_nodes=40 incl. root)", hits.len());
    for hit in hits.iter().take(10) {
        println!(
            "  d{} {}  [{}]",
            hit.depth(),
            hit.smiles(),
            show_path(&hit.path)
        );
    }
}
