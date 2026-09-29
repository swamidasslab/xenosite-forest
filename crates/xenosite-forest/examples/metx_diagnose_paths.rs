use std::time::{Duration, Instant};
use xenosite_forest::{FindPathConfig, PathCounters, find_path_partial, phase_one};

fn run(label: &str, r: &str, p: &str) {
    let set = phase_one();
    let mut c = PathCounters::default();
    let cfg = FindPathConfig {
        max_paths: 1,
        max_nodes: 200,
        use_atom_diff: true,
        lazy_closer: true,
        timeout: Some(Duration::from_secs_f64(2.0)),
        ..FindPathConfig::default()
    };
    let t0 = Instant::now();
    let out = find_path_partial(r, p, &set, &mut c, cfg, None, |_| true).unwrap();
    println!(
        "\n{label}\n  hit={} bill={} nodes={} edits={} secs={:.2} timed_out={} exact={} partials={}",
        !out.exact.is_empty(),
        c.billed(),
        c.nodes,
        c.mol_edits,
        t0.elapsed().as_secs_f64(),
        c.timed_out,
        out.exact.len(),
        out.partials.len()
    );
    if let Some(h) = out.exact.first() {
        println!("  path smiles={} steps={}", h.smiles, h.steps.len());
        for s in &h.steps {
            println!("    {} @{} → {}", s.pattern_name, s.site, s.product);
        }
    }
    if let Some(p) = out.partials.first() {
        println!(
            "  closest cost={} cats={:?} smiles={}",
            p.residual.cost, p.residual.categories, p.smiles
        );
    }
}

fn main() {
    run("N-OH-IQ", "Cn1c(=N)[nH]c2c3cccnc3ccc21", "Cn1c(NO)nc2c3cccnc3ccc21");
    run("tacrine→7-OH", "N=c1c2c([nH]c3ccccc13)CCCC2", "Nc1c2c(nc3c(O)cccc13)CCCC2");
    run("CTRL caffeine→theobromine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Cn1cnc2c(=O)[nH]c(=O)n(C)c12");
    run("thalidomide arene oxide", "O=C1N=C(O)CC[C@@H]1N1C(=O)c2ccccc2C1=O", "O=C1N=C(O)CC[C@@H]1N1C(=O)C2=CC3OC3C=C2C1=O");
    run("leflunomide N-O open", "Cc1oncc1C(=O)Nc1ccc(C(F)(F)F)cc1", r"C/C(O)=C(/C=N)C(=O)Nc1ccc(C(F)(F)F)cc1");
    // What if tacrine product kept imine tautomer with OH?
    run("tacrine→7-OH imine-tautomer guess", "N=c1c2c([nH]c3ccccc13)CCCC2", "N=c1c2c([nH]c3c(O)cccc13)CCCC2");
}
