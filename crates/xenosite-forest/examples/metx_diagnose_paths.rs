//! Path search on MetX misses with call-site stereo strip (no lib change).
use chematic::chem::remove_stereo;
use std::time::{Duration, Instant};
use xenosite_forest::{FindPathConfig, PathCounters, as_forest_mol, find_path_partial, phase_one};

fn nostereo(s: &str) -> String {
    let fm = as_forest_mol(s).unwrap();
    let stripped = remove_stereo(fm.mol());
    as_forest_mol(stripped).unwrap().csmi().as_ref().to_string()
}

fn run(label: &str, r: &str, p: &str) {
    let r = nostereo(r);
    let p = nostereo(p);
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
    let out = find_path_partial(&r, &p, &set, &mut c, cfg, None, |_| true).unwrap();
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
        println!("  path steps={}", h.steps.len());
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
    run(
        "BIOTID00138 S-ox thrash",
        "CC(=O)c1cc2ccccc2s1",
        "CC(O)C1CC2CCCCC2S1=O",
    );
    run(
        "BIOTID00198 chromenone sat",
        "CCCCOc1ccc2ccc(=O)oc2c1",
        "CC(O)CCOC1CCC2CCC(O)OC2C1",
    );
    run(
        "BIOTID00055 amphetamine 'ArOH'",
        "COc1ccccc1C[C@@H](C)N",
        "COC1CCC(O)CC1CC(C)N",
    );
    run(
        "BIOTID00966 tacrine→7-OH",
        "N=c1c2c([nH]c3ccccc13)CCCC2",
        "Nc1c2c(nc3c(O)cccc13)CCCC2",
    );
    run(
        "BIOTID00966 imine-tautomer guess",
        "N=c1c2c([nH]c3ccccc13)CCCC2",
        "N=c1c2c([nH]c3c(O)cccc13)CCCC2",
    );
    run(
        "BIOTID01270 leflunomide N-O",
        "Cc1oncc1C(=O)Nc1ccc(C(F)(F)F)cc1",
        r"C/C(O)=C(/C=N)C(=O)Nc1ccc(C(F)(F)F)cc1",
    );
    run(
        "BIOTID00154 dichlorobenzonitrile epox",
        "N#Cc1c(Cl)cccc1Cl",
        "N#CC1C(Cl)=CC=C2OC21Cl",
    );
    run(
        "BIOTID01090 vicriviroc N-dealk",
        "COC[C@@H](c1ccc(C(F)(F)F)cc1)N1CCN(C2(C)CCN(C(=O)c3c(C)ncnc3C)CC2)C[C@@H]1C",
        "Cc1ncnc(C)c1C(=O)N1CCC(C)(N2CCN[C@@H](C)C2)CC1",
    );
    run(
        "BIOTID00184 nitropyrene diol",
        "O=[N+]([O-])c1cc2cccc3ccc4cccc1c4c32",
        "O=[N+]([O-])c1cc2cccc3c2c2c(cccc12)[C@@H](O)[C@@H]3O",
    );
    run(
        "BIOTID01318 phenol→epoxide",
        "CCc1nn(CCCN2CCN(c3ccc(O)c(Cl)c3)CC2)c(=O)n1CCOc1ccccc1",
        "CCc1nn(CCCN2CCN([C@]34C=C(Cl)C(O)=C[C@H]3O4)CC2)c(=O)n1CCOc1ccccc1",
    );
    run(
        "BIOTID00059 TZD open",
        "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(COc1ccc(C[C@@H]3SC(=O)N=C3O)cc1)O2",
        "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(COc1ccc(C[C@@H](C(=O)N=CO)S(=O)O)cc1)O2",
    );
}
