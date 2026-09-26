//! Shell-cost bench on real mid/hard find_path cases.
//!
//! Site cost = Σ normalized |current − target| at site atoms (`site_shell_cost`,
//! δ = 0). Keep iff cost > 0. Path column is the same residual over all heavies.
//!
//! ```text
//! cargo run -p xenosite-forest --example shell_cost_bench --release
//! cargo run -p xenosite-forest --example shell_cost_bench --release -- --hard
//! ```

use std::time::Instant;

use xenosite_forest::{
    FindPathConfig, PathCounters, atom_diff, candidate_could_help_on, find_path_with,
    molecule_shells, pair_could_help, parse_mol, phase_one, site_atoms_with_leave, site_shell_cost,
};

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
    (
        "TBA→enyne aldehyde",
        "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",
        r"C(#C/C=C/C=O)C(C)(C)C",
    ),
    (
        "2-MeO-naph→1,2-NQ",
        "COc1ccc2ccccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
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
        "bisMeO-naph→1,2NQ",
        "COc1ccc2c(OC)cccc2c1",
        "O=C1C(=O)c2ccccc2C=C1",
    ),
    (
        "tetraMeO-biphenyl→tetraOH",
        "COc1ccc(-c2ccc(OC)c(OC)c2)cc1OC",
        "Oc1ccc(-c2ccc(O)c(O)c2)cc1O",
    ),
    (
        "veratrole-allyl→allylQ",
        "COc1ccc(CC=C)cc1OC",
        "O=C1C=CC(=O)C(CC=C)=C1",
    ),
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

fn all_heavy(shells: &xenosite_forest::MoleculeShells) -> Vec<usize> {
    shells.atoms.keys().copied().collect()
}

struct Row {
    name: String,
    hit: bool,
    secs: f64,
    steps: usize,
    edits: usize,
    n_cand: usize,
    gate_keep: usize,
    shell_keep: usize,
    tp: usize,
    fp: usize,
    tn: usize,
    fn_: usize,
    hop0_site_cost: Option<f64>,
    hop0_rank: Option<(usize, usize)>,
    shell_costs: Vec<f64>,
    atom_costs: Vec<usize>,
    shell_mono: bool,
    atom_mono: bool,
}

fn eval_case(name: &str, reactant: &str, target: &str) -> Row {
    let ra = parse_mol(reactant).unwrap();
    let rb = parse_mol(target).unwrap();
    let set = phase_one();
    let cur0 = molecule_shells(&ra);
    let tgt = molecule_shells(&rb);
    let map0 = atom_diff(&ra, &rb).mapping;
    let ad = atom_diff(&ra, &rb);

    let mut scored: Vec<(f64, bool, Vec<usize>)> = Vec::new();
    for c in set.candidates(&ra).collect::<Result<Vec<_>, _>>().unwrap() {
        let mut atoms = site_atoms_cand(&c);
        if c.pattern.effect.cleaves {
            let leave_n = c.pattern.effect.leave_count.map(|n| n as usize);
            atoms = site_atoms_with_leave(&ra, &atoms, leave_n, &ad.cleavage_bonds);
        }
        let gate = candidate_could_help_on(&c, &ad, Some(&ra), Some(&rb));
        let cost = site_shell_cost(&cur0, None, &tgt, &map0, &atoms);
        scored.push((cost, gate, atoms));
    }
    for p in set
        .pair_candidates(&ra)
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
    {
        let gate = pair_could_help(&p, &ad, &ra, &rb);
        let mut atoms = p.plan_site_atoms();
        atoms.sort_unstable();
        atoms.dedup();
        // Close pairs: joint ends already; no leave fragment.
        let cost = site_shell_cost(&cur0, None, &tgt, &map0, &atoms);
        scored.push((cost, gate, atoms));
    }

    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut tn = 0usize;
    let mut fn_ = 0usize;
    let mut gate_keep = 0usize;
    let mut shell_keep = 0usize;
    for &(cost, gate, _) in &scored {
        let shell = cost > 1e-12;
        if gate {
            gate_keep += 1;
        }
        if shell {
            shell_keep += 1;
        }
        match (gate, shell) {
            (true, true) => tp += 1,
            (false, false) => tn += 1,
            (false, true) => fp += 1,
            (true, false) => fn_ += 1,
        }
    }

    let mut counters = PathCounters::default();
    let config = FindPathConfig {
        max_nodes: 800,
        max_paths: 1,
        use_atom_diff: true,
        ..FindPathConfig::default()
    };
    let t0 = Instant::now();
    let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .collect_all()
        .unwrap();
    let secs = t0.elapsed().as_secs_f64();

    let mut shell_costs = Vec::new();
    let mut atom_costs = Vec::new();
    let mut shell_mono = true;
    let mut atom_mono = true;
    let mut hop0_site_cost = None;
    let mut hop0_rank = None;
    let mut steps = 0usize;

    if let Some(hit) = hits.first() {
        steps = hit.steps.len();
        let mut cur = reactant.to_string();
        let mut prev_s: Option<f64> = None;
        let mut prev_a: Option<usize> = None;
        for (i, step) in hit.steps.iter().enumerate() {
            let mol = parse_mol(&cur).unwrap();
            let cur_s = molecule_shells(&mol);
            let map = atom_diff(&mol, &rb).mapping;
            let sh = site_shell_cost(&cur_s, None, &tgt, &map, &all_heavy(&cur_s));
            let at = atom_diff(&mol, &rb).cost();
            shell_costs.push(sh);
            atom_costs.push(at);
            if let Some(ps) = prev_s {
                if sh > ps + 1e-9 {
                    shell_mono = false;
                }
            }
            if let Some(pa) = prev_a {
                if at > pa {
                    atom_mono = false;
                }
            }
            prev_s = Some(sh);
            prev_a = Some(at);

            if i == 0 {
                let mut atoms = if step.site_orbit.is_empty() {
                    vec![step.site]
                } else {
                    step.site_orbit.clone()
                };
                atoms.sort_unstable();
                atoms.dedup();
                // Cleavage hops: fold in leaving-fragment heavies.
                if !step.sides.is_empty() || ad.has_cleavage() {
                    atoms = site_atoms_with_leave(&ra, &atoms, None, &ad.cleavage_bonds);
                }
                let cost = site_shell_cost(&cur0, None, &tgt, &map0, &atoms);
                hop0_site_cost = Some(cost);
                let mut gate_costs: Vec<f64> = scored
                    .iter()
                    .filter(|(_, g, _)| *g)
                    .map(|(c, _, _)| *c)
                    .collect();
                gate_costs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                let n = gate_costs.len();
                // Rank by cost value among gate-kept (1 = highest residual).
                let rank = gate_costs.iter().filter(|&&c| c > cost + 1e-9).count() + 1;
                hop0_rank = Some(if n == 0 { (0, 0) } else { (rank.min(n), n) });
            }
            cur = step.product.clone();
        }
        shell_costs.push(0.0);
        atom_costs.push(0);
    }

    Row {
        name: name.into(),
        hit: !hits.is_empty(),
        secs,
        steps,
        edits: counters.mol_edits,
        n_cand: scored.len(),
        gate_keep,
        shell_keep,
        tp,
        fp,
        tn,
        fn_,
        hop0_site_cost,
        hop0_rank,
        shell_costs,
        atom_costs,
        shell_mono,
        atom_mono,
    }
}

fn print_table(title: &str, rows: &[Row]) {
    println!("\n=== {title} ===\n");
    println!(
        "{:<28} {:>3} {:>5} {:>4} {:>4}/{:<4} {:>3} {:>3} {:>3} {:>3} {:>5} {:>8}  shell_path            atom_path   mono",
        "case", "hit", "s", "ed", "gK", "sK", "TP", "FP", "TN", "FN", "hop0", "rank"
    );
    let mut sum_tp = 0;
    let mut sum_fp = 0;
    let mut sum_tn = 0;
    let mut sum_fn = 0;
    let mut top1 = 0usize;
    let mut known = 0usize;
    for r in rows {
        sum_tp += r.tp;
        sum_fp += r.fp;
        sum_tn += r.tn;
        sum_fn += r.fn_;
        if let Some((rank, _)) = r.hop0_rank {
            if rank > 0 {
                known += 1;
                if rank == 1 {
                    top1 += 1;
                }
            }
        }
        let hit = if r.hit { "ok" } else { "MISS" };
        let hop = r
            .hop0_site_cost
            .map(|c| format!("{c:.2}"))
            .unwrap_or_else(|| "-".into());
        let rank = match r.hop0_rank {
            Some((0, n)) => format!("?/{n}"),
            Some((r0, n)) => format!("{r0}/{n}"),
            None => "-".into(),
        };
        let shell_path = r
            .shell_costs
            .iter()
            .map(|c| format!("{c:.1}"))
            .collect::<Vec<_>>()
            .join("→");
        let atom_path = r
            .atom_costs
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join("→");
        let mono = format!(
            "s{} a{}",
            if r.shell_mono { "↓" } else { "↑" },
            if r.atom_mono { "↓" } else { "↑" }
        );
        println!(
            "{:<28} {:>3} {:>5.2} {:>4} {:>4}/{:<4} {:>3} {:>3} {:>3} {:>3} {:>5} {:>8}  {:>20}  {:>10}  {mono}",
            r.name,
            hit,
            r.secs,
            r.edits,
            r.gate_keep,
            r.shell_keep,
            r.tp,
            r.fp,
            r.tn,
            r.fn_,
            hop,
            rank,
            shell_path,
            atom_path
        );
        let _ = r.steps;
        let _ = r.n_cand;
    }
    let n = sum_tp + sum_fp + sum_tn + sum_fn;
    let prec = if sum_tp + sum_fp == 0 {
        0.0
    } else {
        sum_tp as f64 / (sum_tp + sum_fp) as f64
    };
    let rec = if sum_tp + sum_fn == 0 {
        0.0
    } else {
        sum_tp as f64 / (sum_tp + sum_fn) as f64
    };
    let agree = if n == 0 {
        0.0
    } else {
        (sum_tp + sum_tn) as f64 / n as f64
    };
    println!("\n  decisions={n}  TP={sum_tp} FP={sum_fp} TN={sum_tn} FN={sum_fn}");
    println!("  site_cost>0 vs live gate: precision={prec:.2} recall={rec:.2} agree={agree:.2}");
    println!("  first-hop site top-1 by site_cost among gate-kept: {top1}/{known}");
}

fn main() {
    let hard = std::env::args().any(|a| a == "--hard");
    let cases = if hard { HARD } else { MID };
    let title = if hard { "HARD" } else { "MID" };
    println!(
        "shell_cost_bench — site_shell_cost Σ norm|current−target| (+leave); keep iff cost>0\n"
    );
    let t0 = Instant::now();
    let rows: Vec<Row> = cases.iter().map(|(n, r, t)| eval_case(n, r, t)).collect();
    print_table(title, &rows);
    println!("\n  wall {:.2}s", t0.elapsed().as_secs_f64());
}
