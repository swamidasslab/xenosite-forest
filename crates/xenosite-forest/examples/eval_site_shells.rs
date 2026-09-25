//! Site selection with **n0, n1, and n2** read separately (not summed).
//!
//! Uses [`AtomNeighborhood::oxy_shape`]: O and H must agree in the same shell.
use xenosite_forest::rules::{dealkylation, hydroxylation, oxygen_reduction};
use xenosite_forest::{
    AtomNeighborhood, OxyShellShape, aligned_shells, atom_diff, candidate_could_help_on,
    format_shell, parse_mol,
};

fn lose_c_n1(env: &AtomNeighborhood) -> bool {
    env.shell_get("n1", "C") < 0
}

fn site_atoms_of(c: &xenosite_forest::Candidate) -> Vec<usize> {
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

fn line(env: &AtomNeighborhood) -> String {
    format!(
        "n0{{{}}} n1{{{}}} n2{{{}}} |δ|={} shape={:?}",
        format_shell(&env.n0),
        format_shell(&env.n1),
        format_shell(&env.n2),
        env.abs_delta(),
        env.oxy_shape()
    )
}

fn main() {
    println!("## Per-shell O+H (n0 / n1 / n2) — alcohol vs carbonyl\n");
    for (a, b, want) in [
        ("CC", "CCO", "Alcohol"),
        ("CC", "CC=O", "Carbonyl"),
        ("CCC", "CC(O)C", "Alcohol"),
        ("CCC", "CCC=O", "Carbonyl"),
        ("c1ccccc1", "Oc1ccccc1", "Alcohol"),
    ] {
        let d = aligned_shells(&parse_mol(a).unwrap(), &parse_mol(b).unwrap());
        let attach: Vec<_> = d
            .atoms
            .iter()
            .filter(|(_, e)| e.oxy_shape() != OxyShellShape::None)
            .map(|(&i, e)| format!("r{i} {}", line(e)))
            .collect();
        println!("  {a}→{b} want={want}\n    {}", attach.join("\n    "));
    }

    println!("\n## at_sites keep vs candidate_could_help\n");
    let mut tp = 0usize;
    let mut tn = 0usize;
    let mut fp = 0usize;
    let mut fn_ = 0usize;
    let mut n = 0usize;

    for (a, b, kind) in [
        ("CC", "CCO", "oh"),
        ("CC", "CC=O", "oh"),
        ("CCC", "CC(O)C", "oh"),
        ("CCC", "CCC=O", "oh"),
        ("c1ccccc1", "Oc1ccccc1", "oh"),
        ("COc1ccccc1", "Oc1ccccc1", "cleave"),
        ("COc1ccc(CCN)cc1OC", "NCCc1ccc(O)c(O)c1", "cleave"),
        ("CC=O", "CCO", "or"),
    ] {
        let ra = parse_mol(a).unwrap();
        let rb = parse_mol(b).unwrap();
        let d = aligned_shells(&ra, &rb);
        let ad = atom_diff(&ra, &rb);
        let set = match kind {
            "oh" => hydroxylation(),
            "cleave" => dealkylation(),
            "or" => oxygen_reduction(),
            _ => continue,
        };
        for c in set.candidates(&ra).collect::<Result<Vec<_>, _>>().unwrap() {
            if kind == "cleave" && !c.pattern.effect.cleaves {
                continue;
            }
            let atoms = site_atoms_of(&c);
            let gate = candidate_could_help_on(&c, &ad, Some(&ra), Some(&rb));
            let site = d.at_sites(&atoms);
            let cost = site.cost();
            let shell_keep = match kind {
                "oh" => site.atoms.values().any(|e| {
                    matches!(
                        e.oxy_shape(),
                        OxyShellShape::Alcohol | OxyShellShape::Carbonyl
                    )
                }),
                "cleave" => d.unaligned_reactant > 0 && site.atoms.values().any(lose_c_n1),
                "or" => cost > 0,
                _ => false,
            };
            n += 1;
            let tag = match (gate, shell_keep) {
                (true, true) => {
                    tp += 1;
                    "TP"
                }
                (false, false) => {
                    tn += 1;
                    "TN"
                }
                (false, true) => {
                    fp += 1;
                    "FP"
                }
                (true, false) => {
                    fn_ += 1;
                    "FN"
                }
            };
            if tag != "TN" {
                let detail: Vec<_> = site
                    .atoms
                    .iter()
                    .map(|(&i, e)| format!("r{i}:{}", line(e)))
                    .collect();
                println!(
                    "  {tag} cost={cost} {a}→{b} {}/{:?}\n       {}",
                    c.pattern.name,
                    atoms,
                    detail.join(" | ")
                );
            }
        }
    }

    println!(
        "\n## Scorecard (per-shell n0/n1/n2)\n  n={n} TP={tp} TN={tn} FP={fp} FN={fn_}\n  precision={:.2} recall={:.2} agree={:.2}",
        if tp + fp == 0 {
            0.0
        } else {
            tp as f64 / (tp + fp) as f64
        },
        if tp + fn_ == 0 {
            0.0
        } else {
            tp as f64 / (tp + fn_) as f64
        },
        if n == 0 {
            0.0
        } else {
            (tp + tn) as f64 / n as f64
        }
    );
}
