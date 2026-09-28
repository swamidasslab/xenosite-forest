//! Site vs full align: cost = Σ|δ| over aromatic + n0/n1/n2 (missing = 0).
use xenosite_forest::rules::{dealkylation, hydroxylation, oxygen_reduction};
use xenosite_forest::{
    AtomNeighborhood, aligned_shells, atom_diff, candidate_could_help_on, format_shell, parse_mol,
    shell_l1,
};

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
    let (n0, n1, n2) = (
        shell_l1(&env.n0, &Default::default()),
        shell_l1(&env.n1, &Default::default()),
        shell_l1(&env.n2, &Default::default()),
    );
    format!(
        "n0{{{}}}={n0} n1{{{}}}={n1} n2{{{}}}={n2} Σ={}",
        format_shell(&env.n0),
        format_shell(&env.n1),
        format_shell(&env.n2),
        env.abs_delta()
    )
}

fn main() {
    println!("## Attachment atoms: per-shell |δ| (missing=0)\n");
    for (a, b) in [
        ("CC", "CCO"),
        ("CC", "CC=O"),
        ("CCC", "CC(O)C"),
        ("CCC", "CCC=O"),
        ("COc1ccccc1", "Oc1ccccc1"),
    ] {
        let d = aligned_shells(&parse_mol(a).unwrap(), &parse_mol(b).unwrap());
        println!("  {a}→{b}  full_cost={}", d.without_unchanged().cost());
        for (&i, env) in &d.atoms {
            if env.is_unchanged() {
                continue;
            }
            println!("    r{i}  {}", line(env));
        }
    }

    println!("\n## at_sites(orbit).cost() vs candidate_could_help\n");
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
            // Simplified keep: any nonzero site |δ|.
            let shell_keep = cost > 0;
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
                    "  {tag} site_cost={cost} {a}→{b} {}/{:?}\n       {}",
                    c.pattern.name,
                    atoms,
                    detail.join(" | ")
                );
            }
        }
    }

    println!(
        "\n## Scorecard (site cost = Σ|δ|, missing=0)\n  n={n} TP={tp} TN={tn} FP={fp} FN={fn_}\n  precision={:.2} recall={:.2} agree={:.2}",
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
