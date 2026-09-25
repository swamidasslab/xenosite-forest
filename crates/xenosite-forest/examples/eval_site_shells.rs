//! Site selection with **n0, n1, and n2** read separately (not summed).
//!
//! O lands in n1 at the attachment atom and in n2 at its neighbor. H for
//! alcohol vs carbonyl is the H delta **in the same shell as the O**.
use xenosite_forest::rules::{dealkylation, hydroxylation, oxygen_reduction};
use xenosite_forest::{
    AtomNeighborhood, Shell, aligned_shells, atom_diff, candidate_could_help_on, format_shell,
    parse_mol,
};

fn get(shell: &Shell, el: &str) -> i32 {
    shell.get(el).copied().unwrap_or(0)
}

fn shells(env: &AtomNeighborhood) -> [(&str, &Shell); 3] {
    [("n0", &env.n0), ("n1", &env.n1), ("n2", &env.n2)]
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OxyShape {
    /// Some shell has O:+1 and H:−1 in that same shell.
    Alcohol,
    /// Some shell has O:+1 and H:≤−2 in that same shell.
    Carbonyl,
    /// O:+ in a shell without a clean H signature.
    OtherO,
    None,
}

/// Classify from n0/n1/n2: O and H must agree **within one shell**.
fn oxy_shape(env: &AtomNeighborhood) -> OxyShape {
    let mut best = OxyShape::None;
    for (_, sh) in shells(env) {
        let o = get(sh, "O");
        if o <= 0 {
            continue;
        }
        let h = get(sh, "H");
        let shape = if h == -1 {
            OxyShape::Alcohol
        } else if h <= -2 {
            OxyShape::Carbonyl
        } else {
            OxyShape::OtherO
        };
        // Prefer a clean alcohol/carbonyl over OtherO.
        best = match (best, shape) {
            (OxyShape::None, s) => s,
            (OxyShape::OtherO, s) if s != OxyShape::OtherO => s,
            (b, _) => b,
        };
    }
    best
}

/// Cleavage leave mark: C lost in **n1** (not only n2 bleed to neighbors).
fn lose_c_n1(env: &AtomNeighborhood) -> bool {
    get(&env.n1, "C") < 0
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
        oxy_shape(env)
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
            .filter(|(_, e)| shells(e).iter().any(|(_, sh)| get(sh, "O") > 0))
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
                // O-add: some site/orbit atom has O in any of n0/n1/n2 with
                // a same-shell H signature (alcohol or carbonyl).
                "oh" => site
                    .atoms
                    .values()
                    .any(|e| matches!(oxy_shape(e), OxyShape::Alcohol | OxyShape::Carbonyl)),
                // Cleave: n1 loses C (leave bond atom), not n2-only bleed.
                "cleave" => d.unaligned_reactant > 0 && site.atoms.values().any(lose_c_n1),
                // Reduction: any shell change on the site.
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
