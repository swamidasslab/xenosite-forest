//! Mid/hard: does keep-if-drop (|proj−tgt| decreases) still reach the target?
//!
//! Compares live `find_path` (atom_diff gate) vs a small BFS that only expands
//! edits whose site residual toward the **final** target drops.
use std::collections::{HashSet, VecDeque};

use xenosite_forest::{
    Candidate, FindPathConfig, ForestMol, PathCounters, SiteShellCostOpts, atom_diff, edit_shells,
    find_path_with, molecule_shells, phase_one, site_shell_cost_opts,
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

fn site_atoms_cand(c: &Candidate) -> Vec<usize> {
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

fn residual_drops(
    parent: &ForestMol,
    child: &ForestMol,
    target: &ForestMol,
    atoms: &[usize],
) -> bool {
    let opts = SiteShellCostOpts { dearomatic: false };
    // Map parent→target by tags when possible; else MCS atom_diff.
    let ad = atom_diff(parent.mol(), target.mol());
    let cur = molecule_shells(parent.mol());
    let tgt = molecule_shells(target.mol());
    let before = site_shell_cost_opts(&cur, None, &tgt, &ad.mapping, atoms, opts);
    let edit = edit_shells(parent, child);
    let after = site_shell_cost_opts(&cur, Some(&edit), &tgt, &ad.mapping, atoms, opts);
    before > after + 1e-12
}

fn bfs_keep_if_drop(reactant: &str, target: &str, max_nodes: usize) -> Option<usize> {
    let set = phase_one();
    let start = ForestMol::parse(reactant).unwrap();
    let goal = ForestMol::parse(target).unwrap();
    let goal_key = goal
        .stable_csmi_key()
        .unwrap_or_else(|| goal.csmi().clone());

    let mut seen = HashSet::new();
    seen.insert(
        start
            .stable_csmi_key()
            .unwrap_or_else(|| start.csmi().clone()),
    );
    let mut q: VecDeque<(ForestMol, usize)> = VecDeque::new();
    q.push_back((start, 0));
    let mut expanded = 0usize;

    while let Some((parent, depth)) = q.pop_front() {
        let pk = parent
            .stable_csmi_key()
            .unwrap_or_else(|| parent.csmi().clone());
        if pk == goal_key {
            return Some(depth);
        }
        if expanded >= max_nodes {
            break;
        }
        expanded += 1;

        let mol = parent.mol();
        for c in set.candidates(mol).collect::<Result<Vec<_>, _>>().unwrap() {
            let atoms = site_atoms_cand(&c);
            let Ok(pieces) = c.materialize_mols(mol) else {
                continue;
            };
            for piece in pieces {
                let child = parent.adopt_product(piece);
                if !residual_drops(&parent, &child, &goal, &atoms) {
                    continue;
                }
                let key = child
                    .stable_csmi_key()
                    .unwrap_or_else(|| child.csmi().clone());
                if seen.insert(key) {
                    q.push_back((child, depth + 1));
                }
            }
        }
        for p in set
            .pair_candidates(mol)
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
        {
            let mut atoms = p.plan_site_atoms();
            atoms.sort_unstable();
            atoms.dedup();
            let Ok(pieces) = p.materialize_mols(mol) else {
                continue;
            };
            for piece in pieces {
                let child = parent.adopt_product(piece);
                if !residual_drops(&parent, &child, &goal, &atoms) {
                    continue;
                }
                let key = child
                    .stable_csmi_key()
                    .unwrap_or_else(|| child.csmi().clone());
                if seen.insert(key) {
                    q.push_back((child, depth + 1));
                }
            }
        }
    }
    None
}

fn live_find(reactant: &str, target: &str) -> Option<usize> {
    let set = phase_one();
    let mut counters = PathCounters::default();
    let config = FindPathConfig {
        max_nodes: 2000,
        max_paths: 1,
        use_atom_diff: true,
        ..FindPathConfig::default()
    };
    let hits = find_path_with(reactant, target, &set, &mut counters, config, |_| true)
        .unwrap()
        .collect_all()
        .unwrap();
    hits.first().map(|h| h.steps.len())
}

fn run(label: &str, cases: &[(&str, &str, &str)]) {
    println!("## {label}");
    println!(
        "{:<28} {:>8} {:>10} {:>12}",
        "case", "live", "drop-BFS", "same?"
    );
    for &(name, reactant, target) in cases {
        let live = live_find(reactant, target);
        let drop = bfs_keep_if_drop(reactant, target, 5000);
        let same = match (live, drop) {
            (Some(_), Some(_)) => "yes",
            (None, None) => "both-miss",
            (Some(_), None) => "DROP-MISS",
            (None, Some(_)) => "live-miss",
        };
        println!(
            "{:<28} {:>8} {:>10} {:>12}",
            name,
            live.map(|h| format!("{h}h"))
                .unwrap_or_else(|| "MISS".into()),
            drop.map(|h| format!("{h}h"))
                .unwrap_or_else(|| "MISS".into()),
            same
        );
    }
}

fn main() {
    println!(
        "live = find_path (atom_diff gate)\n\
         drop-BFS = only expand edits that lower Σ|proj−tgt| toward final target\n"
    );
    run("MID", MID);
    run("HARD", HARD);
}
