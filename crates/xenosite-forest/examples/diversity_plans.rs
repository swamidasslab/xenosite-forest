//! Compare diversity off vs −n: bills + emitted step plans on hard multipath.
//!
//! ```text
//! cargo run -p xenosite-forest --example diversity_plans --release
//! ```

use std::collections::BTreeSet;

use xenosite_forest::{FindPathConfig, PathCounters, PathOutcome, find_path_with, phase_one};

const MAX_NODES: usize = 800;
const MAX_PATHS: usize = 10;

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
        "COc1ccc(CC=C)c(OC)c1OC",
        "O=C1C=C(CC=C)C(=O)C(O)=C1",
    ),
    (
        "tetraMeO-naph→polyOH-NQ",
        "COc1cc(OC)c2c(OC)cc(OC)cc2c1",
        "O=C1C=C(O)C(=O)c2c(O)cc(O)cc12",
    ),
];

fn step_sig(hit: &PathOutcome) -> String {
    hit.plan
        .iter()
        .map(|s| {
            let site: Vec<_> = s
                .site
                .iter()
                .map(|a| match a {
                    xenosite_forest::PlanAtom::Label(t) => t.0.to_string(),
                    other => format!("{other:?}"),
                })
                .collect();
            format!("{}@[{}]", s.rule, site.join(","))
        })
        .collect::<Vec<_>>()
        .join(" → ")
}

fn rule_multiset(hit: &PathOutcome) -> String {
    let mut rules: Vec<_> = hit.plan.iter().map(|s| s.rule.as_str()).collect();
    rules.sort_unstable();
    rules.join(",")
}

struct Run {
    label: &'static str,
    diversity: bool,
    bill: usize,
    drop_dup: usize,
    drop_exact: usize,
    drop_skel: usize,
    diversity_repush: usize,
    hits: usize,
    unique_step_sigs: usize,
    unique_rule_bags: usize,
    same_lin_pairs: usize,
    skeleton_pairs: usize,
    plans: Vec<String>,
}

fn run_case(start: &str, target: &str, diversity: bool) -> Run {
    let mut counters = PathCounters::default();
    let hits = find_path_with(
        start,
        target,
        &phase_one(),
        &mut counters,
        FindPathConfig {
            max_paths: MAX_PATHS,
            max_nodes: MAX_NODES,
            diversity,
            ..FindPathConfig::default()
        },
        |_| true,
    )
    .unwrap()
    .collect_all()
    .unwrap();
    let mut step_sigs = BTreeSet::new();
    let mut rule_bags = BTreeSet::new();
    let mut same_lin_pairs = 0usize;
    let mut skeleton_pairs = 0usize;
    for i in 0..hits.len() {
        step_sigs.insert(step_sig(&hits[i]));
        rule_bags.insert(rule_multiset(&hits[i]));
        for j in (i + 1)..hits.len() {
            if hits[i].plan.same_linearizations(&hits[j].plan) {
                same_lin_pairs += 1;
            }
            if hits[i].plan.same_rule_maybe_skeleton(&hits[j].plan) {
                skeleton_pairs += 1;
            }
        }
    }
    let plans: Vec<String> = hits
        .iter()
        .enumerate()
        .map(|(i, h)| format!("  [{i}] {}", step_sig(h)))
        .collect();
    Run {
        label: if diversity { "−n" } else { "off" },
        diversity,
        bill: counters.billed(),
        drop_dup: counters.dropped_duplicate_plan,
        drop_exact: counters.dropped_exact_plan,
        drop_skel: counters.dropped_skeleton_twin,
        diversity_repush: counters.diversity_repush,
        hits: hits.len(),
        unique_step_sigs: step_sigs.len(),
        unique_rule_bags: rule_bags.len(),
        same_lin_pairs,
        skeleton_pairs,
        plans,
    }
}

fn main() {
    println!(
        "diversity_plans  hard  paths={MAX_PATHS}  nodes={MAX_NODES}  \
         penalty=fixed−n vs off\n"
    );
    println!(
        "{:<28} {:>4} {:>6} {:>5} {:>5} {:>5} {:>5} {:>6} {:>5} {:>5} {:>5} {:>7}",
        "case",
        "mode",
        "bill",
        "hits",
        "uniq",
        "bags",
        "linP",
        "skelP",
        "dup",
        "ex",
        "sk",
        "repush"
    );

    let mut tot_off_bill = 0usize;
    let mut tot_on_bill = 0usize;
    let mut tot_off_dup = 0usize;
    let mut tot_on_dup = 0usize;
    let mut tot_off_uniq = 0usize;
    let mut tot_on_uniq = 0usize;
    let mut tot_off_skel_p = 0usize;
    let mut tot_on_skel_p = 0usize;

    for &(name, start, target) in HARD {
        let off = run_case(start, target, false);
        let on = run_case(start, target, true);
        for r in [&off, &on] {
            println!(
                "{:<28} {:>4} {:>6} {:>5} {:>5} {:>5} {:>5} {:>6} {:>5} {:>5} {:>5} {:>7}",
                if r.diversity { "" } else { name },
                r.label,
                r.bill,
                r.hits,
                r.unique_step_sigs,
                r.unique_rule_bags,
                r.same_lin_pairs,
                r.skeleton_pairs,
                r.drop_dup,
                r.drop_exact,
                r.drop_skel,
                r.diversity_repush,
            );
        }
        println!("  plans off:");
        for p in &off.plans {
            println!("{p}");
        }
        println!("  plans −n:");
        for p in &on.plans {
            println!("{p}");
        }
        println!();
        tot_off_bill += off.bill;
        tot_on_bill += on.bill;
        tot_off_dup += off.drop_dup;
        tot_on_dup += on.drop_dup;
        tot_off_uniq += off.unique_step_sigs;
        tot_on_uniq += on.unique_step_sigs;
        tot_off_skel_p += off.skeleton_pairs;
        tot_on_skel_p += on.skeleton_pairs;
    }

    println!("=== totals ===");
    println!(
        "{:<28} {:>4} {:>6} {:>5} {:>5} {:>5} {:>5} {:>6}",
        "", "mode", "bill", "", "uniq", "", "", "skelP"
    );
    println!(
        "{:<28} {:>4} {:>6} {:>5} {:>5} {:>5} {:>5} {:>6}  drop_dup={}",
        "Σ hard", "off", tot_off_bill, "", tot_off_uniq, "", "", tot_off_skel_p, tot_off_dup
    );
    println!(
        "{:<28} {:>4} {:>6} {:>5} {:>5} {:>5} {:>5} {:>6}  drop_dup={}",
        "", "−n", tot_on_bill, "", tot_on_uniq, "", "", tot_on_skel_p, tot_on_dup
    );
    println!(
        "\nuniq = distinct step@site sequences among yielded hits; \
         bags = distinct sorted rule multisets; \
         linP/skelP = pairwise same_linearizations / same_rule_maybe_skeleton among hits; \
         dup/ex/sk = PathCounters yield drops."
    );
}
