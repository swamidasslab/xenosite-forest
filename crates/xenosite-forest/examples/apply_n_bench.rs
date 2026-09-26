//! ApplyN emit bench: distinct products + covering plan path counts.
//!
//! ```text
//! cargo run -p xenosite-forest --example apply_n_bench --release
//! ```

use std::time::Instant;

use xenosite_forest::{
    ApplyN, RuleSet, apply_n_emit_products, hydroxylation, leaf_rule, phase_one,
};

fn dealkylation() -> RuleSet {
    leaf_rule("Dealkylation").expect("Dealkylation leaf")
}

fn fmt_ms(ms: f64) -> String {
    if ms >= 1000.0 {
        format!("{:.2} s", ms / 1000.0)
    } else if ms >= 1.0 {
        format!("{ms:.2} ms")
    } else {
        format!("{:.1} µs", ms * 1000.0)
    }
}

struct Case {
    name: &'static str,
    smiles: &'static str,
    arms: &'static [&'static str],
    count: u16,
    ruleset: fn() -> xenosite_forest::RuleSet,
}

fn main() {
    let cases = [
        Case {
            name: "benzene OH×1",
            smiles: "c1ccccc1",
            arms: &["Hydroxylation"],
            count: 1,
            ruleset: hydroxylation,
        },
        Case {
            name: "benzene OH×2",
            smiles: "c1ccccc1",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "benzene OH×3",
            smiles: "c1ccccc1",
            arms: &["Hydroxylation"],
            count: 3,
            ruleset: hydroxylation,
        },
        Case {
            name: "toluene OH×2",
            smiles: "Cc1ccccc1",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "ethane OH×2",
            smiles: "CC",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "propane OH×2",
            smiles: "CCC",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "anisole dealk×1",
            smiles: "COc1ccccc1",
            arms: &["Dealkylation"],
            count: 1,
            ruleset: dealkylation,
        },
        Case {
            name: "veratrole dealk×1",
            smiles: "COc1ccc(OC)cc1",
            arms: &["Dealkylation"],
            count: 1,
            ruleset: dealkylation,
        },
        Case {
            name: "ethene PhaseOne×1",
            smiles: "C=C",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "EpoxideOpening",
                "Dehydrogenation",
            ],
            count: 1,
            ruleset: phase_one,
        },
        Case {
            name: "ethene PhaseOne×2",
            smiles: "C=C",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "EpoxideOpening",
                "Dehydrogenation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "phenol OH×1",
            smiles: "Oc1ccccc1",
            arms: &["Hydroxylation"],
            count: 1,
            ruleset: hydroxylation,
        },
        Case {
            name: "n-butylbenzene OH×2",
            smiles: "CCCCc1ccccc1",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
    ];

    println!(
        "{:<28} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10}  sample products",
        "case", "elig", "combo", "prod", "plans", "paths", "time"
    );
    println!("{}", "-".repeat(110));

    for case in &cases {
        let set = (case.ruleset)();
        // Sanity: arms exist when they are leaf names.
        for arm in case.arms {
            let _ = leaf_rule(arm);
        }
        let pool = ApplyN::new(case.arms.iter().copied(), case.count);
        let t0 = Instant::now();
        let (products, stats) = apply_n_emit_products(case.smiles, &set, &pool)
            .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        let ms = t0.elapsed().as_secs_f64() * 1000.0;

        let sample: Vec<_> = products
            .iter()
            .take(3)
            .map(|p| p.smiles.as_str())
            .collect();
        let sample = if products.len() > 3 {
            format!("{}, …", sample.join(", "))
        } else {
            sample.join(", ")
        };

        println!(
            "{:<28} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10}  {sample}",
            case.name,
            stats.n_eligible_sites,
            stats.n_combinations,
            stats.n_products,
            stats.n_plans,
            stats.n_covering_linearizations,
            fmt_ms(ms),
        );

        // Spot-check: every product plan reaches its CSMI.
        for p in &products {
            for plan in &p.plans {
                assert!(
                    plan.reaches(case.smiles, &p.smiles).unwrap(),
                    "{}: plan does not reach {}",
                    case.name,
                    p.smiles
                );
            }
        }
    }
}
