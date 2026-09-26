//! ApplyN emit bench: distinct products + covering plan path counts.
//!
//! Covers Hydroxylation multi-k plus a span of Phase I leaves (oxidation,
//! reduction, cleavage, hydrolysis, mixed OR pools).
//!
//! ```text
//! cargo run -p xenosite-forest --example apply_n_bench --release
//! ```

use std::time::Instant;

use xenosite_forest::{
    ApplyN, RuleSet, apply_n_emit_products, hydroxylation, leaf_rule, phase_one,
};

fn dealkylation() -> RuleSet {
    leaf_rule("Dealkylation").expect("Dealkylation")
}
fn n_dealkylation() -> RuleSet {
    leaf_rule("NDealkylation").expect("NDealkylation")
}
fn sulfur_oxidation() -> RuleSet {
    leaf_rule("SulfurOxidation").expect("SulfurOxidation")
}
fn nitrogen_oxidation() -> RuleSet {
    leaf_rule("NitrogenOxidation").expect("NitrogenOxidation")
}
fn dehydrogenation() -> RuleSet {
    leaf_rule("Dehydrogenation").expect("Dehydrogenation")
}
fn epoxidation() -> RuleSet {
    leaf_rule("Epoxidation").expect("Epoxidation")
}
fn epoxide_hydration() -> RuleSet {
    leaf_rule("EpoxideHydration").expect("EpoxideHydration")
}
fn epoxide_opening() -> RuleSet {
    leaf_rule("EpoxideOpening").expect("EpoxideOpening")
}
fn hydrogenation() -> RuleSet {
    leaf_rule("Hydrogenation").expect("Hydrogenation")
}
fn oxidative_dehalogenation() -> RuleSet {
    leaf_rule("OxidativeDehalogenation").expect("OxidativeDehalogenation")
}
fn reductive_dehalogenation() -> RuleSet {
    leaf_rule("ReductiveDehalogenation").expect("ReductiveDehalogenation")
}
fn hydrolysis() -> RuleSet {
    leaf_rule("Hydrolysis").expect("Hydrolysis")
}
fn dehydration() -> RuleSet {
    leaf_rule("Dehydration").expect("Dehydration")
}
fn sulfur_reduction() -> RuleSet {
    leaf_rule("SulfurReduction").expect("SulfurReduction")
}
fn nitrogen_reduction() -> RuleSet {
    leaf_rule("NitrogenReduction").expect("NitrogenReduction")
}
fn oxygen_reduction() -> RuleSet {
    leaf_rule("OxygenReduction").expect("OxygenReduction")
}
fn benzodioxole_reduction() -> RuleSet {
    leaf_rule("BenzodioxoleReduction").expect("BenzodioxoleReduction")
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
    ruleset: fn() -> RuleSet,
}

fn main() {
    let cases = [
        // --- Hydroxylation multi-k (Aut orbits / combo scaling) ---
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
        // --- Cleavage ---
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
            name: "DMA NDealk×1",
            smiles: "CN(C)c1ccccc1",
            arms: &["NDealkylation"],
            count: 1,
            ruleset: n_dealkylation,
        },
        // --- Heteroatom oxidation ---
        Case {
            name: "EtSH SOx×1",
            smiles: "CCS",
            arms: &["SulfurOxidation"],
            count: 1,
            ruleset: sulfur_oxidation,
        },
        Case {
            name: "Me2S SOx×1",
            smiles: "CSC",
            arms: &["SulfurOxidation"],
            count: 1,
            ruleset: sulfur_oxidation,
        },
        Case {
            name: "MeNH2 NOx×1",
            smiles: "CN",
            arms: &["NitrogenOxidation"],
            count: 1,
            ruleset: nitrogen_oxidation,
        },
        Case {
            name: "EtNH2 NOx×1",
            smiles: "CCN",
            arms: &["NitrogenOxidation"],
            count: 1,
            ruleset: nitrogen_oxidation,
        },
        // --- Dehydrogenation / alkene oxygenations ---
        Case {
            name: "EtOH DH×1",
            smiles: "CCO",
            arms: &["Dehydrogenation"],
            count: 1,
            ruleset: dehydrogenation,
        },
        Case {
            name: "iPrOH DH×1",
            smiles: "CC(O)C",
            arms: &["Dehydrogenation"],
            count: 1,
            ruleset: dehydrogenation,
        },
        Case {
            name: "ethene epox×1",
            smiles: "C=C",
            arms: &["Epoxidation"],
            count: 1,
            ruleset: epoxidation,
        },
        Case {
            name: "ethene hyd×1",
            smiles: "C=C",
            arms: &["EpoxideHydration"],
            count: 1,
            ruleset: epoxide_hydration,
        },
        Case {
            name: "oxirane open×1",
            smiles: "C1OC1",
            arms: &["EpoxideOpening"],
            count: 1,
            ruleset: epoxide_opening,
        },
        Case {
            name: "2-butene epox×1",
            smiles: "C/C=C/C",
            arms: &["Epoxidation"],
            count: 1,
            ruleset: epoxidation,
        },
        // --- Reduction / dehalogenation ---
        Case {
            name: "acetylene H2×1",
            smiles: "C#C",
            arms: &["Hydrogenation"],
            count: 1,
            ruleset: hydrogenation,
        },
        Case {
            name: "ClPh oxdehal×1",
            smiles: "Clc1ccccc1",
            arms: &["OxidativeDehalogenation"],
            count: 1,
            ruleset: oxidative_dehalogenation,
        },
        Case {
            name: "ClPh reddehal×1",
            smiles: "Clc1ccccc1",
            arms: &["ReductiveDehalogenation"],
            count: 1,
            ruleset: reductive_dehalogenation,
        },
        Case {
            name: "DMSO Sred×1",
            smiles: "CS(=O)C",
            arms: &["SulfurReduction"],
            count: 1,
            ruleset: sulfur_reduction,
        },
        Case {
            name: "PhNO Nred×1",
            smiles: "O=Nc1ccccc1",
            arms: &["NitrogenReduction"],
            count: 1,
            ruleset: nitrogen_reduction,
        },
        Case {
            name: "anthraquinone Ored×1",
            smiles: "c1ccc2c(c1)C(=O)c1ccccc1C2=O",
            arms: &["OxygenReduction"],
            count: 1,
            ruleset: oxygen_reduction,
        },
        Case {
            name: "benzodioxole red×1",
            smiles: "c1ccc2c(c1)OCO2",
            arms: &["BenzodioxoleReduction"],
            count: 1,
            ruleset: benzodioxole_reduction,
        },
        // --- Hydrolysis / dehydration ---
        Case {
            name: "MeOAc hydro×1",
            smiles: "CC(=O)OC",
            arms: &["Hydrolysis"],
            count: 1,
            ruleset: hydrolysis,
        },
        Case {
            name: "EtOH dehyd×1",
            smiles: "CCO",
            arms: &["Dehydration"],
            count: 1,
            ruleset: dehydration,
        },
        // --- Mixed PhaseOne OR pools (multi-arm / multi-count) ---
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
            name: "anisole PhaseOne×1",
            smiles: "COc1ccccc1",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dealkylation",
                "Dehydrogenation",
            ],
            count: 1,
            ruleset: phase_one,
        },
        Case {
            name: "toluene PhaseOne×2",
            smiles: "Cc1ccccc1",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dehydrogenation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "EtSH PhaseOne×2",
            smiles: "CCS",
            arms: &[
                "Hydroxylation",
                "SulfurOxidation",
                "Dehydrogenation",
                "Dealkylation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "ClPh PhaseOne×1",
            smiles: "Clc1ccccc1",
            arms: &[
                "Hydroxylation",
                "OxidativeDehalogenation",
                "ReductiveDehalogenation",
                "Epoxidation",
            ],
            count: 1,
            ruleset: phase_one,
        },
        Case {
            name: "2-butene PhaseOne×2",
            smiles: "C/C=C/C",
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
            name: "phenol PhaseOne×2",
            smiles: "Oc1ccccc1",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dehydrogenation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "EtOH OH+DH x2",
            smiles: "CCO",
            arms: &["Hydroxylation", "Dehydrogenation"],
            count: 2,
            ruleset: phase_one,
        },
        // --- Harder / larger scaffolds ---
        Case {
            name: "styrene P1 x2",
            smiles: "C=Cc1ccccc1",
            arms: &[
                "Epoxidation",
                "EpoxideHydration",
                "EpoxideOpening",
                "Hydroxylation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "thioanisole SOx x1",
            smiles: "CSc1ccccc1",
            arms: &["SulfurOxidation"],
            count: 1,
            ruleset: sulfur_oxidation,
        },
        Case {
            name: "indole OH x2",
            smiles: "c1ccc2[nH]ccc2c1",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "coumarin OH x2",
            smiles: "O=c1ccc2ccccc2o1",
            arms: &["Hydroxylation"],
            count: 2,
            ruleset: hydroxylation,
        },
        Case {
            name: "eugenol P1 x2",
            smiles: "COc1cc(CC=C)ccc1O",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dealkylation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "trimethoxyPEA dealk x2",
            smiles: "COc1cc(OC)c(OC)c(CCN)c1",
            arms: &["Dealkylation"],
            count: 2,
            ruleset: dealkylation,
        },
        Case {
            name: "dimethoxyNaph dealk x2",
            smiles: "COc1ccc2c(OC)cccc2c1",
            arms: &["Dealkylation"],
            count: 2,
            ruleset: dealkylation,
        },
        Case {
            name: "allyl-veratrole P1 x2",
            smiles: "COc1ccc(CC=C)c(OC)c1",
            arms: &[
                "Epoxidation",
                "EpoxideHydration",
                "Dealkylation",
                "Hydroxylation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "ClPh P1 x2",
            smiles: "Clc1ccccc1",
            arms: &[
                "OxidativeDehalogenation",
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "BuPh P1 x2",
            smiles: "CCCCc1ccccc1",
            arms: &[
                "Hydroxylation",
                "Epoxidation",
                "EpoxideHydration",
                "Dehydrogenation",
            ],
            count: 2,
            ruleset: phase_one,
        },
        Case {
            name: "phenacetin dealk x1",
            smiles: "CCOc1ccc(NC(C)=O)cc1",
            arms: &["Dealkylation"],
            count: 1,
            ruleset: dealkylation,
        },
        Case {
            name: "anthraquinone Ored x1",
            smiles: "c1ccc2c(c1)C(=O)c1ccccc1C2=O",
            arms: &["OxygenReduction"],
            count: 1,
            ruleset: oxygen_reduction,
        },
    ];

    println!(
        "{:<28} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10}  sample products",
        "case", "elig", "combo", "prod", "plans", "paths", "time"
    );
    println!("{}", "-".repeat(110));

    for case in &cases {
        let set = (case.ruleset)();
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

        assert!(
            stats.n_products > 0,
            "{}: expected ≥1 product (elig={} combo={})",
            case.name,
            stats.n_eligible_sites,
            stats.n_combinations
        );

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
