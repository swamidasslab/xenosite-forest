//! Phase I / conjugation reaction rules as [`PatternInfo`] data.
//!
//! Ported from `xenosite.forest.rules`. Catalogs (`phase_one`, `default_ruleset`)
//! are nested [`RuleSet`] namespaces. Pair-endpoint patterns use
//! [`SiteKind::AtomPair`] + [`Edit::PairEndpoint`]; atom/bond SMIRKS metabolize
//! through the generic door.

use std::collections::BTreeMap;

use crate::pattern::{Edit, Effect, PatternInfo, SiteKind};
use crate::ruleset::RuleSet;

fn smirks_row(
    name: &str,
    smirks: &str,
    site_kind: SiteKind,
    site_map: Vec<u16>,
    effect: Effect,
) -> PatternInfo {
    let smarts = smirks.split(">>").next().unwrap_or(smirks).to_string();
    PatternInfo {
        name: name.into(),
        smarts,
        site_kind,
        site_map,
        edit: Edit::Smirks(smirks.into()),
        effect: effect.sealed(),
        possibilities: Vec::new(),
        skip_same_rings: false,
        cleave_side_group: None,
        search_bias: 0,
    }
}

fn endpoint_row(
    name: &str,
    smarts: &str,
    site_map: Vec<u16>,
    effect: Effect,
    pair_edit: &str,
) -> PatternInfo {
    let skip_same_rings = matches!(name, "single_to_double" | "iminium" | "dealkylate");
    PatternInfo {
        name: name.into(),
        smarts: smarts.into(),
        site_kind: SiteKind::AtomPair,
        site_map,
        edit: Edit::PairEndpoint(pair_edit.into()),
        effect: effect.sealed(),
        possibilities: Vec::new(),
        skip_same_rings,
        cleave_side_group: None,
        search_bias: 0,
    }
}

/// Halogen atomic numbers used by dehalogenation / replace_halogen Whens.
const HALIDE_Z: &[(u8, &str)] = &[(9, "F"), (17, "Cl"), (35, "Br"), (53, "I"), (85, "At")];

fn halide_remove_branches(map: u16, base: Effect) -> Vec<Effect> {
    HALIDE_Z
        .iter()
        .map(|&(z, sym)| {
            let mut effect = Effect {
                when: Some(crate::pattern::When::atomic(map, z)),
                partner: Some(sym.into()),
                ..base.clone()
            };
            // Cleaving: partner stays as a product fragment → leave.
            // Non-cleaving: partner is eliminated → junction removes.
            if effect.cleaves {
                effect.leave_formula = BTreeMap::from([(sym.to_string(), 1)]);
                effect.leave_count = Some(1);
            } else {
                effect.removes = Some(sym.into());
            }
            effect.sealed()
        })
        .collect()
}

/// `Hydroxylation` from Python `xenosite.forest.rules`.
pub fn hydroxylation() -> RuleSet {
    RuleSet::new(
        Some("Hydroxylation".into()),
        [
            PatternInfo::hydroxyl("h", "[#6h1:1]"),
            PatternInfo::hydroxyl("h2", "[#6h2,#6h3:1]"),
        ],
    )
}

/// `Dehydrogenation` from Python `xenosite.forest.rules`.
pub fn dehydrogenation() -> RuleSet {
    RuleSet::new(
        Some("Dehydrogenation".into()),
        [
            smirks_row(
                "sulfoxide",
                "[#16v4:1]-[#8H1:2]>>[*:1]=[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("HH".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "alcohol",
                // chematic: bare `h` is inert; digitize like Hydroxylation partition.
                "[#6h1,#6h2,#6h3:1]-[#8H1:2]>>[*:1]=[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("HH".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "amine",
                "[#6h1,#6h2,#6h3:1]-[#7D1H2,#7D2H1:2]>>[*:1]=[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("HH".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "alkyl",
                "[#6h1,#6h2,#6h3:1]-[#6D1H3,#6D2H2,#6D3H1:2]>>[*:1]=[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("HH".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            endpoint_row(
                "phenol_end",
                "[#6:1]-[#8H:2]",
                vec![2],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "single_to_double",
            ),
            endpoint_row(
                "amine_end",
                "[#6:1]-[#7D1H2,#7D2H1:2]",
                vec![2],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "single_to_double",
            ),
            endpoint_row(
                "methide_end",
                "[#6:1]-[#6D1H3,#6D2H2,#6D3H1:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: Some("C".into()),
                    ..Default::default()
                },
                "single_to_double",
            ),
        ],
    )
}

/// `QuinoneFormation` from Python `xenosite.forest.rules`.
pub fn quinone_formation() -> RuleSet {
    RuleSet::new(
        Some("QuinoneFormation".into()),
        [
            endpoint_row(
                "single_to_double",
                "[#6R:1][#8H,#7D1H2,#7D2H1:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    exclusive_partner: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "single_to_double",
            ),
            endpoint_row(
                "methide_end",
                "[#6R:1][#6D1H3,#6D2H2,#6D3H1:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: Some("C".into()),
                    ..Default::default()
                },
                "single_to_double",
            ),
            endpoint_row(
                "add_carbonyl_o",
                "[#6D2H1;R:1]",
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "add_carbonyl_o",
            ),
            endpoint_row(
                "replace_halogen",
                "[#6H0R:1]-[F,Cl,Br,I:2]",
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    exclusive_partner: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "replace_halogen",
            )
            .with_possibilities(halide_remove_branches(
                2,
                Effect {
                    adds: Some("O".into()),
                    dearomatizes: true,
                    exclusive_partner: true,
                    ..Default::default()
                },
            )),
            endpoint_row(
                "iminium",
                "[#6H0R:1][#7D3:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: false,
                    methide: false,
                    exclusive_partner: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: Some("N".into()),
                    ..Default::default()
                },
                "iminium",
            ),
            endpoint_row(
                "dealkylate",
                "[#6R:1][#7,#8:2][#6:3]",
                vec![1],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    exclusive_partner: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "dealkylate",
            ),
        ],
    )
    .with_canonical_plan(crate::canonical_plan::quinone_canonical_plan)
}

/// `Dealkylation` from Python `xenosite.forest.rules`.
pub fn dealkylation() -> RuleSet {
    RuleSet::new(
        Some("Dealkylation".into()),
        [
            smirks_row(
                "methyl_carboxylic",
                "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methyl_carbonyl",
                "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methyl_alcohol",
                "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methylene_carboxylic",
                "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methylene_carbonyl",
                "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methylene_alcohol",
                "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methine_carbonyl",
                "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methine_alcohol",
                "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "quaternary_alcohol",
                "[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "cc_quaternary_alcohol",
                "[#6H0:1][#6:2]>>(O-[*:1].[*:2])",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "cc_alcohol",
                "[#6h1,#6h2,#6h3:1][#6:2]>>(O-[*:1].[*:2])",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "cc_carbonyl",
                "[#6h1,#6h2,#6h3:1][#6:2]>>(O=[*:1].[*:2])",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "hemiaminal",
                "[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `NDealkylation` from Python `xenosite.forest.rules`.
pub fn n_dealkylation() -> RuleSet {
    RuleSet::new(
        Some("NDealkylation".into()),
        [
            smirks_row(
                "methyl_carboxylic",
                "[#6H3:1][#7:2]>>([*:2].[*:1](=O)O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methyl_carbonyl",
                "[#6H3:1][#7:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methyl_alcohol",
                "[#6H3:1][#7:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: Some(1),
                    partner: None,
                    ..Default::default()
                },
            )
            .with_cleave_side_group("Me", "hetero"),
            smirks_row(
                "methylene_carboxylic",
                "[#6H2:1][#7:2]>>([*:2].[*:1](=O)O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methylene_carbonyl",
                "[#6H2:1][#7:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methylene_alcohol",
                "[#6H2:1][#7:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methine_carbonyl",
                "[#6H1:1][#7:2]>>([*:2].[*:1]=O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "methine_alcohol",
                "[#6H1:1][#7:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "quaternary_alcohol",
                "[#6H0:1][#7:2]>>([*:2].[*:1]-O)",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "hemiaminal",
                "[#8H1:3]-[#6:1]-[#7:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `AzoSplitting` from Python `xenosite.forest.rules`.
pub fn azo_splitting() -> RuleSet {
    RuleSet::new(
        Some("AzoSplitting".into()),
        [smirks_row(
            "azo",
            "[#7:1]=,:[#7:2]>>[*:1].[*:2]",
            SiteKind::Bond,
            vec![1, 2],
            Effect {
                adds: None,
                removes: None,
                cleaves: true,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )
        .with_cleave_side_group("azo", "azo")],
    )
}

/// `BenzodioxoleReduction` from Python `xenosite.forest.rules`.
pub fn benzodioxole_reduction() -> RuleSet {
    RuleSet::new(
        Some("BenzodioxoleReduction".into()),
        [smirks_row(
            "dioxole_methylene",
            "[#6R:1]-[#8R:2]-[#6H2R:3]-[#8R:4]-[#6R:5]>>([*:1]-[*:2].[*:3].[*:4]-[*:5])",
            SiteKind::DirectedBond,
            vec![2, 3],
            Effect {
                cleaves: true,
                leave_count: Some(1),
                leave_formula: crate::pattern::leave_ch2(),
                partner: Some("O".into()),
                ..Default::default()
            },
        )],
    )
}

/// `NitroaromaticReduction` from Python `xenosite.forest.rules`.
pub fn nitroaromatic_reduction() -> RuleSet {
    RuleSet::new(
        Some("NitroaromaticReduction".into()),
        [
            smirks_row(
                "nitro_charged",
                "[#8-1:1]-[#7+1:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    cleaves: true,
                    leave_count: Some(1),
                    leave_formula: crate::pattern::leave_o(),
                    partner: Some("N".into()),
                    ..Default::default()
                },
            ),
            smirks_row(
                "nitro_neutral",
                "[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    cleaves: true,
                    leave_count: Some(1),
                    leave_formula: crate::pattern::leave_o(),
                    partner: Some("N".into()),
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `ThiopheneSulfurOxidation` from Python `xenosite.forest.rules`.
pub fn thiophene_sulfur_oxidation() -> RuleSet {
    RuleSet::new(
        Some("ThiopheneSulfurOxidation".into()),
        [smirks_row(
            "thiophene_s_oxide",
            "[#6:2]1=,:[#6:3][#6:4]=,:[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*&H0&+:1]1[O-]",
            SiteKind::Atom,
            vec![1],
            Effect {
                adds: Some("O".into()),
                removes: None,
                cleaves: false,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )],
    )
}

/// `Dephosphorylation` from Python `xenosite.forest.rules`.
pub fn dephosphorylation() -> RuleSet {
    RuleSet::new(
        Some("Dephosphorylation".into()),
        [smirks_row(
            "phosphate_ester",
            "[#8;$([#8][#6]):1][#15:2](=[#8:3])([#8:4])[#8:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]",
            SiteKind::Atom,
            vec![1],
            Effect {
                adds: None,
                removes: None,
                cleaves: true,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )],
    )
}

/// `EpoxideOpening` from Python `xenosite.forest.rules`.
pub fn epoxide_opening() -> RuleSet {
    RuleSet::new(
        Some("EpoxideOpening".into()),
        [
            smirks_row(
                "rearrange",
                "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    // Open chain saturates: C2H4O → C2H6O (+2H).
                    adds: Some("HH".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "hydrate",
                "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)",
                SiteKind::Atom,
                vec![1],
                Effect {
                    // Vicinal diol: +O and the two H that come with opening/OH.
                    adds: Some("OHH".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `Hydrolysis` from Python `xenosite.forest.rules`.
pub fn hydrolysis() -> RuleSet {
    RuleSet::new(
        Some("Hydrolysis".into()),
        [
            smirks_row(
                "add_water",
                "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2](O).[*:3])",
                SiteKind::Bond,
                vec![2, 3],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "cleave",
                "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2].[*:3])",
                SiteKind::Bond,
                vec![2, 3],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `Dehydration` from Python `xenosite.forest.rules`.
pub fn dehydration() -> RuleSet {
    let o_leave = Effect {
        cleaves: true,
        leave_count: Some(1),
        leave_formula: crate::pattern::leave_o(),
        partner: Some("O".into()),
        ..Default::default()
    };
    RuleSet::new(
        Some("Dehydration".into()),
        [
            smirks_row(
                "alcohol",
                "[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                o_leave.clone(),
            ),
            smirks_row(
                "beta_elimination",
                "[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]",
                SiteKind::Atom,
                // Alcohol carbon + adjacent carbon — site set differs from a
                // lone hydroxylation site, so OH→beta-elim is not circular.
                vec![1, 3],
                o_leave.clone(),
            ),
            smirks_row(
                "carbonyl",
                "[#6,#7:1]=[#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                o_leave,
            ),
        ],
    )
}

/// `Hydrogenation` from Python `xenosite.forest.rules`.
///
/// Patterns carry `search_bias = -1` via [`demote_reductive`]: real pathway, but
/// less common and opposite the usual oxidative direction (prone to undo prior
/// edits). Soft demotion on the `find_path` heap only — never dropped
/// (HEURISTICS: not decided).
pub fn hydrogenation() -> RuleSet {
    RuleSet::new(
        Some("Hydrogenation".into()),
        [
            demote_reductive(smirks_row(
                "alkyne",
                "[#6:1]#[#6:2]>>[*:1]=[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: Some("HH".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            )),
            demote_reductive(smirks_row(
                "alkene",
                "[#6:1]=,:[#6:2]>>[*:1]-[*:2]",
                SiteKind::AtomPair,
                vec![1, 2],
                Effect {
                    adds: Some("HH".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            )),
            demote_reductive(endpoint_row(
                "path_end",
                "[*:1]",
                vec![1],
                Effect {
                    adds: Some("H".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
                "keep",
            )),
        ],
    )
}

/// `NitrogenReduction` from Python `xenosite.forest.rules`.
pub fn nitrogen_reduction() -> RuleSet {
    let o_leave = Effect {
        cleaves: true,
        leave_count: Some(1),
        leave_formula: crate::pattern::leave_o(),
        partner: Some("O".into()),
        ..Default::default()
    };
    let oo_leave = Effect {
        cleaves: true,
        leave_count: Some(2),
        leave_formula: crate::pattern::leave_oo(),
        partner: Some("O".into()),
        ..Default::default()
    };
    RuleSet::new(
        Some("NitrogenReduction".into()),
        [
            smirks_row(
                "nitro_charged",
                "[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                o_leave.clone(),
            ),
            smirks_row(
                "nitro_anion",
                "[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                o_leave.clone(),
            ),
            smirks_row(
                "nitro_neutral",
                "[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                o_leave.clone(),
            ),
            smirks_row(
                "nitro_to_amine",
                "[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                oo_leave.clone(),
            ),
            smirks_row(
                "nitro_both",
                "[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                oo_leave.clone(),
            ),
            smirks_row(
                "hydroxylamine",
                "[#7:1]-,:[#8:2]>>([*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                o_leave.clone(),
            ),
            smirks_row(
                "nitroso",
                "[#7D2:1]=[#8:2]>>([*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                o_leave,
            ),
            smirks_row(
                "nitro_both_any",
                "[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                oo_leave,
            ),
        ],
    )
}

/// Soft-demote reductive / counter-oxidative patterns on the find_path heap.
/// Same rationale as Hydrogenation: real but less common toward typical Phase I
/// oxidative targets; prone to undo prior edits (HEURISTICS: not decided).
fn demote_reductive(pattern: PatternInfo) -> PatternInfo {
    pattern.with_search_bias(-1)
}

/// `OxygenReduction` from Python `xenosite.forest.rules`.
///
/// Patterns carry `search_bias = -1` (see [`demote_reductive`]).
pub fn oxygen_reduction() -> RuleSet {
    RuleSet::new(
        Some("OxygenReduction".into()),
        [
            demote_reductive(smirks_row(
                "carbonyl",
                "[#8:1]=[#6,#7:2]>>[*:1]-[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("HH".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            )),
            demote_reductive(smirks_row(
                "peroxide",
                "[#8:1]-[#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            )),
        ],
    )
}

/// `ReductiveDehalogenation` from Python `xenosite.forest.rules`.
pub fn reductive_dehalogenation() -> RuleSet {
    let base = Effect {
        cleaves: true,
        ..Default::default()
    };
    RuleSet::new(
        Some("ReductiveDehalogenation".into()),
        [
            smirks_row(
                "cleave",
                "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![2],
                base.clone(),
            )
            .with_possibilities(halide_remove_branches(1, base.clone())),
            smirks_row(
                "alkene",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]",
                SiteKind::Atom,
                vec![2],
                base.clone(),
            )
            .with_possibilities(halide_remove_branches(1, base)),
        ],
    )
}

/// `SulfurReduction` from Python `xenosite.forest.rules`.
pub fn sulfur_reduction() -> RuleSet {
    RuleSet::new(
        Some("SulfurReduction".into()),
        [
            smirks_row(
                "sulfoxide",
                "[#16:1]=[#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    cleaves: true,
                    leave_count: Some(1),
                    leave_formula: crate::pattern::leave_o(),
                    partner: Some("O".into()),
                    ..Default::default()
                },
            ),
            smirks_row(
                "disulfide",
                "[#16:1]-[#16:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    cleaves: true,
                    partner: Some("S".into()),
                    ..Default::default()
                },
            ),
            smirks_row(
                "thioether",
                "[#16:1]-[#6,#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    cleaves: true,
                    ..Default::default()
                },
            )
            .with_possibilities(vec![
                Effect {
                    when: Some(crate::pattern::When::atomic(2, 6)),
                    cleaves: true,
                    ..Default::default()
                }
                .sealed(),
                Effect {
                    when: Some(crate::pattern::When::atomic(2, 8)),
                    cleaves: true,
                    leave_count: Some(1),
                    leave_formula: crate::pattern::leave_o(),
                    ..Default::default()
                }
                .sealed(),
            ]),
        ],
    )
}

/// `Epoxidation` from Python `xenosite.forest.rules`.
///
/// `dearomatizes` is capability: aromatic C=C/C=N epoxidation clears the ring
/// bit at the site (shell residual needs |Δaromatic|). Resolved false on
/// aliphatic matches via [`PatternInfo::resolve_for_match`].
pub fn epoxidation() -> RuleSet {
    RuleSet::new(
        Some("Epoxidation".into()),
        [smirks_row(
            "epoxide",
            "[#6:1]=[#6,#7:2]>>[*:1]1-[*:2][O]1",
            SiteKind::Bond,
            vec![1, 2],
            Effect {
                adds: Some("O".into()),
                removes: None,
                cleaves: false,
                methide: false,
                dearomatizes: true,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )],
    )
}

/// Alkene → vicinal diol in one hop (epoxidation then hydrolytic opening).
///
/// Metabolize applies the diol SMIRKS; [`canonical_plan`](crate::ruleset::RuleSet::canonical_plan)
/// records `Epoxidation` then `EpoxideOpening` — the matching elementary leaves.
/// Catalog `dearomatizes` is capability (aromatic alkene clears the ring bit).
pub fn epoxide_hydration() -> RuleSet {
    RuleSet::new(
        Some("EpoxideHydration".into()),
        [smirks_row(
            "diol",
            "[#6:1]=[#6,#7:2]>>[*:1](O)[*:2]O",
            SiteKind::Bond,
            vec![1, 2],
            Effect {
                // Alkene → vicinal diol: +2O +2H (each OH carries H).
                adds: Some("OOHH".into()),
                removes: None,
                cleaves: false,
                methide: false,
                dearomatizes: true,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )],
    )
    .with_canonical_plan(crate::canonical_plan::epoxide_hydration_canonical_plan)
    .with_parity_exception(
        "Rust-only PhaseOne leaf (epoxide→diol one-hop); Python PhaseOne omits it",
    )
}

/// `SulfurOxidation` from Python `xenosite.forest.rules`.
pub fn sulfur_oxidation() -> RuleSet {
    RuleSet::new(
        Some("SulfurOxidation".into()),
        [
            smirks_row(
                "zwitterion",
                "[#16;v2,v4:1]>>[*&H0&+:1][O-]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "hydroxy",
                "[#16;v2,v4:1]>>[*:1][O]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    // Net H is substrate-dependent (thioether may gain H; thiol
                    // may not). Keep O-only; formula_check strips H.
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "oxo",
                "[#16;v2,v4:1]>>[*:1]=O",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `NitrogenOxidation` from Python `xenosite.forest.rules`.
pub fn nitrogen_oxidation() -> RuleSet {
    RuleSet::new(
        Some("NitrogenOxidation".into()),
        [
            smirks_row(
                "hydroxylamine",
                "[#7v3h1,#7v3h2:1]>>[*:1]O",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "nitroso",
                "[#7v3H2:1]>>[*:1]=O",
                SiteKind::Atom,
                vec![1],
                Effect {
                    // Primary amine → nitroso: +O and lose both N–H.
                    adds: Some("O".into()),
                    removes: Some("HH".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "n_oxide",
                "[#7v3H0:1]>>[*&H0&+:1][O-]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `OxidativeDehalogenation` from Python `xenosite.forest.rules`.
pub fn oxidative_dehalogenation() -> RuleSet {
    let o_cleave = Effect {
        adds: Some("O".into()),
        cleaves: true,
        ..Default::default()
    };
    let oo_cleave = Effect {
        adds: Some("OO".into()),
        cleaves: true,
        ..Default::default()
    };
    RuleSet::new(
        Some("OxidativeDehalogenation".into()),
        [
            smirks_row(
                "alcohol",
                "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O",
                SiteKind::Atom,
                vec![2],
                o_cleave.clone(),
            )
            .with_possibilities(halide_remove_branches(1, o_cleave.clone())),
            smirks_row(
                "carbonyl",
                "[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O",
                SiteKind::Atom,
                vec![2],
                o_cleave.clone(),
            )
            .with_possibilities(halide_remove_branches(1, o_cleave.clone())),
            smirks_row(
                "carboxylic",
                "[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O",
                SiteKind::Atom,
                vec![2],
                oo_cleave.clone(),
            )
            .with_possibilities(halide_remove_branches(1, oo_cleave.clone())),
            // Halogen migrates onto the adjacent carbon — net formula keeps X.
            smirks_row(
                "rearrange",
                "[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("O".into()),
                    ..Default::default()
                },
            ),
            smirks_row(
                "gem_carboxylic",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]",
                SiteKind::Atom,
                vec![2],
                oo_cleave.clone(),
            )
            .with_possibilities(halide_remove_branches(1, oo_cleave.clone())),
            smirks_row(
                "gem_hydrate",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]",
                SiteKind::Atom,
                vec![2],
                oo_cleave.clone(),
            )
            .with_possibilities(halide_remove_branches(1, oo_cleave)),
        ],
    )
}

/// `Acetylation` from Python `xenosite.forest.rules`.
pub fn acetylation() -> RuleSet {
    RuleSet::new(
        Some("Acetylation".into()),
        [smirks_row(
            "acetyl",
            "[#7h1,#7h2,#8h1,#16h1:1]>>[*:1][#6](=[#8])[#6]",
            SiteKind::Atom,
            vec![1],
            Effect {
                adds: Some("CCO".into()),
                removes: Some("H".into()),
                cleaves: false,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
                ..Default::default()
            },
        )],
    )
}

/// `Sulfation` from Python `xenosite.forest.rules`.
pub fn sulfation() -> RuleSet {
    RuleSet::new(
        Some("Sulfation".into()),
        [
            smirks_row(
                "alcohol",
                "[#6:1][#8H1:2]>>[*:1][*:2]S(=O)(=O)O",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("SOOO".into()),
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "epoxide_methyl_sulfone",
                "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>[*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1",
                SiteKind::Atom,
                vec![4],
                Effect {
                    adds: Some("CSO".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `Glucuronidation` from Python `xenosite.forest.rules`.
pub fn glucuronidation() -> RuleSet {
    RuleSet::new(
        Some("Glucuronidation".into()),
        [
            smirks_row(
                "alcohol",
                "[#8H1:1][#6:2]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2])1",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCOOOOOO".into()),
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "carboxylate",
                "[#8H1,#8-:1][#6:2](=[#8:3])[#6:4]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2](=[#8:3])[*:4])1",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCOOOOOO".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// `Glutathionation` from Python `xenosite.forest.rules`.
pub fn glutathionation() -> RuleSet {
    RuleSet::new(
        Some("Glutathionation".into()),
        [
            smirks_row(
                "epoxide_ch",
                "[#6H1:1]1[#8:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1][*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "epoxide_ch2",
                "[#6H2:1]1[#8:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1][*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "epoxide_c",
                "[#6H0:1]([!#1:4])1[#8:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1]([*:4])[*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "halide",
                "[#6:1][#9,#17,#35,#53:2]>>C(CC(=O)N[C@@H](CS([*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "thiol",
                "[#16h1:1]>>C(CC(=O)N[C@@H](CS([*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "alkene",
                "[#6H2:1]=[#6:2]>>C(CC(=O)N[C@@H](CS([*:1]-[*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "michael",
                "[#6H1:1]=[#6:2][#6:3]=[#8,#7:4]>>C(CC(=O)N[C@@H](CS([*:1][*:2]=[*:3][*:4]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "carbonyl",
                "[#6;H1,H2:1]=[#8:2]>>C(CC(=O)N[C@@H](CS([*:1]([*:2])))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "aziridine_ch",
                "[#6H1:1]1[#7:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1][*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "aziridine_ch2",
                "[#6H2:1]1[#7:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1][*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "aziridine_c",
                "[#6H0:1]([!#1:4])1[#7:2][#6:3]1>>C(CC(=O)N[C@@H](CS([*:1]([*:4])[*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "mesylate",
                "[#6:1][#8:2]S(=O)(=O)>>C(CC(=O)N[C@@H](CS([*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
            smirks_row(
                "isocyanate",
                "[#7:1]=[#6:2]=[#8,#16:3]>>C(CC(=O)N[C@@H](CS([*:2](=[*:3])[*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: Some("CCCCCCCCCCNNNOOOOOOS".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        ],
    )
}

/// Phase I catalog (Python `PhaseOne`).
pub fn phase_one() -> RuleSet {
    RuleSet::compose(
        Some("PhaseOne".into()),
        [
            hydroxylation(),
            epoxidation(),
            epoxide_hydration(),
            sulfur_oxidation(),
            nitrogen_oxidation(),
            dehydrogenation(),
            quinone_formation(),
            dephosphorylation(),
            epoxide_opening(),
            hydrolysis(),
            dehydration(),
            hydrogenation(),
            nitrogen_reduction(),
            oxygen_reduction(),
            reductive_dehalogenation(),
            sulfur_reduction(),
            dealkylation(),
            oxidative_dehalogenation(),
        ],
    )
}

/// Default search ruleset (Python `find_path.default_ruleset`).
pub fn default_ruleset() -> RuleSet {
    RuleSet::compose(
        Some("Default".into()),
        [
            dealkylation(),
            quinone_formation(),
            hydroxylation(),
            dehydrogenation(),
        ],
    )
}

/// Process-wide [`default_ruleset`] for doors that take `&RuleSet`.
///
/// Callers that need a different catalog pass their own `&RuleSet` to the
/// non-`_default` APIs (`find_path`, `bfs`, `product_graph`, …).
pub fn default_ruleset_ref() -> &'static RuleSet {
    use std::sync::OnceLock;
    static DEFAULT: OnceLock<RuleSet> = OnceLock::new();
    DEFAULT.get_or_init(default_ruleset)
}

/// Process-wide [`phase_one`] (full Phase-I nest).
pub fn phase_one_ref() -> &'static RuleSet {
    use std::sync::OnceLock;
    static PHASE_ONE: OnceLock<RuleSet> = OnceLock::new();
    PHASE_ONE.get_or_init(phase_one)
}

/// Every ported leaf rule as one nested catalog.
pub fn all_rules() -> RuleSet {
    RuleSet::compose(
        Some("All".into()),
        [
            hydroxylation(),
            epoxidation(),
            epoxide_hydration(),
            sulfur_oxidation(),
            nitrogen_oxidation(),
            dehydrogenation(),
            quinone_formation(),
            dephosphorylation(),
            epoxide_opening(),
            hydrolysis(),
            dehydration(),
            hydrogenation(),
            nitrogen_reduction(),
            oxygen_reduction(),
            reductive_dehalogenation(),
            sulfur_reduction(),
            dealkylation(),
            oxidative_dehalogenation(),
            n_dealkylation(),
            azo_splitting(),
            benzodioxole_reduction(),
            nitroaromatic_reduction(),
            thiophene_sulfur_oxidation(),
            acetylation(),
            sulfation(),
            glucuronidation(),
            glutathionation(),
        ],
    )
}

/// Leaf names in catalog order.
pub fn catalog_names() -> &'static [&'static str] {
    &[
        "Hydroxylation",
        "Dehydrogenation",
        "QuinoneFormation",
        "Dealkylation",
        "NDealkylation",
        "AzoSplitting",
        "BenzodioxoleReduction",
        "NitroaromaticReduction",
        "ThiopheneSulfurOxidation",
        "Dephosphorylation",
        "EpoxideOpening",
        "Hydrolysis",
        "Dehydration",
        "Hydrogenation",
        "NitrogenReduction",
        "OxygenReduction",
        "ReductiveDehalogenation",
        "SulfurReduction",
        "Epoxidation",
        "EpoxideHydration",
        "SulfurOxidation",
        "NitrogenOxidation",
        "OxidativeDehalogenation",
        "Acetylation",
        "Sulfation",
        "Glucuronidation",
        "Glutathionation",
    ]
}

/// Named leaf [`RuleSet`] for plan replay (elementary apply).
pub fn leaf_rule(name: &str) -> Option<RuleSet> {
    Some(match name {
        "Hydroxylation" => hydroxylation(),
        "Dehydrogenation" => dehydrogenation(),
        "QuinoneFormation" => quinone_formation(),
        "Dealkylation" => dealkylation(),
        "NDealkylation" => n_dealkylation(),
        "AzoSplitting" => azo_splitting(),
        "BenzodioxoleReduction" => benzodioxole_reduction(),
        "NitroaromaticReduction" => nitroaromatic_reduction(),
        "ThiopheneSulfurOxidation" => thiophene_sulfur_oxidation(),
        "Dephosphorylation" => dephosphorylation(),
        "EpoxideOpening" => epoxide_opening(),
        "Hydrolysis" => hydrolysis(),
        "Dehydration" => dehydration(),
        "Hydrogenation" => hydrogenation(),
        "NitrogenReduction" => nitrogen_reduction(),
        "OxygenReduction" => oxygen_reduction(),
        "ReductiveDehalogenation" => reductive_dehalogenation(),
        "SulfurReduction" => sulfur_reduction(),
        "Epoxidation" => epoxidation(),
        "EpoxideHydration" => epoxide_hydration(),
        "SulfurOxidation" => sulfur_oxidation(),
        "NitrogenOxidation" => nitrogen_oxidation(),
        "OxidativeDehalogenation" => oxidative_dehalogenation(),
        "Acetylation" => acetylation(),
        "Sulfation" => sulfation(),
        "Glucuronidation" => glucuronidation(),
        "Glutathionation" => glutathionation(),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use crate::ruleset::{accept_all_rules, accept_all_sites};

    #[test]
    fn phase_one_nests_eighteen_leaf_rules() {
        let set = phase_one();
        assert_eq!(set.members().len(), 18);
        assert_eq!(set.name.as_deref(), Some("PhaseOne"));
    }

    #[test]
    fn all_rules_registers_every_leaf() {
        assert_eq!(all_rules().members().len(), 27);
        assert_eq!(catalog_names().len(), 27);
    }

    #[test]
    fn epoxidation_dearomatizes_capability_resolves_on_aromatic_site() {
        let set = epoxidation();
        let info = &set.patterns()[0];
        assert!(
            info.effect.dearomatizes,
            "catalog capability must declare dearomatizes"
        );

        let benzene = parse_mol("c1ccccc1").unwrap();
        let arom = set
            .candidates(&benzene)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!arom.is_empty(), "{arom:?}");
        assert!(
            arom.iter().all(|c| c.effect().dearomatizes),
            "aromatic site resolves dearomatizes; got {:?}",
            arom.iter()
                .map(|c| (c.site(), c.effect().dearomatizes))
                .collect::<Vec<_>>()
        );

        let ethene = parse_mol("C=C").unwrap();
        let aliph = set
            .candidates(&ethene)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!aliph.is_empty(), "{aliph:?}");
        assert!(
            aliph.iter().all(|c| !c.effect().dearomatizes),
            "aliphatic site resolves false; got {:?}",
            aliph
                .iter()
                .map(|c| (c.site(), c.effect().dearomatizes))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn epoxide_hydration_pattern_info_and_plan() {
        let set = epoxide_hydration();
        assert!(set.has_plan_hook());
        let patterns = set.patterns();
        assert_eq!(patterns.len(), 1);
        let info = patterns[0];
        assert_eq!(info.name, "diol");
        assert_eq!(info.site_kind, SiteKind::Bond);
        assert_eq!(info.site_map, vec![1, 2]);
        assert_eq!(info.effect.adds.as_deref(), Some("OOHH"));
        assert_eq!(info.effect.delta_formula.get("O"), Some(&2));
        assert_eq!(info.effect.delta_formula.get("H"), Some(&2));
        assert!(info.effect.dearomatizes, "catalog capability");
        assert!(!info.effect.cleaves);
        match &info.edit {
            Edit::Smirks(s) => assert!(s.contains(">>"), "{s}"),
            other => panic!("expected Smirks edit, got {other:?}"),
        }

        let mol = parse_mol("C=C").unwrap();
        let emissions = set
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!emissions.is_empty(), "{emissions:?}");
        let want = canon_of("OCCO").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| e.products.iter().any(|p| canon_of(p).unwrap() == want)),
            "want OCCO, got {:?}",
            emissions.iter().map(|e| &e.products).collect::<Vec<_>>()
        );
        let emission = emissions
            .iter()
            .find(|e| e.products.iter().any(|p| canon_of(p).unwrap() == want))
            .unwrap();
        assert_eq!(emission.pattern_name, "diol");
        assert_eq!(emission.leaf_rule(), Some("EpoxideHydration"));
        let names: Vec<_> = emission.plan.iter().map(|s| s.rule.as_str()).collect();
        assert_eq!(names, ["Epoxidation", "EpoxideOpening"]);
        assert_eq!(emission.site_atoms.len(), 2);

        let benzene = parse_mol("c1ccccc1").unwrap();
        let arom = set
            .metabolize(&benzene, accept_all_rules, accept_all_sites, true)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert!(!arom.is_empty(), "{arom:?}");
        let names: Vec<_> = arom[0].plan.iter().map(|s| s.rule.as_str()).collect();
        assert_eq!(names, ["Epoxidation", "EpoxideOpening"]);
    }

    #[test]
    fn reductive_patterns_carry_negative_search_bias() {
        for set in [hydrogenation(), oxygen_reduction()] {
            let patterns = set.patterns();
            assert!(!patterns.is_empty(), "{:?}", set.name);
            assert!(
                patterns.iter().all(|p| p.search_bias == -1),
                "{:?} soft-demotes on the find_path heap; got {:?}",
                set.name,
                patterns
                    .iter()
                    .map(|p| (p.name.as_str(), p.search_bias))
                    .collect::<Vec<_>>()
            );
        }
        assert_eq!(hydroxylation().patterns()[0].search_bias, 0);
    }

    #[test]
    fn patterns_seal_delta_formula_from_adds_removes() {
        let hydroxy = hydroxylation();
        let h_patterns = hydroxy.patterns();
        let h = &h_patterns[0];
        assert_eq!(h.effect.delta_formula.get("O"), Some(&1));
        assert_eq!(h.effect.delta_formula.get("H"), Some(&-1));
        let dh_set = dehydrogenation();
        let dh_patterns = dh_set.patterns();
        let dh = dh_patterns.iter().find(|p| p.name == "sulfoxide").unwrap();
        assert_eq!(dh.effect.delta_formula.get("H"), Some(&-2));
    }

    #[test]
    fn halide_whens_carry_per_branch_delta_formula() {
        let red = reductive_dehalogenation();
        let red_patterns = red.patterns();
        let cleave = red_patterns.iter().find(|p| p.name == "cleave").unwrap();
        assert_eq!(cleave.possibilities.len(), 5);
        let symbols: std::collections::BTreeSet<_> = cleave
            .possibilities
            .iter()
            .map(|e| {
                let (el, n) = e.delta_formula.iter().next().unwrap();
                (el.clone(), *n)
            })
            .collect();
        assert!(symbols.contains(&("Cl".into(), -1)));
        assert!(symbols.contains(&("F".into(), -1)));
        let ox = oxidative_dehalogenation();
        let ox_patterns = ox.patterns();
        let alcohol = ox_patterns.iter().find(|p| p.name == "alcohol").unwrap();
        let cl = alcohol
            .possibilities
            .iter()
            .find(|e| e.when.as_ref().and_then(|w| w.z) == Some(17))
            .unwrap();
        assert_eq!(cl.delta_formula.get("O"), Some(&1));
        assert_eq!(cl.delta_formula.get("Cl"), Some(&-1));
    }

    #[test]
    fn cleavage_variants_include_leave_and_junction() {
        let dealk = dealkylation();
        let patterns = dealk.patterns();
        let methyl = patterns
            .iter()
            .find(|p| p.name == "methyl_carboxylic")
            .unwrap();
        assert_eq!(methyl.effect.leave_formula.get("C"), Some(&1));
        assert_eq!(methyl.effect.leave_formula.get("H"), Some(&3));
        assert_eq!(methyl.effect.delta_formula.get("C"), Some(&-1));
        assert_eq!(methyl.effect.delta_formula.get("H"), Some(&-3));
        assert_eq!(methyl.effect.delta_formula.get("O"), Some(&2));

        let dioxole = benzodioxole_reduction().patterns()[0].clone();
        assert_eq!(dioxole.effect.delta_formula.get("C"), Some(&-1));
        assert_eq!(dioxole.effect.delta_formula.get("H"), Some(&-2));

        let nitro = nitroaromatic_reduction().patterns()[0].clone();
        assert_eq!(nitro.effect.delta_formula.get("O"), Some(&-1));
    }

    #[test]
    fn each_leaf_has_at_least_one_pattern() {
        for set in [
            hydroxylation(),
            dehydrogenation(),
            quinone_formation(),
            dealkylation(),
            n_dealkylation(),
            azo_splitting(),
            benzodioxole_reduction(),
            nitroaromatic_reduction(),
            thiophene_sulfur_oxidation(),
            dephosphorylation(),
            epoxide_opening(),
            hydrolysis(),
            dehydration(),
            hydrogenation(),
            nitrogen_reduction(),
            oxygen_reduction(),
            reductive_dehalogenation(),
            sulfur_reduction(),
            epoxidation(),
            sulfur_oxidation(),
            nitrogen_oxidation(),
            oxidative_dehalogenation(),
            acetylation(),
            sulfation(),
            glucuronidation(),
            glutathionation(),
        ] {
            assert!(!set.patterns().is_empty(), "{:?}", set.name);
        }
    }

    #[test]
    fn dealkylation_anisole_emits_phenol() {
        let mol = parse_mol("COc1ccccc1").unwrap();
        let emissions = dealkylation()
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let phenol = canon_of("Oc1ccccc1").unwrap();
        assert!(
            emissions
                .iter()
                .any(|e| { e.products.iter().any(|p| canon_of(p).unwrap() == phenol) }),
            "{emissions:?}"
        );
    }

    #[test]
    fn hydroxylation_still_matches_built_in_door() {
        let mol = parse_mol("CC").unwrap();
        let emissions = hydroxylation()
            .metabolize(&mol, accept_all_rules, accept_all_sites, true)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(emissions.len(), 1);
        assert_eq!(
            canon_of(&emissions[0].products[0]).unwrap(),
            canon_of("CCO").unwrap()
        );
    }

    #[test]
    fn glutathionation_has_gsh_product_side() {
        assert!(
            glutathionation()
                .patterns()
                .iter()
                .any(|p| p.name == "epoxide_ch")
        );
        assert!(
            glutathionation()
                .patterns()
                .iter()
                .any(|p| p.smarts.contains("[#6H1:1]"))
        );
    }
}
