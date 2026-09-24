//! Phase I / conjugation reaction rules as [`PatternInfo`] data.
//!
//! Ported from `xenosite.forest.rules`. Catalogs (`phase_one`, `default_ruleset`)
//! are nested [`RuleSet`] namespaces. Pair-endpoint patterns use
//! [`SiteKind::AtomPair`] + [`Edit::PairEndpoint`]; atom/bond SMIRKS metabolize
//! through the generic door.

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
        effect,
        skip_same_rings: false,
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
        effect,
        skip_same_rings,
    }
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
                    partner: None,
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
                "[#6R:1][#8H,#7D1H2,#7D2H1,#6D1H3,#6D2H2,#6D3H1:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("H".into()),
                    cleaves: false,
                    methide: true,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
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
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                },
                "replace_halogen",
            ),
            endpoint_row(
                "iminium",
                "[#6H0R:1][#7D3:2]",
                vec![1],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
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
                    dearomatizes: true,
                    leave_count: None,
                    partner: None,
                },
                "dealkylate",
            ),
        ],
    )
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
                },
            ),
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
                },
            ),
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
                },
            ),
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
                },
            ),
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
                },
            ),
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
                },
            ),
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
            },
        )],
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
                adds: None,
                removes: None,
                cleaves: true,
                methide: false,
                dearomatizes: false,
                leave_count: None,
                partner: None,
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
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_neutral",
                "[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
                SiteKind::DirectedBond,
                vec![1, 2],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
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
                    adds: Some("".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "hydrate",
                "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)",
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
                },
            ),
        ],
    )
}

/// `Dehydration` from Python `xenosite.forest.rules`.
pub fn dehydration() -> RuleSet {
    RuleSet::new(
        Some("Dehydration".into()),
        [
            smirks_row(
                "alcohol",
                "[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("OH".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "beta_elimination",
                "[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("OH".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "carbonyl",
                "[#6,#7:1]=[#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
        ],
    )
}

/// `Hydrogenation` from Python `xenosite.forest.rules`.
pub fn hydrogenation() -> RuleSet {
    RuleSet::new(
        Some("Hydrogenation".into()),
        [
            smirks_row(
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
                },
            ),
            smirks_row(
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
                },
            ),
            endpoint_row(
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
                },
                "keep",
            ),
        ],
    )
}

/// `NitrogenReduction` from Python `xenosite.forest.rules`.
pub fn nitrogen_reduction() -> RuleSet {
    RuleSet::new(
        Some("NitrogenReduction".into()),
        [
            smirks_row(
                "nitro_charged",
                "[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_anion",
                "[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_neutral",
                "[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_to_amine",
                "[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("OO".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_both",
                "[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("OO".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "hydroxylamine",
                "[#7:1]-,:[#8:2]>>([*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitroso",
                "[#7D2:1]=[#8:2]>>([*:1].[*:2])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "nitro_both_any",
                "[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("OO".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
        ],
    )
}

/// `OxygenReduction` from Python `xenosite.forest.rules`.
pub fn oxygen_reduction() -> RuleSet {
    RuleSet::new(
        Some("OxygenReduction".into()),
        [
            smirks_row(
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
                },
            ),
            smirks_row(
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
                },
            ),
        ],
    )
}

/// `ReductiveDehalogenation` from Python `xenosite.forest.rules`.
pub fn reductive_dehalogenation() -> RuleSet {
    RuleSet::new(
        Some("ReductiveDehalogenation".into()),
        [
            smirks_row(
                "cleave",
                "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "alkene",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: None,
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
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
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "disulfide",
                "[#16:1]-[#16:2]>>[*:1].[*:2]",
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
                },
            ),
            smirks_row(
                "thioether",
                "[#16:1]-[#6,#8:2]>>[*:1].[*:2]",
                SiteKind::Atom,
                vec![1],
                Effect {
                    adds: None,
                    removes: Some("O".into()),
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
        ],
    )
}

/// `Epoxidation` from Python `xenosite.forest.rules`.
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
                dearomatizes: false,
                leave_count: None,
                partner: None,
            },
        )],
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
                },
            ),
            smirks_row(
                "hydroxy",
                "[#16;v2,v4:1]>>[*:1][O]",
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
                },
            ),
            smirks_row(
                "nitroso",
                "[#7v3H2:1]>>[*:1]=O",
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
                },
            ),
        ],
    )
}

/// `OxidativeDehalogenation` from Python `xenosite.forest.rules`.
pub fn oxidative_dehalogenation() -> RuleSet {
    RuleSet::new(
        Some("OxidativeDehalogenation".into()),
        [
            smirks_row(
                "alcohol",
                "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "carbonyl",
                "[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "carboxylic",
                "[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "rearrange",
                "[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("O".into()),
                    removes: None,
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "gem_carboxylic",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
            smirks_row(
                "gem_hydrate",
                "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]",
                SiteKind::Atom,
                vec![2],
                Effect {
                    adds: Some("OO".into()),
                    removes: None,
                    cleaves: true,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                },
            ),
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
                },
            ),
            smirks_row(
                "epoxide_methyl_sulfone",
                "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>[*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1",
                SiteKind::Atom,
                vec![4],
                Effect {
                    adds: Some("CSO".into()),
                    removes: Some("O".into()),
                    cleaves: false,
                    methide: false,
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
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

/// Every ported leaf rule as one nested catalog.
pub fn all_rules() -> RuleSet {
    RuleSet::compose(
        Some("All".into()),
        [
            hydroxylation(),
            epoxidation(),
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
        "SulfurOxidation",
        "NitrogenOxidation",
        "OxidativeDehalogenation",
        "Acetylation",
        "Sulfation",
        "Glucuronidation",
        "Glutathionation",
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol::{canon_of, parse_mol};
    use crate::ruleset::{accept_all_rules, accept_all_sites};

    #[test]
    fn phase_one_nests_seventeen_leaf_rules() {
        let set = phase_one();
        assert_eq!(set.members().len(), 17);
        assert_eq!(set.name.as_deref(), Some("PhaseOne"));
    }

    #[test]
    fn all_rules_registers_every_leaf() {
        assert_eq!(all_rules().members().len(), 26);
        assert_eq!(catalog_names().len(), 26);
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
