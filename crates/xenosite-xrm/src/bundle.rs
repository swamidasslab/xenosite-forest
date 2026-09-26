//! Annotation bundles: SKOS concepts name things; reactions combine them.
//!
//! Site-localized display names are generated from templates + [`SiteRef`],
//! not stored as ontology terms.

use crate::skos::Thesaurus;
use crate::term::{SiteRef, Term};
use serde::{Deserialize, Serialize};

/// Spine roots used to facet a flat term list into an annotation bundle.
pub mod spines {
    pub const METABOLISM_PHASE: &str = "xrm:1000000";
    pub const CHEMICAL_TRANSFORMATION: &str = "xrm:1100000";
    pub const PHASE1_FAMILY: &str = "xrm:1200000";
    pub const PHASE2_FAMILY: &str = "xrm:1300000";
    pub const MEDCHEM_LIABILITY: &str = "xrm:1400000";
    pub const REACTIVE_METABOLITE: &str = "xrm:1500000";
    pub const SITE_TYPE: &str = "xrm:1600000";
    pub const STRUCTURAL_DELTA: &str = "xrm:1700000";
    pub const PRODUCT_STATUS: &str = "xrm:1800000";
    pub const RULE_PROVENANCE: &str = "xrm:1900000";
    pub const EVIDENCE: &str = "xrm:2000000";
    pub const BIOLOGICAL_CONTEXT: &str = "xrm:2100000";
    pub const LEAVING_GROUP: &str = "xrm:2200000";
    pub const PHASE_I: &str = "xrm:0000001";
    pub const PHASE_II: &str = "xrm:0000002";
}

/// One site-scoped (or molecule-level) annotation bundle.
///
/// Combinatorial labels are not concepts — they live in [`site_label`](Self::site_label).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnotationBundle {
    /// Most specific chemical transformation / reaction-type term.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transformation: Option<String>,
    /// Phase I / Phase II (or metabolism-phase leaf).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase: Option<String>,
    /// Rainbow / phase I family class (e.g. stable oxygenation).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub phase1_family: Vec<String>,
    /// Phase II conjugation family leaves.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub phase2_family: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub site_environment: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub structural_delta: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub medchem_interpretation: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub product_class: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub product_status: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub evidence_type: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub biological_context: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub leaving_group: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub generated_by_rule: Vec<String>,
    /// Dynamic localized display name (not an ontology concept).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub site_label: Option<String>,
    #[serde(default, skip_serializing_if = "SiteRef::is_empty")]
    pub site: SiteRef,
    /// Operational `xmet:*` property tips present in this bundle.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub xmet_properties: Vec<String>,
}

impl AnnotationBundle {
    /// Build one bundle per distinct [`SiteRef`] (empty site = molecule-level).
    pub fn from_terms(thesaurus: &Thesaurus, terms: &[Term]) -> Vec<Self> {
        let mut by_site: std::collections::BTreeMap<SiteRef, Vec<&Term>> =
            std::collections::BTreeMap::new();
        for t in terms {
            by_site.entry(t.site.clone()).or_default().push(t);
        }
        by_site
            .into_iter()
            .map(|(site, group)| Self::from_site_terms(thesaurus, site, &group))
            .collect()
    }

    fn from_site_terms(thesaurus: &Thesaurus, site: SiteRef, terms: &[&Term]) -> Self {
        let mut b = Self {
            site: site.clone(),
            ..Default::default()
        };

        let under = |spine: &str| -> Vec<&Term> {
            terms
                .iter()
                .copied()
                .filter(|t| under_spine(thesaurus, t.id.as_str(), spine))
                .collect()
        };

        // Prefer leaves (is_leaf) then deeper terms for transformation / phase.
        b.transformation = pick_primary(&under(spines::CHEMICAL_TRANSFORMATION))
            .or_else(|| pick_primary(&under(spines::PHASE1_FAMILY)))
            .or_else(|| pick_primary(&under(spines::PHASE2_FAMILY)))
            .map(|t| t.id.as_str().to_string());

        b.phase = terms
            .iter()
            .copied()
            .find(|t| t.id.as_str() == spines::PHASE_I || t.id.as_str() == spines::PHASE_II)
            .or_else(|| pick_primary(&under(spines::METABOLISM_PHASE)))
            .map(|t| t.id.as_str().to_string());

        b.phase1_family = ids_leaves(&under(spines::PHASE1_FAMILY));
        b.phase2_family = ids_leaves(&under(spines::PHASE2_FAMILY));
        b.site_environment = ids_leaves(&under(spines::SITE_TYPE));
        b.structural_delta = ids_leaves(&under(spines::STRUCTURAL_DELTA));
        b.medchem_interpretation = ids_leaves(&under(spines::MEDCHEM_LIABILITY));
        b.product_class = ids_leaves(&under(spines::REACTIVE_METABOLITE));
        b.product_status = ids_leaves(&under(spines::PRODUCT_STATUS));
        b.evidence_type = ids_leaves(&under(spines::EVIDENCE));
        b.biological_context = ids_leaves(&under(spines::BIOLOGICAL_CONTEXT));
        b.leaving_group = ids_leaves(&under(spines::LEAVING_GROUP));
        b.generated_by_rule = ids_leaves(&under(spines::RULE_PROVENANCE));

        b.site_label = site_label(thesaurus, &b);
        b.xmet_properties = xmet_tips(&b);
        b
    }
}

fn under_spine(thesaurus: &Thesaurus, id: &str, spine: &str) -> bool {
    if id == spine {
        return true;
    }
    thesaurus
        .ancestors(id)
        .iter()
        .any(|a| a.as_str() == spine)
}

fn pick_primary<'a>(terms: &[&'a Term]) -> Option<&'a Term> {
    terms
        .iter()
        .copied()
        .filter(|t| t.specificity.is_leaf)
        .max_by_key(|t| t.specificity.depth)
        .or_else(|| {
            terms
                .iter()
                .copied()
                .max_by_key(|t| t.specificity.depth)
        })
}

fn ids_leaves(terms: &[&Term]) -> Vec<String> {
    let mut leaves: Vec<&Term> = terms
        .iter()
        .copied()
        .filter(|t| t.specificity.is_leaf)
        .collect();
    if leaves.is_empty() {
        leaves = terms.to_vec();
    }
    leaves.sort_by(|a, b| {
        b.specificity
            .depth
            .cmp(&a.specificity.depth)
            .then_with(|| a.id.as_str().cmp(b.id.as_str()))
    });
    let mut out = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for t in leaves {
        if seen.insert(t.id.as_str().to_string()) {
            out.push(t.id.as_str().to_string());
        }
    }
    out
}

/// Generate a localized display label from the bundle (not an ontology term).
fn site_label(thesaurus: &Thesaurus, b: &AnnotationBundle) -> Option<String> {
    let transform_id = b.transformation.as_deref()?;
    let concept = thesaurus.get(transform_id)?;
    let base = concept.pref_label.as_str();
    if b.site.is_empty() {
        return Some(base.to_string());
    }
    if !b.site.map_nums.is_empty() {
        let maps = b
            .site
            .map_nums
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(",");
        return Some(format!("{base} @{maps}"));
    }
    if !b.site.reactant_atoms.is_empty() {
        // Element-index style when maps absent: a{i} → site index phrase.
        let atoms = b
            .site
            .reactant_atoms
            .iter()
            .map(|a| format!("a{a}"))
            .collect::<Vec<_>>()
            .join(",");
        return Some(format!("{base} at {atoms}"));
    }
    Some(base.to_string())
}

fn xmet_tips(b: &AnnotationBundle) -> Vec<String> {
    let mut p = Vec::new();
    if b.phase.is_some() {
        p.push("xmet:hasPhase".into());
    }
    if b.transformation.is_some() {
        p.push("xmet:hasTransformation".into());
    }
    if !b.product_class.is_empty() {
        p.push("xmet:hasProductClass".into());
    }
    if !b.site_environment.is_empty() {
        p.push("xmet:hasSiteEnvironment".into());
    }
    if !b.structural_delta.is_empty() {
        p.push("xmet:hasStructuralDelta".into());
    }
    if !b.medchem_interpretation.is_empty() {
        p.push("xmet:hasMedChemInterpretation".into());
    }
    if !b.evidence_type.is_empty() {
        p.push("xmet:hasEvidenceType".into());
    }
    if !b.biological_context.is_empty() {
        p.push("xmet:hasBiologicalContext".into());
    }
    if !b.leaving_group.is_empty() {
        p.push("xmet:hasLeavingGroup".into());
    }
    if !b.generated_by_rule.is_empty() {
        p.push("xmet:generatedByRule".into());
    }
    if !b.site.is_empty() {
        p.push("xmet:localizesToSite".into());
    }
    p
}
