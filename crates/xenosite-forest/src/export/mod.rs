//! Serde export views for explicit `.to_dict()` / JSON — not the primary API.

use serde::Serialize;

use crate::canonical_plan::{PlanAtom, Step};
use crate::find_path::{PathCounters, PathOutcome};
use crate::pattern::Emission;
use crate::pattern::PatternInfo;
use crate::random_path::{RandomPathOutcome, RandomPathStep};

#[derive(Serialize)]
pub struct PlanStepView<'a> {
    pub rule: &'a str,
    pub site: Vec<String>,
}

#[derive(Serialize)]
pub struct PathOutcomeView<'a> {
    pub smiles: &'a str,
    pub steps: Vec<PlanStepView<'a>>,
}

#[derive(Serialize)]
pub struct PathCountersView {
    pub nodes: usize,
    pub mol_edits: usize,
    pub expansions: usize,
    pub billed: usize,
    pub dropped_duplicate_plan: usize,
    pub dropped_exact_plan: usize,
    pub dropped_skeleton_twin: usize,
    pub signal_contained_plan: usize,
    pub diversity_repush: usize,
    pub unstable_csmi_key: usize,
    pub timed_out: bool,
}

fn plan_atom_label(atom: &PlanAtom) -> String {
    match atom {
        PlanAtom::Index(i) => i.to_string(),
        PlanAtom::WillAdd { element, at } => format!("WillAdd({element}@{at})"),
        PlanAtom::AddedBy { rule, anchors } => format!("AddedBy({rule:?},{anchors:?})"),
    }
}

impl<'a> From<&'a Step> for PlanStepView<'a> {
    fn from(step: &'a Step) -> Self {
        Self {
            rule: step.rule.as_str(),
            site: step.site.iter().map(plan_atom_label).collect(),
        }
    }
}

impl<'a> From<&'a PathOutcome> for PathOutcomeView<'a> {
    fn from(hit: &'a PathOutcome) -> Self {
        Self {
            smiles: hit.smiles.as_str(),
            steps: hit.plan.steps().iter().map(PlanStepView::from).collect(),
        }
    }
}

#[derive(Serialize)]
pub struct EmissionView<'a> {
    pub pattern_name: &'a str,
    pub site: usize,
    pub products: Vec<String>,
    pub rule_path: Vec<Option<String>>,
}

#[derive(Serialize)]
pub struct RandomPathStepView<'a> {
    pub rule: &'a str,
    pub pattern: &'a str,
    pub site: Vec<usize>,
    pub products: Vec<String>,
    pub chosen: usize,
}

#[derive(Serialize)]
pub struct PatternInfoView<'a> {
    pub name: &'a str,
    pub smarts: &'a str,
    pub cleaves: bool,
    pub search_bias: i8,
}

#[derive(Serialize)]
pub struct RandomPathOutcomeView<'a> {
    pub smiles: &'a str,
    pub path: Vec<String>,
    pub steps: Vec<RandomPathStepView<'a>>,
    pub patterns: Vec<PatternInfoView<'a>>,
}

impl<'a> From<&'a Emission> for EmissionView<'a> {
    fn from(e: &'a Emission) -> Self {
        Self {
            pattern_name: e.pattern_name.as_str(),
            site: e.site,
            products: e.product_csmis(),
            rule_path: e.rule_path.clone(),
        }
    }
}

impl<'a> From<&'a RandomPathStep> for RandomPathStepView<'a> {
    fn from(s: &'a RandomPathStep) -> Self {
        Self {
            rule: s.rule.as_str(),
            pattern: s.pattern_name.as_str(),
            site: s.site.clone(),
            products: s.products.clone(),
            chosen: s.chosen,
        }
    }
}

impl<'a> From<&'a PatternInfo> for PatternInfoView<'a> {
    fn from(p: &'a PatternInfo) -> Self {
        Self {
            name: p.name.as_str(),
            smarts: p.smarts.as_str(),
            cleaves: p.effect.cleaves,
            search_bias: p.search_bias,
        }
    }
}

impl<'a> From<&'a RandomPathOutcome> for RandomPathOutcomeView<'a> {
    fn from(o: &'a RandomPathOutcome) -> Self {
        Self {
            smiles: o.smiles.as_str(),
            path: o.path.clone(),
            steps: o.steps.iter().map(RandomPathStepView::from).collect(),
            patterns: o.patterns.iter().map(PatternInfoView::from).collect(),
        }
    }
}

impl From<&PathCounters> for PathCountersView {
    fn from(c: &PathCounters) -> Self {
        Self {
            nodes: c.nodes,
            mol_edits: c.mol_edits,
            expansions: c.expansions,
            billed: c.billed(),
            dropped_duplicate_plan: c.dropped_duplicate_plan,
            dropped_exact_plan: c.dropped_exact_plan,
            dropped_skeleton_twin: c.dropped_skeleton_twin,
            signal_contained_plan: c.signal_contained_plan,
            diversity_repush: c.diversity_repush,
            unstable_csmi_key: c.unstable_csmi_key,
            timed_out: c.timed_out,
        }
    }
}
