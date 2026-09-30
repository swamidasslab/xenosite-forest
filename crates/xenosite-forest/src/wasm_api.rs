//! Browser / Node WASM exports. Enabled with `--features wasm`.
//!
//! Parity with [`crate::python_api`]: same classes and factories as the
//! PyO3 module (including `find_path` timeout and `random_path`). JS does
//! not get Python metabolize filter callbacks — compose patterns in data,
//! then `metabolize`.

use std::collections::HashMap;
use std::time::Duration;

use js_sys::{Array, Reflect};
use serde::Serialize;
use serde_json::{Value, json};
use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::wasm_bindgen;

use crate::find_path::{FindPathConfig, HeapScoreMode, PathCounters, find_path_with};
use crate::forest::Formula;
use crate::forest_mol::ForestMol as Held;
use crate::pathway::PathwayOptions;
use crate::pattern::{Edit, Effect, PatternInfo};
use crate::random_path::{random_path as random_path_rs, random_path_with};
use crate::rules::{
    dealkylation as dealkylation_rs, default_ruleset as default_ruleset_rs,
    dehydrogenation as dehydrogenation_rs, epoxidation as epoxidation_rs,
    epoxide_opening as epoxide_opening_rs, hydrolysis as hydrolysis_rs,
    hydroxylation as hydroxylation_rs, n_dealkylation as n_dealkylation_rs,
    phase_one as phase_one_rs, quinone_formation as quinone_formation_rs,
};
use crate::ruleset::{RuleSet, accept_all_rules, accept_all_sites};
use crate::{canon_smiles, hydroxylate, parse_mol};

fn js_err(err: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&err.to_string())
}

fn reflect_err(err: JsValue) -> JsValue {
    err.as_string()
        .map(JsValue::from)
        .unwrap_or_else(|| JsValue::from_str("Reflect failed"))
}

fn to_js(value: &Value) -> Result<JsValue, JsValue> {
    value
        .serialize(&serde_wasm_bindgen::Serializer::json_compatible())
        .map_err(js_err)
}

fn parse_edit(edit: &str) -> Edit {
    if edit.eq_ignore_ascii_case("hydroxyl") {
        Edit::Hydroxyl
    } else {
        Edit::Smirks(edit.to_string())
    }
}

fn edit_label(edit: &Edit) -> String {
    match edit {
        Edit::Hydroxyl => "hydroxyl".into(),
        Edit::Smirks(smirks) => smirks.clone(),
        Edit::PairEndpoint(name) => format!("pair:{name}"),
    }
}

fn plan_atom_label(atom: &crate::PlanAtom) -> String {
    match atom {
        crate::PlanAtom::Index(i) => i.to_string(),
        other => format!("{other:?}"),
    }
}

fn counters_json(counters: &PathCounters) -> Value {
    json!({
        "nodes": counters.nodes,
        "mol_edits": counters.mol_edits,
        "expansions": counters.expansions,
        "billed": counters.billed(),
        "dropped_duplicate_plan": counters.dropped_duplicate_plan,
        "dropped_exact_plan": counters.dropped_exact_plan,
        "dropped_skeleton_twin": counters.dropped_skeleton_twin,
        "diversity_repush": counters.diversity_repush,
        "unstable_csmi_key": counters.unstable_csmi_key,
        "timed_out": counters.timed_out,
    })
}

fn opt_usize(options: &JsValue, snake: &str, camel: &str) -> Result<Option<usize>, JsValue> {
    if options.is_undefined() || options.is_null() {
        return Ok(None);
    }
    for key in [snake, camel] {
        let v = Reflect::get(options, &JsValue::from_str(key)).map_err(reflect_err)?;
        if v.is_undefined() || v.is_null() {
            continue;
        }
        let n = v
            .as_f64()
            .ok_or_else(|| js_err(format!("{key} must be a number")))?;
        return Ok(Some(n as usize));
    }
    Ok(None)
}

fn opt_bool(options: &JsValue, snake: &str, camel: &str) -> Result<Option<bool>, JsValue> {
    if options.is_undefined() || options.is_null() {
        return Ok(None);
    }
    for key in [snake, camel] {
        let v = Reflect::get(options, &JsValue::from_str(key)).map_err(reflect_err)?;
        if v.is_undefined() || v.is_null() {
            continue;
        }
        return Ok(Some(v.as_bool().unwrap_or(false)));
    }
    Ok(None)
}

fn opt_f64(options: &JsValue, snake: &str, camel: &str) -> Result<Option<f64>, JsValue> {
    if options.is_undefined() || options.is_null() {
        return Ok(None);
    }
    for key in [snake, camel] {
        let v = Reflect::get(options, &JsValue::from_str(key)).map_err(reflect_err)?;
        if v.is_undefined() || v.is_null() {
            continue;
        }
        let n = v
            .as_f64()
            .ok_or_else(|| js_err(format!("{key} must be a number")))?;
        return Ok(Some(n));
    }
    Ok(None)
}

fn opt_string(options: &JsValue, key: &str) -> Result<Option<String>, JsValue> {
    if options.is_undefined() || options.is_null() {
        return Ok(None);
    }
    let v = Reflect::get(options, &JsValue::from_str(key)).map_err(reflect_err)?;
    if v.is_undefined() || v.is_null() {
        return Ok(None);
    }
    Ok(v.as_string())
}

/// Thin helper kept for demos; not part of the Python public stub.
#[wasm_bindgen(js_name = forest_canon_smiles)]
pub fn forest_canon_smiles(smiles: &str) -> Result<String, JsValue> {
    let mol = parse_mol(smiles).map_err(js_err)?;
    Ok(canon_smiles(&mol))
}

/// Thin helper kept for demos; not part of the Python public stub.
#[wasm_bindgen(js_name = forest_hydroxylate)]
pub fn forest_hydroxylate(smiles: &str) -> Result<String, JsValue> {
    let mol = Held::parse(smiles).map_err(js_err)?;
    let products = hydroxylate(&mol).map_err(js_err)?;
    Ok(products.join("\n"))
}

/// JS wrap of [`Formula`].
#[wasm_bindgen(js_name = Formula)]
pub struct JsFormula {
    counts: HashMap<String, i32>,
    charge: i32,
}

#[wasm_bindgen(js_class = Formula)]
impl JsFormula {
    #[wasm_bindgen(getter)]
    pub fn charge(&self) -> i32 {
        self.charge
    }

    /// Element → count map (plain object).
    #[wasm_bindgen(getter)]
    pub fn counts(&self) -> Result<JsValue, JsValue> {
        to_js(&json!(self.counts))
    }
}

impl From<&Formula> for JsFormula {
    fn from(formula: &Formula) -> Self {
        Self {
            counts: formula.counts.clone().into_iter().collect(),
            charge: formula.charge,
        }
    }
}

/// JS class wrapping [`ForestMol`].
#[wasm_bindgen(js_name = ForestMol)]
pub struct JsForestMol {
    inner: Held,
}

#[wasm_bindgen(js_class = ForestMol)]
impl JsForestMol {
    #[wasm_bindgen(constructor)]
    pub fn new(smiles: &str) -> Result<JsForestMol, JsValue> {
        Ok(Self {
            inner: Held::parse(smiles).map_err(js_err)?,
        })
    }

    #[wasm_bindgen(getter)]
    pub fn csmi(&self) -> String {
        self.inner.csmi().to_string()
    }

    /// Fail-closed dedup key. May be `undefined` when Chematic has no stable key.
    #[wasm_bindgen(getter, js_name = stable_csmi_key)]
    pub fn stable_csmi_key(&self) -> Option<String> {
        self.inner.stable_csmi_key().map(|s| s.to_string())
    }

    #[wasm_bindgen(getter)]
    pub fn formula(&self) -> JsFormula {
        JsFormula::from(self.inner.formula().as_ref())
    }

    #[wasm_bindgen]
    pub fn clear_structure(&self) {
        self.inner.clear_structure();
    }

    #[wasm_bindgen]
    pub fn copy(&self) -> JsForestMol {
        Self {
            inner: self.inner.copy_mol(),
        }
    }

    #[wasm_bindgen(js_name = edit_copy)]
    pub fn edit_copy(&self) -> JsForestMol {
        Self {
            inner: self.inner.edit_copy(),
        }
    }

    #[wasm_bindgen(js_name = smarts_matches)]
    pub fn smarts_matches(&self, smarts: &str) -> Result<JsValue, JsValue> {
        let hits = self.inner.smarts_matches(smarts).map_err(js_err)?;
        let rows: Vec<HashMap<String, usize>> = hits
            .iter()
            .map(|mapped| mapped.iter().map(|(&k, &v)| (k.to_string(), v)).collect())
            .collect();
        to_js(&json!(rows))
    }
}

/// JS wrap of [`PatternInfo`].
#[wasm_bindgen(js_name = PatternInfo)]
#[derive(Clone)]
pub struct JsPatternInfo {
    inner: PatternInfo,
}

#[wasm_bindgen(js_class = PatternInfo)]
impl JsPatternInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: String,
        smarts: String,
        edit: String,
        adds: Option<String>,
        removes: Option<String>,
        cleaves: Option<bool>,
        methide: Option<bool>,
    ) -> JsPatternInfo {
        Self {
            inner: PatternInfo::new(
                name,
                smarts,
                parse_edit(&edit),
                Effect {
                    adds,
                    removes,
                    cleaves: cleaves.unwrap_or(false),
                    methide: methide.unwrap_or(false),
                    dearomatizes: false,
                    leave_count: None,
                    partner: None,
                    ..Default::default()
                },
            ),
        }
    }

    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn smarts(&self) -> String {
        self.inner.smarts.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn edit(&self) -> String {
        edit_label(&self.inner.edit)
    }

    #[wasm_bindgen(getter)]
    pub fn adds(&self) -> Option<String> {
        self.inner.effect.adds.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn removes(&self) -> Option<String> {
        self.inner.effect.removes.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn cleaves(&self) -> bool {
        self.inner.effect.cleaves
    }

    #[wasm_bindgen(getter)]
    pub fn methide(&self) -> bool {
        self.inner.effect.methide
    }
}

/// JS wrap of [`crate::BoundPattern`].
#[wasm_bindgen(js_name = BoundPattern)]
#[derive(Clone)]
pub struct JsBoundPattern {
    inner: crate::BoundPattern,
}

#[wasm_bindgen(js_class = BoundPattern)]
impl JsBoundPattern {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.inner.name().to_string()
    }

    #[wasm_bindgen(getter, js_name = ruleName)]
    pub fn rule_name(&self) -> Option<String> {
        self.inner.rule_name().map(str::to_string)
    }

    #[wasm_bindgen(getter)]
    pub fn curie(&self) -> String {
        self.inner.curie()
    }

    #[wasm_bindgen(getter)]
    pub fn iri(&self) -> String {
        self.inner.iri()
    }

    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    #[wasm_bindgen(js_name = asRuleSet)]
    pub fn as_ruleset(&self) -> JsRuleSet {
        JsRuleSet {
            inner: self.inner.as_ruleset(),
        }
    }

    #[wasm_bindgen]
    pub fn metabolize(&self, mol: &JsForestMol) -> Result<JsValue, JsValue> {
        let forest = mol.inner.copy_mol();
        let emissions = self
            .inner
            .metabolize_default(&forest, true)
            .collect::<Result<Vec<_>, _>>()
            .map_err(js_err)?;
        let rows: Vec<Value> = emissions
            .into_iter()
            .map(|e| {
                json!({
                    "pattern_name": e.pattern_name,
                    "site": e.site,
                    "products": e.product_csmis(),
                    "rule_path": e.rule_path,
                })
            })
            .collect();
        to_js(&json!(rows))
    }
}

/// JS wrap of [`RuleSet`].
#[wasm_bindgen(js_name = RuleSet)]
#[derive(Clone)]
pub struct JsRuleSet {
    inner: RuleSet,
}

#[wasm_bindgen(js_class = RuleSet)]
impl JsRuleSet {
    #[wasm_bindgen(constructor)]
    pub fn new() -> JsRuleSet {
        Self {
            inner: RuleSet::new(None, []),
        }
    }

    #[wasm_bindgen(js_name = from_patterns)]
    pub fn from_patterns(patterns: Vec<JsPatternInfo>, name: Option<String>) -> JsRuleSet {
        Self {
            inner: RuleSet::new(name, patterns.into_iter().map(|p| p.inner)),
        }
    }

    #[wasm_bindgen]
    pub fn hydroxylation() -> JsRuleSet {
        Self {
            inner: hydroxylation_rs(),
        }
    }

    #[wasm_bindgen(js_name = o_dealkylation)]
    pub fn o_dealkylation() -> JsRuleSet {
        Self {
            inner: crate::ruleset::o_dealkylation(),
        }
    }

    #[wasm_bindgen]
    pub fn compose(sets: Vec<JsRuleSet>, name: Option<String>) -> JsRuleSet {
        Self {
            inner: RuleSet::compose(name, sets.into_iter().map(|set| set.inner.clone())),
        }
    }

    #[wasm_bindgen(getter)]
    pub fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.inner.members().len()
    }

    #[wasm_bindgen]
    pub fn patterns(&self) -> Vec<JsPatternInfo> {
        self.inner
            .patterns()
            .into_iter()
            .cloned()
            .map(|inner| JsPatternInfo { inner })
            .collect()
    }

    /// Convenience: push a hydroxyl pattern (demo helper).
    #[wasm_bindgen(js_name = add_hydroxyl)]
    pub fn add_hydroxyl(&mut self, name: &str, smarts: &str) {
        self.inner.push(PatternInfo::hydroxyl(name, smarts));
    }

    /// Run owned members. Returns rows
    /// `{pattern_name, site, products, rule_path}` (no Python filter callbacks).
    pub fn metabolize(&self, mol: &JsForestMol) -> Result<JsValue, JsValue> {
        let forest = mol.inner.copy_mol();
        let emissions = self
            .inner
            .metabolize(&forest, accept_all_rules, accept_all_sites, true)
            .collect::<Result<Vec<_>, _>>()
            .map_err(js_err)?;
        // Explicit CSMI downgrade at the JS string-row boundary.
        let rows: Vec<Value> = emissions
            .into_iter()
            .map(|e| {
                json!({
                    "pattern_name": e.pattern_name,
                    "site": e.site,
                    "products": e.product_csmis(),
                    "rule_path": e.rule_path,
                })
            })
            .collect();
        to_js(&json!(rows))
    }

    /// Catalog → child `RuleSet`; leaf → `BoundPattern`. Name or index (as string int).
    #[wasm_bindgen(js_name = get)]
    pub fn get(&self, key: &str) -> Result<JsValue, JsValue> {
        if let Ok(mut index) = key.parse::<isize>() {
            let len = self.inner.len() as isize;
            if index < 0 {
                index += len;
            }
            if index < 0 || index >= len {
                return Err(js_err("RuleSet index out of range"));
            }
            let index = index as usize;
            if self.inner.is_catalog() {
                let child = self
                    .inner
                    .get(index)
                    .ok_or_else(|| js_err("RuleSet index out of range"))?;
                return Ok(JsValue::from(JsRuleSet { inner: child }));
            }
            let bp = self
                .inner
                .bound_pattern_at(index)
                .ok_or_else(|| js_err("RuleSet index out of range"))?;
            return Ok(JsValue::from(JsBoundPattern { inner: bp }));
        }
        if self.inner.is_catalog() {
            if let Some(child) = self.inner.get_str(key) {
                return Ok(JsValue::from(JsRuleSet { inner: child }));
            }
        }
        if let Some(bp) = self.inner.bound_pattern(key) {
            return Ok(JsValue::from(JsBoundPattern { inner: bp }));
        }
        Err(js_err(format!("unknown RuleSet member: {key}")))
    }

    #[wasm_bindgen(js_name = contains)]
    pub fn contains(&self, name: &str) -> bool {
        self.inner.contains_name(name)
    }
}

impl Default for JsRuleSet {
    fn default() -> Self {
        Self::new()
    }
}

fn wrap_ruleset(inner: RuleSet) -> JsRuleSet {
    JsRuleSet { inner }
}

#[wasm_bindgen]
pub fn phase_one() -> JsRuleSet {
    wrap_ruleset(phase_one_rs())
}

#[wasm_bindgen]
pub fn epoxidation() -> JsRuleSet {
    wrap_ruleset(epoxidation_rs())
}

#[wasm_bindgen]
pub fn quinone_formation() -> JsRuleSet {
    wrap_ruleset(quinone_formation_rs())
}

#[wasm_bindgen]
pub fn epoxide_opening() -> JsRuleSet {
    wrap_ruleset(epoxide_opening_rs())
}

#[wasm_bindgen]
pub fn n_dealkylation() -> JsRuleSet {
    wrap_ruleset(n_dealkylation_rs())
}

#[wasm_bindgen]
pub fn hydroxylation() -> JsRuleSet {
    wrap_ruleset(hydroxylation_rs())
}

#[wasm_bindgen]
pub fn dehydrogenation() -> JsRuleSet {
    wrap_ruleset(dehydrogenation_rs())
}

#[wasm_bindgen]
pub fn dealkylation() -> JsRuleSet {
    wrap_ruleset(dealkylation_rs())
}

#[wasm_bindgen]
pub fn hydrolysis() -> JsRuleSet {
    wrap_ruleset(hydrolysis_rs())
}

#[wasm_bindgen]
pub fn default_ruleset() -> JsRuleSet {
    wrap_ruleset(default_ruleset_rs())
}

/// Decompressed Forest↔XMET SSSOM TSV text.
#[wasm_bindgen(js_name = forest_xmet_sssom)]
pub fn forest_xmet_sssom_js() -> String {
    crate::mapping::forest_xmet_sssom().to_string()
}

/// Resolve an `xf:` CURIE / Forest IRI to a `RuleSet` or `BoundPattern`.
#[wasm_bindgen(js_name = resolve)]
pub fn resolve_js(id: &str) -> Result<JsValue, JsValue> {
    match crate::mapping::resolve(id).map_err(js_err)? {
        crate::mapping::Resolved::Rule(inner) => Ok(JsValue::from(JsRuleSet { inner })),
        crate::mapping::Resolved::Pattern(inner) => Ok(JsValue::from(JsBoundPattern { inner })),
    }
}

#[wasm_bindgen(js_name = expand_iri)]
pub fn expand_iri_js(curie_or_iri: &str) -> String {
    crate::mapping::expand_iri(curie_or_iri)
}

#[wasm_bindgen(js_name = to_curie)]
pub fn to_curie_js(iri: &str) -> String {
    crate::mapping::to_curie(iri)
}

fn run_find_path(
    reactant: &str,
    target: &str,
    max_paths: Option<usize>,
    max_nodes: Option<usize>,
    use_atom_diff: Option<bool>,
    lazy_closer: Option<bool>,
    diversity: Option<bool>,
    drop_skeleton_twins: Option<bool>,
    score: Option<String>,
    timeout: Option<f64>,
    normalize_tautomer: Option<bool>,
    invert_target_tautomer: Option<bool>,
) -> Result<JsValue, JsValue> {
    let score_label = score.as_deref().unwrap_or("log-neg-pc");
    let heap_score = HeapScoreMode::from_label(score_label).ok_or_else(|| {
        js_err(format!(
            "unknown score {score_label:?}; try log-neg-pc, soft, add-both, …"
        ))
    })?;
    let timeout = match timeout {
        None => None,
        Some(secs) if secs.is_finite() && secs >= 0.0 => Some(Duration::from_secs_f64(secs)),
        Some(secs) => {
            return Err(js_err(format!(
                "timeout must be a non-negative finite number of seconds; got {secs}"
            )));
        }
    };
    let config = FindPathConfig {
        max_paths: max_paths.unwrap_or(1),
        max_nodes: max_nodes.unwrap_or(800),
        use_atom_diff: use_atom_diff.unwrap_or(true),
        lazy_closer: lazy_closer.unwrap_or(false),
        heap_score,
        drop_skeleton_twins: drop_skeleton_twins.unwrap_or(true),
        diversity: diversity.unwrap_or(false),
        timeout,
        normalize_tautomer: normalize_tautomer.unwrap_or(false),
        invert_target_tautomer: invert_target_tautomer.unwrap_or(false),
        ..FindPathConfig::default()
    };
    let rules = default_ruleset_rs();
    let mut counters = PathCounters::default();
    let hits = find_path_with(reactant, target, &rules, &mut counters, config, |_| true)
        .map_err(js_err)?
        .collect_all()
        .map_err(js_err)?;

    let hit_rows: Vec<Value> = hits
        .into_iter()
        .map(|hit| {
            let steps: Vec<Value> = hit
                .plan
                .iter()
                .map(|step| {
                    json!({
                        "rule": step.rule,
                        "site": step.site.iter().map(plan_atom_label).collect::<Vec<_>>(),
                    })
                })
                .collect();
            json!({
                "smiles": hit.smiles,
                "steps": steps,
            })
        })
        .collect();

    let out = Array::new();
    out.push(&to_js(&json!(hit_rows))?);
    out.push(&to_js(&counters_json(&counters))?);
    Ok(out.into())
}

/// Chematic ``find_path`` (PhaseOne). Returns ``[hits, counters]``.
///
/// Prefer [`find_path_with_options`] from TypeScript.
#[wasm_bindgen(js_name = find_path)]
#[allow(clippy::too_many_arguments)]
pub fn find_path(
    reactant: &str,
    target: &str,
    max_paths: Option<usize>,
    max_nodes: Option<usize>,
    use_atom_diff: Option<bool>,
    lazy_closer: Option<bool>,
    diversity: Option<bool>,
    drop_skeleton_twins: Option<bool>,
    score: Option<String>,
    timeout: Option<f64>,
) -> Result<JsValue, JsValue> {
    run_find_path(
        reactant,
        target,
        max_paths,
        max_nodes,
        use_atom_diff,
        lazy_closer,
        diversity,
        drop_skeleton_twins,
        score,
        timeout,
        None,
        None,
    )
}

/// JS-friendly `find_path(reactant, target, { maxPaths, timeout, … })`.
#[wasm_bindgen(js_name = find_path_with_options)]
pub fn find_path_with_options(
    reactant: &str,
    target: &str,
    options: JsValue,
) -> Result<JsValue, JsValue> {
    run_find_path(
        reactant,
        target,
        opt_usize(&options, "max_paths", "maxPaths")?,
        opt_usize(&options, "max_nodes", "maxNodes")?,
        opt_bool(&options, "use_atom_diff", "useAtomDiff")?,
        opt_bool(&options, "lazy_closer", "lazyCloser")?,
        opt_bool(&options, "diversity", "diversity")?,
        opt_bool(&options, "drop_skeleton_twins", "dropSkeletonTwins")?,
        opt_string(&options, "score")?,
        opt_f64(&options, "timeout", "timeout")?,
        opt_bool(&options, "normalize_tautomer", "normalizeTautomer")?,
        opt_bool(&options, "invert_target_tautomer", "invertTargetTautomer")?,
    )
}

/// Chematic tautomer pick → tagged ForestMol. Returns `{csmi, changed}`.
#[wasm_bindgen(js_name = normalize_tautomer)]
pub fn normalize_tautomer_js(smiles: &str) -> Result<JsValue, JsValue> {
    let out = crate::normalize_tautomer(smiles).map_err(|e| js_err(e.to_string()))?;
    to_js(&json!({
        "csmi": out.mol.csmi().as_ref(),
        "changed": out.changed,
    }))
}

/// Seeded random walk. Returns `{smiles, path, steps, patterns}`.
#[wasm_bindgen(js_name = random_path)]
pub fn random_path(
    reactant: &str,
    seed: u64,
    max_steps: Option<usize>,
    ruleset: Option<JsRuleSet>,
    skip_multicomponent: Option<bool>,
    skip_seen: Option<bool>,
) -> Result<JsValue, JsValue> {
    let owned;
    let rules = match &ruleset {
        Some(rs) => &rs.inner,
        None => {
            owned = phase_one_rs();
            &owned
        }
    };
    let options = PathwayOptions {
        skip_multicomponent: skip_multicomponent.unwrap_or(false),
        skip_seen: skip_seen.unwrap_or(false),
    };
    let max_steps = max_steps.unwrap_or(1);
    let outcome = if options == PathwayOptions::default() {
        random_path_rs(reactant, seed, rules, max_steps)
    } else {
        random_path_with(reactant, seed, rules, max_steps, options)
    }
    .map_err(js_err)?;

    let steps: Vec<Value> = outcome
        .steps
        .iter()
        .map(|step| {
            json!({
                "rule": step.rule,
                "pattern": step.pattern_name,
                "site": step.site,
                "products": step.products,
                "chosen": step.chosen,
            })
        })
        .collect();
    let patterns: Vec<Value> = outcome
        .patterns
        .iter()
        .map(|p| {
            json!({
                "name": p.name,
                "smarts": p.smarts,
                "cleaves": p.effect.cleaves,
                "search_bias": p.search_bias,
            })
        })
        .collect();
    to_js(&json!({
        "smiles": outcome.smiles,
        "path": outcome.path,
        "steps": steps,
        "patterns": patterns,
    }))
}

/// JS-friendly `random_path(reactant, seed, { maxSteps, skipSeen, … })`.
#[wasm_bindgen(js_name = random_path_with_options)]
pub fn random_path_with_options(
    reactant: &str,
    seed: u64,
    options: JsValue,
) -> Result<JsValue, JsValue> {
    let max_steps = opt_usize(&options, "max_steps", "maxSteps")?;
    let skip_multicomponent = opt_bool(&options, "skip_multicomponent", "skipMulticomponent")?;
    let skip_seen = opt_bool(&options, "skip_seen", "skipSeen")?;
    // Optional ruleset: if present as a wasm class instance, leave None here —
    // callers that need a custom ruleset use `random_path` with an explicit set.
    random_path(
        reactant,
        seed,
        max_steps,
        None,
        skip_multicomponent,
        skip_seen,
    )
}
