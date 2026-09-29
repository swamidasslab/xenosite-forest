//! Forest ↔ XMET SSSOM coverage gates (resolve + inventory).

use std::collections::BTreeSet;

use xenosite_forest::{
    Resolved, catalog_names, forest_xmet_sssom, leaf_rule, parse_forest_xmet_sssom, resolve,
};

fn object_first_segment(object_id: &str) -> &str {
    let local = object_id.strip_prefix("xf:").unwrap_or(object_id);
    local.split('/').next().unwrap_or(local)
}

#[test]
fn embed_round_trip() {
    let text = forest_xmet_sssom();
    assert!(text.contains("subject_id\t"));
    assert!(text.contains("xf:Tautomerization/tautomer_h"));
    assert!(text.contains("xf:Tautomerization/path_partner"));
    assert!(!text.contains("forest.rule:"));
    assert!(!text.contains("forest.pattern:"));
    // Short-code ruleset objects must not ship.
    for row in parse_forest_xmet_sssom() {
        let local = row
            .object_id
            .strip_prefix("xf:")
            .unwrap_or(row.object_id.as_str());
        assert!(
            !matches!(
                local,
                "CJ" | "SO" | "UO" | "DH" | "HD" | "RD" | "QF" | "TT" | "BA"
            ),
            "short-code object still present: {}",
            row.object_id
        );
    }
}

#[test]
fn all_sssom_objects_resolve() {
    let rows = parse_forest_xmet_sssom();
    assert!(!rows.is_empty());
    let mut failures = Vec::new();
    for row in &rows {
        let id = &row.object_id;
        match resolve(id) {
            Ok(Resolved::Rule(r)) => {
                if id.contains('/') {
                    failures.push(format!("{id}: expected BoundPattern, got Rule {:?}", r.name));
                }
            }
            Ok(Resolved::Pattern(bp)) => {
                if !id.contains('/') {
                    failures.push(format!("{id}: expected Rule, got BoundPattern {}", bp.name()));
                } else {
                    let expect = id.strip_prefix("xf:").unwrap_or(id);
                    assert_eq!(bp.curie(), format!("xf:{expect}"));
                }
            }
            Err(e) => failures.push(format!("{id}: {e}")),
        }
    }
    assert!(
        failures.is_empty(),
        "SSSOM objects failed to resolve ({}):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn patterns_nest_under_rule() {
    let rows = parse_forest_xmet_sssom();
    let objects: BTreeSet<_> = rows.iter().map(|r| r.object_id.clone()).collect();
    let mut failures = Vec::new();
    for id in &objects {
        let Some(rest) = id.strip_prefix("xf:") else {
            continue;
        };
        let Some((rule, pat)) = rest.split_once('/') else {
            continue;
        };
        let parent = format!("xf:{rule}");
        if !objects.contains(&parent) {
            // Parent may still resolve as a registered leaf without its own row.
            if resolve(&parent).is_err() {
                failures.push(format!("{id}: missing parent {parent}"));
                continue;
            }
        }
        match resolve(&parent) {
            Ok(Resolved::Rule(r)) => {
                if r.bound_pattern(pat).is_none() {
                    failures.push(format!("{id}: parent resolves but [{pat}] missing"));
                }
            }
            other => failures.push(format!("{id}: parent {parent} -> {other:?}")),
        }
        // path ≡ index
        let via_path = resolve(id);
        let via_index = resolve(&parent).ok().and_then(|r| match r {
            Resolved::Rule(rule) => rule.bound_pattern(pat).map(Resolved::Pattern),
            Resolved::Pattern(_) => None,
        });
        match (via_path, via_index) {
            (Ok(Resolved::Pattern(a)), Some(Resolved::Pattern(b))) => {
                assert_eq!(a.name(), b.name());
                assert_eq!(a.rule_name(), b.rule_name());
            }
            (a, b) => failures.push(format!("{id}: path≡index mismatch {a:?} vs {b:?}")),
        }
    }
    assert!(
        failures.is_empty(),
        "pattern nesting failures ({}):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn root_tops_registered() {
    let rows = parse_forest_xmet_sssom();
    let mut firsts = BTreeSet::new();
    for row in &rows {
        firsts.insert(object_first_segment(&row.object_id).to_string());
    }
    let registered: BTreeSet<_> = xenosite_forest::ROOT_CATALOGS
        .iter()
        .chain(xenosite_forest::LEAF_CTORS.iter())
        .map(|(n, _)| (*n).to_string())
        .collect();
    let missing: Vec<_> = firsts
        .into_iter()
        .filter(|first| !registered.contains(first))
        .collect();
    assert!(
        missing.is_empty(),
        "SSSOM first segments not in root registry: {missing:?}"
    );
}

#[test]
fn inventory_fully_covered_in_sssom() {
    let rows = parse_forest_xmet_sssom();
    let objects: BTreeSet<_> = rows.iter().map(|r| r.object_id.clone()).collect();
    let mut missing = Vec::new();
    for name in catalog_names() {
        let id = format!("xf:{name}");
        if !objects.contains(&id) {
            missing.push(id);
        }
        let leaf = leaf_rule(name).expect("leaf");
        for pattern in leaf.patterns() {
            let id = format!("xf:{name}/{}", pattern.name);
            if !objects.contains(&id) {
                missing.push(id);
            }
        }
    }
    assert!(
        missing.is_empty(),
        "inventory not in SSSOM ({}):\n{}",
        missing.len(),
        missing.join("\n")
    );
}
