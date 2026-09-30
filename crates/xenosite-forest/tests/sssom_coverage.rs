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

/// Every Forest object IRI/CURIE in the SSSOM must resolve.
///
/// This is the gate for **dangling xf: references** in the mapping file
/// (object_id as CURIE or absolute IRI). Inventory coverage the other way
/// (catalog → SSSOM) is [`inventory_fully_covered_in_sssom`].
///
/// Nested catalogs (`xf:PhaseOne/…`, `xf:Reactivity/…`) resolve to Rule
/// sets for intermediate segments and BoundPattern only for the leaf
/// pattern tip. [`BoundPattern::curie`] stays `xf:<Leaf>/<pattern>` (short
/// leaf path), not the full catalog CURIE.
#[test]
fn all_sssom_objects_resolve() {
    let rows = parse_forest_xmet_sssom();
    assert!(!rows.is_empty());
    let mut failures = Vec::new();
    for row in &rows {
        let id = &row.object_id;
        // CURIE form (as stored) and expanded absolute IRI must both resolve.
        let expanded = xenosite_forest::expand_iri(id);
        for form in [id.as_str(), expanded.as_str()] {
            match resolve(form) {
                Ok(Resolved::Rule(_)) => {}
                Ok(Resolved::Pattern(bp)) => {
                    let local = id.strip_prefix("xf:").unwrap_or(id);
                    let tip = local.rsplit('/').next().unwrap_or(local);
                    if bp.name() != tip {
                        failures.push(format!(
                            "{form}: pattern tip mismatch {} vs {tip}",
                            bp.name()
                        ));
                    }
                }
                Err(e) => failures.push(format!("{form}: {e}")),
            }
        }
    }
    assert!(
        failures.is_empty(),
        "SSSOM object IRIs failed to resolve ({}):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn unknown_sssom_object_iri_does_not_resolve() {
    assert!(
        resolve("xf:NotARealLeaf/not_a_pattern").is_err(),
        "unknown xf: IRI must fail resolve"
    );
    assert!(
        resolve("https://w3id.org/xenosite/forest/NotARealLeaf").is_err(),
        "unknown absolute Forest IRI must fail resolve"
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
        // Only pattern tips (resolve → BoundPattern) need a leaf parent.
        let Ok(Resolved::Pattern(bp)) = resolve(id) else {
            continue;
        };
        let Some((parent_path, pat)) = rest.rsplit_once('/') else {
            failures.push(format!("{id}: pattern without parent path"));
            continue;
        };
        let parent = format!("xf:{parent_path}");
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
                assert_eq!(a.name(), bp.name());
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
#[ignore = "xfail: new Effect SMARTS-split PatternInfo names not yet in xmet-forest.sssom.tsv — do not remove until next version bump (map rows then)"]
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
