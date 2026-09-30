//! Curate MetX hard-case bins under `tests/data/metx/` (Rust, nostereo call-site).
//!
//! ```text
//! cargo run -p xenosite-forest --example metx_curate_bins --release -- \
//!   artifacts/metx_hard_misses_nostereo.tsv \
//!   artifacts/metx_phase1_pairs.tsv \
//!   artifacts/metx_hard_misses.tsv \
//!   tests/data/metx
//! ```
//!
//! Bins:
//! - `thrash_rules_cover` — thrash/framing but PhaseOne should cover
//! - `thrash_rules_gap` — thrash; DB saturation; rules cannot cover
//! - `quiet_rules_gap` — no thrash; chemistry/rule gap
//! - `near_miss_significant_progress` — goal missed; large residual drop;
//!   `path_to_closest` recorded from a fresh `find_path_partial`
//!
//! Each TSV starts with `# commit=` / `# scan=` comment lines.
//! Provenance: MetXBioDB Phase I pairs (`artifacts/metx_phase1_pairs.tsv`).

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use chematic::chem::remove_stereo;
use xenosite_forest::{FindPathConfig, PathCounters, as_forest_mol, find_path_partial, phase_one};

const SCAN: &str =
    "metx_hard_scan_nostereo (taut=0, nostereo call-site, max_nodes=200, timeout=1.5s)";

const THRASH_GAP: &[&str] = &[
    "BIOTID00138",
    "BIOTID00198",
    "BIOTID00197",
    "BIOTID00199",
    "BIOTID00055",
];

const STEREO_FIXED: &[(&str, &str)] = &[
    (
        "BIOTID01090",
        "stereo_match: N-dealk one-step under nostereo",
    ),
    (
        "BIOTID00184",
        "stereo_match: dihydrodiol one-step under nostereo",
    ),
    (
        "BIOTID01318",
        "stereo_match: phenol epoxide one-step under nostereo",
    ),
];

const QUIET_GAP_CORE: &[(&str, &str)] = &[
    (
        "BIOTID01270",
        "missing_isoxazole_NO_cleavage (+2H ring open)",
    ),
    (
        "BIOTID00154",
        "epoxidation_regio: PhaseOne emits other arene oxides; MetX isomer absent",
    ),
    ("BIOTID00059", "tzd_ring_open+S_ox multi-atom rewrite"),
    ("BIOTID00058", "tzd_ring_open variant"),
    ("BIOTID00882", "N-dealkylation large fragment Δha≪0"),
    ("BIOTID01486", "TCE→chloral hydrate; not simple epoxide"),
];

#[derive(Clone, Default)]
struct MissRow {
    biot_id: String,
    bill: usize,
    residual_cost: usize,
    root_cost: usize,
    categories: String,
    closest: String,
    reactant_smi: String,
    product_smi: String,
    substrate: String,
    product: String,
    reaction_type: String,
}

#[derive(Clone)]
struct BinRow {
    biot_id: String,
    bin: String,
    reason: String,
    commit: String,
    scan: String,
    bill: String,
    residual_cost: String,
    root_cost: String,
    cost_drop: String,
    categories: String,
    path_to_closest: String,
    closest_smi: String,
    reaction_type: String,
    reactant_smi: String,
    product_smi: String,
    substrate: String,
    product: String,
}

fn git_commit() -> String {
    Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
        .unwrap_or_else(|| "UNKNOWN".into())
}

fn nostereo_smi(s: &str) -> Option<String> {
    let fm = as_forest_mol(s).ok()?;
    let stripped = remove_stereo(fm.mol());
    Some(as_forest_mol(stripped).ok()?.csmi().as_ref().to_string())
}

fn load_tsv_map(path: &Path) -> BTreeMap<String, BTreeMap<String, String>> {
    let f = File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut lines = BufReader::new(f).lines();
    let header = lines
        .next()
        .expect("empty tsv")
        .unwrap_or_else(|e| panic!("{e}"));
    // skip comment-only files' leading comments if any
    let header = header
        .strip_prefix('#')
        .map(|_| {
            lines
                .by_ref()
                .map(|l| l.unwrap())
                .find(|l| !l.starts_with('#') && !l.is_empty())
                .expect("no header after comments")
        })
        .unwrap_or(header);
    let cols: Vec<&str> = header.split('\t').collect();
    let mut out = BTreeMap::new();
    for line in lines {
        let line = line.unwrap();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < cols.len() {
            continue;
        }
        let mut row = BTreeMap::new();
        for (i, name) in cols.iter().enumerate() {
            row.insert((*name).to_string(), c[i].to_string());
        }
        if let Some(id) = row.get("biot_id") {
            out.insert(id.clone(), row);
        }
    }
    out
}

fn miss_from_map(m: &BTreeMap<String, String>) -> MissRow {
    MissRow {
        biot_id: m.get("biot_id").cloned().unwrap_or_default(),
        bill: m.get("bill").and_then(|s| s.parse().ok()).unwrap_or(0),
        residual_cost: m
            .get("residual_cost")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0),
        root_cost: m.get("root_cost").and_then(|s| s.parse().ok()).unwrap_or(0),
        categories: m.get("categories").cloned().unwrap_or_default(),
        closest: m.get("closest").cloned().unwrap_or_default(),
        reactant_smi: m.get("reactant_smi").cloned().unwrap_or_default(),
        product_smi: m.get("product_smi").cloned().unwrap_or_default(),
        substrate: m.get("substrate").cloned().unwrap_or_default(),
        product: m.get("product").cloned().unwrap_or_default(),
        reaction_type: m.get("reaction_type").cloned().unwrap_or_default(),
    }
}

fn path_format(steps: &[xenosite_forest::PathStep]) -> String {
    if steps.is_empty() {
        return "(no steps / root)".into();
    }
    steps
        .iter()
        .map(|s| format!("{}@{}", s.pattern_name, s.site))
        .collect::<Vec<_>>()
        .join(" > ")
}

fn rerun_partial(r_smi: &str, p_smi: &str) -> Option<(String, String, usize, String, usize)> {
    let r = nostereo_smi(r_smi)?;
    let p = nostereo_smi(p_smi)?;
    let set = phase_one();
    let mut c = PathCounters::default();
    let cfg = FindPathConfig {
        max_paths: 1,
        max_nodes: 200,
        use_atom_diff: true,
        lazy_closer: true,
        timeout: Some(Duration::from_secs_f64(2.0)),
        normalize_tautomer: false,
        ..FindPathConfig::default()
    };
    let out = find_path_partial(r.as_str(), p.as_str(), &set, &mut c, cfg, None, |_| true).ok()?;
    let pp = out.partials.first()?;
    Some((
        path_format(&pp.steps),
        pp.smiles.clone(),
        pp.residual.cost,
        pp.residual.categories.join(","),
        c.billed(),
    ))
}

fn write_bin(path: &Path, commit: &str, rows: &[BinRow]) {
    let mut f = File::create(path).unwrap_or_else(|e| panic!("create {}: {e}", path.display()));
    let short = if commit.len() >= 12 {
        &commit[..12]
    } else {
        commit
    };
    writeln!(f, "# commit={commit}").unwrap();
    writeln!(f, "# commit_short={short}").unwrap();
    writeln!(f, "# scan={SCAN}").unwrap();
    writeln!(
        f,
        "# provenance=MetXBioDB Phase I (artifacts/metx_phase1_pairs.tsv)"
    )
    .unwrap();
    writeln!(
        f,
        "# source=artifacts/metx_hard_misses_nostereo.tsv (+ phase1 / prior miss for stereo-fixed)"
    )
    .unwrap();
    writeln!(
        f,
        "biot_id\tbin\treason\tcommit\tscan\tbill\tresidual_cost\troot_cost\tcost_drop\tcategories\tpath_to_closest\tclosest_smi\treaction_type\treactant_smi\tproduct_smi\tsubstrate\tproduct"
    )
    .unwrap();
    for r in rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.biot_id,
            r.bin,
            r.reason.replace('\t', " "),
            r.commit,
            r.scan.replace('\t', " "),
            r.bill,
            r.residual_cost,
            r.root_cost,
            r.cost_drop,
            r.categories.replace('\t', " "),
            r.path_to_closest.replace('\t', " "),
            r.closest_smi.replace('\t', " "),
            r.reaction_type.replace('\t', " "),
            r.reactant_smi,
            r.product_smi,
            r.substrate.replace('\t', " "),
            r.product.replace('\t', " "),
        )
        .unwrap();
    }
    eprintln!("wrote {} n={}", path.display(), rows.len());
}

fn base_row(m: &MissRow, bin: &str, reason: &str, commit: &str) -> BinRow {
    let drop = m.root_cost.saturating_sub(m.residual_cost);
    BinRow {
        biot_id: m.biot_id.clone(),
        bin: bin.into(),
        reason: reason.into(),
        commit: commit.into(),
        scan: SCAN.into(),
        bill: m.bill.to_string(),
        residual_cost: m.residual_cost.to_string(),
        root_cost: m.root_cost.to_string(),
        cost_drop: drop.to_string(),
        categories: m.categories.clone(),
        path_to_closest: String::new(),
        closest_smi: m.closest.clone(),
        reaction_type: m.reaction_type.clone(),
        reactant_smi: m.reactant_smi.clone(),
        product_smi: m.product_smi.clone(),
        substrate: m.substrate.clone(),
        product: m.product.clone(),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let miss_tsv = PathBuf::from(
        args.get(1)
            .map(String::as_str)
            .unwrap_or("artifacts/metx_hard_misses_nostereo.tsv"),
    );
    let pairs_tsv = PathBuf::from(
        args.get(2)
            .map(String::as_str)
            .unwrap_or("artifacts/metx_phase1_pairs.tsv"),
    );
    let prior_miss_tsv = PathBuf::from(
        args.get(3)
            .map(String::as_str)
            .unwrap_or("artifacts/metx_hard_misses.tsv"),
    );
    let out_dir = PathBuf::from(args.get(4).map(String::as_str).unwrap_or("tests/data/metx"));

    let commit = git_commit();
    eprintln!("commit={commit}");
    eprintln!("miss_tsv={}", miss_tsv.display());

    let misses_raw = load_tsv_map(&miss_tsv);
    let pairs_raw = load_tsv_map(&pairs_tsv);
    let prior_raw = if prior_miss_tsv.exists() {
        load_tsv_map(&prior_miss_tsv)
    } else {
        BTreeMap::new()
    };

    let misses: BTreeMap<String, MissRow> = misses_raw
        .iter()
        .map(|(k, v)| (k.clone(), miss_from_map(v)))
        .collect();

    // --- thrash_rules_cover ---
    let mut thrash_cover = Vec::new();
    if let Some(m) = misses.get("BIOTID00966") {
        let mut row = base_row(
            m,
            "thrash_rules_cover",
            "tautomer_framing: PhaseOne ArOH emits imine-7-OH; amine product CSMI misses",
            &commit,
        );
        if let Some((path, closest, rcost, cats, _)) =
            rerun_partial(&m.reactant_smi, &m.product_smi)
        {
            row.path_to_closest = path;
            row.closest_smi = closest;
            row.residual_cost = rcost.to_string();
            if !cats.is_empty() {
                row.categories = cats;
            }
        }
        thrash_cover.push(row);
    }
    for (bid, reason) in STEREO_FIXED {
        let Some(p) = pairs_raw.get(*bid) else {
            eprintln!("skip stereo-fixed {bid}: not in pairs");
            continue;
        };
        let prior = prior_raw.get(*bid);
        let reason = if let Some(pr) = prior {
            format!(
                "{reason}; prior_stereo_miss bill={}",
                pr.get("bill").map(String::as_str).unwrap_or("?")
            )
        } else {
            (*reason).to_string()
        };
        thrash_cover.push(BinRow {
            biot_id: (*bid).into(),
            bin: "thrash_rules_cover".into(),
            reason,
            commit: commit.clone(),
            scan: SCAN.into(),
            bill: prior
                .and_then(|p| p.get("bill"))
                .cloned()
                .unwrap_or_else(|| "0".into()),
            residual_cost: prior
                .and_then(|p| p.get("residual_cost"))
                .cloned()
                .unwrap_or_else(|| "0".into()),
            root_cost: prior
                .and_then(|p| p.get("root_cost"))
                .cloned()
                .unwrap_or_else(|| "0".into()),
            cost_drop: String::new(),
            categories: prior
                .and_then(|p| p.get("categories"))
                .cloned()
                .unwrap_or_default(),
            path_to_closest: "(exact hit under nostereo)".into(),
            closest_smi: String::new(),
            reaction_type: p.get("reaction_type").cloned().unwrap_or_default(),
            reactant_smi: p.get("reactant_smi").cloned().unwrap_or_default(),
            product_smi: p.get("product_smi").cloned().unwrap_or_default(),
            substrate: p.get("substrate").cloned().unwrap_or_default(),
            product: p.get("product").cloned().unwrap_or_default(),
        });
    }

    // --- thrash_rules_gap ---
    let mut thrash_gap = Vec::new();
    for bid in THRASH_GAP {
        let Some(m) = misses.get(*bid) else {
            continue;
        };
        thrash_gap.push(base_row(
            m,
            "thrash_rules_gap",
            "db_saturation_product; label≠chemistry; PhaseOne cannot saturate",
            &commit,
        ));
    }

    // --- quiet_rules_gap ---
    let mut quiet_gap = Vec::new();
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for (bid, reason) in QUIET_GAP_CORE {
        let Some(m) = misses.get(*bid) else {
            continue;
        };
        quiet_gap.push(base_row(m, "quiet_rules_gap", reason, &commit));
        seen.insert((*bid).into());
    }
    for (bid, m) in &misses {
        if seen.contains(bid) || THRASH_GAP.contains(&bid.as_str()) {
            continue;
        }
        if m.bill >= 200 {
            continue;
        }
        let rxn = m.reaction_type.to_lowercase();
        let reason = if rxn.contains("n-o")
            || rxn.contains("ring open")
            || rxn.contains("thiazolidinedione")
        {
            Some(format!(
                "heterocycle_open / labeled open: {}",
                m.reaction_type
            ))
        } else if rxn.contains("n-hydrox") || rxn.contains("n-oxidation") {
            Some(format!("N-ox/N-OH gap?: {}", m.reaction_type))
        } else if m.root_cost.saturating_sub(m.residual_cost) == 0
            && m.categories.contains("cleavage")
            && m.root_cost >= 40
        {
            // skip generic large MCS fails
            None
        } else {
            None
        };
        // large fragment from ha via atom counts on nostereo mols
        let reason = reason.or_else(|| {
            let rh = as_forest_mol(&m.reactant_smi)
                .ok()
                .map(|f| f.heavy_atom_count())?;
            let ph = as_forest_mol(&m.product_smi)
                .ok()
                .map(|f| f.heavy_atom_count())?;
            let dha = ph as i32 - rh as i32;
            if dha <= -8 {
                Some(format!("large_fragment Δha={dha}"))
            } else {
                None
            }
        });
        if let Some(reason) = reason {
            quiet_gap.push(base_row(m, "quiet_rules_gap", &reason, &commit));
            seen.insert(bid.clone());
        }
    }

    // --- near_miss_significant_progress ---
    let mut near = Vec::new();
    for (bid, m) in &misses {
        if THRASH_GAP.contains(&bid.as_str()) {
            continue;
        }
        if m.residual_cost == 0 || m.closest.is_empty() {
            continue;
        }
        let drop = m.root_cost.saturating_sub(m.residual_cost);
        let candidate = (m.residual_cost <= 2 && drop >= 3) || (m.residual_cost <= 3 && drop >= 7);
        if !candidate {
            continue;
        }
        eprintln!(
            "near-miss pathing {bid} root={}→{} drop={drop}",
            m.root_cost, m.residual_cost
        );
        let Some((path, closest, rcost, cats, bill)) =
            rerun_partial(&m.reactant_smi, &m.product_smi)
        else {
            eprintln!("  no partial");
            continue;
        };
        // keep only verified close + real path
        if rcost > 3 || path == "(no steps / root)" {
            eprintln!("  drop re-run rcost={rcost} path={path}");
            continue;
        }
        let mut row = base_row(
            m,
            "near_miss_significant_progress",
            &format!(
                "root_cost {}→{rcost} (scan_drop={drop}); path to closest noted; goal not reached",
                m.root_cost
            ),
            &commit,
        );
        row.path_to_closest = path;
        row.closest_smi = closest;
        row.residual_cost = rcost.to_string();
        row.bill = bill.to_string();
        row.cost_drop = m.root_cost.saturating_sub(rcost).to_string();
        if !cats.is_empty() {
            row.categories = cats;
        }
        near.push(row);
    }
    near.sort_by(|a, b| {
        b.cost_drop
            .parse::<usize>()
            .unwrap_or(0)
            .cmp(&a.cost_drop.parse::<usize>().unwrap_or(0))
    });

    // attach near paths onto overlapping quiet rows
    let near_by: BTreeMap<_, _> = near
        .iter()
        .map(|r| (r.biot_id.clone(), r.clone()))
        .collect();
    for row in &mut quiet_gap {
        if let Some(n) = near_by.get(&row.biot_id) {
            if row.path_to_closest.is_empty() {
                row.path_to_closest = n.path_to_closest.clone();
                row.closest_smi = n.closest_smi.clone();
            }
        }
    }

    std::fs::create_dir_all(&out_dir).unwrap();
    write_bin(
        &out_dir.join("metx_thrash_rules_cover.tsv"),
        &commit,
        &thrash_cover,
    );
    write_bin(
        &out_dir.join("metx_thrash_rules_gap.tsv"),
        &commit,
        &thrash_gap,
    );
    write_bin(
        &out_dir.join("metx_quiet_rules_gap.tsv"),
        &commit,
        &quiet_gap,
    );
    write_bin(&out_dir.join("metx_near_miss_progress.tsv"), &commit, &near);

    let mut combined = Vec::new();
    combined.extend(thrash_cover);
    combined.extend(thrash_gap);
    combined.extend(quiet_gap);
    combined.extend(near);
    write_bin(&out_dir.join("metx_hard_case_bins.tsv"), &commit, &combined);

    eprintln!("done bins: cover / gap / quiet / near / combined");
}
