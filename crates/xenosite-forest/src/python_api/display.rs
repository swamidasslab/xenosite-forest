//! Truncated text formatters for notebook / REPL display.

use crate::canonical_plan::Deps;
use crate::forest_mol::ForestMol;
use crate::metabolic_network::MetabolicNetwork;
use crate::pattern::Emission;
use crate::ruleset::{RuleMember, RuleSet};

/// Soft cap on lines for `__str__` / `__repr__` trees.
pub const MAX_LINES: usize = 64;
/// Soft cap on characters.
pub const MAX_CHARS: usize = 3500;

/// Truncate a multi-line display string; append an ellipsis line if cut.
pub fn truncate_display(text: &str) -> String {
    let mut out = String::new();
    let mut lines = 0usize;
    for line in text.lines() {
        if lines >= MAX_LINES || out.len() + line.len() + 1 > MAX_CHARS {
            out.push_str("… (truncated)");
            return out;
        }
        if !out.is_empty() {
            out.push('\n');
        }
        out.push_str(line);
        lines += 1;
    }
    out
}

fn indent(level: usize) -> String {
    "  ".repeat(level)
}

/// Tab-indented RuleSet / pattern hierarchy.
pub fn format_ruleset(set: &RuleSet) -> String {
    let mut lines = Vec::new();
    format_ruleset_into(set, 0, &mut lines);
    truncate_display(&lines.join("\n"))
}

fn format_ruleset_into(set: &RuleSet, level: usize, lines: &mut Vec<String>) {
    let pad = indent(level);
    let title = match &set.name {
        Some(n) => format!("{pad}{n}"),
        None => format!("{pad}(unnamed)"),
    };
    let n_pat = set.patterns().len();
    let n_mem = set.members().len();
    lines.push(format!("{title}  [{n_mem} members, {n_pat} patterns]"));
    for member in set.members() {
        match member {
            RuleMember::Pattern(p) => {
                let smarts = if p.smarts.len() > 48 {
                    format!("{}…", &p.smarts[..47])
                } else {
                    p.smarts.clone()
                };
                lines.push(format!("{}  {}  {}", indent(level), p.name, smarts));
            }
            RuleMember::Set(child) => format_ruleset_into(child, level + 1, lines),
        }
    }
}

pub fn format_forest_mol(mol: &ForestMol) -> String {
    let formula = mol.formula();
    let n = mol.mol().atom_count();
    let mut survivors = 0usize;
    let mut born = 0usize;
    let mut untagged = 0usize;
    let tags: Vec<String> = (0..n)
        .map(|i| match mol.tag_of(i) {
            Some(t) => {
                let raw = t.get();
                if raw < mol.stamp_end() {
                    survivors += 1;
                } else {
                    born += 1;
                }
                raw.to_string()
            }
            None => {
                untagged += 1;
                "-".into()
            }
        })
        .collect();
    let tag_s = if tags.len() > 24 {
        format!("{}… ({} atoms)", tags[..24].join(","), tags.len())
    } else {
        tags.join(",")
    };
    let parts = [
        format!("ForestMol {}", mol.csmi()),
        format!(
            "  formula charge={} atoms={}",
            formula.charge, n
        ),
        format!(
            "  trace stamp_end={} survivors={} born={} untagged={}",
            mol.stamp_end(),
            survivors,
            born,
            untagged
        ),
        format!("  tags [{tag_s}]"),
    ];
    truncate_display(&parts.join("\n"))
}

pub fn format_network(net: &MetabolicNetwork) -> String {
    let root = net.root_csmi.as_deref().unwrap_or("(none)");
    let mut target_csmis: Vec<String> = net
        .targets
        .iter()
        .filter_map(|&i| net.nodes.get(i).map(|n| n.csmi.clone()))
        .collect();
    target_csmis.sort();
    let reached = target_csmis.iter().any(|t| net.reaches(t));
    let targets = if target_csmis.is_empty() {
        "(none marked)".into()
    } else if target_csmis.len() <= 4 {
        target_csmis.join(", ")
    } else {
        format!(
            "{}, … ({} targets)",
            target_csmis[..3].join(", "),
            target_csmis.len()
        )
    };
    truncate_display(&format!(
        "MetabolicNetwork  nodes={} edges={}\n  root: {}\n  targets: {}\n  target_reached: {}",
        net.n_nodes(),
        net.n_edges(),
        root,
        targets,
        if target_csmis.is_empty() {
            "n/a"
        } else if reached {
            "yes"
        } else {
            "no"
        }
    ))
}

pub fn format_step_plan(plan: &Deps) -> String {
    let n_lin = plan.n_linearizations();
    let mut lines = vec![format!(
        "StepPlan  {} steps, ~{} linearization(s)",
        plan.steps().len(),
        n_lin
    )];
    for (i, step) in plan.steps().iter().enumerate() {
        let site: Vec<String> = step.site.iter().map(|a| format!("{a:?}")).collect();
        let site_s = if site.len() > 6 {
            format!("[{}, …]", site[..6].join(", "))
        } else {
            format!("[{}]", site.join(", "))
        };
        lines.push(format!("  {i}: {} @ {}", step.rule, site_s));
    }
    let precedes = plan.precedes();
    if !precedes.is_empty() {
        let edge_s: Vec<String> = precedes
            .iter()
            .take(12)
            .map(|&(a, b)| format!("{a}→{b}"))
            .collect();
        let more = if precedes.len() > 12 { " …" } else { "" };
        lines.push(format!("  precedes: {}{more}", edge_s.join(", ")));
    }
    truncate_display(&lines.join("\n"))
}

pub fn format_emission(em: &Emission) -> String {
    let products = em.product_csmis();
    let prod = if products.len() <= 3 {
        products.join(" | ")
    } else {
        format!(
            "{} | … ({} products)",
            products[..2].join(" | "),
            products.len()
        )
    };
    let path = em
        .rule_path
        .iter()
        .map(|p| p.as_deref().unwrap_or("?"))
        .collect::<Vec<_>>()
        .join("/");
    truncate_display(&format!(
        "Emission  {}  site={}\n  path: {}\n  products: {}",
        em.pattern_name, em.site, path, prod
    ))
}

pub fn format_path_outcome(smiles: &str, n_steps: usize, n_plan_steps: usize, n_lin: usize) -> String {
    truncate_display(&format!(
        "PathOutcome  {smiles}\n  walk_hops={n_steps}  plan_steps={n_plan_steps}  ~{n_lin} linearization(s)\n  (display: ForestMol tag trace + StepPlan; use .hops() for walk trail)"
    ))
}
