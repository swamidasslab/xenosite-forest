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
    for (i, line) in text.lines().enumerate() {
        if i >= MAX_LINES || out.len() + line.len() + 1 > MAX_CHARS {
            out.push_str("… (truncated)");
            return out;
        }
        if !out.is_empty() {
            out.push('\n');
        }
        out.push_str(line);
    }
    out
}

fn indent(level: usize) -> String {
    "  ".repeat(level)
}

/// Indented RuleSet / pattern hierarchy (full tree — catalogs are browsed whole).
pub fn format_ruleset(set: &RuleSet) -> String {
    let mut lines = Vec::new();
    format_ruleset_into(set, 0, &mut lines);
    lines.join("\n")
}

fn format_ruleset_into(set: &RuleSet, level: usize, lines: &mut Vec<String>) {
    let pad = indent(level);
    let n_pat = set.patterns().len();
    let n_mem = set.members().len();
    match &set.name {
        Some(n) => {
            lines.push(format!("{pad}{n}  [{n_mem} members, {n_pat} patterns]"));
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
        // Unnamed composite (e.g. phase_one()): print named members at this level.
        None => {
            for member in set.members() {
                match member {
                    RuleMember::Pattern(p) => {
                        let smarts = if p.smarts.len() > 48 {
                            format!("{}…", &p.smarts[..47])
                        } else {
                            p.smarts.clone()
                        };
                        lines.push(format!("{pad}{}  {}", p.name, smarts));
                    }
                    RuleMember::Set(child) => format_ruleset_into(child, level, lines),
                }
            }
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
        format!("  formula charge={} atoms={}", formula.charge, n),
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

/// Soft cap: expand linearization orders in text/HTML only when small.
pub const MAX_LIN_EXPAND: usize = 4;

pub fn format_step_plan(plan: &Deps) -> String {
    let n_lin = plan.n_linearizations();
    let mut lines = vec![format!(
        "StepPlan  {} steps, ~{} linearization(s)",
        plan.steps().len(),
        n_lin
    )];
    for (i, step) in plan.steps().iter().enumerate() {
        lines.push(format!("  {i}: {}", format_plan_step(step)));
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
    if n_lin > 0 && n_lin <= MAX_LIN_EXPAND {
        lines.push("  linearizations:".into());
        for (li, lin) in plan.linearizations().iter().enumerate() {
            let order: Vec<String> = lin
                .steps
                .iter()
                .map(|s| format!("{}@{}", s.rule, format_site_compact(&s.site)))
                .collect();
            lines.push(format!("    [{li}] {}", order.join(" → ")));
        }
    } else if n_lin > MAX_LIN_EXPAND {
        lines.push(format!(
            "  (~{n_lin} linearizations; not expanded — use .linearizations())"
        ));
    }
    let maybe = plan.maybe();
    if !maybe.is_empty() {
        lines.push(format!("  maybe: {} bag(s)", maybe.entries.len()));
        for (i, e) in maybe.entries.iter().take(8).enumerate() {
            let site: Vec<String> = e.site.iter().map(|a| a.to_string()).collect();
            let opens: Vec<String> = e
                .opens
                .iter()
                .take(4)
                .map(|o| {
                    format!(
                        "{{{}}}",
                        o.iter()
                            .map(|a| a.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect();
            let open_s = if e.opens.is_empty() {
                String::new()
            } else if e.opens.len() > 4 {
                format!(" opens=[{},…]", opens.join(","))
            } else {
                format!(" opens=[{}]", opens.join(","))
            };
            lines.push(format!(
                "    [{i}] site={{{}}} side={}{open_s}",
                site.join(","),
                e.side
            ));
        }
        if maybe.entries.len() > 8 {
            lines.push(format!("    … (+{} bags)", maybe.entries.len() - 8));
        }
    }
    truncate_display(&lines.join("\n"))
}

fn format_site_compact(site: &[crate::canonical_plan::PlanAtom]) -> String {
    let parts: Vec<String> = site.iter().take(6).map(|a| format!("{a:?}")).collect();
    if site.len() > 6 {
        format!("[{},…]", parts.join(","))
    } else {
        format!("[{}]", parts.join(","))
    }
}

fn format_plan_step(step: &crate::canonical_plan::Step) -> String {
    format!("{} @ {}", step.rule, format_site_compact(&step.site))
}

pub fn format_path_counters(c: &crate::find_path::PathCounters) -> String {
    truncate_display(&format!(
        "PathCounters  billed={} nodes={} edits={} expansions={}\n  dropped_dup={} (exact={} skeleton={}) signal_contained={}\n  diversity_repush={} timed_out={}",
        c.billed(),
        c.nodes,
        c.mol_edits,
        c.expansions,
        c.dropped_duplicate_plan,
        c.dropped_exact_plan,
        c.dropped_skeleton_twin,
        c.signal_contained_plan,
        c.diversity_repush,
        c.timed_out
    ))
}

pub fn format_partial_outcome(
    smiles: &str,
    residual_cost: usize,
    n_plan_steps: usize,
    n_lin: usize,
) -> String {
    truncate_display(&format!(
        "PartialOutcome  {smiles}\n  residual_cost={residual_cost}  plan_steps={n_plan_steps}  ~{n_lin} linearization(s)\n  (display: ForestMol tag trace + residual; use .to_dict() for full residual)"
    ))
}

pub fn format_pattern_info(p: &crate::pattern::PatternInfo) -> String {
    let smarts = if p.smarts.len() > 64 {
        format!("{}…", &p.smarts[..63])
    } else {
        p.smarts.clone()
    };
    truncate_display(&format!(
        "PatternInfo  {}\n  smarts: {}\n  adds={:?} removes={:?} cleaves={} methide={}",
        p.name, smarts, p.effect.adds, p.effect.removes, p.effect.cleaves, p.effect.methide
    ))
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
        .rev() // outer→leaf for display (PhaseOne/QuinoneFormation)
        .filter_map(|p| p.as_deref())
        .filter(|p| *p != "Default")
        .collect::<Vec<_>>()
        .join("/");
    truncate_display(&format!(
        "Emission  {}  site={} on {}\n  path: {}\n  products: {}",
        em.pattern_name,
        em.site,
        em.reactant.csmi().as_ref(),
        path,
        prod
    ))
}

pub fn format_path_outcome(
    smiles: &str,
    n_steps: usize,
    n_plan_steps: usize,
    n_lin: usize,
) -> String {
    truncate_display(&format!(
        "PathOutcome  {smiles}\n  walk_hops={n_steps}  plan_steps={n_plan_steps}  ~{n_lin} linearization(s)\n  (display: reactant→product hop trail with SOM on each reactant; use .hops() / .plan)"
    ))
}
