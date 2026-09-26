//! Dump proposed terms for a fixed panel. Re-run when changing ontology/assignments
//! so reviewers can inject feedback between rounds.
//!
//! ```text
//! cargo run -p xenosite-xrm --example sample_terms
//! cargo run -p xenosite-xrm --example sample_terms -- --write
//! ```

use std::env;
use std::fs;
use std::path::PathBuf;
use xenosite_xrm::{Namer, DEFAULT_MANIFEST};

struct Case {
    label: &'static str,
    reactant: &'static str,
    product: &'static str,
    tags: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        label: "ethane → ethanol",
        reactant: "CC",
        product: "CCO",
        tags: &[],
    },
    Case {
        label: "benzene → phenol",
        reactant: "c1ccccc1",
        product: "Oc1ccccc1",
        tags: &[],
    },
    Case {
        label: "ethene → oxirane",
        reactant: "C=C",
        product: "C1CO1",
        tags: &[],
    },
    Case {
        label: "ethanol → acetaldehyde",
        reactant: "CCO",
        product: "CC=O",
        tags: &[],
    },
    Case {
        label: "tag: glucuronidation",
        reactant: "CCO",
        product: "CCO",
        tags: &["forest.rule:Glucuronidation"],
    },
    Case {
        label: "tag: GSH Michael",
        reactant: "C=CC=O",
        product: "C=CC=O",
        tags: &["forest.pattern:Glutathionation/michael"],
    },
    Case {
        label: "tag: N-dealkylation",
        reactant: "CCN(C)C",
        product: "CCNC",
        tags: &["forest.rule:NDealkylation"],
    },
    Case {
        label: "tag: quinone formation",
        reactant: "c1ccccc1",
        product: "O=C1C=CC(=O)C=C1",
        tags: &["forest.rule:QuinoneFormation"],
    },
];

fn main() {
    let write = env::args().any(|a| a == "--write");
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).expect("manifest");
    let mut out = String::from("# Proposed terms — feedback round\n\n");
    out.push_str("Regenerate: `cargo run -p xenosite-xrm --example sample_terms -- --write`\n\n");
    out.push_str("Comment inline or open an issue noting the case label + wanted change.\n\n");

    for case in CASES {
        let terms = namer
            .name_smiles(case.reactant, case.product, case.tags)
            .unwrap_or_else(|e| panic!("{}: {e}", case.label));
        let block = format_case(case, &namer.format_sample_lines(&terms));
        print!("{block}");
        out.push_str(&block);
    }

    if write {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/samples/round-001.md");
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, &out).unwrap();
        eprintln!("wrote {}", path.display());
    }
}

fn format_case(case: &Case, lines: &[String]) -> String {
    let mut s = format!("## {}\n", case.label);
    s.push_str(&format!("- reactant: `{}`\n", case.reactant));
    s.push_str(&format!("- product: `{}`\n", case.product));
    if case.tags.is_empty() {
        s.push_str("- tags: _(none)_\n");
    } else {
        s.push_str(&format!("- tags: `{}`\n", case.tags.join("`, `")));
    }
    s.push_str("- terms:\n");
    if lines.is_empty() {
        s.push_str("  - _(none)_\n");
    } else {
        for line in lines {
            s.push_str(&format!("  - {line}\n"));
        }
    }
    s.push('\n');
    s
}
