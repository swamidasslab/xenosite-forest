//! Dump proposed terms for feedback rounds.
//!
//! ```text
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
        label: "ethane → ethanol (aliphatic hydroxylation)",
        reactant: "CC",
        product: "CCO",
        tags: &[],
    },
    Case {
        label: "benzene → phenol (aromatic hydroxylation)",
        reactant: "c1ccccc1",
        product: "Oc1ccccc1",
        tags: &[],
    },
    Case {
        label: "chem:para-hydroxylation",
        reactant: "CCc1ccccc1",
        product: "CCc1ccc(O)cc1",
        tags: &["chem:para-hydroxylation", "chem:aromatic-hydroxylation"],
    },
    Case {
        label: "chem:benzylic-hydroxylation",
        reactant: "CCc1ccccc1",
        product: "CC(O)c1ccccc1",
        tags: &["chem:benzylic-hydroxylation"],
    },
    Case {
        label: "ethene → oxirane",
        reactant: "C=C",
        product: "C1CO1",
        tags: &[],
    },
    Case {
        label: "chem:arene-oxide + NIH-shift facets",
        reactant: "c1ccccc1",
        product: "Oc1ccccc1",
        tags: &["chem:arene-oxide", "chem:NIH-shift"],
    },
    Case {
        label: "ethanol → acetaldehyde",
        reactant: "CCO",
        product: "CC=O",
        tags: &[],
    },
    Case {
        label: "chem:N-demethylation",
        reactant: "CN(C)C",
        product: "CNC",
        tags: &["chem:N-demethylation"],
    },
    Case {
        label: "chem:oxidative-deamination",
        reactant: "CCN",
        product: "CC=O",
        tags: &["chem:oxidative-deamination"],
    },
    Case {
        label: "chem:O-demethylation",
        reactant: "COc1ccccc1",
        product: "Oc1ccccc1",
        tags: &["chem:O-demethylation"],
    },
    Case {
        label: "chem:acyl-glucuronidation",
        reactant: "CC(=O)O",
        product: "CC(=O)O",
        tags: &["chem:acyl-glucuronidation"],
    },
    Case {
        label: "chem:phenolic-glucuronidation",
        reactant: "Oc1ccccc1",
        product: "Oc1ccccc1",
        tags: &["chem:phenolic-glucuronidation"],
    },
    Case {
        label: "chem:GSH-Michael + conjugate-addition facet",
        reactant: "C=CC=O",
        product: "C=CC=O",
        tags: &["chem:GSH-Michael"],
    },
    Case {
        label: "chem:quinone-formation (dearomatization + bioactivation)",
        reactant: "c1ccccc1",
        product: "O=C1C=CC(=O)C=C1",
        tags: &["chem:quinone-formation"],
    },
    Case {
        label: "chem:quinone-imine + one-step quinone formation",
        reactant: "CC(=O)Nc1ccc(O)cc1",
        product: "CC(=O)N=C1C=CC(=O)C=C1",
        tags: &["chem:quinone-imine", "chem:one-step-quinone-formation"],
    },
    Case {
        label: "chem:two-step-quinone-formation",
        reactant: "c1ccccc1",
        product: "O=C1C=CC(=O)C=C1",
        tags: &["chem:two-step-quinone-formation"],
    },
    Case {
        label: "chem:imine-methide",
        reactant: "Nc1ccc(C)cc1",
        product: "N=C1C=CC(=C)C=C1",
        tags: &["chem:imine-methide"],
    },
    Case {
        label: "chem:dearomatization alone",
        reactant: "c1ccccc1",
        product: "C1=CC=CC=C1",
        tags: &["chem:dearomatization"],
    },
    Case {
        label: "chem:nitroaromatic-reduction",
        reactant: "O=[N+]([O-])c1ccccc1",
        product: "Nc1ccccc1",
        tags: &["chem:nitroaromatic-reduction"],
    },
    Case {
        label: "chem:cyanide-hydrolysis",
        reactant: "CC#N",
        product: "CC(=O)O",
        tags: &["chem:cyanide-hydrolysis"],
    },
    Case {
        label: "chem:carbonyl-reduction",
        reactant: "CC(=O)C",
        product: "CC(O)C",
        tags: &["chem:carbonyl-reduction"],
    },
    Case {
        label: "chem:glycine-conjugation",
        reactant: "c1ccccc1C(=O)O",
        product: "c1ccccc1C(=O)O",
        tags: &["chem:glycine-conjugation"],
    },
    Case {
        label: "chem:tautomerization",
        reactant: "CC(=O)C",
        product: "CC(O)=C",
        tags: &["chem:tautomerization"],
    },
    Case {
        label: "chem:regio-ambiguity + aromatic hydroxylation",
        reactant: "CCc1ccccc1",
        product: "CCc1ccc(O)cc1",
        tags: &[
            "chem:aromatic-hydroxylation",
            "chem:regio-ambiguity",
            "chem:competing-type",
        ],
    },
    Case {
        label: "chem:pathway-depth-ambiguity + two-step quinone",
        reactant: "c1ccccc1",
        product: "O=C1C=CC(=O)C=C1",
        tags: &[
            "chem:two-step-quinone-formation",
            "chem:pathway-depth-ambiguity",
            "chem:intermediate-underspecified",
        ],
    },
    Case {
        label: "chem:mapping-underspecified",
        reactant: "CC",
        product: "CCO",
        tags: &["chem:hydroxylation", "chem:mapping-underspecified"],
    },
];

fn main() {
    let write = env::args().any(|a| a == "--write");
    let namer = Namer::from_manifest(DEFAULT_MANIFEST).expect("manifest");
    let mut out = String::from("# Proposed terms — feedback round 004\n\n");
    out.push_str("ChatGPT-aligned spines: metabolism phase, chemical transformation, ");
    out.push_str("Rainbow phase I family, phase II conjugation family, medchem liability, ");
    out.push_str("reactive metabolite family, site type, structural delta, product status, ");
    out.push_str("rule provenance, evidence, biological context; plus ambiguity and ");
    out.push_str("Metabolic Forest map (alias). Site-localized names are templates, not concepts.\n\n");
    out.push_str("See `data/ontology/ANNOTATION.md` for bundle/link design.\n\n");
    out.push_str("Regenerate: `cargo run -p xenosite-xrm --example sample_terms -- --write`\n\n");

    for case in CASES {
        let terms = namer
            .name_smiles(case.reactant, case.product, case.tags)
            .unwrap_or_else(|e| panic!("{}: {e}", case.label));
        let block = format_case(case, &namer.format_sample_lines(&terms));
        print!("{block}");
        out.push_str(&block);
    }

    if write {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/samples/round-004.md");
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
