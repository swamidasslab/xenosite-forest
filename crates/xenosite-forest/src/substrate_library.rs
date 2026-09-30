//! Shared substrate SMILES for PatternInfo / When coverage scans.
//!
//! **SoT:** [`coverage_substrates.txt`](../../../tests/data/coverage_substrates.txt)
//! (sections `[library]` / `[pattern]`). Seeded from the historical native
//! union (phase1 SMARTS probes, formula-hint probes, conjugation probes,
//! rare OR/When fillers). Catalog Effect / When coverage reads
//! [`coverage_candidates`] — not per-leaf
//! [`crate::rules::LEAF_EXAMPLE_SUBSTRATES`] (those are the short site_kind /
//! `_example_substrates` list).
//!
//! Expand the data file's `[pattern]` section when a PatternInfo or When arm
//! stays mute. Do not bury probes in test modules.

use std::sync::OnceLock;

const RAW: &str = include_str!("../../../tests/data/coverage_substrates.txt");

fn parse_sections() -> (Vec<&'static str>, Vec<&'static str>) {
    let mut library = Vec::new();
    let mut pattern = Vec::new();
    let mut section: Option<&str> = None;
    for line in RAW.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Only known headers — SMILES may contain [...] atom specs.
        if line == "[library]" || line == "[pattern]" {
            section = Some(&line[1..line.len() - 1]);
            continue;
        }
        match section {
            Some("library") => library.push(line),
            Some("pattern") => pattern.push(line),
            other => panic!(
                "coverage_substrates.txt: SMILES {line:?} outside [library]/[pattern] (section={other:?})"
            ),
        }
    }
    assert!(
        !library.is_empty(),
        "coverage_substrates.txt: empty [library] section"
    );
    assert!(
        !pattern.is_empty(),
        "coverage_substrates.txt: empty [pattern] section"
    );
    (library, pattern)
}

fn sections() -> &'static (Vec<&'static str>, Vec<&'static str>) {
    static SECTIONS: OnceLock<(Vec<&'static str>, Vec<&'static str>)> = OnceLock::new();
    SECTIONS.get_or_init(parse_sections)
}

/// Phase-I / quinone / conjugation probes + drug-like crashers.
pub fn substrate_library() -> &'static [&'static str] {
    sections().0.as_slice()
}

/// Extras beyond [`substrate_library`] for rare When / OR branches.
pub fn pattern_substrates() -> &'static [&'static str] {
    sections().1.as_slice()
}

/// Ordered unique union — catalog Effect / When coverage pool.
pub fn coverage_candidates() -> &'static [&'static str] {
    static CANDIDATES: OnceLock<Vec<&'static str>> = OnceLock::new();
    CANDIDATES
        .get_or_init(|| {
            let mut out =
                Vec::with_capacity(substrate_library().len() + pattern_substrates().len());
            for &s in substrate_library()
                .iter()
                .chain(pattern_substrates().iter())
            {
                if !out.contains(&s) {
                    out.push(s);
                }
            }
            out
        })
        .as_slice()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ForestMol;

    #[test]
    fn coverage_file_sections_parse_and_union() {
        let lib = substrate_library();
        let pat = pattern_substrates();
        let cand = coverage_candidates();
        assert_eq!(lib.len(), 54, "library section size (seeded native SoT)");
        assert!(
            pat.len() >= 76,
            "pattern extras size (seeded native + mute-arm fillers), got {}",
            pat.len()
        );
        assert_eq!(
            cand.len(),
            {
                let mut s = std::collections::BTreeSet::new();
                s.extend(lib.iter().copied());
                s.extend(pat.iter().copied());
                s.len()
            },
            "candidates must be unique union of library + pattern"
        );
        assert!(
            cand.len() >= 130,
            "coverage pool should stay at least the native union size, got {}",
            cand.len()
        );
        // Library entries appear first, in order.
        assert_eq!(&cand[..lib.len()], lib);
    }

    #[test]
    fn every_coverage_candidate_parses() {
        let mut bad = Vec::new();
        for &smi in coverage_candidates() {
            if ForestMol::parse(smi).is_err() {
                bad.push(smi);
            }
        }
        assert!(
            bad.is_empty(),
            "coverage_substrates.txt SMILES failed to parse:\n  {}",
            bad.join("\n  ")
        );
    }
}
