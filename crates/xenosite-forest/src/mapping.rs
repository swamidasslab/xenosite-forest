//! Forest ↔ XMET SSSOM embed, CURIE/IRI helpers, and `resolve`.
//!
//! The living TSV lives at repo-root [`mappings/xmet-forest.sssom.tsv`].
//! `build.rs` gzips it into `OUT_DIR`; this module decompresses on demand.

use std::io::Read;
use std::sync::OnceLock;

use flate2::read::GzDecoder;

use crate::bound_pattern::BoundPattern;
use crate::rules::resolve_root;
use crate::ruleset::RuleSet;

/// Embedded gzip produced by `build.rs` from `mappings/xmet-forest.sssom.tsv`.
const SSSOM_GZ: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/xmet-forest.sssom.tsv.gz"));

static SSSOM_TEXT: OnceLock<String> = OnceLock::new();

/// Scheme IRI for `xf:` CURIEs.
pub const XF_PREFIX: &str = "https://w3id.org/xenosite/forest/";
/// Scheme IRI for `xmet:` CURIEs (documentation / mapping subjects).
pub const XMET_PREFIX: &str = "https://xenosite.org/ontology/xmet#";

/// Raw gzip bytes of the shipped Forest↔XMET SSSOM (as embedded).
pub fn forest_xmet_sssom_gz() -> &'static [u8] {
    SSSOM_GZ
}

/// Decompressed SSSOM TSV text (UTF-8). Lazily inflated once.
pub fn forest_xmet_sssom() -> &'static str {
    SSSOM_TEXT.get_or_init(|| {
        let mut decoder = GzDecoder::new(SSSOM_GZ);
        let mut text = String::new();
        decoder
            .read_to_string(&mut text)
            .expect("decompress embedded xmet-forest.sssom.tsv.gz");
        text
    })
}

/// Expand a CURIE (`xf:…` / `xmet:…`) or pass through an absolute IRI.
pub fn expand_iri(curie_or_iri: &str) -> String {
    if let Some(rest) = curie_or_iri.strip_prefix("xf:") {
        format!("{XF_PREFIX}{rest}")
    } else if let Some(rest) = curie_or_iri.strip_prefix("xmet:") {
        format!("{XMET_PREFIX}{rest}")
    } else {
        curie_or_iri.to_string()
    }
}

/// Compact an absolute `xf` / `xmet` IRI to a CURIE when possible.
pub fn to_curie(iri: &str) -> String {
    if let Some(rest) = iri.strip_prefix(XF_PREFIX) {
        format!("xf:{rest}")
    } else if let Some(rest) = iri.strip_prefix(XMET_PREFIX) {
        format!("xmet:{rest}")
    } else {
        iri.to_string()
    }
}

/// Normalize CURIE or IRI to `xf:` path segments (no prefix).
pub fn xf_path_segments(id: &str) -> Result<Vec<String>, String> {
    let path = if let Some(rest) = id.strip_prefix(XF_PREFIX) {
        rest.to_string()
    } else if let Some(rest) = id.strip_prefix("xf:") {
        rest.to_string()
    } else if id.starts_with("xmet:") {
        return Err("xmet: CURIEs are not Forest objects; use xf: path IDs".into());
    } else if id.starts_with("https://") || id.starts_with("http://") {
        return Err(format!("not an xf: Forest IRI: {id}"));
    } else {
        // Bare path "Tautomerization/tautomer_h"
        id.to_string()
    };
    let segs: Vec<String> = path
        .split('/')
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    if segs.is_empty() {
        return Err("empty xf: path".into());
    }
    Ok(segs)
}

/// One row of the embedded SSSOM (subset of columns we care about).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SssomRow {
    pub subject_id: String,
    pub predicate_id: String,
    pub object_id: String,
    pub mapping_justification: String,
    pub subject_label: String,
    pub object_label: String,
    pub always_with: String,
}

/// Parse the embedded SSSOM into rows (skips `#` comments and the header).
pub fn parse_forest_xmet_sssom() -> Vec<SssomRow> {
    let text = forest_xmet_sssom();
    let mut rows = Vec::new();
    let mut header: Option<Vec<&str>> = None;
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if header.is_none() {
            if parts.first().copied() == Some("subject_id") {
                header = Some(parts);
            }
            continue;
        }
        let hdr = header.as_ref().unwrap();
        let get = |name: &str| -> String {
            hdr.iter()
                .position(|h| *h == name)
                .and_then(|i| parts.get(i).copied())
                .unwrap_or("")
                .to_string()
        };
        rows.push(SssomRow {
            subject_id: get("subject_id"),
            predicate_id: get("predicate_id"),
            object_id: get("object_id"),
            mapping_justification: get("mapping_justification"),
            subject_label: get("subject_label"),
            object_label: get("object_label"),
            always_with: get("always_with"),
        });
    }
    rows
}

/// Resolved Forest object: a catalog/leaf rule or a pattern bound to its rule.
#[derive(Clone, Debug)]
#[allow(clippy::large_enum_variant)]
pub enum Resolved {
    Rule(RuleSet),
    Pattern(BoundPattern),
}

impl Resolved {
    pub fn as_rule(&self) -> Option<&RuleSet> {
        match self {
            Self::Rule(r) => Some(r),
            Self::Pattern(bp) => Some(bp.rule()),
        }
    }
}

/// Resolve an `xf:` CURIE / Forest IRI by walking public path segments.
///
/// Always returns Rust product-door objects (never native/legacy).
pub fn resolve(id: &str) -> Result<Resolved, String> {
    let segs = xf_path_segments(id)?;
    let mut current = resolve_root(&segs[0])
        .ok_or_else(|| format!("unknown Forest root segment: {}", segs[0]))?;

    if segs.len() == 1 {
        return Ok(Resolved::Rule(current));
    }

    for (i, seg) in segs[1..].iter().enumerate() {
        let is_last = i + 2 == segs.len();
        if current.is_catalog() {
            current = current.get_str(seg).ok_or_else(|| {
                format!(
                    "no member {seg:?} under {}",
                    current.name.as_deref().unwrap_or("?")
                )
            })?;
            if is_last {
                return Ok(Resolved::Rule(current));
            }
        } else if is_last {
            let bp = current.bound_pattern(seg).ok_or_else(|| {
                format!(
                    "no pattern {seg:?} under {}",
                    current.name.as_deref().unwrap_or("?")
                )
            })?;
            return Ok(Resolved::Pattern(bp));
        } else {
            return Err(format!(
                "unexpected extra segment {seg:?} under leaf {}",
                current.name.as_deref().unwrap_or("?")
            ));
        }
    }
    Ok(Resolved::Rule(current))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embed_round_trip_starts_with_comment_or_header() {
        let text = forest_xmet_sssom();
        assert!(
            text.contains("subject_id") && text.contains("xf:Tautomerization"),
            "embedded SSSOM missing expected content"
        );
        assert!(text.contains("xf:Tautomerization/tautomer_h"));
        assert!(text.contains("xf:Tautomerization/path_partner"));
    }

    #[test]
    fn expand_and_curie_round_trip() {
        assert_eq!(
            expand_iri("xf:Tautomerization/tautomer_h"),
            format!("{XF_PREFIX}Tautomerization/tautomer_h")
        );
        assert_eq!(
            to_curie(&format!("{XF_PREFIX}Tautomerization")),
            "xf:Tautomerization"
        );
    }

    #[test]
    fn resolve_leaf_and_pattern() {
        let rule = resolve("xf:Tautomerization").expect("leaf");
        assert!(matches!(rule, Resolved::Rule(_)));
        let bp = resolve("xf:Tautomerization/tautomer_h").expect("pattern");
        match bp {
            Resolved::Pattern(p) => {
                assert_eq!(p.name(), "tautomer_h");
                assert_eq!(p.rule_name(), Some("Tautomerization"));
            }
            Resolved::Rule(_) => panic!("expected BoundPattern"),
        }
    }
}
