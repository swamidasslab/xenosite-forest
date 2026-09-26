//! Strict boundary: this crate must not depend on xenosite-forest.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn cargo_toml_has_no_xenosite_forest_dependency() {
    let toml = fs::read_to_string(crate_dir().join("Cargo.toml")).unwrap();
    // Only flag dependency keys, not prose in description/repository URL.
    for line in toml.lines() {
        let t = line.trim();
        if t.starts_with('#') {
            continue;
        }
        let dep_line = t.starts_with("xenosite-forest")
            || t.starts_with("xenosite_forest")
            || t.contains("xenosite-forest =")
            || t.contains("xenosite_forest =");
        assert!(
            !dep_line,
            "xenosite-xrm must not depend on xenosite-forest: {t}"
        );
    }
}

#[test]
fn rust_sources_do_not_import_forest() {
    let src = crate_dir().join("src");
    for entry in fs::read_dir(&src).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|e| e.to_str()) != Some("rs") {
            continue;
        }
        let text = fs::read_to_string(&path).unwrap();
        for (i, line) in text.lines().enumerate() {
            let t = line.trim();
            if t.starts_with("//") {
                continue;
            }
            let imports = t.contains("use xenosite_forest")
                || t.contains("extern crate xenosite_forest")
                || t.contains("xenosite_forest::");
            assert!(
                !imports,
                "{}:{} imports forest: {t}",
                path.display(),
                i + 1
            );
        }
    }
}

#[test]
fn cargo_metadata_excludes_forest_from_xrm_deps() {
    let output = Command::new("cargo")
        .args(["metadata", "--format-version", "1", "--manifest-path"])
        .arg(crate_dir().join("Cargo.toml"))
        .output()
        .expect("cargo metadata");
    assert!(
        output.status.success(),
        "cargo metadata failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let meta: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let packages = meta["packages"].as_array().unwrap();
    let xrm = packages
        .iter()
        .find(|p| p["name"] == "xenosite-xrm")
        .expect("xenosite-xrm package");
    for dep in xrm["dependencies"].as_array().unwrap() {
        let name = dep["name"].as_str().unwrap_or("");
        assert_ne!(name, "xenosite-forest", "forbidden dependency");
    }
}
