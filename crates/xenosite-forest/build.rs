//! Gzip `mappings/xmet-forest.sssom.tsv` into OUT_DIR for `include_bytes!`.
//!
//! The uncompressed TSV at repo-root `mappings/` is the revision-controlled
//! SoT. Do not commit a `.gz` — this script builds it each compile.

use std::env;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let tsv = manifest_dir.join("../../mappings/xmet-forest.sssom.tsv");

    println!("cargo:rerun-if-changed={}", tsv.display());

    let raw = fs::read(&tsv).unwrap_or_else(|err| {
        panic!(
            "missing Forest↔XMET SSSOM at {} ({err}); expected repo-root mappings/xmet-forest.sssom.tsv",
            tsv.display()
        )
    });

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let gz_path = out_dir.join("xmet-forest.sssom.tsv.gz");
    let compressed = gzip_bytes(&raw);
    fs::write(&gz_path, &compressed).expect("write SSSOM gzip to OUT_DIR");
}

fn gzip_bytes(raw: &[u8]) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::GzEncoder;

    let mut encoder = GzEncoder::new(Vec::new(), Compression::best());
    encoder.write_all(raw).expect("gzip SSSOM");
    encoder.finish().expect("finish SSSOM gzip")
}
