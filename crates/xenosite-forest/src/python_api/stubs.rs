//! `.pyi` generation for `xenosite.forest._rust` (pyo3-stub-gen).
//!
//! Signatures come from `#[gen_stub_*]` on the pyclasses / pyfunctions;
//! `TypedDict`s come from [`crate::export`] `typed_dict!` structs. Output is
//! one file, `src/xenosite/forest/_rust.pyi`, checked by `make check-stubs`.
//!
//! `stubs` feature only (dev tool): never compiled into wheels.

use std::path::{Path, PathBuf};

use pyo3_stub_gen::generate::{ClassDef, MemberDef};
use pyo3_stub_gen::{StubInfo, TypeInfo};

use crate::export::TypedDictInfo;

/// Workspace `pyproject.toml` (maturin `module-name` / `python-source`).
fn pyproject_toml() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../pyproject.toml")
}

/// Path of the generated stub, relative to the workspace root.
pub const STUB_PATH: &str = "src/xenosite/forest/_rust.pyi";

/// Absolute path of the generated stub.
pub fn stub_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(STUB_PATH)
}

fn leak(s: String) -> &'static str {
    Box::leak(s.into_boxed_str())
}

fn resolved(mut t: TypeInfo, module: &str) -> TypeInfo {
    t.resolve_default_module(module);
    t
}

fn typed_dict_class(info: &TypedDictInfo, module: &str) -> ClassDef {
    ClassDef {
        name: info.name,
        module: None,
        doc: leak(info.doc.trim().to_string()),
        attrs: info
            .fields
            .iter()
            .map(|f| MemberDef {
                name: f.name,
                r#type: resolved((f.r#type)(), module),
                doc: leak(f.doc.trim().to_string()),
                default: None,
                deprecated: None,
            })
            .collect(),
        getter_setters: Default::default(),
        methods: Default::default(),
        bases: vec![TypeInfo::with_module("typing.TypedDict", "typing".into())],
        classes: Vec::new(),
        match_args: None,
        // TypedDict cannot be ``@typing.final``.
        subclass: true,
    }
}

/// Gather pyclass / pyfunction stubs plus `typed_dict!` classes.
pub fn stub_info() -> pyo3_stub_gen::Result<StubInfo> {
    let mut info = StubInfo::from_pyproject_toml(pyproject_toml())?;
    let module_name = info.default_module_name.clone();
    let module = info
        .modules
        .get_mut(&module_name)
        .ok_or_else(|| std::io::Error::other(format!("no stub module {module_name:?}")))?;
    for td in inventory::iter::<TypedDictInfo> {
        module
            .class
            .insert((td.type_id)(), typed_dict_class(td, &module_name));
    }
    Ok(info)
}

/// Rendered `_rust.pyi` text.
pub fn render_stub() -> pyo3_stub_gen::Result<String> {
    let info = stub_info()?;
    let module = &info.modules[&info.default_module_name];
    Ok(module.format_with_config(info.config.use_type_statement))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Drift gate: the checked-in stub must match the Rust annotations.
    #[test]
    fn checked_in_stub_is_current() {
        let want = render_stub().expect("render stub");
        let have = std::fs::read_to_string(stub_path()).unwrap_or_default();
        assert!(
            have == want,
            "{STUB_PATH} is stale — run `make stubs` and commit the result"
        );
    }
}
