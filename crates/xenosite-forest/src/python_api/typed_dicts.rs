//! Runtime side of [`crate::export`] `typed_dict!` structs.
//!
//! Defines each one as a `typing.TypedDict` on `_rust` (the generated stub
//! lists them in `__all__`), and converts [`AttrValue`] to / from Python.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::export::{AttrValue, TypedDictInfo};

/// Define each `typed_dict!` struct as a runtime `typing.TypedDict` on `m`.
///
/// So `from xenosite.forest._rust import HopDict` works outside
/// `TYPE_CHECKING` too. Field types are `typing.Any` at runtime; the precise
/// types live in the stub.
pub fn add_typed_dicts(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let typing = py.import("typing")?;
    let typed_dict = typing.getattr("TypedDict")?;
    let any = typing.getattr("Any")?;
    for td in inventory::iter::<TypedDictInfo> {
        let fields = pyo3::types::PyDict::new(py);
        for f in td.fields {
            fields.set_item(f.name, &any)?;
        }
        let cls = typed_dict.call1((td.name, fields))?;
        cls.setattr("__module__", m.name()?)?;
        cls.setattr("__doc__", td.doc.trim())?;
        m.add(td.name, cls)?;
    }
    Ok(())
}

impl<'py> IntoPyObject<'py> for AttrValue {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(match self {
            Self::Bool(b) => b.into_pyobject(py)?.to_owned().into_any(),
            Self::Int(i) => i.into_pyobject(py)?.into_any(),
            Self::Str(s) => s.into_pyobject(py)?.into_any(),
        })
    }
}

impl<'py> FromPyObject<'py> for AttrValue {
    fn extract_bound(value: &Bound<'py, PyAny>) -> PyResult<Self> {
        // bool before int: Python ``bool`` is an ``int`` subclass.
        if let Ok(b) = value.extract::<bool>() {
            return Ok(Self::Bool(b));
        }
        if let Ok(i) = value.extract::<i64>() {
            return Ok(Self::Int(i));
        }
        if let Ok(s) = value.extract::<String>() {
            return Ok(Self::Str(s));
        }
        Err(PyValueError::new_err(
            "attr value must be bool, int, or str",
        ))
    }
}
