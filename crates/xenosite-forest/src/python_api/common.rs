//! Shared PyO3 helpers for the product API.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub const NETWORK_LOCK_POISONED: &str = "MetabolicNetwork lock poisoned";

pub fn py_err(err: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(err.to_string())
}

pub fn network_lock_err() -> PyErr {
    PyValueError::new_err(NETWORK_LOCK_POISONED)
}
