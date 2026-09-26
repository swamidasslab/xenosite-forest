//! Error types for config load and naming.

use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
    #[error("config: {0}")]
    Config(String),
    #[error("chemistry: {0}")]
    Chemistry(String),
    #[error("unknown concept: {0}")]
    UnknownConcept(String),
}
