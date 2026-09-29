# Archived Metabolic Forest (pre-swap) — READ-ONLY

Temporary archive of the previous `xenosite.forest` implementation.

**Locked:** do not modify files in this tree (except this README). Historical
reference only for parity / H2H / StepPlan apply until those call sites move.

The live public package is `xenosite.forest` (Rust). The RDKit 0.7-era engine
`xenosite.forest.native` is separately **feature-frozen** (see
`docs/forest/NATIVE.md`). Import the archive as `xenosite.forest.legacy`.
Neither is a supported surface for new product features.

**CI:** archived sources and `tests/forest/legacy/` are excluded from pytest
collection, coverage, ruff, and pyright. Do not add them back to CI gates.
