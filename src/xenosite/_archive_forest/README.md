# Archived Metabolic Forest (pre-POC swap) — READ-ONLY

Temporary archive of the pre-refactor `xenosite.forest` implementation.

**Locked:** do not modify files in this tree (except this README). Historical
reference only for parity / H2H / StepPlan apply until those call sites move.

The live public package is `xenosite.forest` (promoted from `refactor_poc`).
Import as `xenosite._archive_forest`. Not a supported public API.

**CI:** archived sources and `tests/_archive_forest/` are excluded from pytest
collection, coverage, ruff, and pyright. Do not add them back to CI gates.
