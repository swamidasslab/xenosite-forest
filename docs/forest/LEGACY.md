# Legacy forest (frozen 0.6.x)

`xenosite.forest.legacy` is the archived Metabolic Forest from the `v0.6.1`
era. It is **not** evolved toward 0.7 / Rust APIs. Importing the package emits
a `DeprecationWarning`. Prefer `xenosite.forest` (Rust).

Native may still import `xenosite.forest.legacy.step_plan` (`StepPlan` / `Deps`)
until those types are re-homed.

## Freeze vs `v0.6.1`

Compared `v0.6.1:src/xenosite/forest` to `src/xenosite/forest/legacy` after the
re-home (promote / rename from `_archive_forest`).

| Path | Status | Notes |
|------|--------|-------|
| `base.py`, `bfs.py`, `guided_path.py`, `net.py`, `path_context.py`, `phaseone.py`, `rules.py`, `rulesets.py`, `step_plan.py`, `trace.py`, `unstable.py`, `utils.py` | approved | Byte-identical to `v0.6.1` (package path only). |
| `README.md` | approved | Archive note; not in `v0.6.1`. |
| `__init__.py` | approved | Docstring rewrite; drop `__version__` re-export (version lives on `xenosite.forest`); add `DeprecationWarning` on import; RDKit `ImportError` guard; import-order only for `FormulaHint`. No chemistry or API changes. |
| Absolute imports (`xenosite.forest.*` → `xenosite.forest.legacy.*`) in `net.py`, `base.py`, `guided_path.py`, `phaseone.py`, `rules.py`, `rulesets.py` | approved | Required by the package re-home; chemistry unchanged. |

Do not silently “fix” or modernize modules under `legacy/`. Record any further
intentional drift here with `Status: approved`, `not approved`, or `not decided`.
