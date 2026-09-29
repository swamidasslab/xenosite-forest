# Migrating to 0.8.0 — package layout

Breaking import and install changes from **0.7.x** → **0.8.0**. Chemistry notes
for the earlier 0.7 rewrite remain in [MIGRATING_0.7.md](MIGRATING_0.7.md).

## Feature policy (read this first)

**New features go only to the Rust door** (`crates/xenosite-forest` + public
`xenosite.forest` wrapper). Do **not** add them to `xenosite.forest.native` or
`xenosite.forest.legacy`. Native is a frozen RDKit reference for parity;
legacy is a frozen 0.6.x archive. Details: [NATIVE.md](NATIVE.md),
[LEGACY.md](LEGACY.md).

## Layout

| Import | Role |
|--------|------|
| `xenosite.forest` | **Product API:** Rust chematic door (`find_path`, allowlisted rulesets). New work lands here. |
| `xenosite.forest.native` | RDKit Python engine. **Reference / validation only — feature-frozen.** Optional `[rdkit]` extra. Emits `UserWarning` on import. |
| `xenosite.forest.legacy` | Frozen 0.6.x archive. Optional `[rdkit]` extra. Emits `DeprecationWarning` on import. See [LEGACY.md](LEGACY.md). |

API divergence across the public stub / native / legacy is expected. Do not force
shared signatures, and do not back-port Rust APIs into native.

## Public stub allowlist

```python
from xenosite.forest import (
    find_path,
    available,
    PhaseOne,
    Epoxidation,
    QuinoneFormation,
    EpoxideOpening,
    NDealkylation,
)
```

Other catalog factories (`hydroxylation`, `dealkylation`, …) remain on the
package for advanced use but are **not** in root `__all__`.

## Install

```bash
# Recommended (Rust door; no RDKit)
pip install xenosite-forest

# RDKit reference + legacy archive
pip install 'xenosite-forest[rdkit]'
```

## Import map (0.7 → 0.8)

| Was (0.7) | Is (0.8) |
|-----------|----------|
| `xenosite.forest.find_path` (RDKit) | `xenosite.forest.native.find_path` |
| `xenosite.forest.find_path_rust` / `native_available` | **Removed.** Use `xenosite.forest.find_path` / `available` |
| `xenosite.forest.rules` / `rulesets` / `bfs` / … | `xenosite.forest.native.…` |
| `xenosite._archive_forest` | `xenosite.forest.legacy` |
| `StepPlan` on the root stub | `xenosite.forest.legacy.step_plan` (native may still import it) |

## Removed in 0.8

- Module `xenosite.forest.find_path_rust` and names `find_path_rust` / `native_available`.
- Nested `xenosite.forest.rust` wrapper package (API lives on `xenosite.forest`).
- Root re-exports of RDKit `bfs` / `rules` / `StepPlan` / wide ruleset catalog.
