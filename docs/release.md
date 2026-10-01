# Release

One package: ``xenosite-forest`` is a maturin mixed wheel (Python under
``src/xenosite`` + Rust extension ``xenosite.forest._rust``). The published
version is the ``vX.Y.Z`` git tag. Release CI rewrites
``crates/xenosite-forest/Cargo.toml`` (plus ``Cargo.lock`` and
``js/package.json``) from that tag before building; you do not bump Cargo by
hand.

Do not run ``towncrier build`` locally for a real release (Protect main blocks
direct pushes; the tag workflow compiles fragments, or open a changelog PR).

## During development

User-facing changes need a towncrier fragment in `changelog.d/` in the same unit of work:

```bash
uv run towncrier create --no-edit -c "Short description." added.md
```

Types: `added` | `changed` | `fixed` | `removed` | `deprecated` | `security`. Skip fragments for internal tests, refactors, and tooling.

Local editable install (builds the extension):

```bash
uv sync --extra rdkit --extra network --group dev
```

Between releases, ``Cargo.toml`` may still show the last cut until the next
tag’s housekeeping commit lands; that is fine for development.

## Cut a release

1. Merge to ``main`` with fragments committed.
2. Tag and push:

```bash
git tag -a v0.8.1 -m "0.8.1"
git push origin v0.8.1
```

Pushing `v*` runs [`.github/workflows/release.yml`](../.github/workflows/release.yml):
[`scripts/set_release_version.py`](../scripts/set_release_version.py) sets
package versions from the tag, then tests, platform wheels (manylinux /
musllinux / macOS / Windows, x86_64 + aarch64), sdist, GitHub Release, PyPI
upload, JS publish, and a housekeeping commit on the default branch
(`CHANGELOG.md` + Cargo/JS versions). That push needs repo secret
`RELEASE_PUSH_TOKEN` (see [`.github/PUBLISH.md`](../.github/PUBLISH.md)).

Wheels are built **only** on that tag workflow — not on PR or push to
`main`. CI profile:

| Trigger | What runs |
|---|---|
| PR | Lint only (ruff / pyright / rustfmt / clippy / WASM compile) |
| `main` | Single Python 3.12: `make test-python-smoke` + `make test-rust-smoke` + JS |
| `v*` tag | `make test-python` on 3.11–3.14 (no legacy); publish blocked on fail |

Local defaults: `make test` / `make test-python` exclude frozen legacy and
native↔legacy `test_parity.py`. Smoke is opt-in (`make test-python-smoke` /
`make test-rust-smoke` / `make test-smoke`). Legacy archive + native↔legacy
pairs: `make test-python-legacy` (not CI). `test-rust-smoke` omits the
`ms1_apply_fuzz` integration binary (full suite on tags / `make test-rust`
still runs it). CI non-test gates (rustfmt, clippy, WASM compile, ruff,
pyright) are `make check` — run that before pushing; `make test*` does not.

### Main smoke contents

**Include**

- Full Rust crate: `make test-rust` and `make test-rust-python`
- Python product door: `tests/forest/rust/`
- Thin native↔Rust parity: `test_rule_parity_pairs.py`,
  `test_rule_parity_fuzz.py` with ``XENOSITE_PARITY_FULL=0`` /
  ``pytest --parity-focused`` (CoverIntent corpus only), and
  `test_rust_parity_coverage.py` (coverage-substrate product sets)
- JS wrapper: `npm test` (tsc + API smoke)

**Exclude (never on CI; opt-in locally)**

- Entire `tests/forest/legacy/` and native↔legacy `test_parity.py`
  (`make test-python-legacy`)
- Native-only RDKit suites on main smoke only (goldens, Hypothesis fuzz, CLI,
  ported archives, `test_canonical_plan`, …) — still on tag `make test-python`
- Full rule×mol cartesian native↔Rust parity on main smoke
  (`XENOSITE_PARITY_FULL=1`); tag / local `make test-python` keep the full
  cartesian

## Forest ↔ XMET SSSOM snapshot

Living SoT: [`mappings/xmet-forest.sssom.tsv`](../mappings/xmet-forest.sssom.tsv)
(uncompressed only — never commit `.gz`). See
[`docs/forest/MAPPING.md`](forest/MAPPING.md).

On each `vX.Y.Z` cut, after the tag is green:

1. Copy the TSV into `xenosite-xmet` as
   `data/mappings/forest/xmet-forest.sssom.vX.Y.Z.tsv`.
2. Point that repo’s living `data/mappings/xmet-forest.sssom.tsv` at the same
   content.
3. Confirm xmet forest-pattern coverage still passes on the snapshot.
