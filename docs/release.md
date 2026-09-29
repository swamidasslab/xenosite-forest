# Release

One package: ``xenosite-forest`` is a maturin mixed wheel (Python under
``src/xenosite`` + Rust extension ``xenosite.forest._rust``). The wheel version is
``crates/xenosite-forest/Cargo.toml`` ``package.version``. Tag releases as
``vX.Y.Z`` with that Cargo version set to ``X.Y.Z`` on the tagged commit.

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
./scripts/vendor-chematic.sh
uv sync --extra rdkit --extra network --group dev
```

## Cut a release

1. Set ``version = "X.Y.Z"`` in ``crates/xenosite-forest/Cargo.toml``.
2. Merge to ``main`` with fragments committed.
3. Tag and push:

```bash
git tag -a v0.8.1 -m "0.8.1"
git push origin v0.8.1
```

Pushing `v*` runs [`.github/workflows/release.yml`](../.github/workflows/release.yml):
tests, platform wheels (manylinux / musllinux / macOS / Windows, x86_64 +
aarch64), sdist, GitHub Release, PyPI upload, and towncrier compiling
`CHANGELOG.md` onto the default branch when the branch ruleset allows.
