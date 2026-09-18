# Release

Package version comes from git tags ([hatch-vcs](https://github.com/ofek/hatch-vcs)), not a static `version` in `pyproject.toml`. Do not bump a version field. Do not run `towncrier build` locally for a real release (hatch-vcs would write a `.devN` version into `CHANGELOG.md`).

## During development

User-facing changes need a towncrier fragment in `changelog.d/` in the same unit of work:

```bash
uv run towncrier create --no-edit -c "Short description." added.md
```

Types: `added` | `changed` | `fixed` | `removed` | `deprecated` | `security`. Skip fragments for internal tests, refactors, and tooling.

## Cut a release

On `main`, with fragments committed:

```bash
git tag -a v0.5.0 -m "0.5.0"
git push origin v0.5.0
```

Pushing `v*` runs [`.github/workflows/release.yml`](../.github/workflows/release.yml): tests, a clean `X.Y.Z` version, `uv build`, GitHub Release, PyPI upload, and towncrier compiling `CHANGELOG.md` onto the default branch (new commit; the tag is not moved).
