# Publishing xenosite-forest packages

## Python (PyPI)

Tag `vX.Y.Z` matching `crates/xenosite-forest/Cargo.toml`. Workflow
[`release.yml`](workflows/release.yml) runs tests, builds manylinux / musllinux /
macOS / Windows wheels (x86_64 + aarch64) + sdist, then publishes to PyPI via
Trusted Publishing (OIDC, environment `pypi`).

## JavaScript (GitHub Packages)

Same `v*` tag also publishes **`@swamidasslab/forest`** to GitHub Packages
(npm). Source `js/package.json` keeps `"name": "@xenosite/forest"`; the release
job rewrites the scope at publish time (GH Packages requires `@OWNER/...`).

Auth is `GITHUB_TOKEN` (`packages: write`). No `NPM_TOKEN`.

```bash
# .npmrc
@swamidasslab:registry=https://npm.pkg.github.com
//npm.pkg.github.com/:_authToken=${GITHUB_TOKEN}

npm install @swamidasslab/forest
```

Package URL: https://github.com/swamidasslab/xenosite-forest/pkgs/npm/forest

Local check:

```bash
./scripts/build_wasm.sh
(cd js && npm install && npm test)
```
