# Publishing xenosite-forest packages

## Python (PyPI)

Tag `vX.Y.Z`. Workflow [`release.yml`](workflows/release.yml) sets
`crates/xenosite-forest/Cargo.toml` from the tag, runs tests, builds manylinux /
musllinux / macOS / Windows wheels (x86_64 + aarch64) + sdist, then publishes to
PyPI via Trusted Publishing (OIDC, environment `pypi`).

## Release housekeeping push (`RELEASE_PUSH_TOKEN`)

Protect main blocks `GITHUB_TOKEN` from pushing the post-tag version /
CHANGELOG commit. Use a PAT for an account on the ruleset bypass list
(currently `swamidass`).

1. Create a token (pick one):
   - **Classic:** https://github.com/settings/tokens/new — scope `repo`,
     note e.g. `xenosite-forest release push`, expiry you are willing to rotate.
   - **Fine-grained:** https://github.com/settings/personal-access-tokens/new
     — Resource owner your user; only repository `swamidasslab/xenosite-forest`;
     Permissions → Repository → Contents: Read and write.
2. Store it as a repo secret (paste when prompted; do not commit the token):

```bash
gh secret set RELEASE_PUSH_TOKEN --repo swamidasslab/xenosite-forest
```

3. Confirm:

```bash
gh secret list --repo swamidasslab/xenosite-forest | grep RELEASE_PUSH_TOKEN
```

Rotate by creating a new PAT and re-running `gh secret set`.

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
