# `@xenosite/forest` (JS / WASM)

Metabolic Forest for Node and the browser. Same public surface as the Python
stub: `findPath`, `randomPath`, `PhaseOne`, leaf factories, `ForestMol` /
`RuleSet`.

## Install (GitHub Packages)

Published as **`@swamidasslab/forest`** (GitHub requires `@OWNER/...`).

```bash
# .npmrc
@swamidasslab:registry=https://npm.pkg.github.com
//npm.pkg.github.com/:_authToken=${GITHUB_TOKEN}

npm install @swamidasslab/forest
```

## Usage

```ts
import { init, findPath, PhaseOne, ForestMol } from "@swamidasslab/forest";

await init();
const [hits, counters] = findPath("CC", "CCO", { maxPaths: 1 });
const rows = PhaseOne().metabolize(new ForestMol("c1ccccc1"));
```

## Develop

```bash
../scripts/build_wasm.sh
npm install
npm test
```
