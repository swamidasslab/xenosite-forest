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
import { init, ForestMol, resolve } from "@swamidasslab/forest";

await init();

const rules = resolve("xf:PhaseOne"); // RuleSet (or BoundPattern for a leaf CURIE)
const smi = "c1ccccc1OC";

/** @returns {Array<[string[], number, Array<string|null>]>} */
function products(rules, reactant) {
  const rows = rules.metabolize(new ForestMol(reactant));
  // rows: { pattern_name, site, products: string[], rule_path: (string|null)[] }[]
  return rows.map((emit) => {
    const path = [...emit.rule_path].reverse();
    path.push(emit.pattern_name);
    return [emit.products, emit.site, path];
  });
}

for (const row of products(rules, smi)) {
  console.log(row);
}
```

## Develop

```bash
../scripts/build_wasm.sh
npm install
npm test
```
