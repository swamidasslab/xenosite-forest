/**
 * Node smoke: ForestMol, findPath, PhaseOne metabolize, randomPath.
 * Run: `npx tsx src/api.smoke.ts` (after wasm build).
 */
import {
  ForestMol,
  PhaseOne,
  ensureInit,
  findPath,
  randomPath,
} from "./index.js";

if (typeof document !== "undefined") {
  throw new Error("api smoke expects Node (server-side) without document");
}

await ensureInit();

const mol = new ForestMol("CCO");
if (!mol.csmi || typeof mol.csmi !== "string") {
  throw new Error("ForestMol.csmi missing");
}
if (mol.formula.charge !== 0) {
  throw new Error("ForestMol.formula.charge");
}

const [hits, counters] = findPath("CC", "CCO", { maxPaths: 1 });
if (!Array.isArray(hits) || hits.length < 1) {
  throw new Error(`findPath expected ≥1 hit, got ${hits?.length}`);
}
if (!hits[0]!.steps?.length) {
  throw new Error("findPath hit missing steps");
}
if (hits[0]!.steps[0]!.rule !== "Hydroxylation") {
  throw new Error(`expected Hydroxylation, got ${hits[0]!.steps[0]!.rule}`);
}
if (typeof counters.billed !== "number" || counters.billed < 1) {
  throw new Error(`counters.billed=${counters.billed}`);
}
if (counters.timed_out !== false) {
  throw new Error("unexpected timed_out");
}

const rules = PhaseOne();
if (rules.len < 1) {
  throw new Error("PhaseOne empty");
}
const benzene = new ForestMol("c1ccccc1");
const rows = rules.metabolize(benzene) as Array<{
  pattern_name: string;
  products: string[];
}>;
if (!Array.isArray(rows) || rows.length < 1) {
  throw new Error("PhaseOne.metabolize returned no rows");
}
if (!rows[0]!.products?.length) {
  throw new Error("metabolize row missing products");
}

const walk = randomPath("c1ccccc1", 1, { maxSteps: 1 });
if (!walk.smiles || !Array.isArray(walk.path) || walk.path.length < 1) {
  throw new Error("randomPath missing smiles/path");
}
if (!Array.isArray(walk.steps)) {
  throw new Error("randomPath missing steps");
}

console.log("ok: ForestMol, findPath, PhaseOne.metabolize, randomPath");
