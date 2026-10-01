/**
 * Node smoke: ForestMol, findPath, PhaseOne / resolve metabolize, randomPath.
 * Run: `npx tsx src/api.smoke.ts` (after wasm build).
 */
import {
  BoundPattern,
  ForestMol,
  PhaseOne,
  Reactivity,
  RuleSet,
  ensureInit,
  expandIri,
  findPath,
  randomPath,
  resolve,
  toCurie,
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

const viaResolve = resolve("xf:PhaseOne");
if (!(viaResolve instanceof RuleSet)) {
  throw new Error("resolve(xf:PhaseOne) expected RuleSet");
}
const anisole = new ForestMol("c1ccccc1OC");
const resolveRows = viaResolve.metabolize(anisole) as Array<{
  pattern_name: string;
  site: number;
  products: string[];
  rule_path: Array<string | null>;
}>;
if (!Array.isArray(resolveRows) || resolveRows.length < 1) {
  throw new Error("resolve(xf:PhaseOne).metabolize returned no rows");
}
const first = resolveRows[0]!;
if (!first.products?.length || typeof first.site !== "number") {
  throw new Error("resolve metabolize row missing products/site");
}
if (!Array.isArray(first.rule_path)) {
  throw new Error("resolve metabolize row missing rule_path");
}

const leaf = resolve("xf:PhaseOne/StableOxygenation/Hydroxylation/h");
if (!(leaf instanceof BoundPattern)) {
  throw new Error("resolve(…/h) expected BoundPattern");
}
if (!leaf.curie.includes("Hydroxylation")) {
  throw new Error(`unexpected BoundPattern.curie=${leaf.curie}`);
}
const leafRows = leaf.metabolize(benzene) as Array<{ products: string[] }>;
if (!Array.isArray(leafRows) || leafRows.length < 1) {
  throw new Error("BoundPattern.metabolize returned no rows");
}

const iri = expandIri("xf:PhaseOne");
if (!iri.includes("PhaseOne") || toCurie(iri) !== "xf:PhaseOne") {
  throw new Error(`expandIri/toCurie roundtrip failed: ${iri}`);
}

const react = Reactivity();
const epoxide = new ForestMol("C1OC1c1ccccc1");
const conj = react.metabolize(epoxide) as Array<{
  pattern_name: string;
  products: string[];
}>;
if (!Array.isArray(conj) || conj.length < 1) {
  throw new Error("Reactivity.metabolize returned no rows");
}
const cx = conj.flatMap((r) => r.products ?? []).find((p) => p.includes("*"));
if (!cx) {
  throw new Error("Reactivity products missing *");
}
if (!cx.includes("|")) {
  throw new Error(`expected CXSMILES block, got ${cx}`);
}
if (!["GSH", "Protein", "DNA", "Cyanide"].some((lab) => cx.includes(lab))) {
  throw new Error(`expected conjugate atomLabel in ${cx}`);
}

const walk = randomPath("c1ccccc1", 1, { maxSteps: 1 });
if (!walk.smiles || !Array.isArray(walk.path) || walk.path.length < 1) {
  throw new Error("randomPath missing smiles/path");
}
if (!Array.isArray(walk.steps)) {
  throw new Error("randomPath missing steps");
}

console.log(
  "ok: ForestMol, findPath, PhaseOne/resolve.metabolize, Reactivity CX stars, randomPath"
);
