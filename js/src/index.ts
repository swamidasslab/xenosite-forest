/**
 * `@xenosite/forest` — Metabolic Forest for the web / Node (WASM).
 *
 * Thin wrapper mirroring the Python public stub (`find_path`, `random_path`,
 * PhaseOne / leaf factories). Call {@link init} once before use.
 */

import { initNative, requireReady } from "./native.js";
import {
  BoundPattern,
  ForestMol,
  Formula,
  PatternInfo,
  RuleSet,
  dealkylation as wasmDealkylation,
  default_ruleset as wasmDefaultRuleset,
  dehydrogenation as wasmDehydrogenation,
  epoxidation as wasmEpoxidation,
  epoxide_opening as wasmEpoxideOpening,
  expand_iri as wasmExpandIri,
  find_path_with_options as wasmFindPathWithOptions,
  forest_canon_smiles as wasmForestCanonSmiles,
  forest_hydroxylate as wasmForestHydroxylate,
  forest_xmet_sssom as wasmForestXmetSssom,
  hydrolysis as wasmHydrolysis,
  hydroxylation as wasmHydroxylation,
  n_dealkylation as wasmNDealkylation,
  phase_one as wasmPhaseOne,
  quinone_formation as wasmQuinoneFormation,
  random_path as wasmRandomPath,
  random_path_with_options as wasmRandomPathWithOptions,
  reactivity as wasmReactivity,
  resolve as wasmResolve,
  to_curie as wasmToCurie,
} from "./wasm/xenosite_forest.js";

export { BoundPattern, ForestMol, Formula, PatternInfo, RuleSet };
export { initNative as init, isNativeReady } from "./native.js";

export type FindPathOptions = {
  maxPaths?: number;
  maxNodes?: number;
  useAtomDiff?: boolean;
  lazyCloser?: boolean;
  diversity?: boolean;
  dropSkeletonTwins?: boolean;
  score?: string;
  timeout?: number;
  // snake_case aliases (Python parity)
  max_paths?: number;
  max_nodes?: number;
  use_atom_diff?: boolean;
  lazy_closer?: boolean;
  drop_skeleton_twins?: boolean;
};

export type FindPathHit = {
  smiles: string;
  steps: Array<{ rule: string; site: string[] }>;
};

export type FindPathCounters = {
  nodes: number;
  mol_edits: number;
  expansions: number;
  billed: number;
  dropped_duplicate_plan: number;
  dropped_exact_plan: number;
  dropped_skeleton_twin: number;
  diversity_repush: number;
  unstable_csmi_key: number;
  timed_out: boolean;
};

export type RandomPathOptions = {
  maxSteps?: number;
  skipMulticomponent?: boolean;
  skipSeen?: boolean;
  max_steps?: number;
  skip_multicomponent?: boolean;
  skip_seen?: boolean;
  ruleset?: RuleSet;
};

export type RandomPathOutcome = {
  smiles: string;
  path: string[];
  steps: Array<{
    rule: string;
    pattern: string;
    site: number[];
    products: string[];
    chosen: number;
  }>;
  patterns: Array<{
    name: string;
    smarts: string;
    cleaves: boolean;
    search_bias: number;
  }>;
};

/** Idempotent wasm init (Node reads the `.wasm` from disk). */
export async function ensureInit(): Promise<void> {
  await initNative();
}

export function findPath(
  reactant: string,
  target: string,
  options: FindPathOptions = {}
): [FindPathHit[], FindPathCounters] {
  requireReady();
  const bag = {
    maxPaths: options.maxPaths ?? options.max_paths,
    maxNodes: options.maxNodes ?? options.max_nodes,
    useAtomDiff: options.useAtomDiff ?? options.use_atom_diff,
    lazyCloser: options.lazyCloser ?? options.lazy_closer,
    diversity: options.diversity,
    dropSkeletonTwins: options.dropSkeletonTwins ?? options.drop_skeleton_twins,
    score: options.score,
    timeout: options.timeout,
  };
  const out = wasmFindPathWithOptions(reactant, target, bag) as [
    FindPathHit[],
    FindPathCounters,
  ];
  return out;
}

export function randomPath(
  reactant: string,
  seed: number | bigint,
  options: RandomPathOptions = {}
): RandomPathOutcome {
  requireReady();
  const seedBig = typeof seed === "bigint" ? seed : BigInt(seed);
  if (options.ruleset) {
    return wasmRandomPath(
      reactant,
      seedBig,
      options.maxSteps ?? options.max_steps ?? null,
      options.ruleset,
      options.skipMulticomponent ?? options.skip_multicomponent ?? null,
      options.skipSeen ?? options.skip_seen ?? null
    ) as RandomPathOutcome;
  }
  return wasmRandomPathWithOptions(reactant, seedBig, {
    maxSteps: options.maxSteps ?? options.max_steps,
    skipMulticomponent: options.skipMulticomponent ?? options.skip_multicomponent,
    skipSeen: options.skipSeen ?? options.skip_seen,
  }) as RandomPathOutcome;
}

export function PhaseOne(): RuleSet {
  requireReady();
  return wasmPhaseOne();
}

export function Reactivity(): RuleSet {
  requireReady();
  return wasmReactivity();
}

export function Epoxidation(): RuleSet {
  requireReady();
  return wasmEpoxidation();
}

export function QuinoneFormation(): RuleSet {
  requireReady();
  return wasmQuinoneFormation();
}

export function EpoxideOpening(): RuleSet {
  requireReady();
  return wasmEpoxideOpening();
}

export function NDealkylation(): RuleSet {
  requireReady();
  return wasmNDealkylation();
}

export function Hydroxylation(): RuleSet {
  requireReady();
  return wasmHydroxylation();
}

export function Dehydrogenation(): RuleSet {
  requireReady();
  return wasmDehydrogenation();
}

export function Dealkylation(): RuleSet {
  requireReady();
  return wasmDealkylation();
}

export function Hydrolysis(): RuleSet {
  requireReady();
  return wasmHydrolysis();
}

export function DefaultRuleset(): RuleSet {
  requireReady();
  return wasmDefaultRuleset();
}

export function forestCanonSmiles(smiles: string): string {
  requireReady();
  return wasmForestCanonSmiles(smiles);
}

export function forestHydroxylate(smiles: string): string {
  requireReady();
  return wasmForestHydroxylate(smiles);
}

/** Resolve an `xf:` CURIE / Forest IRI to a {@link RuleSet} or {@link BoundPattern}. */
export function resolve(id: string): RuleSet | BoundPattern {
  requireReady();
  return wasmResolve(id) as RuleSet | BoundPattern;
}

/** Expand `xf:` / `xmet:` CURIEs (absolute IRIs pass through). */
export function expandIri(curieOrIri: string): string {
  requireReady();
  return wasmExpandIri(curieOrIri);
}

/** Compact an absolute `xf` / `xmet` IRI to a CURIE when possible. */
export function toCurie(iri: string): string {
  requireReady();
  return wasmToCurie(iri);
}

/** Decompressed Forest↔XMET SSSOM TSV text. */
export function forestXmetSssom(): string {
  requireReady();
  return wasmForestXmetSssom();
}
