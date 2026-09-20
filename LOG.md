# Lab log

## 2026-09-19

- Typed `span` as closed `Span` (`_SpanCore` total=True + optional partner/symbol/h/site_aromatic). Bare vs disagreeing branches are `T | tuple[T, ...]`. No type-level `object` left under `refactor_poc` / `tests/refactor_poc` (English “object” in docstrings only).
- Replaced other `object` sites with `ProductInfo`, `TraceInfo`, `PatternInfo`, `SiteInfo`, `EditCounters`, `SitesOn`, `SiteSignature`, `EffectField`, `phase1: None`. `GetMolFrags` frags args are `list[int] | None` / `list[list[int]] | None`.
- pyright poc modules: 0 errors. `pytest tests/refactor_poc -q`: 674 passed. No `Any`.

- Moved composite phase-I plans off `find_path` onto `ReactionRule.canonical_plan` (default identity; `QuinoneFormation` → prep + dehydrogenation via `canonical_plan.py`). Search only calls the hook. Halogen ends read `partner` and emit `OxidativeDehalogenation`. `as_deps` keys `AddedRef` on the full prep site.
- Parameterized plan-apply checks in `tests/refactor_poc/test_canonical_plan.py` over `QUICK_SUBSTRATES` / `SUBSTRATE_LIBRARY`: identity re-metabolize at the plan site; composite forest `Linearization.apply` on a plain mol. `pytest tests/refactor_poc -q`: 674 passed.

- `FilterSites` second arg is now `SiteInfo = SmartsSiteInfo | PairSiteInfo` (closed, total=True) in `records.py`. Shared core: site, rule, options. SMARTS/resonance adds rxn_num, pattern; pair adds ends, end_atoms, end_maps, path_ends.
- Dropped unused ad-hoc keys that closed schemas reject: pair `info["path"]`, and `merge_effects` keys `partners` / `aromatic` (never read). Per-end detail stays on `ends`.
- pyright poc modules: 0 errors. `pytest tests/refactor_poc -q`: 126 passed.

## 2026-09-19

- refactor_poc suite: 126 passed. No skips/xfails. Default-set Phase I rules already have chemical coverage via correspondence/parity; added contract tests for `ensure_forest` identity / `get_forest` read-only, Addition↔PatternInfo link, and FutureSite `AtomRef` nesting. Left FilterSites as `dict[str, object]`.

## 2026-09-19

- Expanded ruff and pyright to every Python module under `src/xenosite/refactor_poc/` and every test under `tests/refactor_poc/`. Still excluding `test.ipynb`, `src/xenosite/forest/`, and `examples/`. Cleared to 0 ruff findings and 0 pyright errors. `ReactionRule.__iter__` / `RuleSet.__iter__` return `Iterator`; `rules.RuleSet` keeps a TYPE_CHECKING re-export beside the runtime `__getattr__`.

## 2026-09-19

- Typing pass on refactor_poc: `mol = ensure_forest(mol)` installs in place and returns the same object as `ForestMol`; `get_forest` only reads a `ForestMol` and returns its `Forest`.
- Site is indices only (`int | tuple[int, ...] | frozenset[int] | frozenset[frozenset[int]]`). FutureSite mirrors that nesting; the top-level single deferred atom is `AtomRef` only (not bare int). Nested leaves are `int | AtomRef`. `AnySite = Site | FutureSite`. `_flat_ints` takes `Site` only.
- Focused `pyright: ignore[reportIncompatibleVariableOverride]` kept on `_forest` overrides in `rdkit_api.py`. Pyright on the poc modules listed in `pyproject.toml`: 0 errors.

## 2026-09-19

- PatternInfo now carries a distinguishing ``name``. Where a short role name is not ready, the rule numbers the pattern ("1", "2", …) on that object at init; that numbering is a disfavored interim, not the final names.
- RuleSet always appends itself on the product chain. Its reported label is the initialization name; an unnamed set (including ``name=""``) stays on the chain and contributes no label. The addition links to the PatternInfo that fired; the pattern is not a chain entry.

## 2026-09-19

- RuleSet was flattening nested sets, so a product chain kept only the outer set and the leaf. It now calls metabolize on each contained rule and appends itself on that product. PhaseOneRS hydroxylation of ethane is Hydroxylation, then StableOxygenation (`SO`), then PhaseOneRS.

## 2026-09-19

- Hydrogenation is a ResonancePairRule. Endpoint `[#6;!$(*#[#6]):1]` with edit `keep`. The path flip writes 2-butene (`CC=CC`) from `C=CC=C`. The even path that wrote 1,2-butadiene is not a pair, so `C=C=CC` stays only on the old side. Ethene and ethyne stay `CC` and `C=C`.
- Ethenediol pair draft `O=[CH2][CH2]=O` failed `SANITIZE_PROPERTIES`. Clearing explicit H and sanitizing again yields glyoxal `O=CC=O`, which round-trips. Ethanol dehydrogenation is still the one-bond SMARTS (`CC=O`, `C=CO`).

## 2026-09-19

- Ruff is a dev dependency (`ruff>=0.16.8`). Select is E, W, F, and I; line length 100; target Python 3.10. Bugbear stays off: B905 wants `zip(..., strict=)`, and `strict=True` can change search. Forest, `phaseone.py`, and `test.ipynb` are excluded until the phase-one class rewrite. `_BOND` keys `1.0`, `2.0`, and `3.0` hashed the same as `1`, `2`, and `3` and stored the same bond types, so dropping the duplicates does not change `GetBondTypeAsDouble()` lookup.

## 2026-09-19

- Moving formal charge, and hydrogen on a neutral carbon, with a kekulé bond flip repairs 1,4-dinitrobenzene: the second nitro parses as `O=Nc1ccc([N+](=O)[O-])cc1` plus `O`. Pyridazine still yields `C=CC=CN`. Both Kekulé forms parse and neither nitrogen is charged, so the dropped atom is the reaction, not a charge left on the wrong atom.

## 2026-09-19

- Forest states are classes on the local RDKit mol. `NoTracingMol` means `atom_trace` is absent and does not say whether `_forest` exists. `NoForestMol` is a `NoTracingMol` with `_forest` none, because a missing forest has no trace. `ForestNoTracingMol` is both `NoTracingMol` and `ForestMol`. `ForestTracingMol` is a `ForestMol` only. Constructors that do not copy `_forest` return `NoForestMol`. `ensure_forest` and `ensure_tracing` install in place and do not reset depth.

## 2026-09-19

- Glutathionation calls ConjugationRule. Styrene oxide with `as_star` false matches the old library as `NC(CCC(=O)NC(CSC(CO)c1ccccc1)C(=O)NCC(=O)O)C(=O)O`. A filter skips the carbon-chlorine site on epichlorohydrin by reading `partner`. Ring carbons are separate patterns so a second opening is not one embedding of the three ring atoms.

## 2026-09-19

- Glucuronidation calls ConjugationRule. Phenol with `as_star` false matches the old library as `O=C(O)C1OC(Oc2ccccc2)C(O)C(O)C1O`. A filter skips the benzylic alcohol on 4-hydroxybenzyl alcohol by reading `partner_h`. The old acid SMARTS rewrites `=N`, `=P`, and `=S` as oxygen; that product is not that group's glucuronide, recorded in `DIVERGENCES.md`. The carbonyl pattern matches `=[#8]` on an OH or an anion, so a carboxylate still forms the acyl glucuronide.

## 2026-09-19

- RDKit tag `Release_2026_03_5` (commit `de8add1e32ff6d3c4e4e406f64b703b662dff1d6`), the installed 2026.03.5. `SubstructMatch` (`Code/GraphMol/Substruct/SubstructMatch.cpp`) calls `boost::vf2_all` (`Code/GraphMol/Substruct/vf2.hpp`), which builds `VF2SubState` with `sortNodes=false`. `NextPair` starts at query atom 0 and then walks every molecule atom. `IsFeasiblePair` calls the atom query (`AtomLabelFunctor` → `atomCompat` → `QueryAtom::Match`) and returns false before unmapped neighbors are walked. `AtomIsotope` is `makeAtomIsotopeQuery` / `queryAtomIsotope` (`Code/GraphMol/QueryOps.cpp`, `QueryOps.h`): `getIsotope()` compared in `EqualityQuery::Match` (`Code/Query/EqualityQuery.h`); `AndQuery::Match` stops at the first failing child. `RunReactants` does not have a second matcher: `getReactantMatchesToTemplate` (`Code/GraphMol/ChemReactions/ReactionRunner.cpp`) calls that same `SubstructMatch` with `uniquify=false`.
- C400 and C8000 chains, one `GetSubstructMatches`, query built once. `[8001#6]` costs the same with the isotope on atom 0, on the last atom, or on no atom (C8000: 0.11 ms). The cost grows with the atom count. There is no isotope index, and a hit does not end the scan. Isotope atom first (`[8001#6]-[#6]`) is 0.11 ms on C8000; a plain carbon first (`[#6]-[8001#6]`) is 0.44 ms. When every carbon is isotope 8001 those two orders are the same (0.13 ms), and `[8001#6]` matches plain `[#6]` (0.05 vs 0.05 ms, 400 hits).
- C400, N isotope-limited `[8001#6]` matches against one `[#6]` match plus a Python index filter. Crossover is between N=1 (0.009 vs 0.055 ms) and N=8 (0.066 vs 0.056 ms). Reparsing the SMARTS each call adds about 0.002 ms. `RunReactants` of `[#6:1]>>[*:1]` is 259 ms (400 products); 32 limited runs with the reaction object reused are 10.5 ms, and rebuilding the SMARTS string each time is the same within noise. No reaction crossover through N=32.
- Already known, one line each. On small molecules, tagged matches cost more than one full match past a small N. Dealkylation needs the map assignment, not the site set. Acetylation's isotope sticks to nitrogen only.
- Next: write the isotope-bearing atom first if the limited match stays. Do not expect RDKit to jump to the stamped atom. Precomputing the SMARTS string is not the match cost.

## 2026-09-19

- Sulfation calls ConjugationRule. Phenol with `as_star` false matches the old library as `O=S(=O)(O)Oc1ccccc1`. A filter skips the benzylic alcohol on 4-hydroxybenzyl alcohol by reading `partner_h`. The old arene-oxide SMARTS writes a methyl sulfone and drops the epoxide oxygen; that product is not a sulfate, recorded in `DIVERGENCES.md`.

- Acetylation calls ConjugationRule. The acetyl SMARTS stay on that class. Phenol with `as_star` false matches the old library as `CC(=O)Oc1ccccc1`. A filter skips the nitrogen on 4-aminophenol by reading `symbol`. On 2-mercaptoethanol both matches still emit the thioacetate; that parent reaction was not changed.

- ConjugationRule collapses the acetyl adduct to a star by default. Those SMARTS stay on this class until Acetylation, Sulfation, Glucuronidation, and Glutathionation are their own classes. A filter reads `symbol`. Star labels are not a mode on every rule.

- Kekulé parents: match once with `=,:`; one parent per assignment of the conjugated system that contains the hit; other systems stay aromatic. Helpers take a dict and do not touch `_forest`; the rule stores that dict.

- NDealkylation is the nitrogen rows of dealkylation. A methyl is `leave_count` 1 because that carbon is the whole leaving piece; a larger alkyl stays `None`. `breaks_ring` is set from the cleaved bond whenever the effect already cleaves, including Dealkylation, not from a rule-name branch. Trimethylamine fragment SMILES match the old library. The macrocycle and both dialdehyde searches find their targets.

- Easy pairs against the old library, canonical fragment SMILES. Ethane hydroxylation, anisole dealkylation, and quinone formation on benzene and on phenol match. Butylbenzene hydroxylation matches except three carbonyls from the old `[#6h2:1]>>[*:1]=O` pattern; recorded in `src/xenosite/refactor_poc/DIVERGENCES.md`. Aromatic SMARTS run on every kekulé form. A repeated site is the same map roles and the same incident bond orders, so ortho and para both emit. `bfs` / `dfs` enumerate one ruleset with the filters the set already forwards, a depth cap, and no atom diff.

- Local RDKit types for the calls the poc uses, signatures taken from Boost error text, and `rdkitutil` is under pyright. Not a full stub package. No `rdkit-stubs` dependency.

- `PatternInfo.site_map` is one atom-map number or a tuple of them. Pair rules already passed `(1, 2)`. `Structure.mcs_matches` and `Structure.mcs_targets` stay `dict[str, McsResult]`, keyed by the target canonical SMILES. `total` stays false. Pyright on records, rules, find_path, and rulesets: 0 errors. Search and reaction chemistry unchanged.

- Project rule `.cursor/rules/data-not-branches.mdc`: generic algorithms read data. Methide is offered by rule data, limited to one site by `PatternInfo`, and excluded by `filter_sites`. Exceptions need a written reason. A less-likely mark is future filter data. That schema is not chosen, so no field was added.

- Atom trace additions are normalized. A new atom's `added_by` is `R1`, `R2`, ... The site, rule hierarchy (including the ruleset), effect, name, phase1, and depth live once on `atom_trace["additions"]`. `formula` is the current heavy-atom counts plus H and formal charge. `delta_formula[id]` is the change at that transform. `metabolize` works on a copy, so the caller's mol is unchanged. Search changes consult `src/xenosite/refactor_poc/HEURISTICS.md` first.
- QuinoneFormation is not a step. The plan is a Hydroxylation for each end that still needs oxygen, then one Dehydrogenation. An end that already has O/N uses that neighbor; no `end_maps` field. Bare hydroxylation is refused when every oxygen-needing site atom needs a carbonyl and loses aromaticity. An alkyl (methide) partner runs only when that exocyclic C-C bond is higher in the target.
- Cleavage filter keeps the bridge between a mapped atom and a leaving atom. Bonds inside the leaving group are not sites. Terbinafine → `CC(C)(C)C#CC=CC=O`: new `mol_edits=3` (old `site_applies=0`, ceiling 40). Acetate → catechol: `mol_edits=12` (old `site_applies=0`, ceiling 80). Phenol/naphthalene → methane and benzene → CF4 miss under 40 edits. Benzene → `[Fe]` is `mol_edits=0`. NDealkylation hard pairs (macrocycle 120, dialdehyde 80/35) skipped; that class is not in the catalog yet.
- Adversarial read of `find_path.py` only: no second search and no `guided_path` ranker to cut. The loose "site touches any leaving atom" check was the extra edits; the bridge test replaced it.
- `RuleSet` lives in `refactor_poc/rulesets.py` and is a rule: running the set runs each child, and `filter_rules` / `filter_sites` still see that child's pattern. `PhaseOne` lists the reaction classes that exist in `rules.py`, including `QuinoneFormation`. Forest-only rules (NDealkylation, Tautomerization, AzoSplitting, BenzodioxoleReduction, NitroaromaticReduction, ThiopheneSulfurOxidation, conjugations) are omitted, not invented. `from rules import RuleSet` still resolves.
- refactor_poc `find_path` calls one `RuleSet`. The set is a rule: running it runs each child, and `filter_rules` / `filter_sites` still see the child pattern. Not a bare list of classes. `phaseone.py` stays dropped.
- refactor_poc `find_path`: one MCS (element match, `CompareAny` bonds) plus the local change at each atom. Not `PathContext`. `mol_edits` is an accepted-site `RunReactants` or a pair-rule kekulé overlay; a refused site is `sites_skipped`. `billed = mol_edits + nodes`. When the target is smaller, non-cleaving patterns fail `filter_rules`. `Deps` edges come from `added_by`, not SMARTS replay.
- refactor_poc rules: `ResonanceRule` / `ResonancePairRule` live in `rules.py` (no AtomTracker). Hydroxylation, then Dehydrogenation, then QuinoneFormation share them. `PatternInfo.possibilities` lists each SMARTS branch; `resolve_effect` narrows that from the matched atoms; `span` is a bare value when every branch agrees and a tuple when it does not.
- Hydroxylation `[#6h2:1]>>[*:1]O` adds OH and removes one H. The `>>[*:1]=O` / `removes: HH` form was wrong.

## 2026-09-18

- AtomTracker.tags → ``mol._forest["atom_trace"]``; stable ``atomLabel`` per heavy atom (stamped once); ``copy_mol`` / ``carry_forest`` preserve forest across RDKit copies; ``install_product_forest`` records ``atom_refs`` from metabolize (shared with Step.apply). Resonance not shared onto reaction products. Profile (`tests/test_atom_trace_profile.py -n0 -s`): **before** wall=0.638s literal_eval_cum=0.139s deepcopy_cum=0.131s tags_cum=0.280s → **after** wall=0.358s literal_eval_cum=0 deepcopy_cum=0 tags_cum=0.006s (benzene PhaseOneRS depth-1→2, products=165).
- AtomTracker.tags before baseline (`uv run pytest tests/test_atom_trace_profile.py -q -n0 -s`): benzene PhaseOneRS depth-1→depth-2 expand — wall=0.638s, products=165, parents=7, `literal_eval_cum`=0.139s, `deepcopy_cum`=0.131s, `tags_cum`=0.280s (cProfile cumtime). Characterization suite locked in `tests/test_atom_trace_char.py` (green on string-prop storage).

## 2026-09-17

- Guided search: a non-molecule product, an unknown `enumerate_for_path` kind, an unknown `search` mode, or a rule that raises is a broken pipeline and now raises. Depth, budget, `could_help`, `child_may_reach`, terminal/star, and seen-duplicates stay ordinary misses. A linearization whose `apply` does not fire is also a miss (caught), not an assertion.

- Dropped unused library helpers (kept only where production calls them): `_as_atoms`, `StepPlan.apply_linearization`, `Deps.canonical_precedes` / `.transitive_closure_masks`, `same_possible_orders`, `deps_from_plan`, `_count_topo_sorts` wrapper, `PathContext.must_create` / `must_remove_or_transform`. Fuzz still uses local `_deps_from_plan` in tests.
- Deps stores transitive reduction on construct (stable precedes/JSON/keys); small graphs so always canonicalize. Closure equality is O(M+N) if needed; `same_linearizations` still aligns Step identity first.
- Chem-disagree miss (`CCCCCCCCCCCC1CCC2OCCCC2C1`→`O=CCCCO`): not a hard-drop — priority poison. CompareAny ring chem-diff (+neighbors) sorted first; winning path is alkyl peels `{10,11}` then `{6,7}`. Fix: chem-disagree stays frontier/never-drop only; priority uses embedding disagreement + geometric frontier only.
- phenol→BQ `sanitize_dropped=35` autopsy: 9 QuinoneFormation (aromatic drafts cleaned after RunReactants) + 26 Dehydrogenation (plan apply ran full resonance fanout then filtered by site). Short circuits: (1) QF `enumerate_for_path` emits phase1 plans from pair matches without materializing/clean; skip plans whose Hydroxylation count exceeds `oxygen_deficit`; (2) `metabolites_from_sites` → `include_sites` into DH/QF/Hydrogenation *before* `apply_modifications`/`clean`; (3) sanitize keep-Hs skipped when strict valence fails, still try `reset_hs=True` (hard-reject broke QF). Result: phenol→BQ **35→2** sanitize drops, billed 8→2; APAP→NAPQI still 1 hop.
- Large-mol A/B (full vs no-chem-disagree vs may_reach-only; |R|≈20–28): theory that chem-diff/frontier helps *a lot* on large mols is only partly true. **Priority+safe-drop** vs may_reach alone: crowded dialdehyde 21 vs 38 billed (0.55×); C16/C22 alkyl-THF ~0.90×; many ND cases already 2 billed (tied). **Chem-disagree alone** rarely changed billed vs geometry-only; once hurt under budget (`C12…bicycle→O=CCCCO`: full miss / no_chem hit). Added `PathSearchCounters.sanitize_dropped` (via `clean` / path scope); phenol→BQ guided run recorded 35 — cleavage benches often 0.
- Failure: CompareAny MCS maps ring atoms onto open-chain T; ring bonds looked deep-interior and were hard-dropped from Required (`CC1CCC2OCCCC2C1`→`O=CCCCC1CCCC1`: chem-diff sites absent from `sites_toward`). Simple THF→`O=CCCCO` has |T|≥|R| so `sites_toward` never runs — not the live miss. Fix: `mcs_chem_disagree_bonds` (mapped bond missing or order/aromaticity differs vs T) counts as frontier; never deep-interior. N-ring-open was already saved by endocyclic special-case; O-ring-open was not.
- Hard-dropping non-frontier cleavage sites from Required is gated, not blanket: (1) deep interior of shared full-size MCS core (core ≥ max(3, |T|//2)); (2) deep exterior only when `|largest emb| + 2 ≥ |T|` (single-cleave size). Else prioritize disagreement → frontier → rest. Bis-benzyl dual embeddings → empty shared core → no interior drops (N–CH2 still first via priority). `unreachable_new_elements` fast-fails metals; Dealk/Hydrolysis declare `elements_may_add={O}` so aldehyde O is not false-impossible. Suite: 44 pass on prune/hard/maybe/guided/classic-counters.
- Cleavage MCS prune hard cases (`tests/test_cleavage_site_prune.py`): (1) union-of-embeddings as one cons → 0 sites on bis-benzyl→PhCHO; need OR per placement. (2) smaller secondary-only (Ph remainder) → miss naphthyl–benzyl→naphthaldehyde; keep full-size + smaller. (3) N-dealk ring-open with benzyl-only MCS → ring∩cons empty; never prune endocyclic N–C. Non-N ring-open only consults full-size embeddings (else TBA reopens naphthalene fanout).
- Guided budget now charges `billed = linearizations_applied + site_applies` (not `rule_expansions`). One Dealk `enumerate_for_path` can hide ~40 plan applies; those now consume the cap. Classic `RuleSet.find_path` takes the same `PathSearchCounters` (dict still mirrors full field set): 1 `site_applies` per frontier `metabolize` preserves prior classic caps; `sites_considered` / `mol_edits` filled for guided-comparable benches.
- `PathOutcome` = Required `plan` + `MaybeFilter` (composition, not StepPlan subclass); `find_path` yields it (5-tuple unpack still works). `CleavageSide.opens` records ring-open sites; bifurcation bags span open+cut (macrocycle → amino-ketone). Distinct Required step-site routes to the same T are kept (seen/dedup by steps key). Defaults: `depth=None`, `max_expansions=200`. Impossible-T tests must abort with empty frontier and `budget_exhausted=False` (miss-after-budget is not a pass).
- Maybe semantics: declarative `CleavageSide` bags only — not traced prefixes. Preceding N-dealkylation on the same nitrogen (different site, e.g. demethyl before TBA-forming dealk) passes `outcome.allows`. Phase-B prefix search removed.
- Guided size-gate: when T is smaller, expand **cleavage peers** only — every `is_cleavage()` rule plus `cleave_alone()` opt-ins (OxidativeDehalogenation, AzoSplitting). TBA under PhaseOneQF ≈ 8 expansions among peers.
- Guided default ruleset → `PhaseOneQF` (QF then PhaseOneRS). Callers need not pick QF vs Phase I: benzene→BQ stays 1 expansion; TBF→TBA still via Dealkylation at Phase I cost (QF inert). Named ruleset registered as `PhaseOneQF`.
- Classic `find_path` gained optional `max_expansions` + mutable `counters` (`rule_expansions` = frontier `metabolize` calls, `budget_exhausted`). Bench vs guided (cap 50): TBF→TBA needs E-stereo for classic; QF+PhaseOne = Phase I Dealkylation for TBA (QF inert); QF alone misses TBA. Tight caps: Phase I benzene→BQ needs ~39 classic expansions (budget 20 exhausts); toluene→o-QM misses on default DH (no methide) even with budget 50.
- Cleavage + formula heuristics: ADD_O on a large parent vs small fragment target is incompatible, so hydroxylation is pruned until after the split. Guided expands cleavage peers alone when T is smaller; `and_cleave_plan` turns cleave-first walks into `And(cleave…, Seq(rest))` so ops before/after the cleave remain plan members. Tests: long phenyl ester → catechol; anisole acetate; biphenyl ester half (`tests/test_guided_path_cleavage.py`).
- Formula hints on SMARTS options are **polymorphic**: `FormulaHint`, `FormulaAny` / `(CLEAVE, ADD_O)`, `FormulaMatch(possible, resolve)`, or `callable(mol, match)`. Pre-match skips use the possibility set (sound); post-match can refine. Dealkylation uses `CLEAVE_OR_ADD_O`.
- Per-SMARTS ``FormulaHint`` on Phase I rules (`formula_effects`); guided `enumerate_for_path` passes `toward_target` to skip incompatible SMARTS before `RunReactants`. Fixed wrong rule-level gates: EpoxideOpening was cleavage-only (now NEUTRAL|ADD_O); OxygenReduction/Dehydration align with hints. Covering tests `tests/test_formula_smarts_hints.py` (55).
- Guided fuzz: prefer ``max_expansions`` (+ ``counters.budget_exhausted``) over wall-clock timeouts so misses surface expansion cost for heuristics. `StepPlan.contains` / `by="rule"` for recipe∈plan without full linearization blow-up. Toluene QF → o-QM recovered via methide pathway.
- Bench classic `find_path` vs `find_path` (median of 3, first path; **never classic BFS at depth ≥ 3**):
  - EtOH→acetaldehyde (d=2): both find; ~same expansions (1); classic faster (~0.5× time — guided overhead on tiny graphs).
  - APAP→NAPQI: guided d=3 vs classic d=1 (DH-only); both find; expansions 1; classic faster.
  - Terbinafine→TBF-A (`CC(C)(C)C#CC=CC=O`, ND, d=2): guided finds (1 expansion); classic BFS reports no hit — target stereo-unspecified vs Forest E emit (`can_smi` strict). Guided strips stereo. Expansions 19→1 when counting classic’s full miss search.
  - Benzene→BQ PhaseOneRS: formula-aware ``could_help`` (OH only while O-deficient; DH when formulas match; skip epoxidation on quinoid targets; cleavage only if T smaller) cut expansions **379 → 15**. QF still 1. Count expansions only when ``could_help`` passes. Rebuild ``PathContext`` every frontier mol.
- Guided `find_path`: `PathContext` (FindMCS `BondCompare.CompareAny`), rule hooks, batched `metabolize` + `include_sites`/`exclude_sites`, phase1 substitution, Required-then-Maybe, stereo-agnostic endpoint match. Counters for optimization claims. Classic `find_path` unchanged.
- StepPlan algebra: `StepPlan` = Seq; `And`/`Or` subclasses. One plan per `phase1_steps` (QF uses top-level `Or`). Linearizations still lazy over nested expr.
- TBF-A gold = `CC(C)(C)C#CC=CC=O` (6,6-dimethylhept-2-en-4-ynal); not the C10 `…CC=O` homolog. Stereo not required for match.
- Forest `Dehydrogenation` missed quinoid products because query SMARTS required `#6h`. Added `#6H0`–`[#8H]` / `#6H0`–NH / `#6H0`–`[#7D3]` (+`addPlus1`) query SMARTS. Hydroquinone→benzoquinone and APAP→NAPQI now emit from Forest DH; quinone StepPlan full `apply` works for those cases. Full suite: 916 passed, 15 xfailed (no new failures).
- Public stamp API: `StepPlan.try_from_mol` / `from_mol` / `attach_to_mol` (dropped `has_on_mol`). Docs cover rule ask, metabolize/ruleset stamp, resolve + `linearizations().apply` / `drop_last`. Ready to ship 0.4.0.
- `AtomRef` + `Step`/`Linearization.apply`: deferred sites (new O after hydroxylation) resolve via `mol._forest["atom_refs"]` (private schema documented next to `_forest_state`); not AtomTracker tags. Quinone DH ends use `added_by` refs. `toward=` / `drop_last=` keep fragments for omitted suffix sites. Creation index stays in `Step.apply` for now; optional rule-side recording deferred (TODO).
- Phase1-equivalent steps (`Step` / `StepPlan`): uniform `phase1_steps(mol, site)`; Phase I + NDealkylation degenerate singletons; QuinoneFormation SMARTS→prep layers + final DH; `attach_phase1_steps` stamps via StepPlan API. `RenumberAtoms` during tag/align drops mol props — copy props across renumber so stamps survive.
- Fuzz uses library apply only. Product identity vs quinone checked via prep+`QuinoneFormation` when Forest DH apply is empty (now often non-empty after quinoid DH SMARTS).
- Mol-scoped resonance cache on `mol._forest["resonance"]` (lazy pull-through joined forms per conjugated/aromatic mode; shared `bfs_all_pairs`). Private `_resonance_cache_disabled()` for parity tests. `EditMol.standardize` propagates `_forest`.
- Full product parity cache on vs off: APAP 69, naph_styryl 228, multi_conj 400 (matched). Form counts matched.
- Full.metabolites call counts (on → off): `_resfrags` 2→7, `join_fragments` APAP 4→14 / naph 10→35 / multi 14→49, `bfs` 1→4.
- Wall (median): Full multi_conj ~1.72s on vs ~1.77s off (~1.02×); resonance-only ruleset ~1.06–1.09×. Repeated `resonance_structures`×5 on same mol ~4.1–4.2×. Full wall dominated by pair-path/SMARTS work after forms exist.
- Issue #3 reports `calcImplicitValence` crash on PhaseOneRS BFS; reproduces on forest 0.1.0 only (already fixed in 0.2.0+).
- Full BFS depth=2 failed with RDKit `Range Error` after acetylation→dehydrogenation: `_remove_dummy_atoms` deleted conjugation `*` adducts. Fix: mark FragmentOnBonds dummies and only remove those.
- BFS `expand_star_conjugates=False` by default: star conjugates are terminal products; opt in to expand them further.
- Added `dfs` pathway search with optional ``shuffle_rng``; fuzz randomly samples Full DFS two-step pathways and persists the Hypothesis example DB (``.hypothesis``, cached on CI).

## 2026-09-13

- Topo emission: `metabolize` topsite keyed on suffix'd rule labels and site ranks only, so UI-equivalent rows leaked and some distinct product SMILES were dropped. Fix: identity = normalized pathway + sorted topo-rank multiset + product SMILES (same as xenosite UI). SMARTS-loop early skip left as TODO.
- N-dealkylation structure bug: `UO.Dealkylation` SMARTS also cleave O/S/C; UI filtered by nitrogen site. Fix: `NDealkylation(Dealkylation)` post-filters to N-containing sites, ruleset `ND`; short-circuit when reactant has no N. Leave `Dealkylation` on `UO` unchanged.

## 2026-09-09

- Framing: reactivity conjugation is high-sensitivity / low-specificity structure enumeration (if a site were positive, what adduct forms) — not a likelihood model.
- Added aldehyde thiohemiacetal; aziridine; sulfonate ester (`C–OSO2`); isocyanate/isothiocyanate. Negatives: azetidine/pyrrolidine/pyridine vs aziridine; sulfonamide/sulfone/sulfonic acid vs sulfonate; nitrile/amide/CO2/carbodiimide vs isocyanate.
- Earlier: Michael + F/Cl/Br/I for quinones / BnBr / BnI; flipped β-enone skip test. Legacy SMARTS were only epoxide, C–Cl, thiol, terminal `CH2=`.
