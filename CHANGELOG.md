# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This project uses [towncrier](https://towncrier.readthedocs.io/) for the next
release; fragments live in [`changelog.d/`](changelog.d/).

<!-- towncrier release notes start -->

## [0.5.2](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.5.2) - 2026-09-18

### Fixed

- SMARTS `metabolites` kekulizes the caller's mol again and drops its resonance cache. `standardize` no longer keeps a kekule prototype.


## [0.5.1](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.5.1) - 2026-09-18

### Added

- Conjugation path policy: one unlabeled `*` peer via `is_redundant`; star products terminal.
- Optional `include_sites` / `exclude_sites` on `metabolize` (default: no filter).
- Optional `max_expansions` / `counters` on classic `RuleSet.find_path` (BFS/DFS): accepts `PathSearchCounters` (or a dict mirror); caps billed work so classic and guided counters are comparable.
- Optional extra `predict` (`xenosite-forest[predict]` → `xenosite-predict`); not installed on CI.
- Package `find_path` (alias `find_path_guided`): MCS-guided search; default `PhaseOneQF`, `depth=None`, `max_expansions=200`; yields `PathOutcome` (Required `plan` + cleavage-side `maybe`). Unstable API.
- Per-SMARTS `formula_hint` / `FormulaHint` on reaction options; `metabolize(..., toward_target=)` skips incompatible SMARTS. Dehydrogenation methide is opt-in (not default Phase I).
- `Deps` (`StepPlan` specialization): flat steps + precedes; linearizations are topological sorts. Guided emission uses `Deps` (not Kahn→And/Seq).
- `NDealkylation.is_redundant` when `Dealkylation` peers are present; duplicate `Epoxidation` instances dropped.
- `StepPlan.n_linearizations()`: count total orders without enumerating (`Deps` uses subset DP, `n≤20`).

### Changed

- Attached plan JSON is an expression tree (`op`: `seq`/`and`/`or`/`step`/`deps`); legacy flat `steps`+`precedes` still loads.
- Cleavage site prune / priority under guided search (MCS frontier, safe-drop interior/exterior, chem-disagree as frontier-only). Early abort when T introduces unreachable elements.
- Guided `max_expansions` caps `PathSearchCounters.billed()` (`linearizations_applied` + `site_applies`), not `rule_expansions`.
- Package `__all__` trimmed: emission helpers (`and_cleave_plan`, `CleavageSide`, `MaybeFilter`, PathContext utils, `@unstable`) are submodule-only. Docs keep `find_path` short.
- QuinoneFormation guided enum can emit phase1 plans without materializing/clean; `include_sites` restricts resonance fanout before `clean`.
- Store atom traces on mol._forest with stable `_forestLabel` (not CX `atomLabel`); metabolize records atom_refs; copy_mol preserves forest.
- `@unstable` / `UnstableWarning` on public guided surface: `find_path`, `PathOutcome`, `PathSearchCounters`, `Deps`, `StepPlan.n_linearizations`.
- `AtomRef.added_by` is `(rule, site)` with frame `depth` for those site idxs; origin refs also carry `depth` (created atoms lack depth-0); `atom_trace` keeps live `records` plus a `removed` event list (cleavage siblings are dropped, not recorded).
- `QuinoneFormation.phase1_steps` / Phase I emitters return a single `StepPlan`. Prefer `plan.branches()` / `plan.linearizations()`.


## [0.4.0] - 2026-09-17

### Added

- Public `AtomRef` / `Step` / `StepPlan` / `Linearization`: reactant-stable sites (origin or created-by), partial orders, `Step.apply` / `Linearization.apply` with `mol._forest["atom_refs"]` creation index and opt-in `resolve_site`.
- `Linearization.apply(..., toward=, drop_last=)` for fragment retention and prep-only replay.
- `ReactionRule.phase1_steps(mol, site)`: Phase I and `NDealkylation` return a degenerate singleton plan; `QuinoneFormation` returns multi-step prep + final dehydrogenation (`AtomRef` ends for new O); other rules raise `NotImplementedError`.
- Opt-in `attach_phase1_steps=True` on `metabolize` / rule `metabolites` / rulesets. Attach and read plans with `StepPlan.attach_to_mol` / `from_mol` / `try_from_mol` (storage is private; do not read mol props by name).
- `Dehydrogenation` query SMARTS for quinoid ends (`#6H0` bonded to OH / NH / tertiary N) so Forest DH can form hydroquinone→benzoquinone and APAP→NAPQI; quinone `StepPlan` full Forest apply works for those cases.
- Hypothesis fuzz via library `Linearization.apply` (prep agreement; full equality when all steps fire).

## [0.3.1] - 2026-09-17

### Changed

- Resonance rules share lazily cached joined resonance forms and `bfs_all_pairs` on `mol._forest` (conjugated vs aromatic modes), avoiding repeated `ResonanceMolSupplier` work within a ruleset. Private `_resonance_cache_disabled()` for parity tests.

## [0.3.0] - 2026-09-17

### Fixed

- Resonance reassembly no longer strips conjugation ``*`` adducts when removing FragmentOnBonds dummies, so `bfs(..., ruleset="Full", depth=2, expand_star_conjugates=True)` no longer raises RDKit `Range Error` after acetylation→dehydrogenation.

### Added

- `dfs(...)` / `find_path(..., search="dfs")` for depth-first pathway search; CLI `--search dfs`.
- Optional `max_paths` on `bfs` / `dfs` (default unlimited) and `shuffle_rng` on both to randomize metabolite order.
- `expand_star_conjugates=False` (default) on BFS/DFS / `find_path`: star (`*`) conjugate adducts are emitted but not metabolized further. Opt in with `expand_star_conjugates=True` or CLI `--expand-star-conjugates`.
- Hypothesis fuzz: randomly sample Full DFS two-step pathways (``shuffle_rng``); persist ``.hypothesis`` example DB on CI.
- Regression lock for GitHub issue #3 parent SMILES (PhaseOneRS valence crash was already fixed since 0.2.0 / 0.2.2).

## [0.2.8] - 2026-09-13

### Fixed

- `metabolize` topological dedup now matches the XenoSite UI identity key (normalized pathway + sorted topological ranks + product SMILES), so equivalent sites collapse without dropping distinct products.

## [0.2.7] - 2026-09-13

### Added

- `NDealkylation` rule and `ND` ruleset: Dealkylation products filtered to nitrogen-containing formation sites (for the N-dealkylation XenoSite model). Short-circuits when the molecule has no nitrogen. Broad `Dealkylation` on `UO` is unchanged.

## [0.2.6] - 2026-09-09

### Fixed

- `Glutathionation` reactivity conjugation (high sensitivity / low specificity): Michael acceptors, C–Br/I/F, aldehyde thiohemiacetals, aziridines, sulfonate esters, and isocyanates/isothiocyanates, in addition to epoxide, C–Cl, thiol, and terminal alkene.

## [0.2.5] - 2026-09-09

### Fixed

- `Glucuronidation` site of metabolism is a single atom: the oxygen that receives GlcA (acid and phenol/alcohol SMARTS), aligned with high `ugt` atom scores.

### Added

- Conjugation regressions that require a single-atom SOM for every UGT / GSH / Protein / DNA / Cyanide SMARTS probe, plus UGT oxygen site checks.

## [0.2.4] - 2026-09-08

### Fixed

- Doctest collection via importlib so `xenosite` stays a PEP 420 namespace package.
- RDKit enumeration drift regressions vs forest 0.2.3 / py2 baselines.
- `clean()` no longer drops whole product sets when one fragment is invalid; invalid metabolites log at DEBUG.

## [0.2.3] - 2026-09-03

### Added

- Conjugation star adducts by default; optional `star_label` / `as_star=False`.
- `Glutathionation(include_thiol=False)` and `load_ruleset("GlutathionationNoThiol")` for DNA / cyanide heads.
- Star-only labels `Protein`, `DNA`, and `Cyanide`.

## [0.2.2] - 2026-08-29

### Fixed

- RDKit valence / property-cache handling around `RunReactants` and SMILES (refresh before react; skip unsanitizable reactants).

[0.4.0]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.4.0
[0.3.1]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.3.1
[0.3.0]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.3.0
[0.2.8]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.8
[0.2.7]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.7
[0.2.6]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.6
[0.2.5]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.5
[0.2.4]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.4
[0.2.3]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.3
[0.2.2]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.2
