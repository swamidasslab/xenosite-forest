# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

[0.3.1]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.3.1
[0.3.0]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.3.0
[0.2.8]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.8
[0.2.7]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.7
[0.2.6]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.6
[0.2.5]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.5
[0.2.4]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.4
[0.2.3]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.3
[0.2.2]: https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.2.2
