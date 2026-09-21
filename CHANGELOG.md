# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This project uses [towncrier](https://towncrier.readthedocs.io/) for the next
release; fragments live in [`changelog.d/`](changelog.d/).

<!-- towncrier release notes start -->

## [0.7.1](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.7.1) - 2026-09-21

### Added

- Declare and CI-test Python 3.14 support (requires-python already >=3.11). ([#20](https://github.com/swamidasslab/xenosite-forest/issues/20))

### Changed

- Post-3.10 cleanup: pin pyright ``pythonVersion = "3.11"``; enable ruff ``UP``
  (pyupgrade) with ``UP031`` ignored; ``match`` for ``Site`` shape checks in forest
  rules. Audit found no remaining ``typing_extensions`` / ``tomli`` / ``Self`` TypeVar
  shims in live forest. ([#21](https://github.com/swamidasslab/xenosite-forest/issues/21))
- Emitted ``site``: directed_bond tuple→frozenset; under ``canonical_emitted_sites``, lex-orbit remap (not unique-edit / CSMI dedup). ``discovered_site`` is always the raw discovery site. Dropped ``_site_tuple`` and unused internal stamp helpers. ([#22](https://github.com/swamidasslab/xenosite-forest/issues/22))

### Fixed

- CI builds pynauty with portable ``CC=-march=x86-64`` so uv-cached wheels do not SIGILL on GitHub runners (``pdobsan/pynauty#49``). pynauty is not bundled into forest wheels. ([#23](https://github.com/swamidasslab/xenosite-forest/issues/23))


## [0.7.0](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.7.0) - 2026-09-20

### Deprecated

- ``AtomTracker`` is deprecated; prefer ``mol.xf`` / ``mol.xf.tracing``
  (see ``docs/forest/XF.md``). ([#16](https://github.com/swamidasslab/xenosite-forest/issues/16))

### Added

- ``xenosite.forest`` typing is at 100% coverage under the package pyright
  gate (informative for callers). ([#18](https://github.com/swamidasslab/xenosite-forest/issues/18))

### Changed

- Metabolic Forest rewrite in ``xenosite.forest``. Previous forest (0.6.x API)
  is the archive directory ``src/xenosite/_archive_forest/`` on GitHub.
  ``metabolize`` yields ``(list[Mol], info)`` — see ``docs/forest/MIGRATING_0.7.md``. ([#15](https://github.com/swamidasslab/xenosite-forest/issues/15))
- **Dependencies:** ``pynauty`` is now required (no longer optional). Needed
  for the correct unique-edit / pair-orbit nauty backend (forest unique-edit
  authority). RDKit/SMILES pair-orbit is not shipped as a product fallback. ([#17](https://github.com/swamidasslab/xenosite-forest/issues/17))
- **Breaking:** Python **3.11+** is required (3.10 is no longer supported). ([#19](https://github.com/swamidasslab/xenosite-forest/issues/19))


## [Unreleased]

Migration details: [`docs/forest/MIGRATING_0.7.md`](docs/forest/MIGRATING_0.7.md).
Previous forest (0.6.x API): archive directory
[`src/xenosite/_archive_forest/`](src/xenosite/_archive_forest/) on GitHub.

### Changed

- Metabolic Forest lives in ``xenosite.forest``. Design notes under
  ``docs/forest/``.
- **Breaking:** ``metabolize`` yields ``(products, info)`` where
  ``products`` is always a ``list[Mol]`` (one element for non-cleavage;
  siblings for cleavage). Product SMILES via ``product.xf.csmi`` (no
  ``info["csmi"]``). See the migration guide.
- Reaction SMIRKS surfaces use ``Smirks`` / ``SmirksReactionRule`` naming
  (match-only SMARTS stay ``Smarts``).
- Rules declare ``site_kind`` (``atom`` / ``bond`` / ``directed_bond`` /
  ``atom_pair``). Public bond sites stay frozensets; directed orientation
  is on ``discovered_site``.
- Cross-rule same-product overlaps under a RuleSet log at INFO.
  Within-rule unique-edit misses still raise ``SiteDeduplicationWarning``.
- ``find_path`` closers distinguish Hydrogenation (adds H) from
  Dehydrogenation (removes H).
- **Dependencies:** ``pynauty`` is now required (no longer optional). Needed
  for the correct unique-edit / pair-orbit nauty backend (forest unique-edit
  authority). RDKit/SMILES pair-orbit is not shipped as a product fallback.
- **Breaking:** Python **3.11+** is required (3.10 is no longer supported).

### Deprecated

- ``AtomTracker`` — prefer ``mol.xf`` / ``mol.xf.tracing``. Migration tutorial:
  [`docs/forest/XF.md`](docs/forest/XF.md#tutorial-replacing-atomtracker-with-molxf).

### Added

- Opt-in ``canonical_emitted_sites`` on metabolize / path search.
- Pair-orbit unique-edit (``swap_group``) for ResonancePair ends.
- ``mol.xf.tracing`` helpers for AtomTracker migration: ``depths()``,
  ``index_at``, ``added_indices``, ``root_map``.
- ``xenosite.forest`` typing is at 100% coverage under the package
  pyright gate (informative for callers).


## [0.6.1](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.6.1) - 2026-09-19

### Fixed

- ``metabolize`` leaves the caller's mol bonds alone. SMARTS kekulize and
  resonance bond-path search run on a copy. Atom maps and tags are still
  written on the input.


## [0.6.0](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.6.0) - 2026-09-19

### Changed

- ``phase1=True`` reports sites as 0-based atom indexes (``1.h`` / ``2.3``
  labels removed).

### Fixed

- Guided search prioritizes dearomatizing rules and sites when a large
  matched region is aromatic on the reactant and not on the target.


## [0.5.3](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.5.3) - 2026-09-18

### Fixed

- ``AtomRef`` for a created atom follows the ``atom_trace`` label stamped
  ``added_by``, with frame ``depth`` for resolution.


## [0.5.2](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.5.2) - 2026-09-18

### Fixed

- SMARTS ``metabolites`` kekulizes the caller's mol and drops its resonance
  cache. ``standardize`` no longer keeps a kekule prototype.


## [0.5.1](https://github.com/swamidasslab/xenosite-forest/releases/tag/v0.5.1) - 2026-09-18

### Added

- Package ``find_path`` (MCS-guided; default ``PhaseOneQF``); yields
  ``PathOutcome``. Unstable API.
- Optional ``include_sites`` / ``exclude_sites`` on ``metabolize``.
- Optional ``max_expansions`` / ``counters`` on classic ``RuleSet.find_path``.
- Per-SMARTS ``formula_hint``; ``metabolize(..., toward_target=)`` skips
  incompatible SMARTS.
- ``Deps`` step plans; ``StepPlan.n_linearizations()``.
- Conjugation path policy: one unlabeled ``*`` peer; star products terminal.
- Optional ``xenosite-forest[predict]`` extra (not on CI).
- ``NDealkylation.is_redundant`` when ``Dealkylation`` peers are present.

### Changed

- Attached plan JSON is an expression tree; legacy flat ``steps``+``precedes``
  still loads.
- Cleavage site prune / priority under guided search; early abort when the
  target introduces unreachable elements.
- Atom traces on ``mol._forest`` with stable labels; ``copy_mol`` preserves
  forest.
- Guided surface marked ``@unstable`` / ``UnstableWarning``.
- Package ``__all__`` trimmed; emission helpers are submodule-only.


## [0.4.0] - 2026-09-17

### Added

- Public ``AtomRef`` / ``Step`` / ``StepPlan`` / ``Linearization``.
- ``ReactionRule.phase1_steps(mol, site)`` and opt-in
  ``attach_phase1_steps=True``.
- ``Dehydrogenation`` query SMARTS for quinoid ends (hydroquinone→benzoquinone,
  APAP→NAPQI).


## [0.3.1] - 2026-09-17

### Changed

- Resonance rules share lazily cached joined forms and pair paths on
  ``mol._forest``.


## [0.3.0] - 2026-09-17

### Fixed

- Resonance reassembly keeps conjugation ``*`` adducts when removing
  FragmentOnBonds dummies (Full BFS depth≥2 after acetylation→dehydrogenation).

### Added

- ``dfs(...)`` / ``find_path(..., search="dfs")``; CLI ``--search dfs``.
- Optional ``max_paths`` and ``shuffle_rng`` on BFS/DFS.
- ``expand_star_conjugates=False`` by default on BFS/DFS / ``find_path``.


## [0.2.8] - 2026-09-13

### Fixed

- ``metabolize`` topological dedup matches the XenoSite UI identity key
  (pathway + topo ranks + product SMILES).


## [0.2.7] - 2026-09-13

### Added

- ``NDealkylation`` rule and ``ND`` ruleset (N-containing formation sites).


## [0.2.6] - 2026-09-09

### Fixed

- ``Glutathionation`` reactivity conjugation: Michael acceptors, C–Br/I/F,
  aldehyde thiohemiacetals, aziridines, sulfonate esters, and
  isocyanates/isothiocyanates, plus existing epoxide / C–Cl / thiol /
  terminal alkene coverage.


## [0.2.5] - 2026-09-09

### Fixed

- ``Glucuronidation`` site of metabolism is the single oxygen that receives
  GlcA (acid and phenol/alcohol SMARTS).


## [0.2.4] - 2026-09-08

### Fixed

- Doctest collection via importlib (PEP 420 namespace package).
- ``clean()`` keeps valid fragments when one product fragment is invalid;
  invalid metabolites log at DEBUG.


## [0.2.3] - 2026-09-03

### Added

- Conjugation star adducts by default; optional ``star_label`` / ``as_star=False``.
- ``Glutathionation(include_thiol=False)`` and
  ``load_ruleset("GlutathionationNoThiol")``.
- Star-only labels ``Protein``, ``DNA``, and ``Cyanide``.


## [0.2.2] - 2026-08-29

### Fixed

- RDKit valence / property-cache handling around ``RunReactants`` and SMILES
  (refresh before react; skip unsanitizable reactants).

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
