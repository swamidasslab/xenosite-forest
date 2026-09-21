# Lab log

## 2026-09-20

- **Epoxidation uses supplier kekulé mols again.** SMARTS is `[#6:1]=[#6,#7:2]`, not `=,:`. `Epoxidation.metabolites` runs that SMARTS on each `ResonanceMolSupplier(..., KEKULE_ALL)` mol. Bond orders are not laid back onto the aromatic parent, and the cached conjugated-component parent is not used. Unique-edit `seen` is shared across those forms; `incident_orders` stay on the aromatic context. Other `ResonanceRule`s still use the cached parent: walking supplier forms dropped anisole and pyridine ring-open products that parity requires. Six neutral goldens are in `tests/forest/test_epoxidation_golden_products.py`. `pytest tests/forest`: 1626 passed, 1 skipped, 2 xfailed. Quinone neutral golden is still absent; pointing its pair parents at the same supplier mols did not restore it.

- **Archive vs live on the current tree (no checkout).** `_fragments` from `tests/forest/test_phase1_correspondence.py` on `metabolites`, and the same fragments from `metabolize`. Sets are identical for every case below; live `xf.csmi` matches. Archive has every neutral golden. Live does not. The earlier “both contain thiazole and coumarin” line on `b45a597` is not reproduced. Log: `artifacts/epox_archive_vs_live.live.log`.
  - Thiazole: archive 7 includes `CSC1=CC2OC23SC(N)=NC3=C1` (no formal charge). Live 8, golden absent: 5 nonzero formal charge (`[NH2+]` / `[S+]`), 3 `[SH]` tautomers.
  - Aryl sulfate: golden `O=S(=O)(O)OC1=CC=C(O)C2OC12` archive only. Live keeps the other two neutral regioisomers and replaces that one with `O=S(=O)(O)[O+]=C1C=CC(=[OH+])C2OC12`.
  - Sulfonamide: golden and both other regioisomers neutral in archive. Live emits only `[S-2](N)([O-])[O-]` analogues (net −4).
  - Nitro: archive epoxides keep `[N+](=O)[O-]` (net 0). Live same skeletons as `N([O-])[O-]` (net −2).
  - Nitroimidazole: both goldens in archive, plus two extra net-0 zwitterions live does not emit. Live’s only product is `CC1N(CCO)C([N+](=O)[O-])C2ON21` (nitro kept, not a golden).
  - Coumarin: pyrone epoxide in archive (7). Live emits nothing.
  - Quinone `CS(=O)(=O)c1ccc(-c2cn3ccccc3n2)cc1` (`QuinoneFormation`): neutral golden in archive (14, no formal charge). Live 9 lacks it and adds `CS(=O)(=O)c1ccc(C2=Nc3cccc[n+]3C2=O)cc1`.

- **Bisect, no fix.** Reproducing call is `Rule().metabolites` then canonical fragment SMILES (same loss at `b45a597` via `metabolize` + `p.xf.csmi`): thiazole 8/8 charged, coumarin empty, quinone golden absent with `[n+]`. Archive still has all three goldens. Epoxide neutrals first disappear at `7ef23ae` (parent `e524683` still has both). That commit still picks kekulé order 2.0; it replaces RDKit `KEKULE_ALL` plus SMARTS `=` with a cached assignment and `=,:`, not a single-bond start. Quinone golden is already missing in the first poc rules (`891a2c2`); `[n+]` first appears at `dd9c8ef` (parent `4ea915a` clean). Seed bond stays double there; the diff scopes parents to aromatic atoms.

- **Epoxidation poc/archive parity does not cover the bioactivation golds.** Live-vs-archive product equality is `tests/forest/test_parity.py` (no Epoxidation; quinone is benzene/phenol only) and `tests/forest/test_phase1_correspondence.py` (Epoxidation is only `CC=C` → `CC1CO1`). Not xfail. On HEAD `b45a597`, single-rule canonical fragments: coumarin sets match and both contain the pyrone epoxide; thiazole live is a subset of archive (4 vs 7) and both contain `CSC1=CC2OC23SC(N)=NC3=C1`. Neither set is charged on those two. Nitroimidazole archive adds two charged products live drops; both still emit the two golden epoxides. A set-equality test on thiazole or nitroimidazole would fail, not pass.

- **0.7.2 hotfix:** ``info["rule"]`` is ``list[ReactionRule]`` (leaf first; RuleSet appends ``self``), matching ``addition["rules"]``.

- **CI 3.14 SIGILL = pynauty `-march=native`, not per-Python uv cache.**
  `astral-sh/setup-uv` keys already include Python (`…-gnu-3.14-<hash>` vs
  `3.11`/`3.12`/`3.13`); Hypothesis cache is also per-py. Do **not** add a
  redundant `cache-suffix: ${{ matrix.python-version }}`. Failures are
  `Fatal Python error: Illegal instruction` in `pynauty.autgrp` /
  `nautywrap…so` on Ubuntu GHA; 3.12 stays green (manylinux wheel). Sdist
  `Makefile.nauty` runs `./configure CFLAGS='-O4 -fPIC'`; nauty then adds
  `-march=native`. Same uv cache hit can pass or SIGILL on different runner
  CPUs. Upstream: https://github.com/pdobsan/pynauty/issues/49 (refs #39); fix PR: https://github.com/pdobsan/pynauty/pull/50.
  Workaround on PR #18: `CC`/`CXX` with `-march=x86-64 -mtune=generic` plus
  `cache-suffix: portable-x86-64` in `test.yml` / `release.yml`.

- **0.7.1** is the dead-code drop plus site-match typing patch.

- **site / discovered_site contract.** Emitted `site`: directed_bond tuple→frozenset only; under `canonical_emitted_sites`, remapped to lex orbit. `discovered_site` always raw (no wrap). Dropped `_site_tuple`.

- **Dead forest shims (coverage).** CI-shaped
  `uv run pytest tests/forest src/xenosite/forest -n auto --cov=xenosite.forest`
  (1616 passed). 0-hit and uncalled: `_top_site` (wrong `topol_equiv` remap;
  canonical path is `canonical_emitted_site`), `stamp_forest_labels` /
  `install_forest` / `forest_trace` / `set_terminal_product` / `_work_copy`
  (inlined `copy_mol`), `reordered_forest_labels` (live helper is
  `rdkitutil._reordered_forest_labels`), `may` / `must`,
  `SmirksReactionRule._get_site` / `_get_product_mappings`,
  `odd_anchor_pairs`, `alternating_path` (plural stays), `_merge_options`,
  `is_canonical_orbit_candidate`, `atom_pair_orbit_pynauty_ordered`,
  `_ensure_smiles_tables`. Left: `ReactionRule.metabolites` base,
  `bfs.main`, AtomTracker facade, xf orbit plumbing, isotope/pynauty
  profile oracles, `cannonicalize_order`, `ordered_bond_pair_orbit`,
  `XfTracing._install`. `_rule_name` stays (trace bookkeeping).

## 2026-09-20

- **FA+force / override experiment (reject; no PR).** Branched from
  `main`@PR#17. Proposed merge: add FA + `force-future-annotations`,
  `reportMissingOverride=error`, `reportDeprecated=warning`,
  `typeCheckingMode=standard`, plus `typing_extensions`/`tomli`.
  Reality on ruff **0.16.8**: `lint.flake8-future-annotations.force-future-annotations`
  is **gone**; successor is `lint.future-annotations` (bool) — only lets
  *other* rules (e.g. TC*) inject `__future__` annotations, does **not**
  force FA100/FA102 on py311+. Correct pyright option is
  `reportImplicitOverride` (`reportMissingOverride` unrecognized).
  Measured on live forest + `tests/forest`: **FA100=0, FA102=0**,
  `ruff check --fix --diff` empty. Only `__init__.py` + `_version.py`
  lack the future import (no annotation surface). Pyright
  `reportImplicitOverride=error`: **9 errors / 2 files** (`rules.py` 6,
  `rulesets.py` 3); `reportDeprecated=warning`: **0**. No tomli need
  (tomllib-only; no TOML read path). **Decision: do not switch.** Keep
  `select=[E,W,F,I,UP]` + UP031 ignore; leave FA off; leave
  `reportImplicitOverride` off (9 is tractable later with
  `typing_extensions.override`, not worth a dep for this alone). No
  mypy; CI already on GHA. Experiment branch deleted; config reverted.

- **Disk: drop multi-version bench venvs.** Removed `.venvs/{py3.11,py3.12,py3.13,py3.14}`
  (~299–302 MiB each). Kept primary `.venv` (~299 MiB). Repo
  `du` **1.51 → 0.34 GiB** (~**1.17 GiB** freed). Added `.venvs/` to
  `.gitignore`. Future version benches: `uv run --python 3.X` without a
  persisted `.venvs/` matrix.

- **typing_extensions / Ruff FA+TID blast (dry-run, no switch).** On top of
  current `select=[E,W,F,I,UP]` + `target-version=py311` / pyright
  `pythonVersion=3.11`, temporarily enabled `FA`+`TID` with
  `banned-api` for `typing.override` / `typing.Buffer` pointing at a
  hypothetical `forest._compat`. Scope: `src/xenosite/forest` +
  `tests/forest`. Counts: **FA100=0, FA102=0, TID251=0**, no other TID
  hits; full proposed select still exit 0 (UP031 remains ignored: 31).
  FA is a no-op at py311 (rules only fire for older `target-version`).
  58/75 forest files already have `from __future__ import annotations`;
  17 omit it — FA does not force them. Zero live imports of
  `typing.override` / `typing.Buffer` / `TypeAliasType` / PEP695 `type`.
  `typing_extensions` is not a direct dep (transitive via pyright etc.;
  deliberately dropped from runtime). Three `# type: ignore[override]`
  only — no `@override` call sites. **Decision: leave UP-as-is.** Skip
  FA; skip `_compat` + re-add of `typing_extensions` until a real
  `@override`/Buffer need; TID251 bans optional later (zero blast) but
  not worth a PR without a destination. Keep
  `reportImplicitOverride` off; no mypy.

- **Tooling toward 3.11+ typing.** Ruff: keep `target-version=py311`; enable
  `UP` (pyupgrade); ignore `UP031` (printf→format flood). `UP040`/`UP046`/
  `UP047` (`type` aliases / PEP 695) stay off via target until floor ≥3.12.
  No FA rules that strip `from __future__ import annotations`. Pyright:
  `pythonVersion = "3.11"` (already). Skip `reportImplicitOverride` while
  floor is 3.11 (`typing.override` is 3.12+; do not re-add
  `typing_extensions` just for `@override`). One UP fix: `Sequence` import
  in `bfs.py` → `collections.abc` (UP035).

- **Post-3.10 antipattern audit (live forest).** Drop commit `7d61dde` already
  removed `typing-extensions` and moved `NotRequired` to stdlib `typing`;
  ruff `target-version=py311`. Re-scan of `src/xenosite/forest` + tests +
  packaging: **zero** hits for `Self` TypeVar / quoted fluid returns,
  `tomllib`/`toml`/`tomli` read paths (tomli only transitive via
  `coverage[toml]` on ≤3.11), `sys.version_info` guards, asyncio /
  `TaskGroup`, `ExceptionGroup`/`except*`/`add_note`, or explicit dict
  merge-copy loops that are just `{**a,**b}`/`|`. Kept all
  `from __future__ import annotations` (still support 3.11–3.13; no
  mass-unquote). Applied: `match` on four `Site` shape helpers in
  `rules.py`; `pythonVersion = "3.11"` in pyright. Not inventing async
  for speed — search is sync RDKit-bound.

- **find_path wallclock 3.11 vs 3.14.** Smoke: `artifacts/find_path_wallclock_py_compare.py` (CASES+LARGER from `tests/forest/bench_find_path_h2h.py`, find_path-only, 40×, 400 hits). uv `.venvs/py3.11` / `.venvs/py3.14`, rdkit 2026.03.5 both. `/usr/bin/time -p` real: 3.11.11 **23.23s**, 3.14.7 **24.28s** (+4.5%); live_sum 22.77s → 23.70s (+4.1%). Quick cProfile on one pass: ~46% xenosite / ~46% rdkit tottime (bindings + C++); interpreter bump does not move this wall much. Log: `artifacts/find_path_wallclock_3.11_vs_3.14.live.log`. TaskGroup/asyncio N/A.

- **Pytest wallclock 3.11–3.14** (CI-matched: `pytest tests/forest src/xenosite/forest -n auto --cov=xenosite.forest --cov-report=term-missing`; uv envs `.venvs/py3.XX`; `/usr/bin/time -p`). Same counts each: 1 failed, 1615 passed, 1 skipped, 2 xfailed. Wallclock real: 3.11.11 **24.38s**, 3.12.9 **26.64s** (+9.3% vs 3.11), 3.13.2 **26.36s** (+8.1%), 3.14.7 **22.04s** (−9.6%). Pytest session: 10.33 / 11.95 / 11.51 / 7.76s. Log: `artifacts/pytest_wallclock_3.11_vs_3.14.live.log`.

- **Python 3.14 support (post v0.7.0).** v0.7.0 test+release CI green on
  main. Local `uv` 3.14.7: deps sync (rdkit 2026.3.5, pynauty, numpy);
  forest suite **1616 passed**, 1 skipped, 2 xfailed. Added classifier,
  README 3.11–3.14, CI matrices in `test.yml` / `release.yml`, towncrier
  `changelog.d/20.added.md`. Not testing freethreaded **3.14t** in CI:
  `rdkit` has no `cp314t` wheels (only `cp314`); `uv sync` fails. Numpy
  installs on 3.14t; pynauty builds from sdist but re-enables the GIL on
  import (`nautywrap`). Revisit when RDKit ships freethreaded wheels.
- **uv vs Homebrew Python.** Prefer **uv-managed** CPython for this repo
  (`uv python install` / `uv sync` / `astral-sh/setup-uv` in CI). Matches
  CI and global-caches guidance; Homebrew `python@3.14` is fine as a
  system interpreter but project envs should stay uv-pinned.

- **`cleared_0_7` mark.** Baseline `v0.6.1`; marked the 28 hard rules-example xfails still open after the cleavage XPASS clear (`b2830ee` `_XFAIL_IDS`), now passing on live forest. Not the quinone `Reaction*` bulk. Still xfail: Glutathionation `epoxide_c` / `aziridine_c` pattern_info. `pytest -m cleared_0_7`: 28 passed.

- **`test_distinct_signatures_distinct_product_sets[dh-aminophenol]` skip.**
  `Nc1ccc(O)cc1` DH emits **1** pair sig → one product `N=C1C=CC(=O)C=C1`;
  skip is correct (not a latent failure). Removing skip would pass
  vacuously. No other param in that test covers DH with ≥2 sigs; QF
  siblings pass. Multi-end DH mols (e.g. `Oc1ccc(N)c(O)c1`,
  `Nc1ccc(O)c(O)c1`) get 2 sigs and the injectivity assert **passes**.
  Still skip-worthy for p-aminophenol; better substrate optional. No code
  change.

- **pynauty required; no SMILES pair-orbit product path.** `pynauty` moved to
  `pyproject.toml` dependencies. `get_pair_orbit_backend()` always `"nauty"`;
  removed `PynautyRecommendedWarning`, env override, and smiles dispatcher.
  Isotope helpers remain for oracle/profile tests only. Docs: PAIR_ORBITS /
  XF internals / MIGRATING / CHANGELOG.
- **AtomTracker deprecated → xf.tracing.** Tutorial in `docs/forest/XF.md`;
  helpers `depths()` / `index_at` / `added_indices` / `root_map`. Pair-orbit /
  `_of_products` stay underscored (not public xf). Scrubbed
  `attach_phase1_steps` from public usage/MIGRATING (archive-only; not a live
  public contract).


- **Drop `info["csmi"]`.** No lazy/eager emission frozenset on ProductInfo
  (zombie risk from retained mols; extra API). Callers use `product.xf.csmi`
  or `frozenset(p.xf.csmi for p in products)`. Internal unique_csmi / check /
  RuleSet still compute emission frozensets from products. `ProductInfo` =
  `SiteInfo`. Docs: MIGRATING / HEURISTICS / XF / usage / README / CHANGELOG.

- **WAE dearomatizes / hydrogenation1.** PatternInfo `dearomatizes` is
  capability; coverage test asserts resolve vs `_site_map_aromatic`, not
  declared==resolved. `filter_rules` all-dearomatizes refuse skips patterns
  that also add H (path_end aliphatic `CC=O`→`CCO` keeps Hydrogenation).
  MeOPhOH bill **180→226** (still ≪900); regression gate `<250`.

- **`metabolize` list-yield restore.** Yield shape is again
  `(list[TracingMol], ProductInfo)` per emission (1-element non-cleavage;
  cleavage siblings together). Dropped `product_index` / `product_count`.
  `unique_csmi` drops duplicate *emissions*, not within-list sibling CSMIs.
  RuleSet / bfs unpack lists; bfs nodes stay one mol.

- **`XF.md`.** Full `docs/forest/XF.md` for `Mol.xf` (public accessor, not
  forest schema). Frames cache-into-`_forest` + link to `records.py`
  TypedDicts as unstable layout. Committed `54b7ded`.

- **H2H post Hydrogenation filter.** `uv run python tests/forest/bench_find_path_h2h.py` → live total **1.133s → 0.264s**. MeOPhOH→HQ **0.993s → 0.206s**, bill **898 → 180**, nd 69 → 28 (same 4-step plan). Eugenol bill 45→19; 2-MeO→1,2-NQ 43→7. `--larger` flat (live 0.269s). Artifacts: `bench_find_path_h2h_3way_post_filter.*`, `bench_find_path_h2h_larger_post_filter.*`. Docs: PERFORMANCE.md current numbers. Confirms dramatic speedup. **No commit.**

- **Larger-mol PERFORMANCE depictions.** Generated SoM-marked SVGs for five `--larger` H2H cases via `uv run python tools/som_depict.py --preset larger --out-dir docs/forest/performance_assets` (find_path SoM + MCS align). Embedded reactant/product pairs under PERFORMANCE.md “Larger mols”; render cmds now point at `som_depict.py` presets (MeOPhOH mid-size assets already present). **No commit.**

- **MeOPhOH bill: Hydrogenation adds H (not DH).** Diagnosis: (1) `filter_sites` dearomatize gate only checks `loses_aromaticity` — allows reductive **Hydrogenation (adds H / saturates)** toward non-aromatic quinone targets that need oxidative **Dehydrogenation (removes H)**. Pre-fix `path_end` also lacked `dearomatizes=True` capability. (2) Retracing: `seen` keys **child CSMI at enqueue**; parents expand once, but new parents still pay `mol_edits` before seen — 56/98 closer edges blocked post-chemistry (mostly Dealkylation). Fix: PatternInfo `path_end` `adds=H`+`dearomatizes=True`; generic adds-H filter (mirror removes-H); alkene SMARTS keeps span dearomatizes false (aliphatic C=C→CC). Bill **898→176**. Tests: `test_hydrogenation_pattern_info.py`. HEURISTICS Status: approved. Artifacts: `artifacts/meophoh_{probe,retrace,ablation}.*`. Residual: Epoxidation + seen-after-edit. **No commit.**

- **Site typing (pyright gate).** Helpers that take site identity now annotate `Site` (not `Collection[int]`): `site_signature`, `bond_rank_key`, `site_orbit`, `canonicalize_smarts_match`, `_remap_site`, `_copy_site`, `_remap_match_via_auto`. Automorphism `atom_map` on `_remap_site` is `Sequence[int]` (tuple index API), matching `remap_mapped_atoms`. Bare `int` remains in the `Site` union; remap/canonicalize use `_flat_ints`. `uv run pyright src/xenosite/forest` → **0 errors**. **No commit.**

- **Larger-mol H2H** (`bench_find_path_h2h.py --larger`). Five HA 17–26 cases from path-outcome / substrate_library. Totals: BFS **4.011s** · DFS **5.165s** · live **0.267s** (hits 2/5 · 0/5 · 5/5). Live wall **0.005–0.135s** (~30×); archive CAP floors ~0.55–1.15s — clearer **wall** spread vs mid-size flat aggregate. Live bills stay small (4–36); bill spread still MeOPhOH-dominated, not size. Artifacts: `bench_find_path_h2h_larger.{out,live.log,tee.log}`. Docs: PERFORMANCE.md “Larger mols”. **No commit.**

- **Smarts→Smirks naming.** Renamed mislabeled SMIRKS surfaces to match `Smarts` / `Smirks` NewTypes: `SmartsReactionRule`→`SmirksReactionRule`, `.smarts`→`.smirks`, `reaction_from_smarts`→`reaction_from_smirks`, `_smarts2rxns`→`_smirks2rxns`, `_isotope_smarts`→`_isotope_smirks`, `_smarts_mapped_bond_order`→`_smirks_mapped_bond_order`, `SmartsSiteInfo`/`SmartsProductInfo`→`SmirksSiteInfo`/`SmirksProductInfo`. Left match-only: `endpoints`, `smarts_matches`, `_reactant_smarts`, `canonicalize_smarts_match`. Docs: HEURISTICS / DIVERGENCES. Pyright: same 5 pre-existing Site/Collection errors; focused rename tests: **318 passed / 2 xfailed**. **No commit.**

- **H2H PERFORMANCE rerun** (after CSMI check/yield + `site_kind`). Command: `uv run python tests/forest/bench_find_path_h2h.py`. Tree: `fec1594` + dirty. Live total **1.133s** vs morning **1.117s** (~+1.4%, flat). Per-case live: eugenol 0.051 (was 0.062), PEA 0.009 (0.015), MeOPhOH 0.993 (0.947; bill 898 vs 785), TBA 0.015 (0.017), 2-MeO→1,2-NQ 0.066 (0.075). Hits unchanged 5/5 live. Artifacts: `bench_find_path_h2h_3way.out`, `bench_find_path_h2h_3way_post_csmi.{out,live.log}`. Docs: PERFORMANCE.md updated. No `_kekule_forms`-drop win in aggregate. **No commit.**

- **Redundant-rules → INFO only.** `_report_redundant_rules_drop` no longer calls `warnings.warn`; overlaps log INFO (`kept_rule` / `dropped_rule` / product CSMI). Removed `RedundantRulesWarning` class and `pyproject.toml` `ignore::RedundantRulesWarning`. `SiteDeduplicationWarning` unchanged (still Warning + WAE error). Docs: REDUNDANT_RULES / HEURISTICS / PAIR_ORBITS / DIVERGENCES. Tests: caplog in `test_ruleset` / `test_csmi_check_vs_yield`. Forest WAE: **1607 passed / 0 failed** (2 skipped, 2 xfailed). **No commit.**

- **Dedup check vs yield (implemented).** Leaf `ReactionRule.metabolize`: **check** always (emission CSMI frozenset + site ranks under `(rule, pattern)` → `SiteDeduplicationWarning` on miss; no drop); **yield** only if `unique_csmi=True` (quiet fragment CSMI drop). Replaced fragment-subset `novel_exists` gate with emission-set equality. `RuleSet.metabolize` forces children `unique_csmi=False` (yield off, check on); parent yield is cross-child CSMI → `RedundantRulesWarning`; nested sets same (no double-drop). pytest `ignore::RedundantRulesWarning` in `pyproject.toml` (see `REDUNDANT_RULES.md`); `SiteDeduplicationWarning` stays error. HEURISTICS Status: **approved**. Tests: `test_csmi_check_vs_yield.py`. Forest WAE: **1607 passed / 0 failed** (2 skipped, 2 xfailed). **No commit.**

- **Dedup check vs yield (design only).** Agree: separate unique-edit-miss **check** from CSMI **yield** drop — check when site unique-edit is on; yield only if `unique_csmi=True`; check-with-yield-off enables RuleSet child-csmi-off bubble-up without losing `SiteDeduplicationWarning`. CSMI caching on `xf.csmi` makes double-read cheap after first post-`of_products` miss. Footguns: no site-dedup toggle today (check follows that layer if one lands); parent must not re-warn same-rule; unequal-rank iso stays quiet. HEURISTICS: Status not decided. **No code change.**

- **`Smarts` / `Smirks` NewTypes.** Match-only vs reaction (`reactant>>product`) brands in `records.py` (same pattern as `TopoGroupId`). Rule data: `SmartsReactionRule.smarts` / `rxns` → `Smirks`; `ResonancePairRule.endpoints` → `Smarts`. Helpers: `smarts_matches` / `_smarts_matches` / `Structure.smarts_matches` → `Smarts`; `reaction_from_smarts` / `run_reactants` / `react_at` / `_isotope_*` / `_smarts2rxns` / `_ndealk` → `Smirks`. Field name `smarts` kept (historical); docstring notes values are SMIRKS. Runtime unchanged (NewType is still `str`). Pyright: **0** new Smarts/Smirks errors; 5 pre-existing Site/Collection noise remain. Smoke: `test_rdkitutil` + `test_ruleset` 13 passed; invariants/pattern_info 203 passed / 2 xfailed. **No commit.**

- **Redundant pairs doc + CsmiDedup audit + uniquify-plan feedback.** Inventory of 6 cross-rule overlaps in `docs/forest/REDUNDANT_RULES.md` (Status: not decided on treating as permanent expected); HEURISTICS links it. Audit of recent CsmiDedup→0 fixes vs data-not-branches: (1) `atom_pair` SMARTS → `bond_rank_key` is PatternInfo/`site_kind` semantics (undirected one-bond on ResonancePair rules; pair paths still `pair_site_signature`) — data-driven, no fix. (2) warn only when full CSMI dup **and** site ranks match keeper is correct unique-edit-miss detector (quiet ≠ drop; unequal ranks = product iso / leaving-group / cleavage siblings) — data-driven, not a silent product filter; no code change. Plan feedback (not implemented): agree site vs CSMI layers + child-csmi-off for cross-rule bubble-up; disagree “any rule-level CSMI = site-dedup error / fail tests” while unequal-rank iso exists; child-off loses within-rule CsmiDedup on RuleSet path unless parent also checks same-rule; prefer kwargs-off over callback; nested PhaseOne RuleSets need outermost-only CSMI. **No commit / no redesign.**

- **SiteDeduplicationWarning → 0 under forest WAE.** Before: ~34–38 WAE fails mostly CsmiDedup (panel census H 34 / NR 6 / Dealk 4). After: **0 CsmiDedup**; WAE **5 failed / 1594 passed** — all `RedundantRulesWarning` (Dehydration↔Hydrolysis ×2, OxygenReduction↔Hydrogenation ×3). No global ignore of CsmiDedup. Fixes: (1) `site_signature` uses undirected `bond_rank_key` for `site_kind in ("bond","atom_pair")` so Hydrogenation alkene/alkyne + DH one-bond SMARTS collapse symmetry (phenol ortho/meta, butadiene ends, DH alkyl). (2) `metabolize` CSMI warn gated: warn only when emission is a *full* CSMI duplicate *and* site ranks match the keeper (true unique-edit miss); quiet isomorphic cleavage siblings (Azo two anilines, ethanol two methanols, gem-Cl, nitro two O), shared leaving groups, and distinct-site product isomorphism (RD vicinal Br→same alkene; epoxide either carbon→same diol). Panel census now 0 CsmiDedup / 11 Redundant. Tests: `test_hydrogenation_unique_edit.py`, `test_csmi_dedup_quiet.py`. Docs: HEURISTICS site_kind + warn semantics. Artifacts: `artifacts/forest_wae_csmi_fix2.live.log`, `artifacts/census_csmi_vs_redundant.out`. **Ask:** should WAE ignore `RedundantRulesWarning` (chem-only previously ignored CsmiDedup only — no existing Redundant ignore)? Prefer leave-as-is until decided. **No commit / no merge / no tag.**

- **RuleSet CSMI → `RedundantRulesWarning` (not dedup bug).** Per-rule `unique_csmi` drops stay `SiteDeduplicationWarning` (`ReactionRule.metabolize` only). `RuleSet.metabolize` keys product CSMI across *different* child rules; later rule dropped with `RedundantRulesWarning` naming both rules + shared CSMI (INFO: `kept_rule` / `dropped_rule`). Same-rule multi-pattern same-CSMI may still both emit. Focused: `test_ruleset_redundant_rules_warning_names_both_rules`, `test_ruleset_redundant_rules_pytest_warns`. **WAE forest** (`tests/forest`, `filterwarnings=error`): **38 failed / 1543 passed** — E-lines ≈ **34 `SiteDeduplicationWarning` / 3 `RedundantRulesWarning`**. Remaining **true** per-rule dedup (do not rebrand): Dealkylation, Hydrogenation, NitrogenReduction, AzoSplitting, OxidativeDehalogenation, Dehydration, Dehydrogenation, ReductiveDehalogenation. Panel census (`artifacts/census_csmi_vs_redundant.{py,out}`): CsmiDedup=44 (H 34 / NR 6 / Dealk 4), Redundant=11 (e.g. OxygenReduction↔Hydrogenation, NitrogenReduction↔Dehydration). RuleSet path never emits `SiteDeduplicationWarning`. Docs: HEURISTICS / PAIR_ORBITS / DIVERGENCES. **No commit / no merge / no tag.**

- **`directed_bond` public site = frozenset; orientation on `discovered_site`.** No `directed_sites` kwarg / no extra field. Unique-edit still keys directed `MapRankKey` before yield. With `canonical_emitted_sites`, `site` is frozenset(lex) and `discovered_site` stays the directed discovery tuple. Trace preserves tuple order for `discovered_site`. Docs: HEURISTICS / PAIR_ORBITS §5b / records. Tests: `test_directed_bond_unique_edit.py`, `test_site_kind.py`. **No merge / no tag.**

- **`directed_bond` emits ordered tuple.** `_site_indexes(..., site_kind=)` returns `tuple` in `site_map` order for `directed_bond` (Dealkylation / NDealkylation / Benzodioxole / Nitroaromatic); undirected `bond` stays `frozenset`. Meta-test asserts container type. Canonical remap preserves tuple vs frozenset. TODO: rerun H2H + update PERFORMANCE.md after CSMI/`site_kind` stabilize (WAE not quiet yet). **No merge / no tag.**

- **NDealkylation → ResonanceRule + directed_bond check.** Unique-edit for `directed_bond` uses ordered `map_rank_key` (vs undirected `bond_rank_key` for `bond`). Anisole Dealkylation: undirected would merge 3 orientation pairs with distinct products (`C=CC(=CC=CO)OC` vs `C=CC=C(C=CO)OC`, …); directed keeps them; true same-orientation duplicates collapse; 0 CSMI on anisole/TMA/aniline ND probes. NDealkylation SMARTS are the N rows of Dealkylation (`site_map=(1,2)`, `partner=N`, methyl `leave_count=1`) — OK. Gap: plain `SmartsReactionRule` matched pyridine methine C–N but emitted 0 (no Kekulé parent); archive ND had ring-open. Reparented to `ResonanceRule` like Dealkylation. Tests: `test_directed_bond_unique_edit.py`, pyridine parity. **No merge / no tag.**

- **Epoxidation CSMI + site taxonomy (`bond` / `directed_bond` / `atom_pair`).** Phenol: only ordered `map_key` split CSMI twins — `((1,4),(2,2))` vs `((1,2),(2,4))`; `incident_orders`/orbit same. Undirected `bond_rank_key` fixes Epoxidation (0 CSMI). Same undirected key on Dealkylation wrongly drops anisole regioisomers (`C=CC(=CC=CO)OC` etc.). Taxonomy: `bond` = undirected ranks (Epoxidation, Azo); `directed_bond` = directed MapRankKey (Dealkylation, NDealk, Benzodioxole, Nitroaromatic); `atom_pair` = ResonancePair only. Meta-test updated. Phenol Hydroxylation still 0 CSMI. Chem-only: **1574 passed**. WAE: **36 failed** (was ~55; no Epoxidation left). Artifacts: `artifacts/epox_directed_{chem_only,wae}.live.log`. **No merge / no v0.7.0 tag.**

- **Ship ResonanceRule SMARTS-implied bond order + drop global `_kekule_forms` on PR #15.** `_smarts_mapped_bond_order` reads reactant template maps 1–2 (`=` / `=,:` → 2.0; `-` / `-,:` / unspecified → 1.0). `_reactant_parent` uses that on aromatic hits; if no Kekulé assignment exists, fall back to aromatic `mol` (charged rings / Reaction94457). Unique-edit `incident_orders` on aromatic parent (not Kekulé work). Reparented + aromatic-matching SMARTS: `ThiopheneSulfurOxidation` (`=,:`), `Dealkylation`, `AzoSplitting` (`=,:`), `NitrogenReduction` (hydroxylamine `-,:`). Plain `SmartsReactionRule.metabolites` loops on input `mol` only (`_kekule_forms` left for tests/other callers). Phenol Hydroxylation alone: **0** CSMI. Chem-only (`ignore SiteDeduplicationWarning`): **1569 passed / 0 failed**. WAE still **55 failed** — residual CSMI mostly Epoxidation (+ Hydrogenation / Dealkylation / OxidativeDehalogenation / …); suite not green under warnings-as-errors. Docs: DIVERGENCES (no global kekule loop; ResonanceRule parent selection). **No merge / no v0.7.0 tag.** Artifacts: `artifacts/resonance_order_chem_only2.live.log`, `artifacts/resonance_order_wae.live.log`, `artifacts/phenol_csmi_after_resonance_order.out`, `artifacts/phenol_hydroxylation_no_kekule.out`, `artifacts/forest_no_kekule_{wae,chem_only}.live.log`.

- **Drop global `_kekule_forms` (pre-ship probe).** After loop removal alone: phenol Hydroxylation **0** CSMI (was 2); WAE **64 failed / 1503 passed**; chem-only **9 chemistry fails** (thio S-ox empty, pyridine dealk ring-open, isoxazole N-red, anisole dealk parity, azo pyridazine filter, …) until ResonanceRule reparent + aromatic SMARTS closed them. Decision: keep loop off; ship with ResonanceRule order fix (entry above).

- **Phenol unique-edit `seen` trace (Hydroxylation alone).** Instrumented `SmartsReactionRule.metabolites` per candidate; reverted. Artifact: `artifacts/phenol_hydroxylation_unique_edit_seen.out`. Kekulé0 keeps 5 sites (`{2,3,4,5,6}`); ranks already match for ortho `{2,6}` (rank 4) and meta `{3,5}` (rank 2), but **`incident_orders` flips 1.0↔2.0** so signatures differ → unique-edit miss. Kekulé1 all SKIP (cross-form hit). CSMI then drops `{5}`/`{6}` (2 warnings). Root: Kekulé bond orders in signature, not rank mismatch.

- **Phenol Hydroxylation alone → CSMI warning YES.** `Hydroxylation().metabolize(Oc1ccccc1)` (not via RuleSet): **2** `SiteDeduplicationWarning` drops, pattern `h` (no pattern named `c`; carbon-OH is `h`/`h2`). Dropped sites `{5}`→`Oc1cccc(O)c1`, `{6}`→`Oc1ccccc1O`; kept 3 (`{2}`,`{3}`,`{4}` → ortho/meta/para). Same 2 Hydroxylation drops inside `PhaseOne` (plus Epoxidation×2, Hydrogenation×2). Probe: `artifacts/probe_phenol_hydroxylation_csmi.py`.

- **CSMI dedup: swap_group theory vs inventory.** `swap_group` is only on ResonancePair/`pair_site_signature` — not on SMARTS `site_signature` (Hydroxylation, most PhaseOne). Phenol `{5}`/`{6}` falsifies swap_group: ranks already match; Kekulé `incident_orders` splits ortho/meta. Warning text now includes **rule name** (identical per rule → once-per-rule under stdlib/`-W once`). Forest suite INFO drop counts (serial, `artifacts/csmi_dedup_rule_inventory_serial.live.log`): Dealkylation 302, NitrogenReduction 140, Hydroxylation 140, Epoxidation 109, Hydrogenation 31, OxidativeDehalogenation 25, Dehydration 24, ReductiveDehalogenation 13, BenzodioxoleReduction 10, NitrogenOxidation 7, NitroaromaticReduction 6, NDealkylation 6, AzoSplitting 5, SulfurReduction 2, Dehydrogenation 2, OxygenReduction 1, Glucuronidation 1 (**17 rules**, ~824 drops). Only DH (+ Hydrogenation as ResonancePair) are pair-signature candidates for swap_group; bulk is non-pair. No signature fix yet.

- **CSMI dedup diagnostics:** `ReactionRule.metabolize` / `RuleSet.metabolize` emit `SiteDeduplicationWarning` every drop (identical generic message; stdlib once-per-message when not error) plus `logger.info` per drop (`substrate`, `rule`, `site`, `pattern`, `product`; INFO off by default). Suite under warnings-as-errors: **87 failed / 1482 passed** — all sampled failures are `SiteDeduplicationWarning` (unique-edit/orbit still misses duplicates). Artifact: `artifacts/csmi_dedup_warn_suite.live.log`. Prefer fixing root causes over globally silencing; do not ignore in pyproject yet.

- **Pytest:** uncaught warnings fail tests (`filterwarnings = error` in pyproject). Archive `UnstableWarning` still ignored; AtomTracker deprecation asserted via `pytest.warns` / local `catch_warnings(ignore)`.

- **`site_kind` + examples (pre-v0.7.0):** `RuleSiteKind` `"atom"` / `"atom_pair"` on every leaf; `_example_substrates` short SMILES; `tests/forest/test_site_kind.py` asserts emit len matches declaration. Expand examples later for all patterns/whens (TODO).
- **Retire `bond_atom` UniqueOrbit:** never an intended Site pattern (POC unique-edit framing from `dd9c8ef` / `d30b6d5` / `ba836e5`). DH/QF/Hydrogenation declare ``site_kind="atom_pair"`` + ``sites_on="atom_pairs"``; one-bond SMARTS ``site_map=(1,2)``. Unique-edit is atom–atom only (`AtomPairOrbitSignature`). Dropped `BondAtom*` types, `unique_orbit`, xf `bond_atom_orbit_key`. Site shape is `site_kind` class data (+ `_example_substrates` meta-test) — not `Generic[SiteT]` (heterogeneous RuleSets; pyright cannot enforce frozenset cardinality). Docs: PAIR_ORBITS / HEURISTICS / DROPPED approved.


- **Promote + v0.7.0:** `refactor_poc` → live `xenosite.forest`; pre-swap forest archived to `src/xenosite/_archive_forest/` (read-only / locked; CI excludes archive tests + lint). Tests → `tests/forest/`. `AtomTracker` facade kept, deprecated (prefer `mol.xf`); StepPlan still imported from archive. Design docs moved to `docs/forest/` (`PERFORMANCE.md`, `PAIR_ORBITS.md`, `HEURISTICS.md`, `DIVERGENCES.md`, `DROPPED.md`, `performance_assets/`); package stays code-focused. Stripped leftover POC/`refactor_poc` wording. Release target **v0.7.0** (hatch-vcs; latest tag `v0.6.1`). Commit: `b323ba6` (promote), docs follow-ups through `1500dd1`.
- Long fuzz post-swap: `HYPOTHESIS_PROFILE=long XENOSITE_FUZZ_EXAMPLES=200` on bfs/guided/phase1/and_cleave/find_path_phase1 fuzz — **18 passed** (~23s). Artifact: `artifacts/long_fuzz_post_swap.out`.
- **PERFORMANCE** three-way H2H (archive BFS/DFS vs live `find_path`): harness `MAX_MOLS=200` (`tests/forest/bench_find_path_h2h.py`); raw `artifacts/bench_find_path_h2h_3way.out`. Wins where both archive modes CAP — eugenol→allyl-Q (4 steps / 0.062s), dimethoxy-PEA→catechol (`&`, 0.015s), MeOPhOH→HQ (4 steps / 0.947s), TBA→aldehyde (22 HA, 0.017s). Retuned naphthalene to reachable **2-MeO→1,2-NQ** (3 steps / 0.075s; BFS early hit, DFS CAP). **2-MeO→1,4-NQ** = no PhaseOne path (frontier empty at nd=95), not a find_path bug. Depictions: xenopict circles on SOM only (find_path plan AtomRef origins); MCS-aligned products; regenerator `docs/forest/performance_assets/_render.py`.
- **ResonancePair / butadiene** (`516d66f`): Hydrogenation false emits unmasked (not caused) by unique-edit at `dd9c8ef`. Reactant `C=CC=C` — invalid `C=C=C=C` (single-first) and `C=C=CC` (even bond-count); valid `C=CCC` / `CC=CC`. Fix: double-first `alternating_paths` + odd bond-count only — lasting shared path chemistry, not butadiene-specific / not PatternInfo. Hard-fail: `test_butadiene_hydrogenation_includes_2_butene`. Do not restore bad dedup.
- **Pair-orbit unique-edit (landed):** two-site key is the unordered pair orbit (own signature field; HEURISTICS approved). Recipes in `graph_isomorphism.py`; unified `site_pair_orbits_{smiles,nauty}` return `atom_atom` / `bond_bond` (atom_bond stays six-family only); `PairGroupId` sequential `0..n-1` via CIP. PatternInfo/When ``swap_group`` (approved): unordered iff shared non-empty group, else ordered by ``name``; default ``swap_group`` = ``name`` (omit redundant annotations). Signatures tagged `ordered: Literal[True, False]`. Bond–atom UniqueOrbit later retired (never a Site; see same-day retire note). Prefer ``swap_group`` over SMARTS narrowing for same-group pairs (narrowing-as-replacement: HEURISTICS not approved). Path-budget after `dd9c8ef`: bill once per unique-edit combo → restored butyl=3 / PhCH2OH=30. Profile note: nauty batch ≪ isotope on drug-like sizes (`artifacts/pair_orbit_profile.out`).
- **Canonical lex-orbit emission:** opt-in `canonical_emitted_sites` via explicit kwargs only (dropped env / `set_canonical_emitted_sites`). Lex-rep cache read from **parent** after product `clear_structure` (`eee738f`); singleton atoms remapped too. Profile after parent-cache: OFF **3.862s**, ON **4.007s**, Δ **+3.7%** (was +7.2%). Artifacts: `artifacts/canonical_emitted_sites_profile.{out,pstats,live.log}`. Fuzz draws the flag with Hypothesis. Docs: PAIR_ORBITS §5, HEURISTICS. PAIR_ORBITS.md rewritten (McKay–Piperno 2014; Sharp 1999; Surge/Laffitte optional).
- **`forest_copy`:** explicit walker in `forest_copy.py` (not deepcopy/pickle hooks). Keys: `immutable` / `cache` (was `structure`) / mutable. Profile (same find_path harness): before wall 4.959s, deepcopy cum **1.637s** (~33%); after wall 4.018s, deepcopy **0**. Artifacts: `artifacts/forest_copy_{before,after}.*`.
- **PatternInfo / chemistry hygiene:** named all 98 PatternInfo; within-rule emit-name uniqueness + when-coverage inventory (~190 possibilities / 135 when-branches). Hydroxylation: `h2` → `[#6h2,#6h3]`; partition `h`/`h2` by H count; GSH epoxide/aziridine → `#6H0`. Dealkylation overlapping C–C alcohol SMARTS → `#6H0` (nevirapine). `unique_csmi` key `(rule, PatternInfo.name|SMARTS, csmi)`. CX `atomLabel` preserved via `_forest["start_labels"]`. Methide always in DH data; refuse via `filter_sites` (DROPPED `pathways=` opt-in). Dropped POC `find_network_paths` / net (DROPPED approved); TautomerRule stub + AtomTracker xf facade. POC CI lint (ruff+pyright); untyped params → error; 0 pyright errors after typing. Coverage snapshot: **84.0%** Cover / **87.5%** statements (`artifacts/poc_coverage.txt`). Ported forest base + hard PathOutcome/bfs fuzz into POC before promote. Did not touch `guided_path.py` / `test.ipynb`.

## 2026-09-19

- Two-site unique-edit collapse stays `Status: not decided` (pre-orbit land). One-site rank collapse stays. Checked both-ends rank signatures on benzene, naphthalene, xylenes, biphenyl, diphenylmethane, two-atom SMARTS rules: no group contained two matches with different products.
- Quinone kekulé parents stay on aromatic atoms of the conjugated component (`aromatic_parent_atoms`); biaryl/exocyclic amide not joined. Neutral aromatic atoms no longer take formal charge when 1.5-order bonds kekulize (`move_charge_with_bonds`). Cleared 76 historical xfails (marked `# progression:`). `pytest tests/refactor_poc`: 1191 passed, 1 skipped, 28 xfailed. `xf.tracing` gained `atom_depths`, `atom_root`, `removed_roots`. Fuzz ports: guided/find_path_phase1/phase1_steps/and_cleave/bfs. Still open then: long-range imine quinones; `CC=O`→`CCO` naming (OxygenReduction vs Hydrogenation); `CCSO` cleavage; GSH Michael/aziridine; sulfation example empty.
- **xf mint-on-read:** `Mol.xf` property mints new `Xf` with strong parent ref (no stored xf/weakref/copy-between-mols). Brands `ForestMol` / `NoForestMol` / `TracingMol` (TYPE_CHECKING). Library on xf; deleted dual free APIs (`get_csmi` / `get_forest` / public structure helpers). Conjugation terminal: `ConjugationRule.is_terminal_rule`; find_path skips expanding terminals. Nested `xf.tracing` (`XfTracing`). `of_products` via `copy_mol`. Dropped `cannonicalize_order` from metabolize/`_finish` hot path; product `csmi` lazy. Profile after xf: wall 1.12s (was 1.98s); deepcopy cum 0.387s (~35%, was ~59%). H2H after xf: both_ok=9/10, poc_only=1 (MeOPhOH→hydroxyQ); forest 4.57s / poc 0.67s; both_ok speed ≈ **17.8×**. Artifacts: `artifacts/poc_find_path_profile_after_xf.out`, `artifacts/bench_find_path_h2h_after_xf.out`.
- Lazy heap on `find_path`: key `(hit_tier, seq)`; rescore on pop. Filter API takes live `ForestTracingMol` first. `order_key` prefers cleave → dearomatize → oxygen; `leave_count` prunes. Equal-cost sideways / sibling cost-sort inflated budgets — HEURISTICS not decided / do-not-revive. Raise `ValueError` on `None` inputs. H2H: poc 10/10 hit, both 9/10; PhCH2OH ed=22; TBA ed=3; hydroxyQ ed=331 (forest EXH).
- Cleavage: react → split connected components → `forest_trace` per fragment (no dotted yields). Cleared 37 XPASS xfails; remaining ~105 were chemical/wrong-rule (quinone ~67, rules ~28, GSH 4, phaseone 5). Ported forest correctness suites with `@pytest.mark.xfail` + `regression`: 1009 passed, 146 xfailed.
- Forest epox / N-dealk `phase1_steps` are identity singletons (not StableOxygenation / UnstableOxygenation look-aheads; reverted dfa17e2). PhCH2OH→quinone: MCS ring-membership mismatch adds bridge → path found (ed=21). TBA: prepend target hits → ed=3. H2H harness prints comparable ed/re/nd/wall; 10/10 poc hit, 9/10 both. Composite plans on `ReactionRule.canonical_plan` (QF → prep + DH); Halogen ends read `partner` → OxidativeDehalogenation.
- Typing: closed `Span` / `SiteInfo`; Site = indices only; FutureSite nests `int | AtomRef`. Expanded ruff+pyright to all POC modules: 0/0. Suite grew 126 → 674 passed during plan/filter typing.
- Hydrogenation as ResonancePairRule: path flip writes `CC=CC` from `C=CC=C`; even path `C=C=CC` not a pair on this side. Ethenediol draft sanitizes to glyoxal.
- Ruff as dev dep (`ruff>=0.16.8`); E/W/F/I; line length 100; target 3.10. Forest/`phaseone.py`/`test.ipynb` excluded then. `_BOND` float keys hashed same as ints — dropping duplicates OK.
- Kekulé charge-with-bonds repairs 1,4-dinitrobenzene parse; pyridazine still yields `C=CC=CN` (dropped atom is the reaction).
- Forest state brands: `NoTracingMol` / `NoForestMol` / `ForestNoTracingMol` / `ForestTracingMol`. `ensure_forest` / `ensure_tracing` install in place.
- ConjugationRule subclasses (GSH / glucuronidation / sulfation / acetylation): filters read `partner` / `partner_h` / `symbol`; star collapse default on ConjugationRule. Divergences recorded for old acid SMARTS / arene-oxide sulfone / etc.
- Kekulé parents: match once with `=,:`; one parent per conjugated-system assignment. NDealkylation = N rows of dealkylation; `leave_count` 1 for methyl; `breaks_ring` from cleaved bond.
- Easy parity vs old library: ethane OH, anisole dealk, benzene/phenol QF match; butylbenzene carbonyls from old `[#6h2:1]>>[*:1]=O` in DIVERGENCES. Local RDKit types under pyright (no stubs package).
- `PatternInfo.site_map` = map number or tuple. Project rule `data-not-branches.mdc` (methide-as-data). Atom-trace additions normalized (`R1`, `R2`, …); formula / delta_formula. QuinoneFormation plan = Hydroxylation per O-needing end + one DH. Cleavage filter keeps bridge only. `RuleSet` in `rulesets.py` is a rule; `PhaseOne` lists classes that exist. `find_path`: one MCS + local atom change; `mol_edits` = accepted RunReactants or pair overlay. ResonanceRule / ResonancePairRule in `rules.py`. Hydroxylation `[#6h2:1]>>[*:1]O` (not `=O` / `removes: HH`).
- RDKit isotope SubstructMatch timing (C400/C8000): no isotope index; put isotope-bearing atom first if limited match stays. Already known: tagged matches cost more past small N on small mols.

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
