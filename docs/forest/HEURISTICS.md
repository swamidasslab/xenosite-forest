# Find-path heuristics and correctness

Consult this file before changing `find_path` behavior, and again when a search miss or a slow pair shows up. Add an entry when a new idea appears. Do not delete entries. Mark one done when the test that names it passes without `xfail`.

Follow `.cursor/rules/data-not-branches.mdc`. Prefer a generic decision that reads `PatternInfo` over a special case in the search.
`leave_count` and `breaks_ring` on `Effect` are that data: the same `cleaves` bit, with the leaving piece named or open, and ring opening on the bond.

## Exceptions

An exception to that rule is allowed only with a good reason. Record it here. Each entry has `Status: approved`, `not approved`, or `not decided`. A branch in the code without an `approved` entry is not an exception. None are approved.

A rule may later carry data that it is less likely, and a filter may read that. The schema is not decided. Do not invent a second field beside the trial ``search_bias`` on `PatternInfo` (soft heap demotion only — see Schema), and do not treat a missing bias as an exception.

## Migration

Status: done. Previous forest: `src/xenosite/_archive_forest/` (archive
directory on GitHub). Live public API is `xenosite.forest`. Prefer `mol.xf`
over `AtomTracker` — accessor API in `docs/forest/XF.md` (caches into
`_forest`; TypedDicts in `records.py` are unstable layout, not the public
contract). Callers: `docs/forest/MIGRATING_0.7.md`.

## Correctness, not yet done

- Do not drop a child only because `atom_diff` cost is not strictly lower. The `closer` check in `find_path` refuses a sideways step. Some real routes do not move the score down on every hop. Tried equal-cost sideways when the site touched cleavage / dearomatization / oxygen; PhCH2OH and multi-oxidation quinones blew up (`mol_edits` 21→68). Status: not decided — needs a tighter predicate than “touches the diff”.
- Do not walk a step plan that is the same steps in another order. Once a `Deps` has been yielded, a later walk that is only a reordering of those steps is not a new path. Rust `find_path` skips a hit when `plan.same_linearizations` matches an already-yielded outcome (`linearization_overlap` counts shared orders).
- **High `n_lin` is not fixed by inventing precedes.** `Deps::bind` only adds edges from `WillAdd` / `AddedBy` (e.g. OH ≺ DH that consumes that oxygen). Independent Dealkylations of distinct methoxys are chemically unordered — 8 steps with one OH≺DH edge → `n_lin = 8!/2 = 20160` is the free DAG, not a miss. Fake precedes between free dealks would lie. **WillAdd bind** matches only prep rules that supply the element (`Hydroxylation` / `OxidativeDehalogenation` for O) — not every earlier step whose site includes the anchor carbon (that overstated Dealk ≺ DH). Status: approved (do not invent order).
- **Multipath redundancy vs free linearizations.** Yield filter already drops reorderings of the *same* step multiset (exact Step identity). Unique-edit already emits one canonical site per class — do **not** widen `same_linearizations` to orbit/`same_site_class` while that holds. Plans that share early dealks but pick **different** (non-equivalent) OH carbons before quinone DH stay as alternate paths. Remapped Index notes on free dealk reorderings are a separate gap (`same_rule_maybe_skeleton` tried). Status: not decided for the remapped-index / dominated-extension levers.
- **Use yielded plans mid-search (novel_site soft-demote).** Removed — see [DROPPED](DROPPED.md). Yield still drops exact / skeleton twin plans; do **not** auto-abort from yield-drop counters.
- **Circular hydration ↔ dehydration.** Block only when site sets match: singletons use unique-edit / top orbits ([`same_site_orbit`]); multi-atom sites (beta-elim `site_map` includes adjacent C) compare by set equality — **allowed if at least one atom differs**. So OH → alcohol dehydration at the same carbon is refused; OH → beta-elim to ethene stays (adjacent C in the site). Counter: `blocked_circular_oxygen`. Test: `hydroxylation_then_dehydration_same_site_yields_ethene` + `alcohol_dehydration_after_hydroxylation_same_site_blocked`. Status: approved (Rust derisk).
- **Ancestor dedup_smi cycle kill.** Refuse any emission whose product `dedup_smi` equals an ancestor frame (`atom_trace["dedup_smi"][:depth]` in Python metabolize; Rust `Walk.ancestors` at find_path enqueue). Kills A→…→A and identity hops. Unstable (`None`) keys do not trigger. Distinct from global `seen` and from site unique-edit / yield CSMI layers. Status: approved.
- **Dehydrogenation heavy-neighbor match.** Pair ends whose effect is dearomatizing removes-H (not cleavage, not adds-O) require that each end atom's heavy neighbors on the **product** (after application) map onto exactly the target atom's heavy neighbors under **the same** MCS placement (not each end's best independent view). Check runs post-keep in `find_path`, not in `filter_sites` / `pair_could_help` before the edit. Parent end atoms map onto the product via forest labels/tags. DH keeps connectivity; mismatched neighbors mean the applied edit did not produce what the target needs. Reads effect shape + ends / mapping / target — not a rule-name branch. Status: approved. Tests: `test_cycle_and_dh_neighbors.py`.

- **Site help under one MCS placement.** `filter_sites` / `candidate_could_help` / `pair_could_help` score the **whole site** against each MCS mapping separately and keep the site if any one placement works. Do not take each atom's best match from the merged top-rank union across placements — that optimistic split disagrees with multi-atom topology (bonds and ResonancePair ends). Single-atom sites are the trivial case of the same rule. Pair-only data (per-end oxygen / methide `partner`) still reads `ends`; the topology gate is not pair-special. Status: approved (Rust derisk + Python).

- **Lift child MCS via extend (bare MCS when not cost 0).** After an edit, tag-lift the parent MCS onto the product and **always extend** where the gap is placeable (unmapped product heavies onto free same-element target atoms bonded to mapped neighbor images — same rules as add-placement / `n_extra`). **Keep the lift only if `field_cost == 0`** (guaranteed — skip MCS). **Otherwise do not Aut-chase or trust a non-zero lift — rematch with fresh bare MCS** (same door as expand [`atom_diff`](../../crates/xenosite-forest/src/atom_diff.rs)) and count `PathCounters::mcs_lift_rematch`. Tag-lift impossible → `mcs_lift_fallback` (expect zero). No double work of Aut-then-MCS. Correctness suite: `lift_mcs_correctness`. Status: approved (Rust derisk).
- **Fresh `atom_diff` alignment is bare MCS.** Placeable grow ([`mcs_extend`](../../crates/xenosite-forest/src/atom_diff.rs)) runs on **lift seeds** only (`extend_mapping_where_possible`); failed non-zero lifts rematch with fresh bare MCS. A process-wide expand opt-in (`USE_MCS_EXTEND` / `find_path_bench --mcs-extend`) and a PatternInfo-gated `FindPathConfig.mcs_extend` were tried and cut — full extend on every expand with PhaseOne+EpoxideHydration regressed hard filter bill **273→449** (veratrole 60→216). Status: **approved** bare expand. Tests: `epoxide_o_extends_onto_target_hydroxyl` / `mcs_align_is_bare_mcs` / `epoxide_toward_diol_drops_with_heavier_n_extra` / `ethene_epoxide_toward_glycol_needs_extend_views`.

- **Or of cleavages that share an identical fragment.** Distinct cleaving steps that discard the **same** fragment CSMI multiset fold into one [`CleavageOr`](../../crates/xenosite-forest/src/cleavage_graph.rs): arms share sorted `fragments` (both sides of the bifurcation are first-class). Fold key is [`CleaveFoldKey`](../../crates/xenosite-forest/src/pattern.rs) = resolved `cleave_side_group` signature + fragments (cross-rule). Choosing an arm **and** a continuation fragment builds linked [`Maybe`](../../crates/xenosite-forest/src/canonical_plan.rs) from the other side(s). `find_path` expand enqueues each `(fold_key, continue_csmi)` once. Not a revival of archive And/Or plan trees ([DROPPED](DROPPED.md)). Status: approved (Rust derisk).
- **Cleavage product graph → Or/Maybe, then non-cleavage across survivors.** Enumerate the cleavage layer (on demand / BFS); both fragment CSMIs are graph nodes. Fold arms by [`CleaveFoldKey`](../../crates/xenosite-forest/src/pattern.rs). **MCS / `atom_diff` still gates** which fragments enqueue. Expand gate: strict cost drop vs parent **and** ha ≥ target; if **both** sides match, **both** stay — no single-winner prune. Shrink (and local-add) diffs **tag-lift + extend**; keep only cost-0 lifts, else **MCS** (`mcs_lift_rematch`); also MCS when lift is impossible (`mcs_lift_fallback`, expect zero). Do **not** add a find_path cleave-only phase — treat cleaving edits like other reactions (`order_key` / filters). Redundancy across cleavage *sides* uses ``cleave_side_group`` data collected at expand — not a second predicate. Status: approved (Rust derisk; free-step leave-group yield key still not decided).
- **General product graph (all edits).** [`product_graph`](../../crates/xenosite-forest/src/product_graph.rs) BFS of PhaseOne metabolites as CSMI nodes with inbound [`ProductHop`](../../crates/xenosite-forest/src/product_graph.rs) edges (rule, pattern, site Tags, added Tags). Cleavage-only graph stays separate (Or fold). Not find_path; not archive NetworkX `MetaboliteNetwork` (DROPPED). Target MCS gate matches cleavage expand. Status: approved (Rust derisk door). Tests: `product_graph` unit; example `product_graph_bench`.
- **Metabolite enumerator (stream).** [`enumerate`](../../crates/xenosite-forest/src/enumerate.rs) / [`bfs`](../../crates/xenosite-forest/src/enumerate.rs) / [`dfs`](../../crates/xenosite-forest/src/enumerate.rs): pull-iterator of metabolites up to `max_depth` with accumulated [`PathInfo`](../../crates/xenosite-forest/src/enumerate.rs) (hops = rule, pattern, site Tags, added Tags). BFS queue vs DFS stack (Python `find_path.bfs` / `dfs` shape). Default [`EnumConfig::unique_csmi`](../../crates/xenosite-forest/src/enumerate.rs) emits each CSMI once; `unique_csmi=false` / [`with_all_paths`](../../crates/xenosite-forest/src/enumerate.rs) streams every path (still drops cycles on the current walk). Optional `max_nodes`. Expand via `product_layer` with no target gate. Not find_path; not the product graph materialization. Status: approved (Rust derisk door). Tests: `enumerate` unit.
- **Yield: remapped free-step twins (keep).** `same_rule_maybe_skeleton` (rule multiset + Maybe side CSMIs + precedes under rule-name align) runs beside exact `same_linearizations` at hit yield. Not orbit-aware — unique-edit canonical sites make that unnecessary. Counters: `PathCounters::dropped_duplicate_plan` (total wasted hits), split into `dropped_exact_plan` vs `dropped_skeleton_twin` (remapped twins only). Dominated-extension yield filter over-collapsed legitimate multipath (dimethylamine→PhCHO); dropped as a yield prune. `PathCounters::signal_contained_plan` still **counts** when a hit would match `dominates_extension_of` — signal only, does not drop. Status: **keep** skeleton twin yield-drop (ablation: off pads `max_paths` with remapped twins); dominated-extension as automatic drop/abort **not approved**.

- **`ApplyN` plan pool (MS1).** [`ApplyN`](../../crates/xenosite-forest/src/canonical_plan.rs) on [`Deps`](../../crates/xenosite-forest/src/canonical_plan.rs): OR of elementary rule-name arms + exact apply `count`. Composable with steps / precedes / [`Maybe`](../../crates/xenosite-forest/src/canonical_plan.rs) (same bag style as Maybe — not archive And/Or trees). Example: long hydroxylation list as arms, `count = N` when the spectrum expects N oxygenations. Pattern-level OR stays on PatternInfo inside a leaf. **Counts / dedup:** eligible sites = union of unique-edit orbits for matching arms; [`ApplyN::n_combinations`](../../crates/xenosite-forest/src/canonical_plan.rs) / [`unordered_site_combinations`](../../crates/xenosite-forest/src/orbits.rs) = Aut-deduped unordered `count`-subsets (benzene OH×2 → 3). **Emit:** [`apply_n_emit_products`](../../crates/xenosite-forest/src/canonical_plan.rs) one materialize per combo (sorted tags); each product carries free [`Deps`] plans whose linearizations cover path orders (`n_covering_linearizations`). Bench: `examples/apply_n_bench.rs`. Plan replay counts: [`Deps::n_linearizations`](../../crates/xenosite-forest/src/canonical_plan.rs) + [`Deps::n_distinct_products`](../../crates/xenosite-forest/src/canonical_plan.rs) / [`Deps::replay_stats`](../../crates/xenosite-forest/src/canonical_plan.rs). Status: approved (Rust derisk door). Tests: `apply_n_benzene_oh2_three_combinations_three_products`, `apply_n_emit_covers_paths_without_permuting_applies`.

- **`find_path_ms1` (m/z target).** Sibling of structure [`find_path`](../../crates/xenosite-forest/src/find_path.rs): hit when formula m/z is within `tol_da` **and** ApplyN pools are exactly filled. Default adduct **[M+H]⁺**. No structure MCS. Closer uses materialized mono mass ([`molecule_mono_mass`](../../crates/xenosite-forest/src/mass.rs): common isotope unless `atom.isotope` labeled). Catalog `delta_formula` soft pre-hint picks the better of full vs heavy-only (hydroxyl `+O −H` vs live `+O`); skipped for `cleaves` (leave bags ≠ fragment mono mass). EpoxideHydration / EpoxideOpening / N-oxidation nitroso bags mass-faithful; hydroxyl junction bags and S-ox hydroxy left soft. **Yield:** same as structure find_path — drop exact `same_linearizations`; drop remapped `same_rule_maybe_skeleton` only when hit CSMI matches (isobar regioisomers must not collapse). Path orders for one product are covering linearizations on the plan, not extra hits. **Property:** apply leaf → product CSMI must appear in MS1 hits at that product's m/z; every emitted hit's final (and last-step) product satisfies tol; path mz-error is monotonic; emitted plans **replay** from the reactant (`Deps::reaches`) — site notes are forest **labels** ([`Tag`](../../crates/xenosite-forest/src/labels.rs) / `Atom.tag`), not chematic indexes, so multi-hop linearizations rebind after earlier edits. Cleavage leave fragments exempt from the mz check (hit species must appear). ApplyN Aut-emit products at the target mass ⊆ MS1 hits under the same pool (covering). Tests: unit apply/hard chains + proptest [`ms1_apply_fuzz`](../../crates/xenosite-forest/tests/ms1_apply_fuzz.rs) (`XENOSITE_FUZZ_EXAMPLES`). Must not break structure find_path. Status: approved (Rust derisk door).

## Plans are Deps (elementary Steps + precedes)

Status: approved (Rust derisk; shape may replace Python `CanonicalStep`).

The plan language is flat [`Deps`](../../crates/xenosite-forest/src/canonical_plan.rs):
elementary `Step`s (rule name + site notes) plus precedes, plus [`Maybe`](../../crates/xenosite-forest/src/canonical_plan.rs)
cleavage bags and [`ApplyN`](../../crates/xenosite-forest/src/canonical_plan.rs) OR-pools
on the same plan (not siblings on the path outcome). There is no
parallel `CanonicalStep` dialect for search.

Site notes are one enum: known **label** ([`Tag`](../../crates/xenosite-forest/src/labels.rs) / `Atom.tag`) | `WillAdd(element @ anchor label)` | `AddedBy(rule, anchor labels)`.
Not chematic atom indexes — indexes shuffle across hops; labels survive so multi-hop plans replay.
Composite hops own a `canonical_plan` hook on the leaf (Python
`ReactionRule.canonical_plan`): it returns elementary `Step`s named after
catalog rules (`Hydroxylation`, `Dehydrogenation`, …). Search asks the leaf;
there is no `PlanKind` enum. Metabolize may still apply a composite leaf in one
hop; the plan is the elementary split for search and replay. QuinoneFormation
expands to Hydroxylation (or OxidativeDehalogenation) then Dehydrogenation.
EpoxideHydration expands to Epoxidation then EpoxideOpening (peers already in
PhaseOne). Status: approved (Rust derisk).

Replay: `Deps::linearizations` → each `Linearization::apply` runs named
elementary rules at resolved sites. Correctness: an accepted product’s plan
must reach that product under some linearization; a `find_path` hit’s plan
must reach the hit. No mid-plan re-metabolize of the composite leaf.

Plan identity: construction stores the transitive reduction via the Python
port of `transitive_closure_masks` / `canonical_dependency_edges` (bitmasks).
`Deps::==` is order-sensitive on steps; `same_linearizations` aligns by step
equality (greedy multiset) then compares reduced precedes — use that to treat
reorderings as the same path (HEURISTICS: do not re-walk a reordering).

`Addition.phase1` stays reserved until expansion schema is named on the
record. Cleaving quinone ends that need a Dealkylation prep are still a gap
(data not yet on the expansion).

## Speed, not yet done

- Put the isotope-bearing atom first in the reaction SMARTS. `SubstructMatch` starts at query atom 0 and does not reorder, so a common atom first is walked before the isotope is tested.
- **Rust door candidates (site–pattern–parent).** `RuleSet::candidates` yields deferred triples; `Candidate::materialize` runs the edit. `PairCandidate` is the ResonancePair analogue (merged `Effect` for filtering; path flip on materialize). `find_path` filters by reading pattern/effect data before materialize (`mol_edits` only counts survivors). Filter closures on `metabolize` stay for Python parity. Status: approved for the derisk door.
- **Pair products split like Python `split_fragments`.** `PairCandidate::materialize_mols` runs chematic `fragments()` after the path flip (same as SMIRKS `apply_smirks_at`). Cleaving pair ends (QF `dealkylate`) must not yield one disconnected `C.quinone` mol: `find_path` reads `n_products == 1 && cleaves` as ring-open and `n_products >= 2` as bifurcation (`PathStep.sides` / `Maybe`). Public doors are polymorphic: [`RuleSet::candidates`] / [`RuleSet::metabolites`] yield edit + ResonancePair under [`Candidate`] / [`Emission`] (materialize via `Candidate::materialize_mols` / `emit`). Pair-only helpers are `pub(crate)` — not underscore prefixes. Status: approved (Rust derisk + Python parity). Tests: `pair_edit::quinone_dealkylate_splits_fragments_like_find_path`, `find_path::quinone_formation_cleavage_records_side_and_maybe`, `enumerate::quinone_formation_cleavage_streams_both_fragments`, `python_suite_port::pair_cleave_splits_like_python_split_fragments`, `ruleset::pair_leaf_door_stamps_rule_path_bare_does_not`.
- Deferred metabolites. `metabolites` should be able to yield a light object instead of a sanitized mol. The object carries the anticipated formula (the same counts as `atom_trace["formula"]` and `delta_formula`) and says whether it is sure materialization can succeed. `metabolize` either returns that object or materializes it, so today's callers still get mols. A later hop materializes only when a filter or the plan needs the real mol. The same object is how a `StepPlan` can check that a later site was `added_by` an earlier step without replaying the edit.
- A richer walk-ranking key than target-hit + enqueue order. Status: approved (default). Default / current best is **`Match(close×improve-both)`** (formula × atom_diff close×improve); SoftStack and other `MatchScoreSpec` recipes remain opt-in (`--score`, `--matrix`). Plain max-heap pop. Trying absolute-cost sort before enqueue preferred locally cheap dead ends on PhCH2OH→quinone — do not revive that. **Tried `PatternInfo.search_bias`:** SoftStack-only soft i8. **Dropped:** alternate DFS/BFS among equal scores. **Tried switch to close-atom** (lowest single-path hard bill): not taken — multipath hard Σbill favors close×improve-both (3994 vs close-both 4224 / improve-both 4555).
- Pattern frequency. Count how often each pattern matches, and how often that match shows up in a real metabolite. Deprioritize or filter patterns that match almost everything and are rarely the one that was used. That number belongs on `PatternInfo` (related trial: `search_bias`), not in a table inside `find_path`.
- Cleaving edits are ordinary reactions in expand. `order_key` already prefers `cleaves` when the target is smaller; atom_diff / leave_count filters gate sites. Do not add a find_path cleave-only phase or another cleavage predicate until a richer walk key has been tried. (Tried `cleavage_first` expand gate — dropped; see DROPPED.md.)

## Already in the search

- Cleavage patterns run first when the target is smaller or a mapped bond is missing. Speed, not a new rule class. `order_key` also prefers dearomatizing spans when `loses_aromaticity`, and oxygen-adding spans when `needs_oxygen`.
- Oxygen is required only at atoms the diff says need it. Speed. Filters also see the live mol as their first argument (`FilterRules` / `FilterSites` take `ForestTracingMol`); a non-cleaving oxygen-only pattern is refused when the live oxygen count already meets the target. Mol is not stuffed into `SiteInfo`.
- Every full-size target placement stays open to filters and cost. `_mappings` / `_merge_views` keep each placement (not one union core; not only the smaller remainder). Tests: `test_multi_mcs_union_would_prune_both_benzyls`, `test_smaller_secondary_embedding_alone_misses_naphthaldehyde`.
- Ring-membership disagreement on a mapped atom adds the bridge into the ring to `cleavage_bonds`. MCS may keep an exocyclic atom on a ring image (PhCH2OH CH2 onto a quinone carbon); the cut is still that bridge. Reads the mapping, not a rule-name branch. Unblocks benzylic dealkylation on PhCH2OH→quinone.
- Target-hit early stop: a child that is already the target ranks via cost/formula 0 on the heap score (not a bool tier), and further site edits on that node stop once `max_paths` hits are queued (TBA).
- Two-site unique-edit identity is a pair-orbit signature from `graph_isomorphism` / `records` (beside the rank key, not in it): same-kind → ``AtomPairOrbitSignature`` / ``BondPairOrbitSignature``; DH ResonancePair → ``AtomPairOrbitSignature`` on the two end atoms. One-site rules still collapse by atom rank. **Nauty six-family API** (`all_site_pair_orbits_nauty`) is the source of truth: `atom_atom` / `bond_bond` each have unordered and ordered families; `atom_bond` ordered ≡ unordered (types already distinguish — do not invent a split). Same-kind unordered = combinations; ordered = permutations. **Ordered vs unordered for pair unique-edit** reads resolved ``swap_group`` (When → PatternInfo → ``name``; absent/None means use ``name``): unordered iff both ends share the same non-empty group; else ordered with ends sorted by ``PatternInfo.name``. Map-rank embeddings are part of the signature (dealkylate methyl vs benzyl). Status: approved. **QF vs DH vs H:** distinct endpoint ``name``s are distinct groups by default — no cross-role unordered; same-name pairs (`add_carbonyl_o`×2, `phenol_end`×2, Hydrogenation `path_end`×2) unordered. DH `phenol_end` / `amine_end` / `methide_end` stay ordered across roles despite shared `edit=single_to_double`. Explicit ``swap_group`` only when grouping differs from ``name``. Dehydrogenation declares ``site_kind = "atom_pair"`` / ``sites_on = "atom_pairs"`` (one-bond emits both endpoints). Unique-edit is atom–atom only; ``swap_group`` orders asymmetric ends. Site shape is class data (``site_kind``), not ``Generic[SiteT]``. Narrowing SMARTS cannot replace ``swap_group`` for these couples (HEURISTICS: not approved). **Opt-in ``canonical_emitted_sites``** (default off; explicit kwargs on metabolize / find_path / bfs / dfs only — no env): ``filter_sites`` sees **discovery**; after accept, chemistry + emitted ``site`` remapped to lex orbit rep (singleton atom and pair kinds); ``discovered_site`` escape hatch when they differ; product forest last-layer restamp only; lex-rep cache read from **parent** after product ``clear_structure``. Filters must not assume ``site`` == discovery. Status: approved. Classical group orbits on pairs are not a novelty claim — see ``PAIR_ORBITS.md`` References. Legacy `site_pair_orbits_nauty` maps to unordered same-kind + atom_bond. Tests: `test_pair_signatures.py`, `test_dh_atom_pair_unique_edit.py`, `test_graph_isomorphism.py`, `test_canonical_emitted_sites.py`, `test_site_kind.py`.
- Two dedup layers (do not conflate): (1) **Site unique-edit** — topological ranks + pair orbit signature; keeps site topology so equivalent embeddings are one edit before `RunReactants`. (2) **Product identity** — split into **check** vs **yield**. Keys are Chematic **fail-closed** `canonical_smiles_stable_key` (Rust `stable_csmi_key` / `ForestMol::stable_csmi_key`), not raw display CSMI — **can return null/`None`** (coupled E/Z, non-idempotent spelling). That is not a bug: do not unwrap or substitute `csmi`. `None` ⇒ do not CSMI-dedup that product (still explore / yield; `PathCounters::unstable_csmi_key`). Check (always while site unique-edit is on): after unique-edit survivors, each emission's frozenset of fragment stable keys + site ranks; a later emission under the same `(rule, PatternInfo.name | SMARTS)` with equal `S` and equal ranks → `SiteDeduplicationWarning` (unique-edit miss). Check does not drop. Yield (`unique_csmi`, default on): drop duplicate `(rule, pattern, emission frozenset)` emissions from the stream — same `S` as check; within one cleavage list, identical sibling CSMIs stay. `unique_csmi=False` yields every emission after site unique-edit; check may still warn. `RuleSet.metabolize` forces children `unique_csmi=False` (yield off, check on) so alternate rules bubble up; the RuleSet's own yield (caller setting) is the cross-child CSMI layer — nested sets do the same, so only each asked-to-uniquify level drops. Cross-rule same emission frozenset → INFO log (kept/dropped rule + shared CSMI), not a Warning; parent does not re-fire same-rule `SiteDeduplicationWarning`. Inventory: ``REDUNDANT_RULES.md``. Status: approved (stable-key fail-closed: Rust derisk, not decided vs Python still on display CSMI).
- **Pass down unique-edit orbits on the site bag** (Rust derisk). Unique-edit still emits one canonical/representative site. The collapsed primary-map orbit rides on `SiteInfo.orbit` / `Candidate.orbit` / `Emission.site_orbit` / `PathStep.site_orbit` so (1) multiplicity is visible (benzene hydroxylation orbit len 6) and (2) plan / path steps can test site-class equivalence via `same_site_orbit` without re-emitting every embedding. Higher-order plan orbits (GraphTree) read these; full plan-automorphism collapse is not decided. Status: approved for passing orbits; plan-orbit enumeration not decided.
- **Step plans carry generator-based site orbits** (Rust derisk). `Step.orbit` is the automorphism orbit of the site atoms (closure under atom+bond generators), filled at emit via ForestMol-cached gens when present. `pair_edit` site-pair signatures reuse one generator run per `pair_candidates` call (`atom_pair_orbit_id_with_gens`). Distinct from unique-edit SMARTS-class orbits on `SiteInfo`. Status: approved.
- Pattern SMARTS that implement the same edit must not nest applicability ranges (H-count nesting that double-emits under distinct names). Partition like Dealkylation `#6H0` vs `#6h`, Hydroxylation `#6h1` vs `#6h2,#6h3`, Glutathionation epoxide/aziridine `#6H0` for the substituted-carbon row. Residual multi-name same-csmi cases that are chemically distinct edits are listed under DIVERGENCES.md. Status: approved for partitions done; residuals documented.
- `find_path` and `ReactionRule.metabolize` / `RuleSet.metabolize` raise `ValueError` on `None` input (invalid parse must not soft-pass as any-path).
- Lazy heap / stale-priority handling on the find_path frontier. Status: approved (shape).

  **Current default — `HeapScoreMode::Match(MatchScoreSpec::log_neg_pc())`.** Heap key is score then `seq` (with opt-in diversity: `(score + −n, score)`). Base score = `Σ −ln(1+child) − ln(1+parent)` over atom+formula. **Diversity** (`FindPathConfig.diversity`, default **off**): count accepted pops per key `(rule, site Tags, added-atom Tags)`; primary += fixed-point `−n`. On pop, recompute; if stale, re-push (`diversity_repush`); else accept and increment. SoftStack ignores. Ablate: `diversity`. Hard p5: off **1266** / −n on **1588**; p10 nearly tied. Status: **approved** log-neg-pc default; diversity opt-in **not decided**. Prior close×improve-both / Product matrix (below) remains the ablation baseline (`--score product-both` / `close×improve-both`).

  **Match-score matrix (filter-only, plain pop; mid + hard).** Axes: combine = close | improve | close×improve; metric = atom | formula | both. SoftStack included as baseline. Cost = MCS HA gaps only; `formula_l1` includes H (no special `h_off`).

  | score | mid Σbill | hard Σbill | hard s |
  | --- | ---: | ---: | ---: |
  | close-atom | 79 | **241** | **0.97** |
  | close-formula | 79 | 250 | 0.93 |
  | close×improve-atom | 79 | 255 | 1.05 |
  | close×improve-both (current best) | 79 | 260 | 0.99 |
  | close-both | 79 | 260 | 0.99 |
  | soft | 79 | 354 | 1.35 |
  | improve-atom | 79 | 354 | 1.35 |
  | close×improve-formula | 79 | 376 | 1.33 |
  | improve-formula | 79 | 474 | 1.91 |
  | improve-both | 79 | 494 | 1.93 |

  Read: mid bills collapsed (all 79; MeOPhOH edits up). Hard single-path ranks **close-atom** lowest bill; close×improve-both ahead of SoftStack. No misses. Multipath (max_paths=5): hard Σbill **close×improve-both 3994** < close-both 4224 < improve-both 4555 — keep close×improve-both as default. Vs prior HA+H cost (`|h_delta≠0|` in `field_cost`): hard close×improve-both was 217 / close-atom 210 — simplification costs ~40–50 hard bill. Switch to close-atom: **not approved** (multipath loss).

  **Current — aromatic |Δ| in atom cost; n_extra weighs more than cleavage.** `field_cost = 3·(|cleaved|+|cleavage_bonds|) + 4·n_extra + |aromatic_delta|`. Unmapped target heavies cost more than a cleavage bond so ethene→epoxide toward a diol drops (8→7) instead of staying lateral at uniform ×3. Per mapped atom, aromatic is 0/1; mismatch joins `aromatic_delta` and costs 1. H is not a cost term (`formula_l1` / on-demand `atom_h_delta`). Status: not decided (n_extra weight trial; shells stay full n0/n1/n2). Test: `epoxide_toward_diol_drops_with_heavier_n_extra`.

  **Trial — `MoleculeShells` / `AlignedShells` (n0/n1/n2).** Every heavy atom stores `aromatic` plus element bags `n0`/`n1`/`n2` that **include H**. Align → deltas + `unaligned_*`. **Gate bench baseline** (best so far on mid/hard `cost>0` vs live gate): legacy `at_sites` Σ|δ| on kept atoms only (no cleaved count) — mid P/R/A 0.07/1.00/0.48, hard 0.05/1.00/0.22. **`|current−target|` (δ=0) as a score was a bug** — removed from the ablation. Cost is always Σ|projected−target| with `projected = current + editδ` (δ = product−reactant). Keep when residual drops vs the no-edit residual (progress baseline only). Optional `|Δaromatic|` ∈ {0,1} per atom (aromatic-bit mismatch; off by default) and `+leave` are still under trial on that projected cost. Status: not decided.

  **Tried — HA+H in both atom cost and `formula_l1`.** `field_cost` included `|h_delta≠0|`. Better hard bills (close×improve-both 217). Status: not approved (H as special cost term).

  **Previous — `HeapScoreMode::SoftStack` (opt-in `--score soft`).** Lexicographic soft key: `search_bias` → `site_progress` → `cost_gain` → `seq`. Kept for comparison. Status: not decided (superseded as default). Tests: `heap_prefers_higher_search_bias_over_seq` / `heap_lack_of_improvement_counters_dfs` / `hop_cost_gain_is_parent_minus_child` / `heap_pops_best_ord_value_only` / `match_product_prefers_joint_improvement_and_closeness` / `match_combine_and_metric_axes`.

  **`delta_formula` is not a heap score** — soft mismatch warn/counter only. **Dropped:** alternate DFS/BFS among equal scores.
- **Adds-H site partners + H-progress score (Rust).** OxygenReduction carbonyl sites the heteroatom (`site_map=[1]`); H change is on the double-bond partner carbon. `scope_could_help` / expand `order_key` / heap `site_progress` extend scope with those partners when a live mol is present so adds-H helps only where the target needs more H (undo would not). Soft only after the gate. Status: not decided (Rust derisk). Tests: `oxygen_reduction_carbonyl_partners_allow_toward_alcohol` / `oxygen_reduction_refused_when_undoing_toward_carbonyl` / `order_key_h_progress_prefers_adds_h_toward_alcohol` / `heap_prefers_higher_site_progress_over_seq`.
- `leave_count` on a cleaving effect is read by `filter_sites`: when it is an int, the smaller fragment across the cleaved bond must have that many heavy atoms. Test: `test_leave_count_one_refuses_a_larger_leaving_fragment`.
- A nitrogen with two double bonds is not a product. `C=[N+]=C` sanitizes (degree 2, valence 4) but it is not an iminium. An iminium is one double bond and two single bonds. `_sanitize_piece` drops that fragment, and the split fails with it. Carbamazepine site `{4, 17}` was iminium on one ring carbon plus dealkylation of the carboxamide on the other, which left `C1=c2ccccc2=[N+]=c2ccccc2=C1` and `NC=O`. Status: approved. Test: `test_carbamazepine_does_not_emit_two_double_nitrogen`.
- A dearomatizing pair edit must change the system that was kekulized. Other aromatic systems may stay. If every aromatic atom of that system is still aromatic and no non-aromatic double or triple bond touches it, the product re-aromatized and is dropped. Alprazolam site `{1, 4}` is that cation: `Cc1nnc2cnc(-c3ccccc3)c3cc(Cl)ccc3[n+]1-2`. A quinone, quinone-imine, or quinodimethane keeps a localized double bond, so it stays. Status: approved. Test: `test_alprazolam_re_aromatized_cation_is_not_a_quinone`.
- Methide is always on the rule data (Dehydrogenation `methide_end`, QuinoneFormation alkyl branch). Opt-in `pathways=("methide",)` is dropped (DROPPED.md). Both ends may be methide: a para-quinodimethane is two alkyl single-to-double ends. The one-side skip was a mistake. There is no separate "no methides" mode. `merge_effects` sets `methide` when either end has it. Search already skips alkyl `partner=="C"` ends unless `_alkyl_bond_raises`. Tests: `test_dh_methide_pathways.py`.
- **Adds-H vs removes-H filters (Hydrogenation ≠ Dehydrogenation).** `filter_sites` already refused removes-H effects unless some site atom has `atom_h_delta < 0` (or intersects `loses_aromaticity`). Symmetric: effects whose `adds` contains H (and that do not also add O / cleave) are refused unless some site atom has `atom_h_delta > 0`. `filter_rules` mirrors this on span.adds. **Hydrogenation adds H** (reduction / saturation); **Dehydrogenation removes H**. Do not conflate. `dearomatizes` alone is not enough toward a quinone target: both H and DH can dearomatize, but only DH/QF match oxidative H loss. Status: approved. Tests: `test_hydrogenation_pattern_info.py`.
- **Hydrogenation `path_end` PatternInfo.** Declares `adds="H"` and `dearomatizes=True` (capability); `merge_effects` resolves dearomatizes against `system_aromatic`. Alkene SMARTS keep span `dearomatizes=False` so aliphatic C=C→CC is not refused by the pattern-level “all dearomatizes” check when `loses_aromaticity` is empty. The same check **skips** patterns whose span also adds H or O (path_end / Epoxidation): those resolve at the site / adds-H gate, so aliphatic carbonyl path reduction (`CC=O`→`CCO`) and ethene→epoxide (`C=C`→`C1OC1`) are not refused. Status: approved.
- **Epoxidation / EpoxideHydration `dearomatizes` capability.** Catalog `Effect.dearomatizes=true` (aromatic epoxidation / diol look-ahead clears the ring bit at the site). [`PatternInfo::resolve_for_match`](../../crates/xenosite-forest/src/pattern.rs) / Python `resolve_effect` resolve false on aliphatic site_map atoms (context mol, not Kekulé form). Site-shell residual reads the resolved bit for `|Δaromatic|`. Same resolve shape as Hydrogenation / pair `merge_dearomatizes`. Python `Epoxidation` span matches Rust. Status: approved (Rust derisk + Python parity). Tests: `epoxidation_dearomatizes_capability_resolves_on_aromatic_site` / `epoxidation_aromatic_residual_drops_with_dearomatic` / `epoxide_hydration_pattern_info_and_plan` / `test_epoxidation_pattern_info.py` / always-on [`pattern_info_catalog`](../../crates/xenosite-forest/src/pattern_info_catalog.rs).
- **Always-on PatternInfo catalog audit.** Part of `cargo test -p xenosite-forest` **and** `cargo test -p xenosite-forest --lib` ([`pattern_info_catalog`](../../crates/xenosite-forest/src/pattern_info_catalog.rs)): every `catalog_names` leaf — including `EpoxideHydration` — unique names, sealed `delta_formula`, site_map/kind shape, edit present; chemistry probe that kept-site aromatic→non-aromatic edits declare `dearomatizes` capability; resolve check on aromatic vs aliphatic probes. Status: approved (Rust derisk).

- **Site-shell leave matches the cleavage gate.** Cleaving residuals expand via [`site_atoms_with_leave`](../../crates/xenosite-forest/src/matched_atom.rs) preferring SMARTS-mapped partners (not MCS `cleavage_bonds` noise). Leave heavies are unmatched debt until cleaved ([`site_shell_cost_leave`](../../crates/xenosite-forest/src/matched_atom.rs)) — no MCS target mate, so phenol-OH→quinone no longer treats the leaving O as a needed map. Status: approved (Rust derisk). Tests: `dehydration_alcohol_residual_drops_with_leave` / `site_atoms_with_leave_includes_methyl`.

- **find_path ablation (`find_path_ablate`).** One-at-a-time vs baseline (eager, **log-neg-pc**, skeleton on, PhaseOne). Mid+hard, paths=1 and paths=5. No misses. Novel-site / `target_hit` heap tiers **removed** (DROPPED). Soft wins hard p5 bill, loses hard p1 — keep SoftStack opt-in.

  | variant | mid p1 Δbill | hard p1 Δbill | mid p5 Δbill | hard p5 Δbill |
  | --- | ---: | ---: | ---: | ---: |
  | lazy | 0 | 0 | −22 | −116 |
  | soft | −7 | **+142** | 0 | **−1140** |
  | close-atom | +7 | +28 | 0 | −830 |
  | with-novel (removed) | 0 | 0 | 0 | **+997** |
  | no-skeleton | 0 | 0 | −117 | −821 |
  | no-epox-hyd | −24 | −26 | −244 | −264 |

  Read: SoftStack hurts hard p1, wins hard p5. Novel demotion burned hard p5 bill (~2×) — **dropped**. Skeleton-off pads `max_paths` with remapped twins — keep skeleton drop on. EpoxideHydration taxes every suite without helping these quinone/dealk cases. Status: not decided for SoftStack default / EpoxideHydration PhaseOne tax.

  **Propose drop / simplify:**
  1. **novel_site soft-demote** — **dropped** (code removed).
  2. **`target_hit` heap / score sentinel** — **dropped**; atom_diff / formula closeness cover hits.
  3. **SoftStack as live path** — keep opt-in only; do not revive as default.
  4. **lazy_closer** — leave default eager; near-neutral.
  5. **skeleton twin yield-drop** — **keep** (on). Watch `dropped_skeleton_twin`.
  6. **EpoxideHydration in PhaseOne** — **keep for now**; record the tax.

  Example: `cargo run -p xenosite-forest --example find_path_ablate --release -- --hard --paths 5`.

## Schema in PatternInfo (approved)

- **``swap_group: str``** on `PatternInfo` (optional When override). Pair unique-edit resolves as When → PatternInfo → ``name``. Same non-empty resolved group on both ends → unordered; unequal → ordered by ``name``. Omit the field when it would equal ``name``; set it only when grouping should differ. Status: approved. Tests: `test_pair_signatures.py`.

- **``cleave_side_group: (leave, keep)``** on `PatternInfo` (Rust: `Option<(String, String)>`). Cleavage analogue of ``swap_group``, but pooled **across rules** at expand / `cleavage_layer` — not within one ResonancePair unique-edit. `None` = ungrouped (fold by fragment multiset only). Equal non-empty strings = swappable sides. Shared leave labels (e.g. ``Me`` on O- and N-methyl dealks) fold distinct rules when fragments match. `find_path` enqueues each `(CleaveFoldKey, continue_csmi)` once. Defaults unequal only when the field is set to unequal labels; omit the field when ungrouped is correct. Status: approved (Rust derisk). Tests: `cleavage_graph` fold / side-sig unit tests. Name: tuple aligned to ``site_map`` (leave = map 1). Free-step yield equivalence for independent same-leave sites still not decided.

- **``search_bias: i8``** on `PatternInfo` (default `0`). Soft heap preference under **`HeapScoreMode::SoftStack`** (not used by default `Match(log-neg-pc)`). Higher pops sooner; negative demotes. **Hydrogenation** and **OxygenReduction** use `-1` (`demote_reductive`). Never drops or aborts. Pair emissions take `min(left, right)`. Status: not decided — schema for “less likely” still open; SoftStack-only score knob. Tests: `heap_prefers_higher_search_bias_over_seq` / `reductive_patterns_carry_negative_search_bias`.

- **``delta_formula: {element: delta}``** on `Effect` / possibility (zeros omitted). Declared net formula change, sealed from junction ``adds`` / ``removes`` bags **minus** ``leave_formula``. ForestMol caches formula as element→count (incl. H). Cleavage: named leave as negative counts, plus O/H at the cut from the bags. When OR arms disagree on the delta (halogen removal), each arm carries its own map under a ``When`` — do not invent a search branch. **Mismatch vs materialized product:** soft warn (`FormulaDeltaMismatchWarning` / `log::warn!`), append structured ``FormulaDeltaMismatch`` to the counters list / suite collector, and increment ``PathCounters.formula_delta_mismatch`` — never drops. Check runs in `ReactionRule.metabolize` **and** `find_path` (metabolites bypass metabolize). **Suite gate:** autouse collector must stay empty for **non-pair** emissions (`tests/conftest.py`); mark intentional mismatch tests ``allow_formula_delta_mismatch``. ResonancePair (`ends` / ``pair=True``) mismatches are still recorded but do not fail the suite until one-placement site scoring is fully trusted for sealed end bags. Prefer the list over catching warnings / parsing message strings. ``removes_partner``: cleaving → ``leave_formula``; non-cleaving → junction ``removes``. Soft-skips: open leave; ring-retained leave (cleaves + leave on a single product that kept those atoms). Status: approved (Rust derisk + Python annotate). Redundancy with existing keys: see below.

### Formula vs existing Effect keys (audit)

| Key | Redundant with ``delta_formula``? | Notes |
| --- | --- | --- |
| ``adds`` / ``removes`` | **Yes** (bag encoding of the junction part) | Positive ≡ adds bag; negative ≡ removes. Keep bags until callers migrate. Note: sealed pattern delta is usually the **net** mol change; hydroxyl is the documented exception (edit stoichiometry ``O:+1,H:-1`` vs sanitized ``O:+1`` — OH restores H). EpoxideHydration ``OOHH``, EpoxideOpening hydrate ``OHH`` / rearrange ``HH``, N-ox nitroso ``O``/``HH`` match sanitized nets (were undercounting H). S-ox hydroxy stays ``O`` — net H is substrate-dependent (thioether vs thiol). |
| ``leave_formula`` | **Yes** (negative part of the sealed map) | Positive counts of the named leave; sealed as negatives. ``Me`` → ``{C:1,H:3}``. Open leaves (`leave_count` None) stay empty. |
| ``leave_count`` | **Partial** | When ``leave_formula`` is set, equals its heavy-atom count. Still useful when leave is open (None) or as a cheap filter without walking the map. |
| ``partner`` / ``partner_h`` | **No** | Site role / identity, not stoichiometry (though When+delta covers partner *removal*). |
| ``needs`` | **No** | Capability / gating (e.g. needs O), not a count. |
| ``cleaves`` / ``breaks_ring`` / ``dearomatizes`` / ``methide`` | **No** | Structural / pathway bits. |
| ``symbol`` / ``h`` / ``site_aromatic`` | **No** | Filled from the matched atom at resolve time. |

## Schema proposals (not yet in PatternInfo)

- **when → name.** Optional: a resolved `When` refines the pattern label so a globbed SMARTS OR (e.g. `[#6h2,#6h3]`) can still report `h2` vs `h3` without a second overlapping pattern. Candidates: `When` grows optional `name=`; or `PatternInfo.name` is a base and the chosen branch supplies a suffix. `_unique_csmi_key` / trace `pattern` would need a defined resolved token (static `PatternInfo.name` vs effect-time name). Prefer data on `When` / `Effect`, not a code branch in metabolize. Status: not decided — do not invent the field until a caller needs branch-specific names after partition. Test hook ready: `emitable_pattern_names` / `optional_when_name` in `tests/forest/pattern_info_inventory.py` already fold an optional `When["name"]` into per-rule uniqueness (`test_emitable_names_unique_within_each_reaction_rule`).

## Dedup check vs yield (approved)

``unique_csmi`` controls **yield** only (suppress duplicate *emissions* whose
fragment-CSMI frozenset already appeared under the same rule/pattern). **Check**
(emission frozenset of product CSMIs + site ranks → ``SiteDeduplicationWarning``
on unique-edit miss) runs whenever site unique-edit is on — including with yield
off, so a RuleSet can force children ``unique_csmi=False`` for cross-rule
bubble-up without losing within-rule miss detection. Site unique-edit has no
public off switch today; if one lands, check-off with it. ``product.xf.csmi``
stays cached after first read. Emission identity for check/yield is
``frozenset(p.xf.csmi for p in products)`` computed internally — not a
``csmi`` field on ``info`` (callers use ``xf.csmi``). Status: approved.

## Site kind taxonomy (approved)

`RuleSiteKind`: ``"atom"`` | ``"bond"`` | ``"directed_bond"`` | ``"atom_pair"``.

- **``atom``** — singleton Site; unique-edit uses directed `MapRankKey` `((mapno, rank), …)`.
- **``bond``** — undirected bond frozenset (two adjacent atoms). Unique-edit first field is sorted site ranks (`bond_rank_key`). Epoxidation, AzoSplitting. ResonancePair **SMARTS** (Hydrogenation alkene/alkyne, Dehydrogenation one-bond) also use `bond_rank_key` via `site_signature` when `site_kind="atom_pair"` — pair-path emissions still use `pair_site_signature`.
- **``directed_bond``** — unique-edit uses ordered map ranks (`MapRankKey` / ``site_map`` order) because map 1 is chemically distinct (Dealkylation / NDealkylation / Benzodioxole / Nitroaromatic). **Public API:** ``info["site"]`` is always a frozenset; orientation is on ``info["discovered_site"]`` as the ordered map-order tuple (same ``Site`` union — no extra field). With ``canonical_emitted_sites``, ``site`` is the frozenset of the lex representative and ``discovered_site`` remains the directed discovery tuple. Coercion is yield-only — unique-edit `seen` never keys on the frozenset.
- **``atom_pair``** — ResonancePair ends only (DH / QF / Hydrogenation / TautomerRule stub). Never on plain SMARTS. Pair unique-edit stays `pair_site_signature` + `swap_group`. One-bond SMARTS on these rules share undirected bond ranks with ``bond`` (see above).

Status: approved. Tests: `test_epoxidation_unique_edit.py`, `test_site_kind.py`, `test_parity.py` (anisole dealk).

## Product-identical distinct sites / Dealk other-side double-emit

Status: **not decided** (emission shape); tracking in [RUST_PYTHON_PARITY.md](RUST_PYTHON_PARITY.md) **C13**.

Automorphism-equivalent sites already emit once with ``orbit`` / ``site_orbit``.
Ester Dealkylation (two directed ``quaternary_alcohol`` sites on the bridging
O, unequal ranks, one product bag) and similar cases are a different class.
``unique_csmi_compliant=False`` on those leaves until data (``product_equiv`` /
identity maps) folds them — only the CSMI-dup test xfails (C11).

## ``unique_csmi_compliant`` (rule data)

Status: **approved**. Class attribute on ``ReactionRule`` (default ``True``).
Yield-layer CSMI dedup runs only when compliant. Goal: all leaves ``True``.
Non-compliant: ``Dealkylation``, ``OxidativeDehalogenation`` (corpus hits).
Tests: ``test_csmi_dedup_soft_failure`` only — hard-fail if compliant + hit;
xfail if non-compliant + hit. Other suites do not xfail on this flag.

## Antipattern: global `_kekule_forms` on SmirksReactionRule

Status: **not approved**. Materializing every Kekulé form inside plain `SmirksReactionRule.metabolites` is combinatorial expansion — that is why `ResonanceRule` exists. Match once on the aromatic parent; react on a cached Kekulé parent selected by SMARTS-implied bond order. Do not reintroduce an all-forms loop or an easy opt-in that restores the tax. `_kekule_forms` may remain for tests / helpers that need the list explicitly.

## Narrowing SMARTS vs ``swap_group`` (ResonancePair)

Status: **not approved** as a replacement for ``swap_group`` / name-default groups. Hydroxylation-style H-count partitions fix *nested same-atom* SMARTS that double-emit under ``unique_csmi``. Pair same-role couples are different: two path ends that both match the same edit (hydroquinone ``phenol_end``×2, benzene ``add_carbonyl_o``×2, diene ``path_end``×2, QF ``dealkylate``×2, …) are real chemistry and stay unordered via resolved group (= ``name``). Narrowing SMARTS so those couples never co-apply would drop pathways. Keep resolved groups + nauty ordered/unordered. Residual: map ranks still distinguish dealkylate embeddings.

## TautomerRule (stub only)

Status: not decided for chemistry / unique-edit; stub lands first.

Archived ``Tautomerization`` extended alternating paths by one H-bearing
neighbor and flipped bonds (net heavy-atom formula and H count unchanged;
path swap, not ``RunReactants``). Live ``TautomerRule`` subclasses
``ResonancePairRule``, is patternless, raises ``NotImplementedError`` from
``metabolites``, and is not in PhaseOne. Orthogonal deferred work: tautomer
*SMARTS matching* via RDKit ``TautomerQuery`` (preferred direction only;
Status: not decided). When the rule is implemented, put ends/effects on
``PatternInfo`` rather than a silent search branch. Docstring on the class
and TODO.md carry the same split. Test: `test_tautomer_rule_stub.py`.
