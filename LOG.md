# Lab log

## 2026-09-17

- `AtomRef` + `Step`/`Linearization.apply`: deferred sites (new O after hydroxylation) resolve via `mol._forest["atom_refs"]` (private schema documented next to `_forest_state`); not AtomTracker tags. Quinone DH ends use `added_by` refs. `toward=` keeps fragments when applying prep prefixes.
- Phase1-equivalent steps (`Step` / `StepPlan`): uniform `phase1_steps(mol, site)`; Phase I + NDealkylation degenerate singletons; QuinoneFormation SMARTS→prep layers + final DH; `attach_phase1_steps` stamps mol prop. `RenumberAtoms` during tag/align drops mol props — copy props across renumber so stamps survive.
- Fuzz uses library apply only. Aromatic Forest `Dehydrogenation` often still cannot fire; prep-order agreement + full match when steps succeed.
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
