# TODO

Parity gap progress gauge (open / in progress / done):
`docs/forest/RUST_PYTHON_PARITY.md` § Work order (also mirrored in agent store
`parity-gaps.md`). Flip Status there; do not delete rows.

- find_path: `seen` blocks re-enqueue by child CSMI but chemistry (`mol_edits`) still runs first — MeOPhOH residual bill≈176 is mostly Dealkylation re-attempts into already-seen products (see PERFORMANCE / meophoh_retrace)
- Expand each rule's `_example_substrates` so examples cover all patterns/whens (meta-test currently only checks site_kind shape)

- Later (optional): hardcode `canonical_emitted_sites=True` on flaky forest tests (topo / guided gold / path hard / unique-edit) if lex site identity helps — fuzz already draws the flip via Hypothesis
- PatternInfo coverage gaps: product-side SMARTS parse; `pin` maps present in reactant; `skip_same_rings` / `edit` checked beyond presence; `breaks_ring` and `partner`/`partner_h` filled by `resolve_effect` (not only declared survival)
- `find_path_to_MS1`: condensed StepPlan-like form of reactions + possible sites to reach an MS1 m/z (z/charge) within tolerance; isotope-aware; all in-tolerance resolutions; alongside PhaseOne + conjugation, rules for common MS changes and MS-relevant input adjustments (e.g. C12→C14, add H+)
- After MS1 is solid/tested: `find_path_to_MS2` same idea plus MS2 fragmentation rules to develop
- Optimize `alternating_paths`: quadratic because each BFS step copies the path; odd rings can leave the search (parity to the other side is fixed by whether the path skips parity), and those decomposed subproblems are then 2-colorable so all-pairs odd-distance BFS per anchor is not required
- xf: hang remaining helpers still free (`resonance_bond_maps`, MCS pairwise, `carry_forest` / split, `_reordered_forest_labels`)
- Deferred: tautomer SMARTS matching / TautomerQuery (preferred direction; Status: not decided) — orthogonal to landed `TautomerRule` stub (`NotImplementedError`; see docs/forest/HEURISTICS)
- Later: SMARTS least-common atoms first
- forest: phase I names still come from the search, not from `atom_trace["additions"][id]["phase1"]` (the field is reserved; the id, rules, effect, and formula delta are stored)
- Guided `and_cleave_plan` SMARTS: optional per-SMARTS flags naming atoms made dirty on match, so match-before/after deps are less leaky (fewer false precedes; less rely on apply-replay)
- Rank ``PathOutcome`` hits by likelihood via ``xenosite-predict`` using cleavage-side bags (without enumerating side reactions) — not this pass
- Guided-search heuristics from fuzz cases that hit ``max_expansions`` / miss under budget
- Use ``sanitize_dropped`` to find remaining rules/SMARTS that still clean-fail after include_sites / QF plan-only enum; further pre-clean rejects if needed
- Short-circuit mol edits that would yield invalid structures before sanitize/clean (edits are costly; reject early)
- Speed up `clean` / `_sanitize_kekulize` (dominant Full wall cost: sanitize + kekulize + SMILES round-trip per product)

## After Rust↔Python leaf parity (execute once parity suite is green)

Do **not** start these while work order #4–12 / parametric product parity is still
red. Tracked also in `docs/forest/RUST_PYTHON_PARITY.md` work order #16–17.

1. **Chematic atom tracking (upstream landed).** Migrate off the vendored
   chematic patch (`vendor/chematic` +
   `patches/chematic-v1.0.21-atom-tag-visit-order.patch` /
   `./scripts/vendor-chematic.sh`) onto the released Chematic atom-tag /
   visit-order API. Pass atom-tracking tests (`AtomTracker` /
   `stamp` / `src_to_new` / `adopt_born` / `write_parse` and successors).
   Remove the vendored submodule and patch once green. See
   `docs/forest/RUST.md` § Atom identity.

2. **Centralize ForestMol cache access and copy/edit.** Every read/write of
   the forest mol cache (structure / ranks / systems / resonance / …) and
   every copy/edit path that must invalidate or preserve cache correctly
   goes through `ForestMol` methods — Python `mol.xf` / Rust `ForestMol`
   alike. Outside accessors need a compelling reason. Temporarily make the
   cache private (or rename) to surface stray mutators/readers, then restore
   a **user-visible** cache surface for Python scripting (not hidden from
   interactive use — just not mutated ad hoc from free functions).
