# Dropped behavior

The old library under `src/xenosite/forest/` has these. This search does not.
Each entry names where it lived, why it was left out, and which old test
depended on it. Do not delete entries. Add one in the commit that leaves
the behavior out.

## PathContext and its MCS bookkeeping

Lived in `src/xenosite/forest/path_context.py`. Multi-embedding rounds,
chemical disagreement, frontier versus deep interior, safe-to-drop
thresholds, and oxygen slack. The atom diff in `find_path.py` is the part
search needs: one pairing, then the local change at each atom.

Old tests: `tests/test_cleavage_site_prune.py`, `tests/test_guided_path_fuzz.py`
(the context object, not the metabolites).

## Formula hint algebra

`FormulaHint`, `FormulaAny`, `FormulaMatch`, and `expand_formula_effects`
in `path_context.py`. Local atom differences cover compatibility. The spec
language does not.

Old tests: `tests/test_formula_smarts_hints.py`.

## and_cleave_plan SMARTS replay

`and_cleave_plan` and the apply-replay correction in `guided_path.py`.
Precedence comes from `added_by` and the resolved `Effect`.

Old tests: `tests/test_and_cleave_plan_fuzz.py`,
`tests/test_guided_path_cleavage.py`.

## And / Or plan trees and plan JSON

`And` / `Or` in `step_plan.py`, plan JSON canonicalization, and exact
linearization matching against the old emission. The poc yields a `Deps`.
`Step` and `Deps` are imported, not reimplemented. Callers do not get the
old JSON shape.

Old tests: `tests/test_step_plan.py`, `tests/test_phase1_steps_fuzz.py`
(the JSON and `And`/`Or` assertions, not the chemistry).

## AtomTracker as a base class, and phase1_steps on each rule

`AtomTracker` in `base.py`, and `phase1_steps` on each rule in
`forest/rules.py`. The poc rules do not subclass `AtomTracker`. A composite
such as quinone formation is a `Deps` this search builds from `added_by`
and the resolved effect, not a method on the rule.

Old tests: `tests/test_phase1_steps_fuzz.py`, `tests/test_atom_trace_char.py`.

## Phase-one reconstruction module

`phaseone.py` (`PhaseOneRS`, `PhaseOneQF`). A ruleset is still a rule: it
holds the classes and running it runs each of them. Phase I is the `Deps`
this search yields, not a second module that rebuilds steps from SMARTS.

Old tests: `tests/test_phaseone.py`.

## Search controls that only bound the old expansion

`max_expansions` as a public mode, depth counters, and BFS versus DFS
metabolite enumeration. A billed counter stays, as a measurement
(`mol_edits + nodes`), not as a second search mode.

Old tests: `tests/test_classic_path_counters.py`, `tests/test_bfs_fuzz.py`,
`tests/test_path_outcome_hard.py` (the budget flag, not the metabolite).

## Conjugate star-label modes and the dehydrogenation methide opt-in

Bare `*` versus a full glucuronide or glutathione structure, and the
dehydrogenation methide pathway opt-in, unless a rule's own test needs
them. If a later rule needs one, take it off this list in that commit
rather than inventing a second API.

Old tests: `tests/test_conjugates.py`, `tests/test_star_expand.py`,
`tests/test_dh_methide_pathways.py`.
