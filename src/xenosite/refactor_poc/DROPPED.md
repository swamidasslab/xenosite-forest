This file lists what the proof of concept left out. A listing is not approval.
Each entry has `Status: approved`, `not approved`, or `not decided`.
`approved` means the cut is accepted for the long term.
`not approved` means the capability comes back. The old shape may still change.
`not decided` means it is absent for now and nobody has accepted that.
Do not delete entries. Change the status in the same commit that decides.

## PathContext and its MCS bookkeeping

Status: not decided.

Lived in `src/xenosite/forest/path_context.py`. Multi-embedding rounds,
chemical disagreement, frontier versus deep interior, safe-to-drop
thresholds, and oxygen slack. The atom diff in `find_path.py` is the part
search needs: one pairing, then the local change at each atom.

Old tests: `tests/test_cleavage_site_prune.py`, `tests/test_guided_path_fuzz.py`
(the context object, not the metabolites).

## Formula hint algebra

Status: not decided.

`FormulaHint`, `FormulaAny`, `FormulaMatch`, and `expand_formula_effects`
in `path_context.py`. Local atom differences cover compatibility. The spec
language does not.

Old tests: `tests/test_formula_smarts_hints.py`.

## and_cleave_plan SMARTS replay

Status: not decided.

`and_cleave_plan` and the apply-replay correction in `guided_path.py`.
Precedence comes from `added_by` and the resolved `Effect`.

Old tests: `tests/test_and_cleave_plan_fuzz.py`,
`tests/test_guided_path_cleavage.py`.

## And / Or plan trees and plan JSON

Status: approved.

`And` / `Or` in `step_plan.py`, plan JSON canonicalization, and exact
linearization matching against the old emission. Dependency step plans
(`Deps`) replace them. `Step` and `Deps` are imported, not reimplemented.
Callers do not get the old JSON shape.

Old tests: `tests/test_step_plan.py`, `tests/test_phase1_steps_fuzz.py`
(the JSON and `And`/`Or` assertions, not the chemistry).

## AtomTracker as a base class, and phase1_steps on each rule

Status: not decided for `AtomTracker` as a base class.
Status: not approved for dropping the composite look-ahead. It does not
have to be a `phase1_steps` method on every rule, and the name need not
stay "phase I". See the phase-one entry.

`AtomTracker` in `base.py`, and `phase1_steps` on each rule in
`forest/rules.py`. The poc rules do not subclass `AtomTracker`.

Old tests: `tests/test_phase1_steps_fuzz.py`, `tests/test_atom_trace_char.py`.

## Phase-one reconstruction module

Status: not approved.

The old module (`phaseone.py`, `PhaseOneRS`, `PhaseOneQF`) is not what the
poc calls. The capability stays. Quinone rules are composite look-aheads.
Epoxidation is a complete subset of stable oxidation. N-dealkylation is a
complete subset of unstable oxidation. That is what avoids rewalking the
same chemistry and what explains which rule fired. A ruleset is still a
rule.

Old tests: `tests/test_phaseone.py`.

## Search controls that only bound the old expansion

Status: not decided.

`max_expansions` as a public mode and the old depth counters. A billed
counter stays, as a measurement (`mol_edits + nodes`). Poc `bfs` and `dfs`
enumerate one ruleset up to a depth cap. They take `filter_rules` and
`filter_sites`, and they do not use the atom diff or the closer drop.
Adding them is not approval to drop enumeration, and it is not a port of
the old budget machinery.

Old tests: `tests/test_classic_path_counters.py`, `tests/test_bfs_fuzz.py`,
`tests/test_path_outcome_hard.py` (the budget flag, not the metabolite).

## Conjugate star-label modes

Status: not approved as a general drop. They belong only to the conjugation rules.

Bare `*` versus a full glucuronide or glutathione structure. Do not build
this as a mode on every rule.

Old tests: `tests/test_conjugates.py`, `tests/test_star_expand.py`.

## Dehydrogenation methide opt-in

Status: approved to drop the opt-in. Methide itself is not dropped.

The rule's data always offers methide. `PatternInfo` resolves at most one
methide site. A caller who wants none uses `filter_sites`. There is no
`pathways=("methide",)` flag.

Old tests: `tests/test_dh_methide_pathways.py`.
