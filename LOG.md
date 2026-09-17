# Lab log

## 2026-09-17

- Issue #3: `bfs(PhaseOneRS)` on `CCC(=O)NCC[C@@H]1CCC2=CC=C3OCCC3=C21` crashed after seven Dehydrogenation products (`calcImplicitValence` / RDKit 2026). Root cause: `join_fragments` resonance copies had empty implicit-H caches; `RunReactants` asserted. Fix: `refresh_mol` at end of `join_fragments` and again before each SMARTS `RunReactants`.

## 2026-09-13

- Topo emission: `metabolize` topsite keyed on suffix'd rule labels and site ranks only, so UI-equivalent rows leaked and some distinct product SMILES were dropped. Fix: identity = normalized pathway + sorted topo-rank multiset + product SMILES (same as xenosite UI). SMARTS-loop early skip left as TODO.
- N-dealkylation structure bug: `UO.Dealkylation` SMARTS also cleave O/S/C; UI filtered by nitrogen site. Fix: `NDealkylation(Dealkylation)` post-filters to N-containing sites, ruleset `ND`; short-circuit when reactant has no N. Leave `Dealkylation` on `UO` unchanged.

## 2026-09-09

- Framing: reactivity conjugation is high-sensitivity / low-specificity structure enumeration (if a site were positive, what adduct forms) — not a likelihood model.
- Added aldehyde thiohemiacetal; aziridine; sulfonate ester (`C–OSO2`); isocyanate/isothiocyanate. Negatives: azetidine/pyrrolidine/pyridine vs aziridine; sulfonamide/sulfone/sulfonic acid vs sulfonate; nitrile/amide/CO2/carbodiimide vs isocyanate.
- Earlier: Michael + F/Cl/Br/I for quinones / BnBr / BnI; flipped β-enone skip test. Legacy SMARTS were only epoxide, C–Cl, thiol, terminal `CH2=`.
