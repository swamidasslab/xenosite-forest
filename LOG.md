# Lab log

## 2026-09-17

- Issue #3 reports `calcImplicitValence` crash on PhaseOneRS BFS; reproduces on forest 0.1.0 only (already fixed in 0.2.0+).
- Full BFS depth=2 failed with RDKit `Range Error` after acetylation→dehydrogenation: `_remove_dummy_atoms` deleted conjugation `*` adducts. Fix: mark FragmentOnBonds dummies and only remove those.
- BFS `expand_star_conjugates=False` by default: star conjugates are terminal products; opt in to expand them further.
- Added `dfs` pathway search; fuzz samples a few Full DFS two-step pathways instead of draining BFS.

## 2026-09-13

- Topo emission: `metabolize` topsite keyed on suffix'd rule labels and site ranks only, so UI-equivalent rows leaked and some distinct product SMILES were dropped. Fix: identity = normalized pathway + sorted topo-rank multiset + product SMILES (same as xenosite UI). SMARTS-loop early skip left as TODO.
- N-dealkylation structure bug: `UO.Dealkylation` SMARTS also cleave O/S/C; UI filtered by nitrogen site. Fix: `NDealkylation(Dealkylation)` post-filters to N-containing sites, ruleset `ND`; short-circuit when reactant has no N. Leave `Dealkylation` on `UO` unchanged.

## 2026-09-09

- Framing: reactivity conjugation is high-sensitivity / low-specificity structure enumeration (if a site were positive, what adduct forms) — not a likelihood model.
- Added aldehyde thiohemiacetal; aziridine; sulfonate ester (`C–OSO2`); isocyanate/isothiocyanate. Negatives: azetidine/pyrrolidine/pyridine vs aziridine; sulfonamide/sulfone/sulfonic acid vs sulfonate; nitrile/amide/CO2/carbodiimide vs isocyanate.
- Earlier: Michael + F/Cl/Br/I for quinones / BnBr / BnI; flipped β-enone skip test. Legacy SMARTS were only epoxide, C–Cl, thiol, terminal `CH2=`.
