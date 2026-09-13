# Lab log

## 2026-09-13

- N-dealkylation structure bug: `UO.Dealkylation` SMARTS also cleave O/S/C; UI filtered by nitrogen site. Fix: `NDealkylation(Dealkylation)` post-filters to N-containing sites, ruleset `ND`; short-circuit when reactant has no N. Leave `Dealkylation` on `UO` unchanged.

## 2026-09-09

- Framing: reactivity conjugation is high-sensitivity / low-specificity structure enumeration (if a site were positive, what adduct forms) — not a likelihood model.
- Added aldehyde thiohemiacetal; aziridine; sulfonate ester (`C–OSO2`); isocyanate/isothiocyanate. Negatives: azetidine/pyrrolidine/pyridine vs aziridine; sulfonamide/sulfone/sulfonic acid vs sulfonate; nitrile/amide/CO2/carbodiimide vs isocyanate.
- Earlier: Michael + F/Cl/Br/I for quinones / BnBr / BnI; flipped β-enone skip test. Legacy SMARTS were only epoxide, C–Cl, thiol, terminal `CH2=`.
