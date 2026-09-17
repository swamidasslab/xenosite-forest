# Lab log

## 2026-09-17

- Mol-scoped resonance cache on `mol._forest["resonance"]` (lazy pull-through joined forms per conjugated/aromatic mode; shared `bfs_all_pairs`). Private `_resonance_cache_disabled()` for parity tests. `EditMol.standardize` propagates `_forest`.
- Full product parity cache on vs off: APAP 69, naph_styryl 228, multi_conj 400 (matched). Form counts matched.
- Full.metabolites call counts (on → off): `_resfrags` 2→7, `join_fragments` APAP 4→14 / naph 10→35 / multi 14→49, `bfs` 1→4.
- Wall (median): Full multi_conj ~1.72s on vs ~1.77s off (~1.02×); resonance-only ruleset ~1.06–1.09×. Repeated `resonance_structures`×5 on same mol ~4.1–4.2×. Full wall dominated by pair-path/SMARTS work after forms exist.

## 2026-09-09

- Framing: reactivity conjugation is high-sensitivity / low-specificity structure enumeration (if a site were positive, what adduct forms) — not a likelihood model.
- Added aldehyde thiohemiacetal; aziridine; sulfonate ester (`C–OSO2`); isocyanate/isothiocyanate. Negatives: azetidine/pyrrolidine/pyridine vs aziridine; sulfonamide/sulfone/sulfonic acid vs sulfonate; nitrile/amide/CO2/carbodiimide vs isocyanate.
- Earlier: Michael + F/Cl/Br/I for quinones / BnBr / BnI; flipped β-enone skip test. Legacy SMARTS were only epoxide, C–Cl, thiol, terminal `CH2=`.
