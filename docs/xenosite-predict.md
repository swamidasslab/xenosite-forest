# Notes from xenosite-predict

Downstream package: sibling [`xenosite-predict`](https://github.com/swamidasslab/xenosite-predict) (`../xenosite-predict`). It enumerates Phase I and conjugate structures with this library.

## RDKit 2026 valence caches

**Folded into this package** (v0.2.2+): `refresh_mol` / `UpdatePropertyCache(strict=False)` before `RunReactants` and SMILES, including `clean` fragments and unique-key SMILES. Failed reactants **continue** to the next SMARTS / resonance copy instead of aborting the rule.

Predict’s temporary shim was:

[`../xenosite-predict/src/xenosite/predict/forest_rdkit.py`](../../xenosite-predict/src/xenosite/predict/forest_rdkit.py)

Regressions live in `tests/test_rdkit_valence.py` (diphenhydramine, ibuprofen, skip-and-continue).

### Examples that crashed on the unpatched wheel, then enumerated

| Name / note | SMILES |
| --- | --- |
| Diphenhydramine | `CN(C)CCOC(c1ccccc1)c1ccccc1` |
| Ibuprofen | `CC(C)Cc1ccc(C(C)C(=O)O)cc1` |
| Atenolol | `CC(C)NCC(O)COc1ccc(CC(N)=O)cc1` |
| Cinnarizine | `C(=Cc1ccccc1)CN1CCN(C(c2ccccc2)c2ccccc2)CC1` |
| Warfarin | `CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O` |
| Trimethoprim | `COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC` |
| Omeprazole-like | `COc1cc2nc(SCc3ccccc3C)[nH]c2cc1OC` |
| Fluconazole | `OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F` |
| Sulfisoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` |
| Phenobarbital | `CCC1(c2ccccc2)C(=O)NC(=O)NC1=O` |
| Chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` |
| Atropine | `CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2` |
| Penicillin G core | `CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N21` |
| Sudoxicam-like | `CN1C(C(=O)Nc2nccs2)=C(O)c2ccccc2S1(=O)=O` |

## Conjugation rules (DNA / cyanide vs GSH)

**Folded into this package** (v0.2.3+): conjugation rules default to bare `*` adducts; optional `star_label` / `as_star=False`; `include_thiol` and `load_ruleset("GlutathionationNoThiol")`. See [usage](usage.md#conjugation-phase-ii).

**UGT SOM** (v0.2.5+): `Glucuronidation` reports a **single-atom** site — the oxygen that receives GlcA (acid OH or phenolic/alcoholic OH). Earlier versions included neighboring mapped atoms in the site.

Predict head → forest:

| Predict head | Forest |
| --- | --- |
| `ugt` | `Glucuronidation(star_label="GlcA")` or `CJ.Glucuronidation` |
| `reactivity.gsh` | `Glutathionation(star_label="GSH")` |
| `reactivity.protein` | `Glutathionation(star_label="Protein")` (star-only) |
| `reactivity.dna` | `Glutathionation(include_thiol=False, star_label="DNA")` or `GlutathionationNoThiol` |
| `reactivity.cyanide` | `Glutathionation(include_thiol=False, star_label="Cyanide")` (label is `Cyanide`, not `CN`) |

Predict’s older adapter ([`conjugates.py`](../../xenosite-predict/src/xenosite/predict/conjugates.py)) can drop once it calls these options. Forest keeps AtomTracker on star products; `mol_to_cxsmiles` copies before stripping props so tracing is unchanged.

### Examples where GSH vs DNA/Cyanide output differs

Thiol substrates match `[#16h1` and must **not** appear as DNA/Cyanide conjugates:

| Note | SMILES |
| --- | --- |
| Aliphatic thiol | `CCS` |
| Thiophenol | `Sc1ccccc1` |
| Cysteine-like thiol | `SC[C@H](N)C(=O)O` |

Electrophiles that **should** still conjugate for GSH, DNA, and Cyanide (no thiol SMARTS):

| Note | SMILES |
| --- | --- |
| Styrene oxide | `c1ccccc1C1OC1` |
| Benzyl chloride | `ClCc1ccccc1` |
| Terminal alkene | `C=CC` |
