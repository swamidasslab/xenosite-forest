# Leaf-rule product parity — remaining gaps (notes)

Snapshot for adjudication. Reproduce with:

```bash
pytest tests/forest/test_rule_parity_fuzz.py -k '<Rule> and <substring of SMILES>' -q
# or compare bags:
python -c "from tests.forest.test_rule_parity_fuzz import *; ..."
```

Identity: RDKit CSMI (C1). Site bags: topological, engine-local (C2).
Arom-model ruled out (C17). Chemical correctness over blind match (C8).

Live recount 2026-09-27: **29 failing / 3 now PASS** (QF Ph-aziridine, PhNMe2, 4-OH-PhNMe2).
Artifact `parity_full_exocyclic.out` still lists 32 (stale on those three).

Columns: **sym** = harness symptom; **chem** = who looks more correct;
**fix** = fixability.

## Summary table

| rule | mol | SMILES | sym | pyS/rsS | pyP/rsP | chem | fix |
|------|-----|--------|-----|---------|---------|------|-----|
| Dealkylation | olsalazine | `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` | rust_extra | 5/8 | 13/22 | rust_ahead | fixable_now |
| Dealkylation | nitrobenzene | `[O-][N+](=O)c1ccccc1` | form_mismatch | 4/4 | 13/13 | rust_ahead | fixable_now |
| Dealkylation | nitrosobenzene | `O=Nc1ccccc1` | rust_extra | 4/4 | 9/13 | rust_ahead | fixable_now |
| Dealkylation | chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` | form_mismatch | 9/9 | 31/31 | rust_ahead | fixable_now |
| Dealkylation | aspirin | `CC(=O)Oc1ccccc1C(=O)O` | site_count | 10/9 | 22/22 | rust_count_better_wrong_reason | needs_decision |
| Dealkylation | sulfamethoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` | rust_extra | 5/6 | 13/15 | rust_ahead | fixable_now |
| Dealkylation | cinnoline | `c1ccc2nnccc2c1` | form_mismatch | 10/10 | 21/21 | rust_ahead | fixable_now |
| Dehydration | olsalazine | `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` | py_extra | 2/2 | 5/4 | rust_ahead | fixable_now |
| Dehydration | chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` | py_extra | 3/3 | 6/5 | rust_ahead | fixable_now |
| Dehydrogenation | olsalazine | `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` | rust_empty | 1/0 | 1/0 | py_ahead | hard |
| Dehydrogenation | benzene-oxide | `C1=CC2OC2C=C1` | site_count | 3/2 | 2/2 | unclear | needs_decision |
| Dehydrogenation | dihydroacridine | `c1ccc2c(c1)Nc1ccccc1C2` | py_empty | 0/1 | 0/1 | rust_ahead | fixable_now |
| Hydrogenation | olsalazine | `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` | py_extra | 23/14 | 23/14 | unclear | hard |
| Hydrogenation | nitrobenzene | `[O-][N+](=O)c1ccccc1` | form_mismatch | 10/7 | 12/8 | both_wrong | needs_decision |
| Hydrogenation | PhNCO | `O=C=Nc1ccccc1` | py_extra | 8/6 | 9/7 | rust_ahead | fixable_now |
| Hydrogenation | PhNCS | `S=C=Nc1ccccc1` | py_extra | 8/6 | 9/7 | rust_ahead | fixable_now |
| Hydrogenation | chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` | form_mismatch | 11/8 | 15/9 | both_wrong | needs_decision |
| Hydrogenation | PhNCN | `N=C=Nc1ccccc1` | py_extra | 8/6 | 9/7 | rust_ahead | fixable_now |
| Hydrogenation | sulfamethoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` | rust_extra | 7/8 | 8/9 | py_ahead | fixable_now |
| Hydrogenation | cinnoline | `c1ccc2nnccc2c1` | form_mismatch | 25/23 | 30/36 | unclear | hard |
| NDealkylation | chloramphenicol | `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` | form_mismatch | 3/3 | 8/8 | rust_ahead | fixable_now |
| NDealkylation | sulfamethoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` | rust_extra | 1/2 | 2/4 | rust_ahead | fixable_now |
| NDealkylation | cinnoline | `c1ccc2nnccc2c1` | form_mismatch | 2/2 | 3/3 | rust_ahead | fixable_now |
| NitrogenReduction | sulfamethoxazole | `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` | form_mismatch | 1/1 | 1/2 | unclear | needs_decision |
| QuinoneFormation | olsalazine | `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` | py_extra | 9/3 | 9/3 | py_ahead | hard |
| QuinoneFormation | Ph-aziridine | `c1ccccc1N1CC1` | PASS | 4/4 | 5/5 | — | pass |
| QuinoneFormation | benzene-oxide | `C1=CC2OC2C=C1` | py_extra | 2/2 | 5/2 | unclear | hard |
| QuinoneFormation | PhNMe2 | `CN(C)c1ccccc1` | PASS | 4/4 | 6/6 | — | pass |
| QuinoneFormation | 4-OH-PhNMe2 | `CN(C)c1ccc(O)cc1` | PASS | 4/4 | 6/6 | — | pass |
| QuinoneFormation | Ph2NMe | `c1ccc(N(C)c2ccccc2)cc1` | py_extra | 4/4 | 9/6 | py_ahead | fixable_now |
| QuinoneFormation | dihydroacridine | `c1ccc2c(c1)Nc1ccccc1C2` | rust_extra | 9/9 | 11/12 | unclear | needs_decision |
| QuinoneFormation | cinnoline | `c1ccc2nnccc2c1` | form_mismatch | 9/7 | 9/14 | unclear | hard |

## Per-case notes (reproduce + context)

### Dealkylation — olsalazine

- **SMILES:** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`
- **Symptom:** `rust_extra` · sites py/rs `5/8` · products py/rs `13/22` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** more cc_* ring-open sites; Py misses embeddings
- **only Rust:** `C=C(C=CC(O)=C(O)C(=O)O)N=Nc1ccc(O)c(C(=O)O)c1`, `C=CC(=CC(C(=O)O)=C(O)O)N=Nc1ccc(O)c(C(=O)O)c1`, `C=CC(O)=C(C=C(O)N=Nc1ccc(O)c(C(=O)O)c1)C(=O)O`, `O=C(O)C(C=C(C=CO)N=Nc1ccc(O)c(C(=O)O)c1)=CO`, `O=C(O)C(C=CN=Nc1ccc(O)c(C(=O)O)c1)=C(O)C=CO`, `O=C(O)C=C(O)C=CC(=CO)N=Nc1ccc(O)c(C(=O)O)c1`, … +3 more

### Dealkylation — nitrobenzene

- **SMILES:** `[O-][N+](=O)c1ccccc1`
- **Symptom:** `form_mismatch` · sites py/rs `4/4` · products py/rs `13/13` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** nitro charge writing [N+](=O)[O-] vs N([O-])[O-]
- **only Python:** `C=C(C=CC=C=O)N([O-])[O-]`, `C=C(C=CC=CO)N([O-])[O-]`, `C=CC=C(C=C=O)N([O-])[O-]`, `C=CC=C(C=CO)N([O-])[O-]`, `C=CC=CC=C(O)N([O-])[O-]`, `O=C=CC=CC=CN([O-])[O-]`, … +1 more
- **only Rust:** `C=C(C=CC=C=O)[N+](=O)[O-]`, `C=C(C=CC=CO)[N+](=O)[O-]`, `C=CC=C(C=C=O)[N+](=O)[O-]`, `C=CC=C(C=CO)[N+](=O)[O-]`, `C=CC=CC=C(O)[N+](=O)[O-]`, `O=C=CC=CC=C[N+](=O)[O-]`, … +1 more

### Dealkylation — nitrosobenzene

- **SMILES:** `O=Nc1ccccc1`
- **Symptom:** `rust_extra` · sites py/rs `4/4` · products py/rs `9/13` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** extra cc_* ring-open; Py empty those sites
- **only Rust:** `C=C(C=CC=C=O)N=O`, `C=C(C=CC=CO)N=O`, `C=CC=C(C=C=O)N=O`, `C=CC=C(C=CO)N=O`

### Dealkylation — chloramphenicol

- **SMILES:** `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl`
- **Symptom:** `form_mismatch` · sites py/rs `9/9` · products py/rs `31/31` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** nitro leave charge form (same as nitrobenzene)
- **only Python:** `C=C(C=CC(=C=O)C(O)CNC(=O)C(Cl)Cl)N([O-])[O-]`, `C=C(C=CC(=C=O)N([O-])[O-])C(O)CNC(=O)C(Cl)Cl`, `C=C(C=CC(=CO)C(O)CNC(=O)C(Cl)Cl)N([O-])[O-]`, `C=C(C=CC(=CO)N([O-])[O-])C(O)CNC(=O)C(Cl)Cl`, `C=CC(=CC=C(O)C(O)CNC(=O)C(Cl)Cl)N([O-])[O-]`, `C=CC(=CC=C(O)N([O-])[O-])C(O)CNC(=O)C(Cl)Cl`, … +5 more
- **only Rust:** `C=C(C=CC(=C=O)C(O)CNC(=O)C(Cl)Cl)[N+](=O)[O-]`, `C=C(C=CC(=C=O)[N+](=O)[O-])C(O)CNC(=O)C(Cl)Cl`, `C=C(C=CC(=CO)C(O)CNC(=O)C(Cl)Cl)[N+](=O)[O-]`, `C=C(C=CC(=CO)[N+](=O)[O-])C(O)CNC(=O)C(Cl)Cl`, `C=CC(=CC=C(O)C(O)CNC(=O)C(Cl)Cl)[N+](=O)[O-]`, `C=CC(=CC=C(O)[N+](=O)[O-])C(O)CNC(=O)C(Cl)Cl`, … +5 more

### Dealkylation — aspirin

- **SMILES:** `CC(=O)Oc1ccccc1C(=O)O`
- **Symptom:** `site_count` · sites py/rs `10/9` · products py/rs `22/22` · `site_kind=directed_bond`
- **Chem / fix:** rust_count_better_wrong_reason / needs_decision
- **Likely cause:** C13: two quaternary_alcohol on ester O → one bag; Py keeps both (unique_csmi_compliant=False); Rust unique_csmi drops one
- **Reproduce sites:** Python `discovered_site` `(1,3)` acetyl C–O and `(4,3)` aryl C–O, both `quaternary_alcohol`, both emit `{CC(=O)O, O=C(O)c1ccccc1O}`. Rust SMARTS hits both; yield `unique_csmi=true` keeps only `(1,3)`.
- **Adjudication note:** one ester hydrolysis → **one** product bag is chemically right. Rust count (9) is closer than Python (10) but for the wrong reason (silent CSMI). Prefer C13 `product_equiv` or SMARTS/`when` partition — do not force Rust to 10.

### Dealkylation — sulfamethoxazole

- **SMILES:** `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1`
- **Symptom:** `rust_extra` · sites py/rs `5/6` · products py/rs `13/15` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** quaternary_alcohol on isox C–N(sulfonamide); Py silent
- **only Rust:** `Cc1cc(O)no1`, `Nc1ccc(S(N)(=O)=O)cc1`

### Dealkylation — cinnoline

- **SMILES:** `c1ccc2nnccc2c1`
- **Symptom:** `form_mismatch` · sites py/rs `10/10` · products py/rs `21/21` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** ring-open keeps N=N (Rust) vs NN hydrazine (Py)
- **only Python:** `C=C1C=CNNC1=CC=C=O`, `C=C1C=CNNC1=CC=CO`, `C=C1NNC=CC1=CC=C=O`, `C=C1NNC=CC1=CC=CO`, `C=CNNc1ccccc1O`, `NNC=Cc1ccccc1O`, … +5 more
- **only Rust:** `C=CN=Nc1ccccc1O`, `C=c1ccnnc1=CC=C=O`, `C=c1ccnnc1=CC=CO`, `C=c1nnccc1=CC=C=O`, `C=c1nnccc1=CC=CO`, `N=NC=Cc1ccccc1O`, … +5 more

### Dehydration — olsalazine

- **SMILES:** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`
- **Symptom:** `py_extra` · sites py/rs `2/2` · products py/rs `5/4` · `site_kind=atom`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** Py ketene-like dearom junk; Rust accept_o_leave
- **only Python:** `O=C=C1CC(N=Nc2ccc(O)c(C(=O)O)c2)CCC1O`

### Dehydration — chloramphenicol

- **SMILES:** `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl`
- **Symptom:** `py_extra` · sites py/rs `3/3` · products py/rs `6/5` · `site_kind=atom`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** Py quinoid junk; Rust refuse
- **only Python:** `O=C(NCC=C1CCC([N+](=O)[O-])CC1)C(Cl)Cl`

### Dehydrogenation — olsalazine

- **SMILES:** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`
- **Symptom:** `rust_empty` · sites py/rs `1/0` · products py/rs `1/0` · `site_kind=atom_pair`
- **Chem / fix:** py_ahead / hard
- **Likely cause:** bis-quinone across azo; Rust rematch fails linked π
- **only Python:** `O=C(O)C1=CC(=NN=C2C=CC(=O)C(C(=O)O)=C2)C=CC1=O`

### Dehydrogenation — benzene-oxide

- **SMILES:** `C1=CC2OC2C=C1`
- **Symptom:** `site_count` · sites py/rs `3/2` · products py/rs `2/2` · `site_kind=atom_pair`
- **Chem / fix:** unclear / needs_decision
- **Likely cause:** C13-like: two alkyl sites → same benzofuran bag
- **Bag:** both engines emit benzofuran `c1ccc2c(c1)O2`; Python counts that bag twice (two alkyl DH sites), Rust once. Same C13 shape as aspirin.

### Dehydrogenation — dihydroacridine

- **SMILES:** `c1ccc2c(c1)Nc1ccccc1C2`
- **Symptom:** `py_empty` · sites py/rs `0/1` · products py/rs `0/1` · `site_kind=atom_pair`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** acridine aromatization; Py empty
- **only Rust:** `c1ccc2nc3ccccc3cc2c1`

### Hydrogenation — olsalazine

- **SMILES:** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`
- **Symptom:** `py_extra` · sites py/rs `23/14` · products py/rs `23/14` · `site_kind=atom_pair`
- **Chem / fix:** unclear / hard
- **Likely cause:** azo-bridge path saturations Py-only
- **only Python:** `O=C(O)C1=C(O)C=CC(=NN=C2C=CC(O)=C(C(=O)O)C2)C1`, `O=C(O)C1=C(O)C=CC(=NNc2ccc(O)c(C(=O)O)c2)C1`, `O=C(O)C1=CC(=NN=C2C=C(C(=O)O)C(O)=CC2)CC=C1O`, `O=C(O)C1=CC(=NN=C2C=CC(O)=C(C(=O)O)C2)C=CC1O`, `O=C(O)C1=CC(=NN=C2C=CC(O)=C(C(=O)O)C2)CC=C1O`, `O=C(O)C1=CC(=NN=C2C=CC(O)C(C(=O)O)=C2)C=CC1O`, … +3 more

### Hydrogenation — nitrobenzene

- **SMILES:** `[O-][N+](=O)c1ccccc1`
- **Symptom:** `form_mismatch` · sites py/rs `10/7` · products py/rs `12/8` · `site_kind=atom_pair`
- **Chem / fix:** both_wrong / needs_decision
- **Likely cause:** Py radicals / N([O-])[O-] vs Rust [NH+](O)
- **only Python:** `[O-]N([O-])C1=CC=CCC1`, `[O-][N+](O)c1ccccc1`, `[O][N+]([O-])=C1C=CC=CC1`, `[O][N+]([O-])=C1C=CCC=C1`, `[O][N+]([O-])c1ccccc1`
- **only Rust:** `[O-][NH+](O)c1ccccc1`

### Hydrogenation — PhNCO

- **SMILES:** `O=C=Nc1ccccc1`
- **Symptom:** `py_extra` · sites py/rs `8/6` · products py/rs `9/7` · `site_kind=atom_pair`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** cumulated×ring O=CN=C1…; Rust skips (exclusive-seed/cumulated)
- **only Python:** `O=CN=C1C=CC=CC1`, `O=CN=C1C=CCC=C1`

### Hydrogenation — PhNCS

- **SMILES:** `S=C=Nc1ccccc1`
- **Symptom:** `py_extra` · sites py/rs `8/6` · products py/rs `9/7` · `site_kind=atom_pair`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** same as PhNCO for S=CN=C1…
- **only Python:** `S=CN=C1C=CC=CC1`, `S=CN=C1C=CCC=C1`

### Hydrogenation — chloramphenicol

- **SMILES:** `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl`
- **Symptom:** `form_mismatch` · sites py/rs `11/8` · products py/rs `15/9` · `site_kind=atom_pair`
- **Chem / fix:** both_wrong / needs_decision
- **Likely cause:** same nitro-H cluster as nitrobenzene
- **only Python:** `O=C(NCC(O)C1=CC=C(N([O-])[O-])CC1)C(Cl)Cl`, `O=C(NCC(O)C1=CCC(N([O-])[O-])C=C1)C(Cl)Cl`, `O=C(NCC(O)C1C=CC(N([O-])[O-])=CC1)C(Cl)Cl`, `O=C(NCC(O)c1ccc([N+]([O-])O)cc1)C(Cl)Cl`, `[O][N+]([O-])=C1C=CC(C(O)CNC(=O)C(Cl)Cl)=CC1`, `[O][N+]([O-])=C1C=CC(C(O)CNC(=O)C(Cl)Cl)C=C1`, … +1 more
- **only Rust:** `O=C(NCC(O)c1ccc([NH+]([O-])O)cc1)C(Cl)Cl`

### Hydrogenation — PhNCN

- **SMILES:** `N=C=Nc1ccccc1`
- **Symptom:** `py_extra` · sites py/rs `8/6` · products py/rs `9/7` · `site_kind=atom_pair`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** same as PhNCO for N=CN=C1…
- **only Python:** `N=CN=C1C=CC=CC1`, `N=CN=C1C=CCC=C1`

### Hydrogenation — sulfamethoxazole

- **SMILES:** `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1`
- **Symptom:** `rust_extra` · sites py/rs `7/8` · products py/rs `8/9` · `site_kind=atom_pair`
- **Chem / fix:** py_ahead / fixable_now
- **Likely cause:** Rust path saturates S=O → N[SH](=O)(O) junk
- **only Rust:** `Cc1cc(N[SH](=O)(O)c2ccc(N)cc2)no1`

### Hydrogenation — cinnoline

- **SMILES:** `c1ccc2nnccc2c1`
- **Symptom:** `form_mismatch` · sites py/rs `25/23` · products py/rs `30/36` · `site_kind=atom_pair`
- **Chem / fix:** unclear / hard
- **Likely cause:** Py keeps NN= dihydros; Rust also emits NN hydrazine paths
- **only Python:** `C1=CC2=NN=CCC2=CC1`, `C1=CCC2=NN=CCC2=C1`
- **only Rust:** `C1=CC2=C(C=CNN2)CC1`, `C1=CC2=C(CC1)NNC=C2`, `C1=CC2=CCNNC2=CC1`, `C1=CC2=CCNNC2C=C1`, `C1=CC2C=CNNC2=CC1`, `C1=CC2NNC=CC2=CC1`, … +2 more

### NDealkylation — chloramphenicol

- **SMILES:** `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl`
- **Symptom:** `form_mismatch` · sites py/rs `3/3` · products py/rs `8/8` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** nitro leave charge form
- **only Python:** `[O-]N[O-]`
- **only Rust:** `O=N[O-]`

### NDealkylation — sulfamethoxazole

- **SMILES:** `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1`
- **Symptom:** `rust_extra` · sites py/rs `1/2` · products py/rs `2/4` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** same C–N sulfonamide cleavage as Dealk
- **only Rust:** `Cc1cc(O)no1`, `Nc1ccc(S(N)(=O)=O)cc1`

### NDealkylation — cinnoline

- **SMILES:** `c1ccc2nnccc2c1`
- **Symptom:** `form_mismatch` · sites py/rs `2/2` · products py/rs `3/3` · `site_kind=directed_bond`
- **Chem / fix:** rust_ahead / fixable_now
- **Likely cause:** ring-open N=N vs NN
- **only Python:** `NNC=Cc1ccccc1O`, `NNc1ccccc1C=C=O`, `NNc1ccccc1C=CO`
- **only Rust:** `N=NC=Cc1ccccc1O`, `N=Nc1ccccc1C=C=O`, `N=Nc1ccccc1C=CO`

### NitrogenReduction — sulfamethoxazole

- **SMILES:** `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1`
- **Symptom:** `form_mismatch` · sites py/rs `1/1` · products py/rs `1/2` · `site_kind=atom`
- **Chem / fix:** unclear / needs_decision
- **Likely cause:** isoxazole N–O open: Py amide-amine vs Rust enol/imine
- **only Python:** `CC(=O)C=C(N)NS(=O)(=O)c1ccc(N)cc1`
- **only Rust:** `CC(O)=CCNS(=O)(=O)c1ccc(N)cc1`, `CC=CC(=N)NS(=O)(=O)c1ccc(N)cc1`

### QuinoneFormation — olsalazine

- **SMILES:** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O`
- **Symptom:** `py_extra` · sites py/rs `9/3` · products py/rs `9/3` · `site_kind=atom_pair`
- **Chem / fix:** py_ahead / hard
- **Likely cause:** cross-azo quinones; Rust rematch thin
- **only Python:** `O=C(O)C1=C(O)C=CC(=NN=C2C=CC(O)=C(C(=O)O)C2=O)C1=O`, `O=C(O)C1=CC(=NN=C2C=C(C(=O)O)C(O)=CC2=O)C(=O)C=C1O`, `O=C(O)C1=CC(=NN=C2C=C(C(=O)O)C(O)=CC2=O)C=CC1=O`, `O=C(O)C1=CC(=NN=C2C=CC(=O)C(C(=O)O)=C2)C=CC1=O`, `O=C(O)C1=CC(=NN=C2C=CC(O)=C(C(=O)O)C2=O)C(=O)C=C1O`, `O=C(O)C1=CC(=NN=C2C=CC(O)=C(C(=O)O)C2=O)C=CC1=O`

### QuinoneFormation — Ph-aziridine

- **SMILES:** `c1ccccc1N1CC1`
- **Symptom:** `PASS` · sites py/rs `4/4` · products py/rs `5/5` · `site_kind=atom_pair`
- **Chem / fix:** — / pass
- **Likely cause:** iminium refuse now green
- **Status:** green on current tree (artifact stale).

### QuinoneFormation — benzene-oxide

- **SMILES:** `C1=CC2OC2C=C1`
- **Symptom:** `py_extra` · sites py/rs `2/2` · products py/rs `5/2` · `site_kind=atom_pair`
- **Chem / fix:** unclear / hard
- **Likely cause:** Py extra oxide/quinone + benzofuran; Rust thinner
- **only Python:** `O=C1C=CC2OC2=C1`, `O=C1C=CC=C2OC12`, `c1ccc2c(c1)O2`

### QuinoneFormation — PhNMe2

- **SMILES:** `CN(C)c1ccccc1`
- **Symptom:** `PASS` · sites py/rs `4/4` · products py/rs `6/6` · `site_kind=atom_pair`
- **Chem / fix:** — / pass
- **Likely cause:** iminium refuse now green
- **Status:** green on current tree (artifact stale).

### QuinoneFormation — 4-OH-PhNMe2

- **SMILES:** `CN(C)c1ccc(O)cc1`
- **Symptom:** `PASS` · sites py/rs `4/4` · products py/rs `6/6` · `site_kind=atom_pair`
- **Chem / fix:** — / pass
- **Likely cause:** iminium refuse now green
- **Status:** green on current tree (artifact stale).

### QuinoneFormation — Ph2NMe

- **SMILES:** `c1ccc(N(C)c2ccccc2)cc1`
- **Symptom:** `py_extra` · sites py/rs `4/4` · products py/rs `9/6` · `site_kind=atom_pair`
- **Chem / fix:** py_ahead / fixable_now
- **Likely cause:** Py dealkylates Ph leave; Rust only Me leave
- **only Python:** `CN=C1C=CC(=O)C=C1`, `CN=C1C=CC=CC1=O`, `c1ccccc1`

### QuinoneFormation — dihydroacridine

- **SMILES:** `c1ccc2c(c1)Nc1ccccc1C2`
- **Symptom:** `rust_extra` · sites py/rs `9/9` · products py/rs `11/12` · `site_kind=atom_pair`
- **Chem / fix:** unclear / needs_decision
- **Likely cause:** Rust also emits acridine under QF (same as DH)
- **only Rust:** `c1ccc2nc3ccccc3cc2c1`

### QuinoneFormation — cinnoline

- **SMILES:** `c1ccc2nnccc2c1`
- **Symptom:** `form_mismatch` · sites py/rs `9/7` · products py/rs `9/14` · `site_kind=atom_pair`
- **Chem / fix:** unclear / hard
- **Likely cause:** Py azo-diones vs Rust hydrazine-diones [nH][nH]
- **only Python:** `O=C1C=CC2=NN=CC(=O)C2=C1`, `O=C1C=NN=C2C(=O)C=CC=C12`
- **only Rust:** `O=c1[nH][nH]c2ccccc2c1=O`, `O=c1cc2c(=O)cccc-2[nH][nH]1`, `O=c1cc2cc[nH][nH]c-2cc1=O`, `O=c1ccc(=O)c2[nH][nH]ccc1=2`, `O=c1ccc2[nH][nH]ccc=2c1=O`, `O=c1ccc2cc(=O)[nH][nH]c-2c1`, … +1 more

## Adjudication queue (unclear / both_wrong / needs_decision)

Focus these when picking chemistry; SMILES in table above.

- **Dealkylation / aspirin** `CC(=O)Oc1ccccc1C(=O)O` — C13: two quaternary_alcohol on ester O → one bag; Py keeps both (unique_csmi_compliant=False); Rust unique_csmi drops one
- **Dehydrogenation / benzene-oxide** `C1=CC2OC2C=C1` — C13-like: two alkyl sites → same benzofuran bag
- **Hydrogenation / olsalazine** `OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O` — azo-bridge path saturations Py-only
- **Hydrogenation / nitrobenzene** `[O-][N+](=O)c1ccccc1` — Py radicals / N([O-])[O-] vs Rust [NH+](O)
- **Hydrogenation / chloramphenicol** `O=C(NCC(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl` — same nitro-H cluster as nitrobenzene
- **Hydrogenation / cinnoline** `c1ccc2nnccc2c1` — Py keeps NN= dihydros; Rust also emits NN hydrazine paths
- **NitrogenReduction / sulfamethoxazole** `Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1` — isoxazole N–O open: Py amide-amine vs Rust enol/imine
- **QuinoneFormation / benzene-oxide** `C1=CC2OC2C=C1` — Py extra oxide/quinone + benzofuran; Rust thinner
- **QuinoneFormation / dihydroacridine** `c1ccc2c(c1)Nc1ccccc1C2` — Rust also emits acridine under QF (same as DH)
- **QuinoneFormation / cinnoline** `c1ccc2nnccc2c1` — Py azo-diones vs Rust hydrazine-diones [nH][nH]

## See also

- [RUST_PYTHON_PARITY.md](RUST_PYTHON_PARITY.md) C8, C11, C13, C16, C17
- [HEURISTICS.md](HEURISTICS.md) directed_bond, unique_csmi_compliant, product-identical sites
