# Divergences

Canonical-SMILES disagreements between live ``xenosite.forest`` and the archived
pre-swap library (``xenosite._archive_forest``). A listing is the record. It is not permission to drop a real
metabolite the old code emits, and it is not permission to hide a product
only the new code emits. Sites may differ when `topol_equiv` puts them in
the same atom class. Order does not matter.

Also give an assessment of whether the divergence is a problem or not, whether it is more or less correct and why, or if it is neutral.

## Hydroxylation of butylbenzene

The old second pattern writes a carbonyl. Both patterns add OH and remove
one H. The `>>[*:1]=O` products are not metabolites.

- Reactant: `c1ccc(CCCC)cc1`
- Only old: `CC(=O)CCc1ccccc1`, `CCC(=O)Cc1ccccc1`, `CCCC(=O)c1ccccc1`
- Only new: none
- Shared: `CC(O)CCc1ccccc1`, `CCC(O)Cc1ccccc1`, `CCCC(O)c1ccccc1`, `CCCCc1ccc(O)cc1`, `CCCCc1cccc(O)c1`, `CCCCc1ccccc1O`, `OCCCCc1ccccc1`

Ethane (`CC`) does not hit that carbonyl. Both libraries yield `CCO` only.

More correct. Not a problem in live forest. The three old-only strings are the `>>[*:1]=O` products, and those are not metabolites. Both patterns here add OH and remove one H, which is the shared set. Ethane does not hit that carbonyl, and both libraries already yield `CCO` only.

## Azo splitting of pyridazine

The old library kekulizes once. That form of pyridazine has an N-N single bond, so the split does not run. Live forest matches aromatic ring N=N with `=,:` on a ResonanceRule and reacts on the Kekulé parent where that bond is double. RDKit's `[#7:1]=,:[#7:2]>>[*:1].[*:2]` then drops one nitrogen. `C=CC=CN` is not a metabolite: pyridazine is C4H4N2. Both Kekulé forms parse, and neither atom is charged, so the dropped nitrogen is the reaction, not a charge left on the wrong atom.

- Reactant: `c1ccnnc1`
- Only old: none
- Only new: `C=CC=CN`
- Shared: none

The same SMARTS on `c1ccc(N=Nc2ccnnc2)cc1` is what both libraries emit, including the ring-opened fragments. Those sets match.

Problem. Live forest is less correct on pyridazine. `C=CC=CN` is not a metabolite: the SMARTS drops one nitrogen, and the rule says both fragments stay. The old library emits nothing only because its one Kekulé form has an N-N single bond. The phenylazo-pyridazine sets already match, so the bad product is this ring writing. Moving formal charge with the bond flip does not change it.

## Nitro reduction of 1,4-dinitrobenzene

A Kekulé copy used to rewrite one nitro as `[N+](=[O-])O` and leave the charge on the oxygen that had become double-bonded. That oxygen has valence 2, so the SMILES did not parse. The charge now follows the bond order. Both libraries yield the nitroso, and the leaving oxygen.

- Reactant: `O=[N+]([O-])c1ccc([N+](=O)[O-])cc1`
- Only old: none
- Only new: none
- Shared: `O=Nc1ccc([N+](=O)[O-])cc1`, `O`

Not a problem. The nitroso is the metabolite the rule describes, and it parses. Nitrobenzene was already that pair.

## Thiophene S-oxidation of benzothiophene

The old library kekulizes once. That form is `C1=CC=C2SC=CC2=C1`, so the thiophene ring is not `C=C-C=C-S` and the oxidation does not run. Live forest matches aromatic thiophene with `=,:` on a ResonanceRule and reacts on the Kekulé parent where maps 1–2 are the S–C single (SMARTS-implied order). That writes the S-oxide. Thiophene itself has no second writing to miss: both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.

- Reactant: `c1ccc2sccc2c1`
- Only old: none
- Only new: `[O-][s+]1ccc2ccccc21`
- Shared: none

Dibenzothiophene (`c1ccc2c(c1)sc1ccccc12`): only new `[O-][s+]1c2ccccc2c2ccccc21`.

More correct. Not a problem in live forest. ResonanceRule picks the Kekulé parent that matches the S–C single implied by the SMARTS, so live forest emits the S-oxide the old one-Kekulé walk misses. Thiophene has no second writing, and both libraries emit `[O-][s+]1cccc1`. Dibenzothiophene is the same miss.

## No global Kekulé loop on SmartsReactionRule

Status: documented (shipped on PR #15; residual CSMI open)

Plain `SmartsReactionRule.metabolites` no longer fans out over `_kekule_forms(mol)`. Matching and unique-edit run on the input (usually aromatic) mol. Rules that need a specific Kekulé bond order for `RunReactants` are `ResonanceRule` subclasses: aromatic-matching SMARTS (`=,:` / `-,:`), then `_reactant_parent` selects the assignment whose maps 1–2 match the SMARTS-implied order (`_smarts_mapped_bond_order`). Unique-edit `incident_orders` stay on the aromatic parent so Kekulé 1.0↔2.0 flips do not split equivalent sites (phenol Hydroxylation: 0 CSMI). Reparented: Dealkylation, AzoSplitting, ThiopheneSulfurOxidation, NitrogenReduction (hydroxylamine).

Chem suite green with `CsmiDedupWarning` ignored. Remaining WAE failures are still mostly Epoxidation CSMI (plus Hydrogenation / Dealkylation / OxidativeDehalogenation / peers) — unique-edit does not yet fully collapse those symmetric sites. Not a parity claim against the archive; residual same-product drops under live `unique_csmi`.

## Arene-oxide methyl sulfone

The old second SMARTS deletes the epoxide oxygen and writes `-S(C)(=O)(=O)`. The methyl carbon is not in the reactant. The product is not a sulfate conjugate.

- Reactant: `C1=CC2OC2C=C1`
- Only old: `CS(=O)(=O)c1ccccc1`
- Only new: none
- Shared: none

The same writing on `Clc1ccc(c(c1)Cl)C1=CC2(Cl)OC2C=C1Cl` is only old `CS(=O)(=O)c1cc(Cl)c(-c2ccc(Cl)cc2Cl)cc1Cl`.

`OC1=CC2OC2C=C1` still shares the sulfate `O=S(=O)(O)OC1=CC2OC2C=C1`. Only old there: `CS(=O)(=O)c1ccc(O)cc1`, `CS(=O)(=O)c1cccc(O)c1`.

More correct. Not a problem in live forest. Those strings are methyl sulfones. A sulfate conjugate keeps the oxygen and adds `S(=O)(=O)O`. The epoxide oxygen has no hydrogen, so it is not a sulfate site.

## Carbonyl heteroatom rewritten as oxygen

The old first SMARTS is `[#8:1][#6:2](=[O,N,P,S:3])[#6:4]` and the product writes `=[#8:3]`. Nitrogen, phosphorus, and sulfur become oxygen. That product is a carboxylic glucuronide, not the glucuronide of the matched group.

- Reactant: `N=C(O)c1ccccc1`
- Only old: `O=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`
- Only new: none
- Shared: `N=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`

The same writing on `S=C(O)c1ccccc1` is only old `O=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1`. Shared: `O=C(O)C1OC(OC(=S)c2ccccc2)C(O)C(O)C1O`.

On `N=C([O-])C` the alcohol pattern does not match. Only old: `CC(=O)OC1OC(C(=O)O)C(O)C(O)C1O`. Shared: none. `S=C([O-])C` and `P=C([O-])C` are the same only-old carboxylate.

More correct. Not a problem in live forest. Those strings replace the heteroatom with oxygen. A glucuronide of `N=C(O)-` or `S=C(O)-` keeps that atom. The pattern here matches `=[#8]` on an OH or an anion, so benzoate still gives the acyl glucuronide. An ester oxygen is not that site: the old product does not sanitize, and this pattern does not match it.

## Hydrogenation of 1,3-butadiene

The class is a ResonancePairRule. The one-bond SMARTS writes 1-butene. The pair flips an alternating path between two carbons and writes 2-butene. That gap is closed. The old walk also writes 1,2-butadiene. That path has an even number of bonds, so it is not a pair here.

- Reactant: `C=CC=C`
- Only old: `C=C=CC`
- Only new: none
- Shared: `C=CCC`, `CC=CC`

Ethene (`C=C`) and ethyne (`C#C`) do not hit that path. Both libraries yield `CC` and `C=C`.

More correct to omit `C=C=CC`. It is still C4H6, the same formula as butadiene, not a hydrogenation. `CC=CC` is C4H8, the 1,4 product, and it is shared.

**Bug note (unmasked by dd9c8ef unique-edit, fixed):** Hydrogenation / ResonancePair path search could start on a *single* bond and apply even bond-count walks, falsely emitting `C=C=C=C` (C4H4, path ends `(1,2)`) and `C=C=CC` (C4H6). Unique-edit dedup had been *masking* those false emits by suppressing them — the chemistry bug is independent of dedup. Fix: double-first paths only (either endpoint may hold the opening double) and odd bond-count walks only. Do not “fix” by restoring bad dedup.

## Dehydrogenation of ethenediol

The pair path writes glyoxal. The draft was `O=[CH2][CH2]=O` and failed sanitization. Clearing those explicit hydrogens and sanitizing again yields `O=CC=O`. That gap is closed. The one-bond SMARTS still writes the ketene.

- Reactant: `OC=CO`
- Only old: none
- Only new: none
- Shared: `O=C=CO`, `O=CC=O`

Ethanol (`CCO`) is one bond, not that path. Both libraries yield `CC=O` and `C=CO`.

Not a problem. The gap is closed. `O=CC=O` is the dehydrogenation of both alcohols (C2H2O2). The shared ketene is the one-end SMARTS.

## Nitrogen reduction of nitrosomethane

The old nitroso SMARTS is `[#7D2:1]=[#8:2]>>([*:1].[*2])`. `[*2]` is not an atom map. The oxygen comes out as `*`. The pattern here is `[*:2]`, so the leaving group is oxygen.

- Reactant: `CN=O`
- Only old: `*`
- Only new: `O`
- Shared: `CN`

Nitromethane (`C[N+](=O)[O-]`) does not hit that string. Both libraries yield `CN`, `CN=O`, and `O`.

More correct. Not a problem in live forest. `*` is not oxygen. The nitroso oxygen is the leaving group, and both libraries already emit the amine.

## Dephosphorylation of methyl phosphate

Substructure matching uniquifies the query. Maps 1, 4, and 5 are the same oxygen class in the old SMARTS, so one match was kept — and RDKit kept a P-OH cleavage (water + methyl phosphite) instead of the ester.

The pattern now requires map 1 to be carbon-bound (`[#8;$([#8][#6]):1]`), so uniquify cannot prefer P-OH.

- Reactant: `COP(=O)(O)O` / `O=P(O)(O)OC`
- Only old (unconstrained SMARTS + uniquify): `CO[PH](=O)O`, `O`
- New / shared: `CO`, `O=[PH](O)O`

More correct. Dephosphorylation takes the phosphate off the carbon.

## Symmetry-equivalent sites

Status: approved

Live forest collapses symmetry-equivalent sites. The unique-edit signature uses atom rank (`topol_equiv`), so equivalent carbons are one edit. The archived library emits every atom index. This is intentional, not a regression to fix.

Neutral. The metabolite structure is the same either way. The reported site is one atom of the class, not the full orbit. Switching the signature to atom index would match the archived library's indices and would change what a unique edit is. That change is not wanted.

This approval is the one-site case. Two-site edits do not collapse both ends to rank. The signature also carries the unordered atom-pair orbit key from `graph_isomorphism` (isotope marks by default, or pynauty when that package imports; see `HEURISTICS.md`). A check of the old both-ends rank keys (benzene, naphthalene, the xylenes, biphenyl, diphenylmethane, diphenyl ether, aniline dimer, triethylamine, hydroquinone, resorcinol, and the SMARTS rules whose site has two atoms) found no signature that grouped two matches and then emitted two different products. Incident bond orders split ortho from para. That check is not a proof that both-ends rank is enough: benzene meta and para share a rank key.

## Product csmi dedup vs archived site+product identity

Status: not decided (RuleSet cross-rule); rule-level key approved

The archived library's `metabolize` with `only_emit_topologically_distinct_sites` keys `(rule name, topo ranks, product SMILES set)` — site topology stays in the product-emit signature. Live forest splits that into two layers (see `HEURISTICS.md`): unique-edit keeps ranks + pair orbit; `unique_csmi` on `ReactionRule` keys `(rule name, PatternInfo.name | SMARTS, product csmi)` and drops site topology. Same product from two sites of one rule is one yield when the pattern token matches; overlapping SMARTS that share a product need partitioned data (e.g. Dealkylation C–C alcohol is `#6H0` vs `#6h`; Hydroxylation is `#6h1` vs `#6h2,#6h3`; Glutathionation epoxide/aziridine substituted-carbon patterns are `#6H0`), not archived-style site+product cross-collapse.

`RuleSet.metabolize` now uses the same `_unique_csmi_key` as `ReactionRule`, so two child rules that share a product structure both emit. The archive keeps a separate seen per child rule as well; the remaining divergence is site topology in the archived emit signature, not cross-rule csmi merging.

## Residual same-product multi-pattern emissions

Status: documented (not a forest parity bug)

A library scan still finds reactant/product pairs that appear under more than one `PatternInfo.name` of the same rule. Those are not H-count nesting of the Hydroxylation kind. Categories:

- **Distinct edits, shared cleavage fragment.** Dealkylation / NDealkylation methyl(ene) alcohol vs carbonyl vs carboxylic oxygenate the carbon side differently; the heteroatom-side fragment can be the same molecule, and each fragment is its own `unique_csmi` yield. Dehydration alcohol vs beta-elimination both leave water. OxidativeDehalogenation alcohol vs carboxylic share the halide piece.
- **Dual nitro writings.** NitroaromaticReduction and NitrogenReduction keep charged and neutral (and multi-O) patterns so both RDKit nitro forms fire; products can coincide.
- **Pair path vs one-bond SMARTS.** Hydrogenation ResonancePairRule emissions (no pattern token) often match the alkene SMARTS product.
- **Conjugation site choice.** Glucuronidation alcohol vs carboxylate on a carboxylic acid can share an acyl glucuronide when both patterns match.

Full non-overlap is not required. Partition when two patterns implement the same edit on overlapping applicability (H-count nesting, or `#6` vs `#6h` / `#6H0`). Leave chemically distinct edits alone.
