# Proposed terms — feedback round 002

Forest-mapped Phase I class parents use full names (stable/unstable oxygenation, dehydrogenation, hydrolysis, reduction). Leaves and facets from Rainbow (21 types), Metabolic Forest (conjugation / quinone / tautomerization), and the quinone-formation paper (species + one-/two-step). No SO/UO/DH/HD/RD as XRM labels.

Regenerate: `cargo run -p xenosite-xrm --example sample_terms -- --write`

## ethane → ethanol (aliphatic hydroxylation)
- reactant: `CC`
- product: `CCO`
- tags: _(none)_
- terms:
  - aliphatic hydroxylation [xrm:0000102] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aliphatic hydroxylation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - hydroxylation [xrm:0000100] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## benzene → phenol (aromatic hydroxylation)
- reactant: `c1ccccc1`
- product: `Oc1ccccc1`
- tags: _(none)_
- terms:
  - aromatic hydroxylation [xrm:0000101] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aromatic hydroxylation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - hydroxylation [xrm:0000100] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:para-hydroxylation
- reactant: `CCc1ccccc1`
- product: `CCc1ccc(O)cc1`
- tags: `chem:para-hydroxylation`, `chem:aromatic-hydroxylation`
- terms:
  - para-hydroxylation [xrm:0000103] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aromatic hydroxylation > para-hydroxylation
  - aromatic hydroxylation [xrm:0000101] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aromatic hydroxylation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - hydroxylation [xrm:0000100] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:benzylic-hydroxylation
- reactant: `CCc1ccccc1`
- product: `CC(O)c1ccccc1`
- tags: `chem:benzylic-hydroxylation`
- terms:
  - aliphatic hydroxylation [xrm:0000102] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aliphatic hydroxylation
  - benzylic hydroxylation [xrm:0000106] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > benzylic hydroxylation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - hydroxylation [xrm:0000100] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## ethene → oxirane
- reactant: `C=C`
- product: `C1CO1`
- tags: _(none)_
- terms:
  - alkene epoxidation [xrm:0000111] ← xenobiotic biotransformation > phase I > stable oxygenation > epoxidation > alkene epoxidation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - epoxidation [xrm:0000110] ← xenobiotic biotransformation > phase I > stable oxygenation > epoxidation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:arene-oxide + NIH-shift facets
- reactant: `c1ccccc1`
- product: `Oc1ccccc1`
- tags: `chem:arene-oxide`, `chem:NIH-shift`
- terms:
  - arene oxide formation [xrm:0000112] ← xenobiotic biotransformation > phase I > stable oxygenation > epoxidation > arene oxide formation
  - aromatic hydroxylation [xrm:0000101] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation > aromatic hydroxylation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - epoxidation [xrm:0000110] ← xenobiotic biotransformation > phase I > stable oxygenation > epoxidation
  - hydroxylation [xrm:0000100] ← xenobiotic biotransformation > phase I > stable oxygenation > hydroxylation
  - NIH shift [xrm:0004003] ← xenobiotic biotransformation > process facet > NIH shift
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - rearomatization [xrm:0004002] ← xenobiotic biotransformation > process facet > rearomatization
  - stable oxygenation [xrm:0000010] ← xenobiotic biotransformation > phase I > stable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## ethanol → acetaldehyde
- reactant: `CCO`
- product: `CC=O`
- tags: _(none)_
- terms:
  - primary alcohol oxidation [xrm:0000311] ← xenobiotic biotransformation > phase I > dehydrogenation > alcohol oxidation > primary alcohol oxidation
  - alcohol oxidation [xrm:0000310] ← xenobiotic biotransformation > phase I > dehydrogenation > alcohol oxidation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:N-demethylation
- reactant: `CN(C)C`
- product: `CNC`
- tags: `chem:N-demethylation`
- terms:
  - N-demethylation [xrm:0000205] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > N-dealkylation > N-demethylation
  - N-dealkylation [xrm:0000201] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > N-dealkylation
  - carbon–heteroatom oxidative cleavage [xrm:0000023] ← xenobiotic biotransformation > phase I > oxidation > carbon–heteroatom oxidative cleavage
  - dealkylation [xrm:0000200] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation
  - carbinolamine cleavage [xrm:0004006] ← xenobiotic biotransformation > process facet > carbinolamine cleavage
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - unstable oxygenation [xrm:0000011] ← xenobiotic biotransformation > phase I > unstable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:oxidative-deamination
- reactant: `CCN`
- product: `CC=O`
- tags: `chem:oxidative-deamination`
- terms:
  - oxidative deamination [xrm:0000208] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > N-dealkylation > oxidative deamination
  - N-dealkylation [xrm:0000201] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > N-dealkylation
  - carbon–heteroatom oxidative cleavage [xrm:0000023] ← xenobiotic biotransformation > phase I > oxidation > carbon–heteroatom oxidative cleavage
  - dealkylation [xrm:0000200] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation
  - carbinolamine cleavage [xrm:0004006] ← xenobiotic biotransformation > process facet > carbinolamine cleavage
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - unstable oxygenation [xrm:0000011] ← xenobiotic biotransformation > phase I > unstable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:O-demethylation
- reactant: `COc1ccccc1`
- product: `Oc1ccccc1`
- tags: `chem:O-demethylation`
- terms:
  - O-demethylation [xrm:0000209] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > O-dealkylation > O-demethylation
  - O-dealkylation [xrm:0000202] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation > O-dealkylation
  - carbon–heteroatom oxidative cleavage [xrm:0000023] ← xenobiotic biotransformation > phase I > oxidation > carbon–heteroatom oxidative cleavage
  - dealkylation [xrm:0000200] ← xenobiotic biotransformation > phase I > unstable oxygenation > dealkylation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - unstable oxygenation [xrm:0000011] ← xenobiotic biotransformation > phase I > unstable oxygenation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:acyl-glucuronidation
- reactant: `CC(=O)O`
- product: `CC(=O)O`
- tags: `chem:acyl-glucuronidation`
- terms:
  - acyl glucuronidation [xrm:0001004] ← xenobiotic biotransformation > phase II > glucuronidation > O-glucuronidation > acyl glucuronidation
  - O-glucuronidation [xrm:0001001] ← xenobiotic biotransformation > phase II > glucuronidation > O-glucuronidation
  - alcohol oxidation [xrm:0000310] ← xenobiotic biotransformation > phase I > dehydrogenation > alcohol oxidation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - acyl migration [xrm:0004017] ← xenobiotic biotransformation > process facet > acyl migration
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - glucuronidation [xrm:0001000] ← xenobiotic biotransformation > phase II > glucuronidation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - phase II [xrm:0000002] ← xenobiotic biotransformation > phase II
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:phenolic-glucuronidation
- reactant: `Oc1ccccc1`
- product: `Oc1ccccc1`
- tags: `chem:phenolic-glucuronidation`
- terms:
  - phenolic glucuronidation [xrm:0001002] ← xenobiotic biotransformation > phase II > glucuronidation > O-glucuronidation > phenolic glucuronidation
  - O-glucuronidation [xrm:0001001] ← xenobiotic biotransformation > phase II > glucuronidation > O-glucuronidation
  - glucuronidation [xrm:0001000] ← xenobiotic biotransformation > phase II > glucuronidation
  - phase II [xrm:0000002] ← xenobiotic biotransformation > phase II
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:GSH-Michael + conjugate-addition facet
- reactant: `C=CC=O`
- product: `C=CC=O`
- tags: `chem:GSH-Michael`
- terms:
  - Michael glutathionation [xrm:0001031] ← xenobiotic biotransformation > phase II > glutathionation > Michael glutathionation
  - conjugate addition [xrm:0004014] ← xenobiotic biotransformation > process facet > conjugate addition
  - glutathionation [xrm:0001030] ← xenobiotic biotransformation > phase II > glutathionation
  - phase II [xrm:0000002] ← xenobiotic biotransformation > phase II
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:quinone-formation (dearomatization + bioactivation)
- reactant: `c1ccccc1`
- product: `O=C1C=CC(=O)C=C1`
- tags: `chem:quinone-formation`
- terms:
  - quinone formation [xrm:0000300] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - bioactivation [xrm:0003000] ← xenobiotic biotransformation > bioactivation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:quinone-imine + one-step quinone formation
- reactant: `CC(=O)Nc1ccc(O)cc1`
- product: `CC(=O)N=C1C=CC(=O)C=C1`
- tags: `chem:quinone-imine`, `chem:one-step-quinone-formation`
- terms:
  - one-step quinone formation [xrm:0000307] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation > one-step quinone formation
  - quinone-imine formation [xrm:0000303] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation > quinone-imine formation
  - alcohol oxidation [xrm:0000310] ← xenobiotic biotransformation > phase I > dehydrogenation > alcohol oxidation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - quinone formation [xrm:0000300] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - bioactivation [xrm:0003000] ← xenobiotic biotransformation > bioactivation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:two-step-quinone-formation
- reactant: `c1ccccc1`
- product: `O=C1C=CC(=O)C=C1`
- tags: `chem:two-step-quinone-formation`
- terms:
  - two-step quinone formation [xrm:0000308] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation > two-step quinone formation
  - quinone formation [xrm:0000300] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - bioactivation [xrm:0003000] ← xenobiotic biotransformation > bioactivation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:imine-methide
- reactant: `Nc1ccc(C)cc1`
- product: `N=C1C=CC(=C)C=C1`
- tags: `chem:imine-methide`
- terms:
  - imine-methide formation [xrm:0000306] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation > imine-methide formation
  - quinone formation [xrm:0000300] ← xenobiotic biotransformation > phase I > dehydrogenation > quinone formation
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - bioactivation [xrm:0003000] ← xenobiotic biotransformation > bioactivation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:dearomatization alone
- reactant: `c1ccccc1`
- product: `C1=CC=CC=C1`
- tags: `chem:dearomatization`
- terms:
  - dearomatization [xrm:0004001] ← xenobiotic biotransformation > process facet > dearomatization
  - process facet [xrm:0004000] ← xenobiotic biotransformation > process facet
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:nitroaromatic-reduction
- reactant: `O=[N+]([O-])c1ccccc1`
- product: `Nc1ccccc1`
- tags: `chem:nitroaromatic-reduction`
- terms:
  - nitroaromatic reduction [xrm:0000513] ← xenobiotic biotransformation > phase I > reduction > nitrogen reduction > nitro reduction > nitroaromatic reduction
  - nitro reduction [xrm:0000511] ← xenobiotic biotransformation > phase I > reduction > nitrogen reduction > nitro reduction
  - nitrogen reduction [xrm:0000510] ← xenobiotic biotransformation > phase I > reduction > nitrogen reduction
  - reduction [xrm:0000014] ← xenobiotic biotransformation > phase I > reduction
  - bioactivation [xrm:0003000] ← xenobiotic biotransformation > bioactivation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:cyanide-hydrolysis
- reactant: `CC#N`
- product: `CC(=O)O`
- tags: `chem:cyanide-hydrolysis`
- terms:
  - cyanide hydrolysis [xrm:0000409] ← xenobiotic biotransformation > phase I > hydrolysis > cyanide hydrolysis
  - hydrolysis [xrm:0000013] ← xenobiotic biotransformation > phase I > hydrolysis
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:carbonyl-reduction
- reactant: `CC(=O)C`
- product: `CC(O)C`
- tags: `chem:carbonyl-reduction`
- terms:
  - carbonyl reduction [xrm:0000524] ← xenobiotic biotransformation > phase I > reduction > oxygen reduction > carbonyl reduction
  - oxygen reduction [xrm:0000520] ← xenobiotic biotransformation > phase I > reduction > oxygen reduction
  - reduction [xrm:0000014] ← xenobiotic biotransformation > phase I > reduction
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:glycine-conjugation
- reactant: `c1ccccc1C(=O)O`
- product: `c1ccccc1C(=O)O`
- tags: `chem:glycine-conjugation`
- terms:
  - alcohol oxidation [xrm:0000310] ← xenobiotic biotransformation > phase I > dehydrogenation > alcohol oxidation
  - carbon oxidation [xrm:0000021] ← xenobiotic biotransformation > phase I > oxidation > carbon oxidation
  - glycine conjugation [xrm:0001051] ← xenobiotic biotransformation > phase II > amino acid conjugation > glycine conjugation
  - amino acid conjugation [xrm:0001050] ← xenobiotic biotransformation > phase II > amino acid conjugation
  - dehydrogenation [xrm:0000012] ← xenobiotic biotransformation > phase I > dehydrogenation
  - oxidation [xrm:0000020] ← xenobiotic biotransformation > phase I > oxidation
  - phase I [xrm:0000001] ← xenobiotic biotransformation > phase I
  - phase II [xrm:0000002] ← xenobiotic biotransformation > phase II
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

## chem:tautomerization
- reactant: `CC(=O)C`
- product: `CC(O)=C`
- tags: `chem:tautomerization`
- terms:
  - tautomerization [xrm:0002000] ← xenobiotic biotransformation > tautomerization
  - xenobiotic biotransformation [xrm:0000000] ← xenobiotic biotransformation

