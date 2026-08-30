"""Define specific reaction rules. """

# Standard Library
import itertools

# Third Party
from rdkit import Chem, rdBase
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import GetMolFrags

# First Party
from .base import (
    AromaticSystems,
    AtomTracker,
    ResonancePairRule,
    ResonanceRule,
    SmartsReactionRule,
    clean,
)
from .utils import canon_smi

# Prevents spammy rdkit messages
rdBase.DisableLog("rdApp.*")


class QuinoneFormation(AromaticSystems, ResonancePairRule):
    """Produces quinones based on the resonance structures of the input molecule.

    All resonance are considered for the input molecule. Consequently, the generated structures
    are not dependent on whether input bonds are marked as aromatic, double/single.

    For example, consider the kekulized and aromatic forms of acetaminophen (APAP).

    >>> QFG = QuinoneFormation()
    >>> apap_arom = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
    >>> apap_kek = Chem.MolFromSmiles('CC(=O)NC1=CC=C(C=C1)O')

    Both input molecules produce the same quinone, NAPQI, for the site [4,7].

    >>> canon_smi(next(QFG.metabolites_from_sites(apap_arom,[4,7], just_smiles=True))) == canon_smi(['CC(=O)N=C1C=CC(=O)C=C1'])
    True
    >>> canon_smi(next(QFG.metabolites_from_sites(apap_kek,[4,7], just_smiles=True))) == canon_smi(['CC(=O)N=C1C=CC(=O)C=C1'])
    True

    This rule considers pairs of ring carbon to be sites of quinone formation.
    >>> site, metabolites = next(QFG.metabolites_from_sites(apap_arom,[4,7]))
    >>> site
    ('QuinoneFormation', frozenset({4, 7}))


    Long range quinones are also produced. For example, the polycyclic aromatic
    hydrocarbon Benzo[a]pyrene has possible quinone structures that are separated across
    four rings, do to this large molecule's entirely aromatic structure.

    >>> benzopyrene = Chem.MolFromSmiles('c1ccc2c(c1)cc3ccc4cccc5c4c3c2cc5')
    >>> canon_smi(next(QFG.metabolites_from_sites(benzopyrene,[0,12], just_smiles=True))) == canon_smi(['O=C1C=CC2=C3C=CC4=CC(=O)C=C5C=CC(=CC2=C1)C3=C54'])
    True

    The loss of aromaticity inherent in quinone formation requires that there be an
    odd number of atoms betweeen each pair of carbons. For example, in the previous
    example, atoms 0 and 12 are separated by 9 bonds. If instead atoms 0 and 11 are
    submitted, no structures are produced, because the rule could not find a valid
    path between these two atoms.

    >>> list(QFG.metabolites_from_sites(benzopyrene,[0,11], just_smiles=True))
    []

    mol = Chem.MolFromSmiles('c1(C(=C)Cc3ccncc3)nc2c(O)c(N(C)C)c(NC)c(N)c2c(C(C#C)C)c(CC)1')


    """

    def __init__(self):
        super(QuinoneFormation, self).__init__(
            query_smarts=[
                ("single2double", "[#6R:1][#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#8H:2]"),
                ("addO", "[#6D2H1:1]"),
                (("replaceHalogenWithO", "single2double"), "[#6H0R:1]-[F,Cl,Br,I:2]"),
                (("single2double", "addPlus1"), "[#6H0R:1][#7D3:2]"),
                (("dealk", "single2double"), "[#6R:1][#7,#8:2][#6:3]"),
            ],
            systems="aromatic",
        )

    def metabolites(self, mol, attach_matches=False, **kwargs):

        # I tried moving this to within each downstream function, but this increased the run time
        # of the test suite so it seems most efficient to precompute the standardized molecule.
        template_mol = self.standardize(mol)

        # If a template_mol could not be successfully constructed due to problems in the input,
        # exit.
        if not template_mol:
            return

        # Preconstruct a dict mapping from each atom index to a list of tuples specifiying the atom
        # indexes of all rings containing that atom.
        rings = self.rings(template_mol)

        matches = self.match_queries(template_mol)
        num = 0

        for res_struct, pair, path in self.resonate_with_pair_paths(
            template_mol, valid_atoms=set(matches)
        ):

            if [self.in_ring_size(rings, a, 6) for a in pair] == [True, True]:
                if len(path) % 2:  # can't form meta quinones on 6-membered rings
                    continue

            for matches1, matches2 in itertools.product(*[matches[x] for x in pair]):

                num += 1

                product = self.apply_modifications(
                    res_struct, [matches1, matches2], rings=rings
                )

                if not product:
                    continue

                self.swap_bonds_along_path(product, path, **kwargs)

                if attach_matches:
                    product.SetProp("matches", str((matches1, matches2)))

                outsite = self.name + "_%d" % num, frozenset(pair)
                products = self.tag_quinone_fragments(product)

                yield outsite, clean(products)

    def tag_quinone_fragments(self, mol):
        frags = list(GetMolFrags(mol, asMols=True, sanitizeFrags=False))
        out = []
        if len(frags) > 1:
            for frag in frags:

                prop_names = AtomTracker.all_atom_prop_names(frag)
                if "dealk-noncarbon" in prop_names:
                    frag.SetProp("Quinone", "False")
                else:
                    frag.SetProp("Quinone", "True")
                out.append(frag)
        elif frags:
            frags[0].SetProp("Quinone", "True")
            out.append(frags[0])
        return out

    def modify(self, emol, mapid2atomidx, modifications, rings):

        if isinstance(modifications, str):
            modifications = [modifications]

        modifications = set(modifications)

        try:
            assert modifications.issubset(self.valid_modifications)

        except AssertionError:
            raise ValueError(
                "Invalid modifications %s. Valid modifications: %s."
                % (
                    " ".join(modifications - self.valid_modifications),
                    " ".join(self.valid_modifications),
                )
            )

        if "single2double" in modifications:
            if rings[mapid2atomidx[1]] == rings[mapid2atomidx[2]]:
                return False
            self.change_bond(emol, mapid2atomidx[1], mapid2atomidx[2], new_bond_type=2)

        if "addO" in modifications:
            self.add_atom(emol, mapid2atomidx[1], 8, 2)

        if "replaceHalogenWithO" in modifications:
            self.replace_atom(emol, mapid2atomidx[2], new_atomic_num=8)

        if "addPlus1" in modifications:
            self.set_charge(emol, mapid2atomidx[2], 1)

        if "dealk" in modifications:
            AtomTracker.set_atom_prop(emol, mapid2atomidx[3], "dealk-noncarbon", 1)
            self.break_bond(emol, mapid2atomidx[2], mapid2atomidx[3])

        return True


class Dehydrogenation(ResonancePairRule):
    """Dehydrogenate

    >>> D = Dehydrogenation()
    >>> mol = Chem.MolFromSmiles('OC=CC=CC=CC=CN')
    >>> canon_smi(next(D.metabolites_from_sites(mol,[9,0],just_smiles=True))) == canon_smi(['N=CC=CC=CC=CC=O'])
    True

    """

    def __init__(self):
        super(Dehydrogenation, self).__init__(
            rxns=[
                "[#16v4:1]-[Oh:2]>>[*:1]=[*:2]",
                "[#6h:1]-[#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#7D3,#8H1:2]>>[*:1]=[*:2]",
            ],
            query_smarts=[
                ("single2double", "[#6h:1][#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#8H:2]"),
                (("single2double", "addPlus1"), "[#6h:1][#7D3:2]"),
            ],
            phase1_sites_on="atom_hydrogen",
            sites_on="atom_pairs",
        )

    def modify(self, emol, mapid2atomidx, modifications):

        if "single2double" in modifications:
            self.change_bond(emol, mapid2atomidx[1], mapid2atomidx[2], new_bond_type=2)

        if "addPlus1" in modifications:
            self.set_charge(emol, mapid2atomidx[2], 1)

        return True

    def site_from_matches(self, matches):
        site = []
        for idx2mapids, query in matches:
            site.append(idx2mapids[2])
        return frozenset(site)

    def metabolites(self, mol, **kwargs):

        for site, metabolites in super(Dehydrogenation, self).metabolites(mol):
            yield site, metabolites

        # I tried moving this to within each downstream function, but this increased the run time
        # of the test suite so it seems most efficient to precompute the standardized molecule.
        template_mol = self.standardize(mol)

        # If a template_mol could not be successfully constructed due to problems in the input,
        # exit.
        if not template_mol:
            return

        matches = self.match_queries(template_mol)

        for res_struct, pair, path in self.resonate_with_pair_paths(
            template_mol, valid_atoms=set(matches)
        ):

            if not set(pair).issubset(set(matches)):
                continue

            for matches1, matches2 in itertools.product(*[matches[x] for x in pair]):

                product = self.apply_modifications(res_struct, [matches1, matches2])

                if not product:
                    continue

                self.swap_bonds_along_path(product, path)
                yield (self.name, self.site_from_matches([matches1, matches2])), clean(
                    product
                )
                # yield (self.name, frozenset(pair)), clean(product)


class Hydrogenation(ResonancePairRule):
    """Swaps single/double bonds for paths within a molecule's resonance structures.

    >>> mol = Chem.MolFromSmiles('CC(=O)N=C1C=CC(=O)C=C1')
    >>> site, metabolites = next(Hydrogenation().metabolites_from_sites(mol, frozenset({8, 3})))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(=O)NC1=CC=C(O)C=C1')
    True

    """

    def __init__(self):
        super(Hydrogenation, self).__init__(
            rxns=["[#6:1]#[#6:2]>>[*:1]=[*:2]", "[#6:1]=[#6:2]>>[*:1]-[*:2]"],
            phase1_sites_on="atoms",
            query_smarts=[("", "[*:1]")],
        )

    def metabolites(self, mol, **kwargs):

        for site, metabolites in super(Hydrogenation, self).metabolites(mol):
            yield site, metabolites

        for res_struct, pair, path in self.resonate_with_pair_paths(mol):

            # For all bonds in bond_path, replace single bonds with double
            # bonds and double bonds with single bonds.
            if len(path) < 2:
                continue

            self.swap_bonds_along_path(res_struct, path)

            yield (self.name, frozenset(pair)), clean(res_struct)


class Tautomerization(ResonanceRule):
    """Swap single and double bonds without changing overall hydrogen count."""

    def __init__(self):

        super(Tautomerization, self).__init__(
            query_smarts=[
                ("", "[#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#8H:1]"),
            ],
        )

    # @staticmethod
    # def hydrogen_count(mol):
    # # Chem.SanitizeMol(mol,catchErrors=True)
    # try:
    # mol_with_hydrogens =  Chem.AddHs(mol)
    # except ValueError:
    # return None
    # return mol_with_hydrogens.GetNumAtoms() - mol_with_hydrogens.GetNumHeavyAtoms()

    @staticmethod
    def no_double_bonds_per_atom(mol):
        return False not in [
            [b.GetBondType() for b in a.GetBonds()].count(Chem.rdchem.BondType.DOUBLE)
            <= 1
            for a in mol.GetAtoms()
        ]

    def metabolites(self, mol, **kwargs):
        # input_hydrogens = self.hydrogen_count(mol)

        matches = set(self.match_queries(mol))

        for res_struct, pair, path in self.resonate_with_pair_paths(mol):

            for idx in pair:
                for endpoint in matches & self.neighbors(mol, idx) - set(path):

                    if idx == path[0]:
                        full_path = [endpoint] + path
                    else:
                        full_path = path + [endpoint]

                    # make a copy
                    metabolite = Mol(res_struct)

                    self.swap_bonds_along_path(metabolite, full_path)

                    # output_hydrogens = self.hydrogen_count(metabolite)
                    if not self.no_double_bonds_per_atom(metabolite):
                        continue

                    # if output_hydrogens is None:
                    # continue

                    # if output_hydrogens != input_hydrogens:
                    # continue

                    yield (self.name, frozenset([full_path[0], full_path[-1]])), [
                        metabolite
                    ]


class Epoxidation(ResonanceRule):
    """Adds epoxides to carbon-carbon and carbon-nitrogen bonds.

    Importantly, this rule automatically considers all resonance
    structures for the input molecule. Consequently, the generated structures
    are not dependent on whether input bonds are marked as aromatic, double,
    or single.

    For example, consider the kekulized and aromatic forms of the same input
    molecule:

    >>> aromatic = Chem.MolFromSmiles('c1c(Cl)cccc1C=C')
    >>> kekulized = Chem.MolFromSmiles('C1=C(Cl)C=CC=C1C=C')

    Both input molecules produce the same structure for the input site [4,5]
    despite that site being specified as a single bond in the kekulized input:

    >>> E = Epoxidation()
    >>> site,metabolites = next(E.metabolites_from_sites(aromatic,[4,5]))
    >>> canon_smi(metabolites[0]) == canon_smi('C=CC1=CC(Cl)=CC2OC12')
    True
    >>> site,metabolites = next(E.metabolites_from_sites(kekulized,[4,5]))
    >>> canon_smi(metabolites[0]) == canon_smi('C=CC1=CC(Cl)=CC2OC12')
    True

    """

    def __init__(self):
        super(Epoxidation, self).__init__(
            rxns="[#6:1]=[#6,#7:2]>>[*:1]1-[*:2][O]1", sites_on="bonds"
        )


class Acetylation(SmartsReactionRule):
    """Adds an acetyl to OH, NH, and SH.

    >>> mol = Chem.MolFromSmiles('C1=CC(=C(C=C1N)C(=O)O)O')
    >>> site, metabolites = next(Acetylation().metabolize(mol))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(=O)NC1=CC(C(=O)O)=C(O)C=C1')
    True
    """

    smarts = "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]"


class Glutathionation(SmartsReactionRule):
    """Adds glutathione to epoxide, halogen, and sulfur motifs.

    >>> G = Glutathionation()
    >>> mol = Chem.MolFromSmiles('c1ccccc1-C1OC1')
    >>> site, metabolites = next(G.metabolize(mol))
    >>> site
    ('Glutathionation', frozenset({6}))
    >>> canon_smi(metabolites, isomericSmiles=False) == canon_smi(['NC(CCC(=O)NC(CSC(CO)C1=CC=CC=C1)C(=O)NCC(=O)O)C(=O)O'], isomericSmiles=False)
    True

    """

    smarts = [
        "[#6:1]1[#8:2][#6:3]1>>\
                C(CC(=O)N[C@@H](CS([*:1][*:3][*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
        "[#6:1][Cl:2]>>C(CC(=O)N[C@@H](CS([*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
        "[#16h1:1]>>C(CC(=O)N[C@@H](CS([*:1]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
        "[#6H2:1]=[#6:2]>>C(CC(=O)N[C@@H](CS([*:1]-[*:2]))C(=O)NCC(=O)O)[C@@H](C(=O)O)N",
    ]
    mapid_site = [1]


class Dealkylation(SmartsReactionRule):
    """Breaks C-C, C-N, C-O, C-S bonds.

    >>> D = Dealkylation()
    >>> mol = Chem.MolFromSmiles('CCO')
    >>> site, metabolites = next(D.metabolize(mol))
    >>> site
    ('Dealkylation', frozenset({0, 1}))
    >>> canon_smi(metabolites) == canon_smi(['CO', 'CO'])
    True

    """

    smarts = [
        "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
        "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
        "[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
        "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)",
        "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
        "[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
        "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)",
        "[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
        "[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)",
        "[#6:1][#6:2]>>(O-[*:1].[*:2])",
        "[#6h:1][#6:2]>>(O-[*:1].[*:2])",
        "[#6h:1][#6:2]>>(O=[*:1].[*:2])",
        "[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])",
    ]
    mapid_site = [1, 2]
    sites_on = "bonds"


class Hydrolysis(SmartsReactionRule):
    """Cleaves the single bond in carbonyls and adds oxygen to the carbon.

    >>> mol = Chem.MolFromSmiles('C1=CC=CC=C1C(=O)OC(C)(C)C')
    >>> site, metabolites = next(Hydrolysis().metabolize(mol))
    >>> site
    ('Hydrolysis', frozenset({8, 6}))
    >>> canon_smi(metabolites) == canon_smi(['O=C(O)C1=CC=CC=C1', 'CC(C)(C)O'])
    True

    """

    smarts = [
        "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2](O).[*:3])",
        "[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2].[*:3])",
    ]
    mapid_site = [2, 3]
    sites_on = "bonds"


class ReductiveDehalogenation(SmartsReactionRule):
    """Cleaves halogen-carbon bonds.

    >>> RD = ReductiveDehalogenation()
    >>> mol = Chem.MolFromSmiles('ClCC(=O)C(NC(=O)c1cc(Cl)c(c(c1)Cl)CO)(CC)C')
    >>> site, metabolites = next(RD.metabolites_from_sites(mol,[0,1]))
    >>> canon_smi(metabolites) == canon_smi(['Cl', 'CCC(C)(NC(=O)C1=CC(Cl)=C(CO)C(Cl)=C1)C(C)=O'])
    True

    """

    smarts = [
        "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]",
        "[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]",
    ]
    sites_on = "bonds"


class OxidativeDehalogenation(SmartsReactionRule):
    """Replaces halogens bound to carbons with hydroxyls, double bonded oxygens, orcarboxylic acids.

    >>> mol = Chem.MolFromSmiles('ClCC(=O)C(NC(=O)c1cc(Cl)c(c(c1)Cl)CO)(CC)C')
    >>> site, metabolites = next(OxidativeDehalogenation().metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['Cl', 'CCC(C)(NC(=O)C1=CC(Cl)=C(CO)C(Cl)=C1)C(=O)CO'])
    True

    """

    smarts = [
        "[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O",
        "[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O",
        "[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O",
        "[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]",
        "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]",
        "[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]",
    ]
    sites_on = "bonds"


class NitrogenOxidation(SmartsReactionRule):
    """Adds a hydroxyl to carbons, nitrogens, and sulfurs.

    Attaches hydroxyls to:

    - Carbons with at least one hydrogen.
    - Nitrogens of valence 3.
    - Sulfurs of valence 2 or 4.

    >>> mol = Chem.MolFromSmiles('CCC(=C)C')
    >>> site, metabolites = next(Hydroxylation().metabolites_from_sites(mol,[frozenset({0})]))
    >>> canon_smi(metabolites[0]) == canon_smi('C=C(C)CCO')
    True

    """

    smarts = ["[#7v3:1]>>[*:1]O", "[#7v3:1]>>[*:1]=O"]
    phase1_sites_on = "atoms"


class SulfurOxidation(SmartsReactionRule):
    """Adds a hydroxyl to carbons, nitrogens, and sulfurs.

    Attaches hydroxyls to:

    - Carbons with at least one hydrogen.
    - Nitrogens of valence 3.
    - Sulfurs of valence 2 or 4.

    >>> mol = Chem.MolFromSmiles('O=C(c1ccc(cc1)C(C(=O)O)C)c2sccc2')
    >>> site, metabolites = next(SulfurOxidation().metabolize(mol))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(C(=O)O)C1=CC=C(C(=O)C2=CC=C[S+]2[O-])C=C1')
    True
    >>> site[1]
    frozenset({14})

    """

    smarts = [
        "[#16;v2,v4:1]>>[*&H0&+:1][O-]",
        "[#16;v2,v4:1]>>[*:1][O]",
        "[#16;v2,v4:1]>>[*:1]=O",
    ]
    phase1_sites_on = "atoms"


class Hydroxylation(SmartsReactionRule):
    """Adds a hydroxyl to carbons, nitrogens, and sulfurs.

    Attaches hydroxyls to:

    - Carbons with at least one hydrogen.
    - Nitrogens of valence 3.
    - Sulfurs of valence 2 or 4.

    >>> mol = Chem.MolFromSmiles('CCC(=C)C')
    >>> site, metabolites = next(Hydroxylation().metabolites_from_sites(mol,[frozenset({0})]))
    >>> canon_smi(metabolites[0]) == canon_smi('C=C(C)CCO')
    True

    """

    smarts = ["[#6h:1]>>[*:1]O", "[#6h2:1]>>[*:1]=O"]
    phase1_sites_on = "atom_hydrogen"


class OxygenReduction(SmartsReactionRule):
    """Converts C=O/N=O bonds to CO/NO.

    >>> mol = Chem.MolFromSmiles('Oc1ccccc1')
    >>> site, metabolites = next(Dehydration().metabolites_from_sites(mol,[frozenset({1,0})]))
    >>> canon_smi(metabolites) == canon_smi(['C1=CC=CC=C1', 'O'])
    True

    """

    smarts = ["[#8:1]=[#6,#7:2]>>[*:1]-[*:2]", "[#8:1]-[#8:2]>>[*:1].[*:2]"]


class Dehydration(SmartsReactionRule):
    """Removes OH groups.

    >>> mol = Chem.MolFromSmiles('Oc1ccccc1')
    >>> site, metabolites = next(Dehydration().metabolites_from_sites(mol,[frozenset({1,0})]))
    >>> canon_smi(metabolites) == canon_smi(['C1=CC=CC=C1', 'O'])
    True

    """

    phase1_sites_on = "bonds"
    sites_on = "bonds"

    smarts = [
        "[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]",
        "[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]",
        "[#6,#7:1]=[#8:2]>>[*:1].[*:2]",
    ]
    mapid_site = [1, 2]


class Dephosphorylation(SmartsReactionRule):
    """Breaks bonds between carbons and phosphates, and hydroxylates the carbon.

    >>> smiles = 'O=C1O[Zn]OC(=O)CN(CCN(C1)Cc1c(cnc(c1O)C)COP(=O)(O)O)Cc1c(cnc(c1O)C)COP(=O)(O)O'
    >>> mol = Chem.MolFromSmiles(smiles)
    >>> site, metabolites = next(Dephosphorylation().metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['CC1=C(O)C(CN2CCN(CC3=C(COP(=O)(O)O)C=NC(C)=C3O)CC(=O)[O][Zn][O]C(=O)C2)=C(CO)C=N1', 'O=[PH](O)O'])
    True

    """

    smarts = "[#8:1][#15:2](=[#8:3])([#8:4])[#8:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]"
    sites_on = "bonds"
    mapid_site = [1, 2]

    # '[#8:1]=[#15:2]([#8:3])([#8:4])[#8:5][#6:6]>>[*:1]=[*:2]([*:3])([*:4])[*:5].O[*:6]',


class BenzodioxoleReduction(SmartsReactionRule):
    """Breaks both carbon-oxygen bonds of the carbon between the two oxygens in the motif C1OCOC1.

    >>> mol = Chem.MolFromSmiles('COc1cc(cc(c1OC)OC)[C@@H]1c2cc3OCOc3cc2C[C@@H]2[C@@H]1C(=O)OC2')
    >>> site, metabolites = next(BenzodioxoleReduction().metabolize(mol))
    >>> canon_smi(metabolites, isomericSmiles=False) == canon_smi(['COC1=CC(C2C3=C(C=C(O)C(O)=C3)CC3COC(=O)C32)=CC(OC)=C1OC', 'C'], isomericSmiles=False)
    True

    """

    smarts = (
        "[#6R:1]-[#8R:2]-[#6H2R:3]-[#8R:4]-[#6R:5]>>([*:1]-[*:2].[*:3].[*:4]-[*:5])"
    )


class SulfurReduction(SmartsReactionRule):
    """
    Converts:
        carbon-carbon double bonds to single bonds
        carbon-oxygen double bonds to single bonds

    Breaks:
        Sulfur-Sulfur single bonds
        Sulfur-Oxygen double bonds

    To do:
        - Sulfur - phosphorus bonds.

    >>> SR = SulfurReduction()
    >>> mol = Chem.MolFromSmiles('c1ccccc1SSc1ccccc1')
    >>> site, metabolites = next(SR.metabolites_from_sites(mol,[7,6]))
    >>> site
    ('SulfurReduction', frozenset({6, 7}))
    >>> canon_smi(metabolites) == canon_smi(['SC1=CC=CC=C1', 'SC1=CC=CC=C1'])
    True

    """

    smarts = [
        "[#16:1]=[#8:2]>>[*:1].[*:2]",
        "[#16:1]-[#16:2]>>[*:1].[*:2]",
        "[#16:1]-[#6,#8:2]>>[*:1].[*:2]",
    ]


class Glucuronidation(SmartsReactionRule):
    """
    >>> mol = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
    >>> site, metabolites = next(Glucuronidation().metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['CC(=O)NC1=CC=C(OC2OC(C(=O)O)C(O)C(O)C2O)C=C1'])
    True

    """

    smarts = [
        "[#6:1][#6:2](=[O,N,P,S:3])[#8:4]>>\
                        O1C(C(=O)O)C(O)C(O)C(O)C([*:4][*:2](=[#8:3])[*:1])1",
        "[#8H1:1][#6:2]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2])1",
    ]


class Sulfation(SmartsReactionRule):
    """
    >>> mol = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
    >>> site, metabolites = next(Sulfation().metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['CC(=O)NC1=CC=C(OS(=O)(=O)O)C=C1'])
    True

    >>> site
    ('Sulfation', frozenset({8, 7}))

    """

    def __init__(self):
        super(Sulfation, self).__init__(
            [
                "[#6:1][#8:2]>>[*:1][*:2]S(=O)(=O)O",
                "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>\
                    [*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1",
            ],
            [1, 2],
        )


class AzoSplitting(SmartsReactionRule):
    """Splits N=N bonds.

    >>> mol = Chem.MolFromSmiles('OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O')
    >>> site,metabolites = next(AzoSplitting().metabolize(mol))
    >>> site
    ('AzoSplitting', frozenset({6, 7}))
    >>> canon_smi(metabolites) == canon_smi(['NC1=CC(C(=O)O)=C(O)C=C1', 'NC1=CC(C(=O)O)=C(O)C=C1'])
    True

    """

    smarts = "[#7:1]=[#7:2]>>[*:1].[*:2]"


class EpoxideOpening(SmartsReactionRule):
    """Breaks one of the carbon-oxygen bonds in epoxides.

    >>> mol = Chem.MolFromSmiles('c1ccccc1C1OC1')
    >>> site,metabolites = next(EpoxideOpening().metabolize(mol))
    >>> site
    ('EpoxideOpening', frozenset({6, 7}))
    >>> canon_smi(metabolites) == canon_smi(['OCCC1=CC=CC=C1'])
    True

    """

    smarts = [
        "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])",
        "[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)",
    ]
    mapid_site = [1, 2]
    sites_on = "bonds"


class NitrogenReduction(SmartsReactionRule):
    """Converts nitro and nitroso groups to nitrogens.

    >>> mol = Chem.MolFromSmiles('CC[C@](c1cc2c3nc4cccc(c4cc3Cn2c(=O)c1CO)N(=O)=O)(C(=O)O)O')
    >>> site, metabolites = next(NitrogenReduction().metabolize(mol))
    >>> canon_smi(metabolites, isomericSmiles=False) == canon_smi(['CCC(O)(C(=O)O)C1=C(CO)C(=O)N2CC3=CC4=C(N=O)C=CC=C4N=C3C2=C1', 'O'], isomericSmiles=False)
    True
    >>> site
    ('NitrogenReduction', frozenset({24, 23}))

    """

    smarts = [
        "[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
        "[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])",
        "[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])",
        "[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])",
        "[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])",
        "[#7:1]-[#8:2]>>([*:1].[*:2])",
        "[#7D2:1]=[#8:2]>>([*:1].[*2])",
        "[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])",
    ]
    mapid_site = [1, 2]


class NitroaromaticReduction(SmartsReactionRule):
    """Converts nitro and nitroso adjacent to a ring carbon to nitrogens.

    >>> mol = Chem.MolFromSmiles('[O-][N+](C1=CC2=C(C=C1)NC(CN=C2C3=CC=CC=C3Cl)=O)=O')
    >>> site, metabolites = next(NitroaromaticReduction().metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['O', 'O=NC1=CC2=C(C=C1)NC(=O)CN=C2C1=CC=CC=C1Cl'])
    True
    >>> site
    ('NitroaromaticReduction', frozenset({0, 1}))

    """

    smarts = [
        "[#8-1:1]-[#7+1:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
        "[#8:1]-[#7:2]([#6R:4])=[#8:3]>>[*:1].[*:2]([*:4])=[*:3]",
    ]
    mapid_site = [1, 2]


class ThiopheneSulfurOxidation(SmartsReactionRule):
    """Adds a hydroxyl to carbons, nitrogens, and sulfurs.

    Attaches hydroxyls to:

    - Carbons with at least one hydrogen.
    - Nitrogens of valence 3.
    - Sulfurs of valence 2 or 4.

    >>> mol = Chem.MolFromSmiles('O=C(c1ccc(cc1)C(C(=O)O)C)c2sccc2')
    >>> site, metabolites = next(ThiopheneSulfurOxidation().metabolize(mol))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(C(=O)O)C1=CC=C(C(=O)C2=CC=C[S+]2[O-])C=C1')
    True
    >>> site[1]
    frozenset({14})

    >>> furan = Chem.MolFromSmiles('O=C(c1ccc(cc1)C(C(=O)O)C)c2occc2')
    >>> list( ThiopheneSulfurOxidation().metabolize(furan))
    []

    """

    smarts = [
        "[#6:2]1=[#6:3][#6:4]=[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*&H0&+:1]1[O-]",
        # '[#6:2]1=[#6:3][#6:4]=[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*:1]1[O]',
        # '[#6:2]1=[#6:3][#6:4]=[#6:5][#16;v2,v4:1]1>>[*:2]1=[*:3][*:4]=[*:5][*:1]1=[O]',
    ]
    phase1_sites_on = "atoms"
    mapid_site = [1]


if __name__ == "__main__":
    import doctest

    doctest.testmod()
