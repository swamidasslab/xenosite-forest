# Standard Library
import csv
import itertools
import re
from collections import OrderedDict, defaultdict

# Third Party
from . import rules as all_rules
from .base import (AtomTracker, ReactionRule, can_smi,
                                        clean)
from rdkit import Chem, rdBase

# Prevents spammy rdkit messages
rdBase.DisableLog('rdApp*')

RULESETS = {}


def load_ruleset(ruleset):

    if isinstance(ruleset, str):
        if ' ' in ruleset:
            RST = load_ruleset(ruleset.split())
        elif '.' not in ruleset:
            if ruleset in RULESETS:
                RST = RULESETS[ruleset]
            else:
                raise ValueError("Valid rulesets names are %s, received %s" %
                                 (" ".join(sorted(RULESETS)), ruleset))
        else:
            ruleset_name, rule_name = ruleset.split('.')
            base_ruleset = RULESETS[ruleset_name]

            rule = [x for x in base_ruleset if x.name == rule_name]

            if not rule:
                raise ValueError(
                    "Valid rule names for %s are %s, instead rule name was \'%s\'"
                    % (ruleset_name, str([x.name for x in base_ruleset]),
                       rule_name))

            RST = RuleSet(rule, name=' '.join([x.name for x in rule]))

    elif isinstance(ruleset, (list, tuple)):
        if len(ruleset) == 1:
            RST = load_ruleset(ruleset[0])
        else:
            RST = RuleSet(
                rules=[load_ruleset(x) for x in ruleset], name='custom')

    elif isinstance(ruleset, RuleSet):
        RST = ruleset

    else:
        raise ValueError(
            "Must submit string, ruleset instance, or list/tuple of either.")

    return RST


def unique_rulenames(ruleset):
    return set([x.name for x in load_ruleset(ruleset)])


class Phase1Site(object):
    def __init__(self,
                 sdf_site_fields=['DH', 'HD', 'RD', 'SO', 'UO'],
                 index_root=1,
                 strict=False,
                 *args,
                 **kwargs):

        self.index_root = index_root
        self.strict = strict
        self.sdf_site_fields = sdf_site_fields

        super(Phase1Site, self).__init__(*args, **kwargs)

    def convert_site_to_phase1_format(self, rxnname, site):
        sites_on = getattr(all_rules, rxnname.split('_')[0])().phase1_sites_on

        site = [x + self.index_root for x in site]
        try:
            site_func = getattr(self, '%s_sites' % sites_on)
        except AttributeError as err:
            err.args = ("Need to add %s_sites function in class %s" %
                        (sites_on, self.name)),
            raise

        return site_func(site)

    def add_topologically_equivalent_sites(self, mol, colors_to_sites):
        """

        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("CCC")
        >>> original_sites = {"field":frozenset([0]), "field2" : frozenset([1])}
        >>> Phase1Site().add_topologically_equivalent_sites(mol, original_sites)
        {'field': frozenset({0, 2}), 'field2': frozenset({1})}

        """
        out = {}
        atom_indexes_to_topological_ids = AtomTracker.topol_equiv(mol)

        for name, sites in list(colors_to_sites.items()):
            topological_ids_of_annotated_sites = [
                atom_indexes_to_topological_ids[idx] for idx in sites
            ]
            out[name] = frozenset([
                a for a, t in list(atom_indexes_to_topological_ids.items())
                if t in topological_ids_of_annotated_sites
            ])
        return out

    def extract_annotated_sites(self, mol):
        colors_to_sites = {
            f: mol.GetProp(f)
            for f in self.sdf_site_fields if mol.HasProp(f)
        }
        colors_to_sites = {
            f: re.split(r'\s|,', v)
            for f, v in list(colors_to_sites.items()) if v
        }
        colors_to_sites = {
            f: frozenset([int(x) - self.index_root for x in v])
            for f, v in list(colors_to_sites.items())
        }

        if colors_to_sites:
            colors_to_sites['all'] = frozenset.union(*list(colors_to_sites.values()))
        else:
            colors_to_sites['all'] = frozenset([])

        return colors_to_sites

    def convert_path_to_phase1_format(self, path):
        outsites = []
        for rule_name, site in path:
            phase1_site = self.convert_site_to_phase1_format(rule_name, site)
            outsites.append((rule_name, phase1_site))
        return outsites

    def atom_hydrogen_sites(self, site):
        return frozenset(["%d.h" % s for s in site])

    def atoms_sites(self, site):
        if self.strict:
            assert len(site) == 1
        return frozenset(["%d.%d" % (s, s) for s in site])

    def bonds_sites(self, site):
        try:
            assert len(site) == 2
        except AssertionError as err:
            if self.strict:
                err.args = ("Site was", site, "needs to have length 2")
                raise
            else:
                return site
        return frozenset(["%d.%d" % tuple(site)])


class RuleSet(Phase1Site, ReactionRule):
    """Collections of rules used to generate metabolic pahtways.

    >>> from rdkit import Chem
    >>> from xenosite.metabolite.rulesets import StableOxygenationRS, DehydrogenationRS

    >>> reactant = Chem.MolFromSmiles("CC")
    >>> product = Chem.MolFromSmiles("OCCO")
    >>> smi, path, rdmols = next(StableOxygenationRS.find_path(reactant, product, depth=2, phase1=True))
    >>> [name for name, site in path]
    ['Hydroxylation', 'Hydroxylation']
    >>> [sorted(site) for name, site in path]
    [['1.h'], ['3.h']]

    >>> reactant = Chem.MolFromSmiles("OC=CC=CC=CC=CN")
    >>> product = Chem.MolFromSmiles('N=CC=CC=CC=CC=O')
    >>> smi, path, rdmols = next(DehydrogenationRS.find_path(reactant, product, phase1=True))
    >>> path[0][0]
    'Dehydrogenation'
    >>> sorted(path[0][1])
    ['1.h', '10.h']

    >>> reactant = Chem.MolFromSmiles("C=C")
    >>> product = Chem.MolFromSmiles('C1OC1')
    >>> smi, path, rdmols = next(StableOxygenationRS.find_path(reactant, product, phase1=True))
    >>> path[0][0], sorted(path[0][1])
    ('Epoxidation', ['1.2'])


   """

    def __init__(self, rules=None, name=None, longname=None, *args, **kwargs):

        super(RuleSet, self).__init__(*args, **kwargs)

        self.name = ''
        if name is None:
            if len(rules) == 1:
                rule = rules[0]
                if isinstance(rule, str):
                    self.name = rule
                elif isinstance(rule, (ReactionRule, RuleSet)):
                    self.name = rule.name
                else:
                    self.name = ''
            else:
                self.name = ''
        else:
            self.name = name

        if longname is not None:
            self.longname = longname
        else:
            self.longname = self.name

        global RULESETS

        self.rulenames = []
        self.rules = []

        # for rule in itertools.chain(*self):
        if rules is not None:
            for rule in rules:

                if isinstance(rule, str):
                    self.rulenames.append(rule)
                    if rule in RULESETS:
                        self.rules.append(RULESETS[rule])
                    elif rule in RULES:
                        self.rules.append(RULES[rule])
                    else:
                        try:
                            self.rules.append(getattr(all_rules, rule))
                        except AttributeError:
                            raise ValueError(
                                "This rule was not found: %s." % rule)

                elif hasattr(rule, 'rules'):
                    self.rules.extend(rule.rules)
                    self.rulenames.extend(rule.rulenames)
                elif hasattr(rule, 'name'):
                    self.rules.append(rule)
                    self.rulenames.append(rule.name)
                else:
                    rule_instance = rule()
                    self.rules.append(rule_instance)
                    self.rulenames.append(rule_instance.name)

        # self.rulenames = [x.name for x in itertools.chain(*self)]

        RULESETS[self.name] = self

        # self.rules = list(itertools.chain(*self))

    def __iter__(self):
        for rule in itertools.chain(*self.rules):
            yield rule

    def add_rule(self, rule):
        if isinstance(rule, RuleSet):
            [self.add_rule(x) for x in rule]
        else:
            rule.ruleset_name = self.name
            self.rules.append(rule)

        self.rulenames = [x.name for x in itertools.chain(*self)]

    def __call__(self, start, end=None, **kwargs):
        return self.find_path(start, end_mol=end, **kwargs)

    def initialize_metabolite_table(self,
                                    openfilehandle,
                                    depth=1,
                                    delimiter='\t'):
        columns = ['Substrate', 'Metabolite']
        while depth > 1:
            columns.append('Intermediate_%d' % (depth - 1))
            depth -= 1
        columns.append('Rules')
        columns.append('Sites')
        columns.append('Molecule')
        columns.append('Model')
        writer = csv.writer(openfilehandle, delimiter=delimiter)
        writer.writerow(columns)
        return writer

    def write_metabolite_table_row(self, rules_and_sites, mols, writer):
        row = []

        rdmol = mols.pop(0)

        if not rdmol.HasProp('_Name'):
            name = 'Molecule%d' % 1
        else:
            name = rdmol.GetProp('_Name')

        row.append(Chem.MolToSmiles(rdmol))
        row.append(Chem.MolToSmiles(mols.pop()))
        for mol in mols:
            row.append(Chem.MolToSmiles(mol))
        rules, sites = list(zip(*rules_and_sites))
        row.append(';'.join(rules))
        row.append(';'.join(map(str, sites)))
        row.append(name)
        row.append(self.longname)

        writer.writerow(row)

    def find_path(self,
                  start,
                  end_mol=None,
                  outmols=True,
                  depth=1,
                  metabolite_table='',
                  openfilehandle=None,
                  **kwargs):
        """Find path between start and end mol.

        Args:
            start: an rdkit.Mol instance
            end_mol: an rdkit.Mol instance. If None, then all metabolites will be generated.
            all_paths: if True, then generate all found paths between start and end (default False).
            outmols: output rdmols and sites in addition to string
            termination_rulenames: A list of rulenames that terminate the search.

        """

        if end_mol:
            desired_endpoint_structures = set(can_smi(rdmol=end_mol))
        else:
            desired_endpoint_structures = set([])

        seen = set()

        if metabolite_table and openfilehandle is None:
            openfilehandle = open(metabolite_table, 'w')

        if metabolite_table or openfilehandle:
            openfilehandle = open(metabolite_table, 'w')
            writer = self.initialize_metabolite_table(
                openfilehandle, depth=depth)

        for smi, sites, mols in self._metabolite_paths_bfs(
            [(start, [], [start])],
                desired_endpoint_structures=desired_endpoint_structures,
                format_output_site=False,
                depth=depth,
                **kwargs):
            
            fsite = sorted(sites)

            fline = str((list(map(Chem.MolToSmiles, mols)), fsite))

            if fline not in seen:
                seen.add(fline)

                if metabolite_table:
                    self.write_metabolite_table_row(sites, mols, writer)

                if outmols:
                    yield smi, sites, mols
                else:
                    yield fline
        if metabolite_table:
            openfilehandle.close()

    def _process_sites(self, sites):
        return eval(sites)

    def _site_match(self, specified_sites, outputted_site, ruleset_name=None):
        if isinstance(specified_sites, frozenset):
            return bool(specified_sites & outputted_site)

        elif isinstance(specified_sites, dict):
            if ruleset_name not in specified_sites:
                return False

            return bool(specified_sites[ruleset_name] & outputted_site)

    def metabolites(self, mol, sites=None, unique=False, strict=True,
                    **kwargs):
        """

        >>> RS = RuleSet([all_rules.Hydroxylation()])
        >>> mol = Chem.MolFromSmiles('CC(C)C')

        >>> [x[0] for x in RS.metabolites(mol)]
        [('Hydroxylation', frozenset({0})), ('Hydroxylation', frozenset({1}))]

        >>> [x[0] for x in RS.metabolites(mol,sites=frozenset([0]))]
        [('Hydroxylation', frozenset({0}))]

        >>> [x[0] for x in RS.metabolites(mol,sites=frozenset([1]))]
        [('Hydroxylation', frozenset({1}))]

        >>> [x[0] for x in RS.metabolites(mol,sites=frozenset([0,1]))]
        [('Hydroxylation', frozenset({0})), ('Hydroxylation', frozenset({1}))]

        >>> RS2 = RuleSet([all_rules.Dealkylation(),all_rules.Hydroxylation()])
        >>> mol = Chem.MolFromSmiles('CC(C)C')

        >>> [x[0] for x in RS2.metabolites(mol,sites={'Hydroxylation':frozenset([0])})]
        [('Hydroxylation', frozenset({0}))]

        >>> [x[0] for x in RS2.metabolites(mol)]
        [('Dealkylation', frozenset({0, 1})), ('Dealkylation', frozenset({0, 1})), ('Dealkylation', frozenset({0, 1})), ('Hydroxylation', frozenset({0})), ('Hydroxylation', frozenset({1}))]


        """
        if unique:
            seen = []

        if sites:
            if isinstance(sites, str):
                sites = self._process_sites(sites=sites)

        for rule in self.rules:

            for (rxnname, site), metabolites in rule.metabolize(
                    mol, strict=strict, **kwargs):

                if sites:
                    if not self._site_match(
                            sites, site, ruleset_name=rule.name):
                        continue

                if unique:
                    fline = '_'.join([
                        '.'.join(map(Chem.MolToSmiles, metabolites)),
                        rule.name,
                        str(site)
                    ])
                    if fline in seen:
                        continue
                    seen.append(fline)

                if site == frozenset([0, 3]):
                    pass
                    #

                yield (rxnname, site), metabolites

    def _prepare_next_mols(self,
                           next_product,
                           path,
                           next_step,
                           product_paths,
                           next_mols=None):

        if next_mols is None:
            next_mols = []

        next_mols.append((next_product, path + [next_step],
                          product_paths + [next_product]))

    def _prep_output(self,
                     new_canonical_product,
                     path,
                     next_step,
                     product_paths,
                     next_product,
                     phase1=False):

        outpath = path + [next_step]
        if phase1:
            outpath = self.convert_path_to_phase1_format(outpath)

        outmols = product_paths + [next_product]

        return new_canonical_product, self.format_site(outpath), outmols

    def _metabolite_paths_bfs(self,
                              resmols_and_paths,
                              desired_endpoint_structures=None,
                              depth=1,
                              current_level=1,
                              all_paths=False,
                              phase1=False,
                              limit_to_phase1_sites_in_sdf=False,
                              termination_rulenames=[],
                              quit_if_not_ended_in_termination_rulenames=True,
                              strict=True,
                              **kwargs):

        if isinstance(depth, str):
            depth = int(depth)

        if limit_to_phase1_sites_in_sdf:
            if depth != 1:
                raise ValueError(
                    "If limit_to_phase1_sites_in_sdf=True, then depth must be 1."
                )

            start = resmols_and_paths[0][0]

            colors_to_sites = self.extract_annotated_sites(start)

            colors_to_sites = self.add_topologically_equivalent_sites(
                start, colors_to_sites)

        # If there are not desired_endpoint_structures,
        # the search will continue until the depth limit
        # is reached.
        if desired_endpoint_structures is None:
            desired_endpoint_structures = []

        next_mols = []

        for resmol, path, product_paths in resmols_and_paths:

            if not resmol:
                continue

            if quit_if_not_ended_in_termination_rulenames and termination_rulenames and path:
                rulename = self.format_site(path[-1])[0]
                if rulename in termination_rulenames:
                    continue

            for next_step, next_products in self.metabolize(
                    resmol, strict=strict, **kwargs):

                if limit_to_phase1_sites_in_sdf:
                    next_rule, next_site = next_step

                    if not next_site & colors_to_sites['all']:
                        continue

                for next_product in clean(next_products):

                    self._prepare_next_mols(next_product, path, next_step,
                                            product_paths, next_mols)

                    new_canonical_product = can_smi(rdmol=next_product)
                    if not new_canonical_product:
                        continue
                    new_canonical_product = new_canonical_product[0]

                    if not desired_endpoint_structures or new_canonical_product in desired_endpoint_structures:

                        if quit_if_not_ended_in_termination_rulenames:
                            if termination_rulenames and self.format_site(
                                    next_step)[0] not in termination_rulenames:
                                continue

                        yield self._prep_output(
                            new_canonical_product,
                            path,
                            next_step,
                            product_paths,
                            next_product,
                            phase1=phase1)

                    if new_canonical_product in desired_endpoint_structures and not all_paths:
                        desired_endpoint_structures.remove(
                            new_canonical_product)

                        if not desired_endpoint_structures: return

        if current_level < depth:

            for pro, pat, propath in self._metabolite_paths_bfs(
                    next_mols,
                    desired_endpoint_structures=desired_endpoint_structures,
                    depth=depth,
                    current_level=current_level + 1,
                    all_paths=all_paths,
                    phase1=phase1,
                    termination_rulenames=termination_rulenames,
                    quit_if_not_ended_in_termination_rulenames=
                    quit_if_not_ended_in_termination_rulenames,
                    **kwargs):

                yield pro, pat, propath


ConjugationRS = RuleSet(name='CJ', longname='Conjugation')
ConjugationRS.add_rule(all_rules.Acetylation())
ConjugationRS.add_rule(all_rules.Glucuronidation())
ConjugationRS.add_rule(all_rules.Glutathionation())
ConjugationRS.add_rule(all_rules.Sulfation())

DehydrogenationRS = RuleSet(name='DH', longname='Dehydrogenation')
DehydrogenationRS.add_rule(all_rules.Dehydrogenation())

HydrolysisRS = RuleSet(name='HD', longname='Hydrolysis')
HydrolysisRS.add_rule(all_rules.Dephosphorylation())
HydrolysisRS.add_rule(all_rules.EpoxideOpening())
HydrolysisRS.add_rule(all_rules.Hydrolysis())
HydrolysisRS.add_rule(all_rules.AzoSplitting())

QuinoneFormationRS = RuleSet(name='QF', longname='Quinone Formation')
QuinoneFormationRS.add_rule(all_rules.QuinoneFormation())

ReductionRS = RuleSet(name='RD', longname='Reduction')
ReductionRS.add_rule(all_rules.BenzodioxoleReduction())
ReductionRS.add_rule(all_rules.Dehydration())
ReductionRS.add_rule(all_rules.Hydrogenation())
ReductionRS.add_rule(all_rules.NitrogenReduction())
ReductionRS.add_rule(all_rules.SulfurReduction())
ReductionRS.add_rule(all_rules.OxygenReduction())
ReductionRS.add_rule(all_rules.ReductiveDehalogenation())

StableOxygenationRS = RuleSet(name='SO', longname='Stable Oxygenation')
StableOxygenationRS.add_rule(all_rules.Epoxidation())
StableOxygenationRS.add_rule(all_rules.Hydroxylation())
StableOxygenationRS.add_rule(all_rules.NitrogenOxidation())
StableOxygenationRS.add_rule(all_rules.SulfurOxidation())

TautomerizationRS = RuleSet(name='TT', longname='Tautomerization')
TautomerizationRS.add_rule(all_rules.Tautomerization())

UnstableOxygenationRS = RuleSet(name='UO', longname='Unstable Oxygenation')
UnstableOxygenationRS.add_rule(all_rules.Dealkylation())
UnstableOxygenationRS.add_rule(all_rules.OxidativeDehalogenation())

Full = RuleSet(
    rules=[
        ConjugationRS, DehydrogenationRS, HydrolysisRS, QuinoneFormationRS,
        ReductionRS, StableOxygenationRS, TautomerizationRS,
        UnstableOxygenationRS
    ],
    name='Full')

Bioactivation = RuleSet(name='BA', longname='BioactivationPathways')
Bioactivation.add_rule(all_rules.QuinoneFormation())
Bioactivation.add_rule(all_rules.Epoxidation())
Bioactivation.add_rule(all_rules.NitroaromaticReduction())
Bioactivation.add_rule(all_rules.ThiopheneSulfurOxidation())

ThiopheneSulfurOxidationRS = RuleSet(
    name='TSO', longname='ThiopheneSulfurOxidation')
ThiopheneSulfurOxidationRS.add_rule(all_rules.ThiopheneSulfurOxidation())

RULES = {r.name: r for r in itertools.chain(*Full)}

# RULESETS = {r.name: r for r in Full.rules}
# RULESETS['Full'] = Full

RULES2RULESETS = {
    rule.name: ruleset.longname
    for ruleset in list(RULESETS.values()) for rule in ruleset.rules
    if ruleset.name != 'Full'
}

RULES2SITETYPE = {
    rule.name: rule.sites_on
    for ruleset in list(RULESETS.values()) for rule in ruleset.rules
    if ruleset.name != 'Full'
}

RULES2PHASEONESITETYPE = {
    rule.name: rule.phase1_sites_on
    for ruleset in list(RULESETS.values()) for rule in ruleset.rules
    if ruleset.name != 'Full'
}

for rule in list(RULESETS.values()):
    RULESETS[rule.longname.replace(' ', '')] = rule

if __name__ == '__main__':
    import doctest
    doctest.testmod()
