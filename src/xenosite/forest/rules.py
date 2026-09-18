"""Define specific reaction rules. """

# Standard Library
import itertools
import json

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
    copy_mol,
)
from .path_context import (
    ADD_O,
    CLEAVE,
    CLEAVE_OR_ADD_O,
    NEUTRAL,
    REMOVE_O,
    FormulaAny,
    FormulaHint,
    FormulaMatch,
)
from .step_plan import Step, StepPlan, AtomRef
from .utils import (
    apply_star_conjugate,
    canon_smi,
    mol_to_cxsmiles,
    unmapped_smiles,
)

# Prevents spammy rdkit messages
rdBase.DisableLog("rdApp.*")

# Labels that are abstract (no real adduct molecule in the SMARTS).
STAR_ONLY_LABELS = frozenset({"Protein", "DNA", "Cyanide"})


class ConjugationRule(SmartsReactionRule):
    """Phase II conjugation with optional star collapse and CX ``atomLabel``.

    By default products are bare ``*`` adducts (the conjugate group is collapsed).
    Pass ``as_star=False`` for the full chemical adduct (GlcA, GSH, acetyl,
    sulfate). Labels ``Protein``, ``DNA``, and ``Cyanide`` are star-only.
    """

    as_star = True
    star_label = None

    def __init__(self, as_star=None, star_label=None, *args, **kwargs):
        if as_star is None:
            as_star = type(self).as_star
        if star_label is None:
            star_label = type(self).star_label
        if star_label in STAR_ONLY_LABELS and as_star is False:
            raise ValueError(
                "%r cannot emit a full conjugate structure; use a star adduct"
                % (star_label,)
            )
        if star_label in STAR_ONLY_LABELS:
            as_star = True
        self.as_star = as_star
        self.star_label = star_label
        super(ConjugationRule, self).__init__(*args, **kwargs)

    def is_terminal_product(self, mol) -> bool:
        from .utils import has_star_conjugate

        return bool(has_star_conjugate(mol))

    def is_redundant(self, peer_rules) -> bool:
        """Keep only the first conjugation rule in the peer list for path search."""
        for peer in peer_rules:
            if isinstance(peer, ConjugationRule):
                return peer is not self
        return False

    def enumerate_for_path(
        self,
        mol,
        ctx,
        expand_phase1_plans=True,
        include_sites=None,
        exclude_sites=None,
        **kwargs,
    ):
        """Emit unlabeled ``*`` adducts only; products are terminal."""
        old_star, old_label = self.as_star, self.star_label
        self.as_star = True
        self.star_label = None
        try:
            yield from SmartsReactionRule.enumerate_for_path(
                self,
                mol,
                ctx,
                expand_phase1_plans=False,
                include_sites=include_sites,
                exclude_sites=exclude_sites,
                **kwargs,
            )
        finally:
            self.as_star = old_star
            self.star_label = old_label

    def metabolites(self, mol, **kwargs):
        for site, products in super(ConjugationRule, self).metabolites(mol, **kwargs):
            if self.as_star:
                products = [
                    apply_star_conjugate(p, label=self.star_label) for p in products if p
                ]
            yield site, products

    def metabolites_from_sites(
        self,
        mol,
        sites,
        just_smiles=False,
        deplete_sites=False,
        tag_atoms=True,
        **kwargs
    ):
        sites = self._cast_sites(sites)

        for site, metabolite in self.metabolize(mol, **kwargs):
            if tuple(site) in sites:
                if just_smiles:
                    if self.star_label:
                        yield [mol_to_cxsmiles(m) for m in metabolite]
                    else:
                        yield list(map(unmapped_smiles, metabolite))
                else:
                    yield tuple(site), metabolite

                if deplete_sites:
                    sites.remove(tuple(site))
                    if not sites:
                        break


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

    @staticmethod
    def _mod_names(modifications):
        if isinstance(modifications, str):
            return frozenset([modifications])
        return frozenset(modifications)

    def _prep_and_dh_end(self, match, mol=None):
        """Return (prep Steps, DH endpoint AtomRef, is_methide) for one match.

        Mapids are GetIdx on ``mol`` (caller). Origin refs stamp that frame's
        depth — created atoms (e.g. OH oxygen) have no depth-0 identity and
        must use ``added_by`` instead.
        """
        from .step_plan import _frame_depth, _origin_refs

        mapids, modifications = match[0], match[1]
        names = self._mod_names(modifications)
        depth = _frame_depth(mol) if mol is not None else 0
        prep = []
        methide = False
        if "addO" in names:
            carbon = mapids[1]
            carbon_refs = _origin_refs([carbon], mol, depth=depth)
            prep.append(Step("Hydroxylation", carbon_refs))
            # added_by site idxs share the origin refs' frame.
            c_ref = next(iter(carbon_refs))
            end = AtomRef(
                added_by=("Hydroxylation", frozenset([c_ref.origin])),
                depth=c_ref.depth,
            )
        elif "replaceHalogenWithO" in names:
            at_refs = _origin_refs([mapids[1], mapids[2]], mol, depth=depth)
            prep.append(Step("OxidativeDehalogenation", at_refs))
            at = frozenset(r.origin for r in at_refs)
            # Homogeneous frame (all projected or all current).
            site_depth = next(iter(at_refs)).depth if at_refs else depth
            end = AtomRef(added_by=("OxidativeDehalogenation", at), depth=site_depth)
        elif "dealk" in names:
            prep.append(
                Step(
                    "Dealkylation",
                    _origin_refs([mapids[2], mapids[3]], mol, depth=depth),
                )
            )
            end = next(iter(_origin_refs([mapids[2]], mol, depth=depth)))
        else:
            # single2double / addPlus1: heteroatom or alkyl already on reactant
            end_idx = mapids[2] if 2 in mapids else mapids[1]
            end = next(iter(_origin_refs([end_idx], mol, depth=depth)))
            if mol is not None:
                try:
                    methide = mol.GetAtomWithIdx(int(end_idx)).GetAtomicNum() == 6
                except Exception:
                    methide = False
        return prep, end, methide

    def _plan_from_match_pair(self, match1, match2, mol=None):
        from .step_plan import And, StepPlan

        prep1, end1, m1 = self._prep_and_dh_end(match1, mol=mol)
        prep2, end2, m2 = self._prep_and_dh_end(match2, mol=mol)
        prep = prep1 + prep2
        # One alkylidene end ⇒ methide DH; both ends never (skip pathway stamp).
        pathways = frozenset()
        if (m1 or m2) and not (m1 and m2):
            pathways = frozenset({Dehydrogenation.PATHWAY_METHIDE})
        final = Step(
            "Dehydrogenation", frozenset([end1, end2]), pathways=pathways
        )
        if not prep:
            return StepPlan((final,))
        if len(prep) == 1:
            return StepPlan([prep[0], final])
        return StepPlan([And(prep), final])

    def phase1_steps(self, mol, site, _matches=None, _match_pair=None, **kwargs):
        """One Phase1-equivalent :class:`StepPlan` for a quinone site.

        Multiple SMARTS match pairs become a top-level ``Or`` inside that
        single plan. Linearizations expand every alternative total order.
        Alkyl quinone-methide ends stamp ``pathways=(\"methide\",)`` on the
        final Dehydrogenation step.
        """
        from .step_plan import StepPlan, or_of_plans

        want = frozenset(self._cast_sites(site)[0][1])
        if _match_pair is not None:
            return self._plan_from_match_pair(*_match_pair, mol=mol)

        template_mol = self.standardize(mol)
        if not template_mol:
            return StepPlan.empty()

        matches = _matches if _matches is not None else self.match_queries(template_mol)
        if len(want) != 2 or not want.issubset(matches):
            return StepPlan.empty()

        a, b = tuple(want)
        plans = []
        seen = set()
        for match1, match2 in itertools.product(matches[a], matches[b]):
            plan = self._plan_from_match_pair(match1, match2, mol=template_mol)
            key = json.dumps(plan.to_json(), sort_keys=True)
            if key in seen:
                continue
            seen.add(key)
            plans.append(plan)
        return or_of_plans(plans)

    def metabolites(self, mol, attach_matches=False, attach_phase1_steps=False, **kwargs):

        # I tried moving this to within each downstream function, but this increased the run time
        # of the test suite so it seems most efficient to precompute the standardized molecule.
        template_mol = self.standardize(mol)

        # If a template_mol could not be successfully constructed due to problems in the input,
        # exit.
        if not template_mol:
            return

        from .base import _include_atom_sets

        want = _include_atom_sets(kwargs.get("include_sites"))

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

            pair_fs = frozenset(pair)
            if want is not None and pair_fs not in want:
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
                products = clean(self.tag_quinone_fragments(product))

                if attach_phase1_steps and products:
                    plan = self._plan_from_match_pair(
                        matches1, matches2, mol=template_mol
                    )
                    for frag in products:
                        if (
                            not frag.HasProp("Quinone")
                            or frag.GetProp("Quinone") == "True"
                        ):
                            plan.attach_to_mol(frag)

                yield outsite, products

    def enumerate_for_path(self, mol, ctx, expand_phase1_plans=True, **kwargs):
        """Emit phase1 plans from aromatic pair matches — no RunReactants/clean."""
        from .path_context import oxygen_deficit
        from .step_plan import pathway_alternatives

        if not self.could_help(mol, ctx.target, ctx):
            return
        if not expand_phase1_plans:
            yield from super(QuinoneFormation, self).enumerate_for_path(
                mol, ctx, expand_phase1_plans=False, **kwargs
            )
            return

        template_mol = self.standardize(mol)
        if not template_mol:
            return
        rings = self.rings(template_mol)
        matches = self.match_queries(template_mol)
        try:
            need_o = max(0, int(oxygen_deficit(mol, ctx.target)))
        except Exception:
            need_o = None
        seen = set()
        for _res, pair, path in self.resonate_with_pair_paths(
            template_mol, valid_atoms=set(matches)
        ):
            if [self.in_ring_size(rings, a, 6) for a in pair] == [True, True]:
                if len(path) % 2:
                    continue
            site = frozenset(pair)
            if site in seen or len(site) != 2:
                continue
            if not site.issubset(matches):
                continue
            seen.add(site)
            try:
                expr = self.phase1_steps(mol, site, _matches=matches)
            except NotImplementedError:
                continue
            if not pathway_alternatives(expr):
                continue
            # Skip plans that add more O than the target needs (e.g. bis-OH when
            # phenol→BQ only needs one).
            if need_o is not None:
                max_oh = 0
                for lin in expr.linearizations():
                    max_oh = max(
                        max_oh,
                        sum(1 for s in lin.steps if s.rule == "Hydroxylation"),
                    )
                if max_oh > need_o:
                    continue
            yield ("plan", expr, site)

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
    """Dehydrogenate.

    Optional uncommon chemotypes are declared as ``query_smarts`` options
    (``pathway``, ``one_sided``, ``recommend``), not separate code branches.
    Default Phase I keeps those pathways off; enable via
    ``Dehydrogenation(pathways={\"methide\"})``, ``metabolize(..., pathways=...)``,
    a :class:`~xenosite.forest.step_plan.Step` with ``pathways=...`` (quinone
    ``phase1_steps`` stamps methide when an end is alkyl), or guided search when
    formula hints + ``recommend`` SMARTS agree.

    >>> D = Dehydrogenation()
    >>> mol = Chem.MolFromSmiles('OC=CC=CC=CC=CN')
    >>> canon_smi(next(D.metabolites_from_sites(mol,[9,0],just_smiles=True))) == canon_smi(['N=CC=CC=CC=CC=O'])
    True

    """

    PATHWAY_METHIDE = "methide"

    def __init__(self, pathways=()):
        super(Dehydrogenation, self).__init__(
            rxns=[
                ("[#16v4:1]-[Oh:2]>>[*:1]=[*:2]", {"formula_hint": NEUTRAL}),
                (
                    "[#6h:1]-[#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#7D3,#8H1:2]>>[*:1]=[*:2]",
                    {"formula_hint": NEUTRAL},
                ),
            ],
            query_smarts=[
                (
                    "single2double",
                    "[#6h:1][#6D1H3,#6D2H2,#6D3H1,#7D2H1,#7D1H2,#8H:2]",
                    {"formula_hint": NEUTRAL},
                ),
                (
                    ("single2double", "addPlus1"),
                    "[#6h:1][#7D3:2]",
                    {"formula_hint": NEUTRAL},
                ),
                # Quinoid / phenol ends: ring C without H bonded to OH or NH
                # (existing patterns require #6h and miss hydroquinone / NAPQI).
                ("single2double", "[#6H0:1][#8H:2]", {"formula_hint": NEUTRAL}),
                (
                    "single2double",
                    "[#6H0:1][#7D2H1,#7D1H2:2]",
                    {"formula_hint": NEUTRAL},
                ),
                (
                    ("single2double", "addPlus1"),
                    "[#6H0:1][#7D3:2]",
                    {"formula_hint": NEUTRAL},
                ),
                # Optional methide: same single2double edit; gated by pathway options.
                (
                    "single2double",
                    "[#6R:1][#6D1H3,#6D2H2,#6D3H1:2]",
                    {
                        "pathway": "methide",
                        "one_sided": True,
                        "recommend": "[#6;R]=[#6D1H2,#6D2H1,#6D3H0]",
                        "formula_hint": NEUTRAL,
                    },
                ),
            ],
            phase1_sites_on="atom_hydrogen",
            sites_on="atom_pairs",
        )
        self.pathways = frozenset(pathways or ())

    @staticmethod
    def _mod_names(modifications):
        if isinstance(modifications, str):
            return frozenset([modifications])
        return frozenset(modifications)

    def _active_pathways(self, pathways=None):
        if pathways is None:
            return self.pathways
        return frozenset(pathways)

    def pathways_toward(self, mol, target) -> frozenset:
        """Instance pathways plus those recommended by query options + formula."""
        return frozenset(self.pathways) | self.pathways_recommended_by_queries(
            mol, target
        )

    def modify(self, emol, mapid2atomidx, modifications):
        names = self._mod_names(modifications)
        if "single2double" in names:
            self.change_bond(emol, mapid2atomidx[1], mapid2atomidx[2], new_bond_type=2)

        if "addPlus1" in names:
            self.set_charge(emol, mapid2atomidx[2], 1)

        return True

    def site_from_matches(self, matches):
        site = []
        for item in matches:
            idx2mapids = item[0]
            site.append(idx2mapids[2])
        return frozenset(site)

    def could_help(self, mol, target, ctx) -> bool:
        """DH when heavy-atom formula already matches (bond / aromatic change)."""
        from .path_context import heavy_formula_equal

        return heavy_formula_equal(mol, target)

    def enumerate_for_path(self, mol, ctx, **kwargs):
        """Opt in query ``pathway`` options when formula + recommend SMARTS agree."""
        pathways = kwargs.pop("pathways", None)
        if pathways is None and ctx is not None:
            pathways = self.pathways_toward(mol, ctx.target)
        kwargs["pathways"] = pathways
        yield from super(Dehydrogenation, self).enumerate_for_path(
            mol, ctx, **kwargs
        )

    def metabolites(self, mol, pathways=None, **kwargs):
        active = self._active_pathways(pathways)

        for site, metabolites in super(Dehydrogenation, self).metabolites(
            mol, **kwargs
        ):
            yield site, metabolites

        # I tried moving this to within each downstream function, but this increased the run time
        # of the test suite so it seems most efficient to precompute the standardized molecule.
        template_mol = self.standardize(mol)

        # If a template_mol could not be successfully constructed due to problems in the input,
        # exit.
        if not template_mol:
            return

        from .base import _include_atom_sets

        want = _include_atom_sets(kwargs.get("include_sites"))

        matches = self.match_queries(template_mol)

        for res_struct, pair, path in self.resonate_with_pair_paths(
            template_mol, valid_atoms=set(matches)
        ):

            if not set(pair).issubset(set(matches)):
                continue

            for matches1, matches2 in itertools.product(*[matches[x] for x in pair]):
                if not self.pair_query_options_allowed(
                    matches1, matches2, active_pathways=active
                ):
                    continue

                site_atoms = self.site_from_matches([matches1, matches2])
                if want is not None and frozenset(site_atoms) not in want:
                    continue

                product = self.apply_modifications(res_struct, [matches1, matches2])

                if not product:
                    continue

                self.swap_bonds_along_path(product, path)
                yield (self.name, site_atoms), clean(product)


class Hydrogenation(ResonancePairRule):
    """Swaps single/double bonds for paths within a molecule's resonance structures.

    >>> mol = Chem.MolFromSmiles('CC(=O)N=C1C=CC(=O)C=C1')
    >>> site, metabolites = next(Hydrogenation().metabolites_from_sites(mol, frozenset({8, 3})))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(=O)NC1=CC=C(O)C=C1')
    True

    """

    def __init__(self):
        super(Hydrogenation, self).__init__(
            rxns=[
                ("[#6:1]#[#6:2]>>[*:1]=[*:2]", {"formula_hint": NEUTRAL}),
                ("[#6:1]=[#6:2]>>[*:1]-[*:2]", {"formula_hint": NEUTRAL}),
            ],
            phase1_sites_on="atoms",
            query_smarts=[("", "[*:1]", {"formula_hint": NEUTRAL})],
        )

    def metabolites(self, mol, **kwargs):

        for site, metabolites in super(Hydrogenation, self).metabolites(mol, **kwargs):
            yield site, metabolites

        # Resonance path swaps are formula-neutral; skip when target formula differs.
        toward = kwargs.get("toward_target")
        if toward is not None and not NEUTRAL.compatible(mol, toward):
            return

        from .base import _include_atom_sets

        want = _include_atom_sets(kwargs.get("include_sites"))

        for res_struct, pair, path in self.resonate_with_pair_paths(mol):

            # For all bonds in bond_path, replace single bonds with double
            # bonds and double bonds with single bonds.
            if len(path) < 2:
                continue

            pair_fs = frozenset(pair)
            if want is not None and pair_fs not in want:
                continue

            self.swap_bonds_along_path(res_struct, path)

            yield (self.name, pair_fs), clean(res_struct)


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
                    metabolite = copy_mol(res_struct)

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
            rxns=[
                (
                    "[#6:1]=[#6,#7:2]>>[*:1]1-[*:2][O]1",
                    {"formula_hint": ADD_O},
                )
            ],
            sites_on="bonds",
        )

    def is_redundant(self, peer_rules) -> bool:
        """Drop duplicate Epoxidation instances (keep the first)."""
        for peer in peer_rules:
            if type(peer) is type(self):
                return peer is not self
        return False

    def could_help(self, mol, target, ctx) -> bool:
        """Skip epoxidation when target is quinoid / non-aromatic with ≥2 O."""
        from .path_context import _aromatic_bond_count, _formula, oxygen_deficit

        if oxygen_deficit(mol, target) <= 0:
            return False
        # Quinone-like targets: prefer Hydroxylation → Dehydrogenation.
        if _aromatic_bond_count(target) == 0 and _formula(target).get("O", 0) >= 2:
            return False
        return True


class Acetylation(ConjugationRule):
    """Adds an acetyl to OH, NH, and SH.

    Defaults to a bare ``*`` adduct. Use ``as_star=False`` for the full acetyl.

    >>> mol = Chem.MolFromSmiles('C1=CC(=C(C=C1N)C(=O)O)O')
    >>> site, metabolites = next(Acetylation(as_star=False).metabolize(mol))
    >>> canon_smi(metabolites[0]) == canon_smi('CC(=O)NC1=CC(C(=O)O)=C(O)C=C1')
    True
    """

    smarts = "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]"


class Glutathionation(ConjugationRule):
    """Adds glutathione to common soft-electrophile motifs.

    Defaults to a bare ``*`` adduct. Use ``as_star=False`` for the full GSH
    peptide. ``star_label`` may be ``GSH``, ``Protein``, ``DNA``, or ``Cyanide``
    (CXSMILES). ``Protein`` / ``DNA`` / ``Cyanide`` are star-only. Set
    ``include_thiol=False`` to drop the substrate-thiol disulfide SMARTS (DNA /
    cyanide reactivity).

    SMARTS are high-sensitivity / low-specificity structure enumerators: if a
    site were positive, these are the adducts that would most likely form.
    Coverage includes epoxides, aziridines, C–halogen (F/Cl/Br/I), thiols,
    terminal alkenes, Michael β-carbons (enones / quinones / quinone-imines),
    aldehyde thiohemiacetals, sulfonate esters, and isocyanates /
    isothiocyanates.

    >>> G = Glutathionation(as_star=False)
    >>> mol = Chem.MolFromSmiles('c1ccccc1-C1OC1')
    >>> site, metabolites = next(G.metabolize(mol))
    >>> site
    ('Glutathionation', frozenset({6}))
    >>> canon_smi(metabolites, isomericSmiles=False) == canon_smi(['NC(CCC(=O)NC(CSC(CO)C1=CC=CC=C1)C(=O)NCC(=O)O)C(=O)O'], isomericSmiles=False)
    True

    """

    _GSH = "C(CC(=O)N[C@@H](CS({attach}))C(=O)NCC(=O)O)[C@@H](C(=O)O)N"
    _ALL_SMARTS = [
        "[#6:1]1[#8:2][#6:3]1>>" + _GSH.format(attach="[*:1][*:3][*:2]"),
        "[#6:1][F,Cl,Br,I:2]>>" + _GSH.format(attach="[*:1]"),
        "[#16h1:1]>>" + _GSH.format(attach="[*:1]"),
        "[#6H2:1]=[#6:2]>>" + _GSH.format(attach="[*:1]-[*:2]"),
        # Michael: β-CH (not CH2) of enone / quinone / quinone-imine → enol adduct.
        "[#6H1:1]=[#6:2][#6:3]=[#8,#7:4]>>"
        + _GSH.format(attach="[*:1][*:2]=[*:3][*:4]"),
        # Aldehyde → thiohemiacetal (formaldehyde is H2; others H1).
        "[#6;H1,H2:1]=[#8:2]>>" + _GSH.format(attach="[*:1]([*:2])"),
        # Aziridine ring opening (epoxide analog).
        "[#6:1]1[#7:2][#6:3]1>>" + _GSH.format(attach="[*:1][*:3][*:2]"),
        # Sulfonate / sulfate ester alkylation (C–OSO2–).
        "[#6:1][#8:2]S(=O)(=O)>>" + _GSH.format(attach="[*:1]"),
        # Isocyanate / isothiocyanate → thiocarbamate / dithiocarbamate.
        "[#7:1]=[#6:2]=[#8,#16:3]>>" + _GSH.format(attach="[*:2](=[*:3])[*:1]"),
    ]
    mapid_site = [1]

    def __init__(self, include_thiol=True, as_star=None, star_label=None, **kwargs):
        smarts = list(self._ALL_SMARTS)
        if not include_thiol:
            smarts = [s for s in smarts if "[#16h1" not in s]
        self.include_thiol = include_thiol
        super(Glutathionation, self).__init__(
            as_star=as_star,
            star_label=star_label,
            rxns=smarts,
            **kwargs
        )


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
        ("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)", {"formula_hint": CLEAVE}),
        ("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)", {"formula_hint": CLEAVE}),
        ("[#6H3:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)", {"formula_hint": CLEAVE}),
        ("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1](=O)O)", {"formula_hint": CLEAVE}),
        ("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)", {"formula_hint": CLEAVE}),
        ("[#6H2:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)", {"formula_hint": CLEAVE}),
        ("[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]=O)", {"formula_hint": CLEAVE}),
        ("[#6H1:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)", {"formula_hint": CLEAVE}),
        ("[#6H0:1][#7,#8H0,#16:2]>>([*:2].[*:1]-O)", {"formula_hint": CLEAVE}),
        ("[#6:1][#6:2]>>(O-[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#6h:1][#6:2]>>(O-[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#6h:1][#6:2]>>(O=[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#8H1:3]-[#6:1]-[#7,#8,#16:2]>>([*:3]=[*:1].[*:2])", {"formula_hint": CLEAVE}),
    ]
    mapid_site = [1, 2]
    sites_on = "bonds"

    def elements_may_add(self):
        # SMARTS emit carbonyl / alcohol O on the carbon fragment (CLEAVE hint alone
        # does not list ADD_O — still introduce oxygen for impossibility checks).
        return frozenset({"O"})


class NDealkylation(Dealkylation):
    """N-dealkylation: Dealkylation products whose formation site includes nitrogen.

    Reuses Dealkylation SMARTS; drops O-/S-/C-only sites. Short-circuits when the
    reactant has no nitrogen atoms.
    """

    def metabolites(self, mol, kekulize=True, **kwargs):
        if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
            return
        for site_key, products in super(NDealkylation, self).metabolites(
            mol, kekulize=kekulize, **kwargs
        ):
            _rule, site = site_key
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in site):
                yield site_key, products

    def sites_toward(self, mol, target, ctx):
        sites = super(NDealkylation, self).sites_toward(mol, target, ctx)
        if sites is None:
            return None
        return [
            s
            for s in sites
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in s)
        ]

    def is_redundant(self, peer_rules) -> bool:
        """Redundant when Phase I Dealkylation is already among peers."""
        for peer in peer_rules:
            if type(peer) is type(self):
                return peer is not self
            if getattr(peer, "name", type(peer).__name__) == "Dealkylation":
                return True
        return False


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
        ("[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2](O).[*:3])", {"formula_hint": CLEAVE}),
        ("[#8,#16:1]=[#6:2]-[#7,#8,#16:3]>>([*:1]=[*:2].[*:3])", {"formula_hint": CLEAVE}),
    ]
    mapid_site = [2, 3]
    sites_on = "bonds"

    def elements_may_add(self):
        return frozenset({"O"})


class ReductiveDehalogenation(SmartsReactionRule):
    """Cleaves halogen-carbon bonds.

    >>> RD = ReductiveDehalogenation()
    >>> mol = Chem.MolFromSmiles('ClCC(=O)C(NC(=O)c1cc(Cl)c(c(c1)Cl)CO)(CC)C')
    >>> site, metabolites = next(RD.metabolites_from_sites(mol,[0,1]))
    >>> canon_smi(metabolites) == canon_smi(['Cl', 'CCC(C)(NC(=O)C1=CC(Cl)=C(CO)C(Cl)=C1)C(C)=O'])
    True

    """

    smarts = [
        ("[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
        ("[#9,#17,#35,#53,#85:1]-[#6:2]-[#6:3]>>[*:1].[*:2]=[*:3]", {"formula_hint": CLEAVE}),
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
        ("[#9,#17,#35,#53,#85:1]-[#6:2]>>[*:1].[*:2]O", {"formula_hint": ADD_O}),
        ("[#9,#17,#35,#53,#85:1]-[#6h1:2]>>[*:1].[*:2]=O", {"formula_hint": ADD_O}),
        ("[#9,#17,#35,#53,#85:1]-[#6H2:2]>>[*:1].[*:2](O)=O", {"formula_hint": ADD_O}),
        ("[#9,#17,#35,#53,#85:1]-[#6:2][#6H1:3]>>[*:2](O)[*:3]-[*:1]", {"formula_hint": ADD_O}),
        ("[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)=O.[*:3]", {"formula_hint": ADD_O}),
        ("[#9,#17,#35,#53,#85:1]-[#6:2]-[#9,#17,#35,#53,#85:3]>>[*:1].[*:2](O)O.[*:3]", {"formula_hint": ADD_O}),
    ]
    sites_on = "bonds"

    def cleave_alone(self) -> bool:
        return True


class NitrogenOxidation(SmartsReactionRule):
    """P450-like N-oxidation of trivalent nitrogen.

    - N-H (amines, pyrrole-like N) -> hydroxylamine / N-OH.
    - Primary amine NH2 -> nitroso.
    - Tertiary amine or pyridine-like N (no H) -> N-oxide.

    Adding OH or =O to an already fully substituted N (pyridine N-OH,
    pyrrole N=O, tertiary amine N-OH) would be pentavalent and is not emitted.

    >>> mol = Chem.MolFromSmiles('CCN')
    >>> site, metabolites = next(NitrogenOxidation().metabolize(mol))
    >>> canon_smi(metabolites[0]) == canon_smi('CCNO')
    True

    """

    smarts = [
        ("[#7v3h:1]>>[*:1]O", {"formula_hint": ADD_O}),
        ("[#7v3H2:1]>>[*:1]=O", {"formula_hint": ADD_O}),
        ("[#7v3H0:1]>>[*&H0&+:1][O-]", {"formula_hint": ADD_O}),
    ]
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
        ("[#16;v2,v4:1]>>[*&H0&+:1][O-]", {"formula_hint": ADD_O}),
        ("[#16;v2,v4:1]>>[*:1][O]", {"formula_hint": ADD_O}),
        ("[#16;v2,v4:1]>>[*:1]=O", {"formula_hint": ADD_O}),
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

    smarts = [
        ("[#6h:1]>>[*:1]O", {"formula_hint": ADD_O}),
        ("[#6h2:1]>>[*:1]=O", {"formula_hint": ADD_O}),
    ]
    phase1_sites_on = "atom_hydrogen"

    def could_help(self, mol, target, ctx) -> bool:
        """Hydroxylation only when the target still needs more oxygen."""
        from .path_context import oxygen_deficit

        return oxygen_deficit(mol, target) > 0

    def sites_toward(self, mol, target, ctx):
        """Conserved MCS atoms, one representative per topological orbit."""
        if ctx is None or not getattr(ctx, "conserved_r_atoms", None):
            return None
        conserved = set(ctx.conserved_r_atoms)
        r_only = set(getattr(ctx, "r_only_atoms", ()) or ())
        if not r_only or len(conserved) == 0:
            return None
        ranks = self.topol_equiv(mol)
        seen_ranks = set()
        out = []
        for a in sorted(conserved):
            rank = ranks.get(a)
            if rank in seen_ranks:
                continue
            seen_ranks.add(rank)
            out.append(frozenset([a]))
        return out or None

    def child_may_reach(self, parent, child, target, ctx) -> bool:
        if ctx is not None:
            try:
                if ctx.distance_proxy(child) > ctx.distance_proxy(parent):
                    return False
            except Exception:
                pass
        return True


class OxygenReduction(SmartsReactionRule):
    """Converts C=O/N=O bonds to CO/NO.

    >>> mol = Chem.MolFromSmiles('Oc1ccccc1')
    >>> site, metabolites = next(Dehydration().metabolites_from_sites(mol,[frozenset({1,0})]))
    >>> canon_smi(metabolites) == canon_smi(['C1=CC=CC=C1', 'O'])
    True

    """

    smarts = [
        ("[#8:1]=[#6,#7:2]>>[*:1]-[*:2]", {"formula_hint": NEUTRAL}),
        ("[#8:1]-[#8:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
    ]


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
        ("[#6,#7:1]-[#8H1:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
        ("[#6:3]-[#6:1]-[#8H1:2]>>[*:3]=[*:1].[*:2]", {"formula_hint": CLEAVE}),
        ("[#6,#7:1]=[#8:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
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

    smarts = [
        (
            "[#8:1][#15:2](=[#8:3])([#8:4])[#8:5]>>[*:1].[*:2](=[*:3])([*:4])[*:5]",
            {"formula_hint": CLEAVE},
        )
    ]
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
        ("[#16:1]=[#8:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
        ("[#16:1]-[#16:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
        ("[#16:1]-[#6,#8:2]>>[*:1].[*:2]", {"formula_hint": CLEAVE}),
    ]


class Glucuronidation(ConjugationRule):
    """Glucuronidation of acids and phenols.

    Defaults to a bare ``*`` adduct. Use ``as_star=False`` for the full GlcA, or
    ``star_label="GlcA"`` for a CXSMILES-labeled star.

    >>> mol = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
    >>> site, metabolites = next(Glucuronidation(as_star=False).metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['CC(=O)NC1=CC=C(OC2OC(C(=O)O)C(O)C(O)C2O)C=C1'])
    True

    """

    # Attachment oxygen is always map :1 so SOM is a single atom.
    smarts = [
        "[#8:1][#6:2](=[O,N,P,S:3])[#6:4]>>\
                        O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2](=[#8:3])[*:4])1",
        "[#8H1:1][#6:2]>>O1C(C(=O)O)C(O)C(O)C(O)C([*:1][*:2])1",
    ]
    mapid_site = [1]


class Sulfation(ConjugationRule):
    """
    >>> mol = Chem.MolFromSmiles('CC(=O)Nc1ccc(O)cc1')
    >>> site, metabolites = next(Sulfation(as_star=False).metabolize(mol))
    >>> canon_smi(metabolites) == canon_smi(['CC(=O)NC1=CC=C(OS(=O)(=O)O)C=C1'])
    True

    >>> site
    ('Sulfation', frozenset({8, 7}))

    """

    def __init__(self, as_star=None, star_label=None, **kwargs):
        super(Sulfation, self).__init__(
            as_star=as_star,
            star_label=star_label,
            rxns=[
                "[#6:1][#8:2]>>[*:1][*:2]S(=O)(=O)O",
                "[#6:1]1=[#6:2][#6:3]2[#8:7][#6:4]2[#6:5]=[#6:6]1>>\
                    [*:1]1=[*:2][*:3]=[*:4](-S(C)(=O)(=O))[*:5]=[*:6]1",
            ],
            mapid_site=[1, 2],
            **kwargs
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

    def cleave_alone(self) -> bool:
        return True


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
        ("[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1])", {"formula_hint": NEUTRAL}),
        ("[#6:1]1[#8:2][#6:3]1>>([*:2][*:3][*:1]O)", {"formula_hint": ADD_O}),
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
        ("[#8:3]=[#7+1:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#8:3]=[#7:1]-[#8-1:2]>>([*:3]=[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#8:3]=[#7:1]-[#8:2]>>([*:3]=[*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#7:1](=[#8:2])-[#8:3]>>([*:1].[*:2].[*:3])", {"formula_hint": CLEAVE}),
        ("[#8:3]=[#7:1]-[#8:2]>>([*:1].[*:2].[*:3])", {"formula_hint": CLEAVE}),
        ("[#7:1]-[#8:2]>>([*:1].[*:2])", {"formula_hint": CLEAVE}),
        ("[#7D2:1]=[#8:2]>>([*:1].[*2])", {"formula_hint": CLEAVE}),
        ("[#7:1](~[#8:2])~[#8:3]>>([*:1].[*:2].[*:3])", {"formula_hint": CLEAVE}),
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


# Phase I rules (and N-dealkylation) expose degenerate phase1_steps singletons.
# Formula heuristics live on SMARTS entry options (``formula_hint``);
# Hydroxylation / Epoxidation / Dehydrogenation keep custom could_help overlays.
for _phase1_rule in (
    Dehydrogenation,
    Dephosphorylation,
    EpoxideOpening,
    Hydrolysis,
    Dehydration,
    Hydrogenation,
    NitrogenReduction,
    OxygenReduction,
    ReductiveDehalogenation,
    SulfurReduction,
    Hydroxylation,
    Epoxidation,
    SulfurOxidation,
    NitrogenOxidation,
    Dealkylation,
    OxidativeDehalogenation,
    NDealkylation,
):
    _phase1_rule.phase1_equivalent = True


if __name__ == "__main__":
    import doctest

    doctest.testmod()
