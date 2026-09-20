"""The four easy pairs, both libraries, canonical fragment SMILES.

Order does not matter. A product only the old library emits is a failure
unless `src/xenosite/refactor_poc/DIVERGENCES.md` records both sets and the
reason. These four are not xfail.
"""

from rdkit import Chem

from xenosite.forest.rules import Acetylation as OldAcetylation
from xenosite.forest.rules import AzoSplitting as OldAzoSplitting
from xenosite.forest.rules import BenzodioxoleReduction as OldBenzodioxoleReduction
from xenosite.forest.rules import ConjugationRule as OldConjugationRule
from xenosite.forest.rules import Dealkylation as OldDealkylation
from xenosite.forest.rules import Glucuronidation as OldGlucuronidation
from xenosite.forest.rules import Glutathionation as OldGlutathionation
from xenosite.forest.rules import Hydroxylation as OldHydroxylation
from xenosite.forest.rules import NDealkylation as OldNDealkylation
from xenosite.forest.rules import NitroaromaticReduction as OldNitroaromaticReduction
from xenosite.forest.rules import QuinoneFormation as OldQuinone
from xenosite.forest.rules import Sulfation as OldSulfation
from xenosite.forest.rules import (
    ThiopheneSulfurOxidation as OldThiopheneSulfurOxidation,
)
from xenosite.refactor_poc.find_path import bfs
from xenosite.refactor_poc.rules import (
    Acetylation,
    AzoSplitting,
    BenzodioxoleReduction,
    ConjugationRule,
    Dealkylation,
    Dehydrogenation,
    Glucuronidation,
    Glutathionation,
    Hydroxylation,
    NDealkylation,
    NitroaromaticReduction,
    QuinoneFormation,
    Sulfation,
    ThiopheneSulfurOxidation,
)
from xenosite.refactor_poc.rulesets import RuleSet

# Old `[#6h2:1]>>[*:1]=O`. Both patterns add OH. See DIVERGENCES.md.
_OLD_BUTYL_KETONES = {
    "CC(=O)CCc1ccccc1",
    "CCC(=O)Cc1ccccc1",
    "CCCC(=O)c1ccccc1",
}

_QUINONE = "O=C1C=CC(=O)C=C1"
_CHAIN = "OCCCCc1ccccc1"


def _canon(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _fragments(mols):
    found = set()
    for mol in mols:
        if mol is None:
            continue
        if isinstance(mol, str):
            parsed = Chem.MolFromSmiles(mol)
            assert parsed is not None, mol
            mol = parsed
        copied = Chem.Mol(mol)
        for atom in copied.GetAtoms():
            atom.SetAtomMapNum(0)
        text = Chem.MolToSmiles(copied)
        parsed = Chem.MolFromSmiles(text)
        assert parsed is not None, text
        pieces = Chem.GetMolFrags(parsed, asMols=True, sanitizeFrags=True)
        for piece in pieces:
            found.add(Chem.MolToSmiles(piece, canonical=True, isomericSmiles=False))
    return found


def _old(rule, smiles):
    pieces = []
    for _site, metabolites in rule.metabolize(Chem.MolFromSmiles(smiles)):
        pieces.extend(metabolites)
    return _fragments(pieces)


def _new(rule, smiles, **kwargs):
    return _fragments(
        product
        for product, _info in rule.metabolize(Chem.MolFromSmiles(smiles), **kwargs)
    )


def test_ethane_hydroxylation_matches_old():
    old = _old(OldHydroxylation(), "CC")
    new = _new(Hydroxylation(), "CC")
    assert new == old == {"CCO"}


def test_butylbenzene_chain_alcohol_matches_old():
    old = _old(OldHydroxylation(), "c1ccc(CCCC)cc1")
    new = _new(Hydroxylation(), "c1ccc(CCCC)cc1")
    assert _canon(_CHAIN) in old
    assert _canon(_CHAIN) in new
    assert old - new == _OLD_BUTYL_KETONES
    assert new - old == set()


def test_anisole_dealkylation_matches_old():
    old = _old(OldDealkylation(), "COc1ccccc1")
    new = _new(Dealkylation(), "COc1ccccc1")
    assert new == old
    assert "Oc1ccccc1" in new


def test_trimethylamine_ndealkylation_matches_old():
    old = _old(OldNDealkylation(), "CN(C)C")
    new = _new(NDealkylation(), "CN(C)C")
    assert "CNC" in new
    assert new == old


def test_filter_skips_named_methyl_and_keeps_open_alkyl():
    """The skip reads ``leave_count``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        count = info["options"].get("leave_count")
        if count is not None:
            refused.append((site, count, info["options"].get("partner")))
            return False
        return True

    products = _new(
        NDealkylation(),
        "CCN(C)C",
        filter_sites=filter_sites,
    )
    assert refused
    assert all(count == 1 and partner == "N" for _site, count, partner in refused)
    assert "CC=O" in products
    assert "C=O" not in products


_OLSALAZINE = "OC(=O)c1cc(/N=N/c2ccc(c(c2)C(=O)O)O)ccc1O"
_AZO_PYRIDAZINE = "c1ccc(N=Nc2ccnnc2)cc1"


def test_olsalazine_azo_splitting_matches_old():
    old = _old(OldAzoSplitting(), _OLSALAZINE)
    new = _new(AzoSplitting(), _OLSALAZINE)
    assert "Nc1ccc(O)c(C(=O)O)c1" in new
    assert new == old


def test_filter_skips_ring_nitrogen_and_keeps_open_azo():
    """The skip reads ``breaks_ring``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        if info["options"].get("breaks_ring"):
            refused.append((site, info["options"].get("partner")))
            return False
        return True

    products = _new(
        AzoSplitting(),
        _AZO_PYRIDAZINE,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(partner == "N" for _site, partner in refused)
    assert "Nc1ccccc1" in products
    assert "C=CC(=CN)N=Nc1ccccc1" not in products
    assert "C=C(C=CN)N=Nc1ccccc1" not in products


_BENZODIOXOLE = "c1ccc2c(c1)OCO2"


def test_benzodioxole_reduction_matches_old():
    old = _old(OldBenzodioxoleReduction(), _BENZODIOXOLE)
    new = _new(BenzodioxoleReduction(), _BENZODIOXOLE)
    assert "Oc1ccccc1O" in new
    assert new == old


def test_filter_skips_named_methylene():
    """The skip reads ``leave_count``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        count = info["options"].get("leave_count")
        if count is not None:
            refused.append(
                (
                    site,
                    count,
                    info["options"].get("partner"),
                    info["options"].get("breaks_ring"),
                )
            )
            return False
        return True

    products = _new(
        BenzodioxoleReduction(),
        _BENZODIOXOLE,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        count == 1 and partner == "O" and ring
        for _site, count, partner, ring in refused
    )
    assert "Oc1ccccc1O" not in products
    assert "C" not in products


_NITROBENZENE = "O=[N+]([O-])c1ccccc1"


def test_nitrobenzene_reduction_matches_old():
    old = _old(OldNitroaromaticReduction(), _NITROBENZENE)
    new = _new(NitroaromaticReduction(), _NITROBENZENE)
    assert "O=Nc1ccccc1" in new
    assert new == old


def test_filter_skips_named_nitro_oxygen():
    """The skip reads ``leave_count``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        count = info["options"].get("leave_count")
        if count is not None:
            refused.append(
                (
                    site,
                    count,
                    info["options"].get("partner"),
                    info["options"].get("breaks_ring"),
                )
            )
            return False
        return True

    products = _new(
        NitroaromaticReduction(),
        _NITROBENZENE,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        count == 1 and partner == "N" and not ring
        for _site, count, partner, ring in refused
    )
    assert "O=Nc1ccccc1" not in products
    assert "O" not in products


_THIOPHENE = "c1ccsc1"


def test_thiophene_sulfur_oxidation_matches_old():
    old = _old(OldThiopheneSulfurOxidation(), _THIOPHENE)
    new = _new(ThiopheneSulfurOxidation(), _THIOPHENE)
    assert "[O-][s+]1cccc1" in new
    assert new == old


def test_filter_skips_sulfur_oxygen_addition():
    """The skip reads ``adds``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        added = info["options"].get("adds")
        if added:
            refused.append(
                (
                    site,
                    added,
                    info["options"].get("symbol"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        ThiopheneSulfurOxidation(),
        _THIOPHENE,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        added == "O" and symbol == "S" and count is None and not cleaved
        for _site, added, symbol, count, cleaved in refused
    )
    assert "[O-][s+]1cccc1" not in products


_CONJUGATE = "[#7,#8,#16;h:1]>>[*:1][#6](=[#8])[#6]"
_PHENOL = "Oc1ccccc1"
_AMINOPHENOL = "Nc1ccc(O)cc1"


def test_phenol_star_conjugate_matches_old():
    """The target contains the star. Canonical SMILES; order does not matter."""

    old = _old(OldConjugationRule(rxns=_CONJUGATE), _PHENOL)
    new = _new(ConjugationRule(), _PHENOL)
    assert "*Oc1ccccc1" in new
    assert all("*" in smiles for smiles in new)
    assert new == old


def test_filter_skips_nitrogen_and_keeps_the_phenol_star():
    """The skip reads ``symbol``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        symbol = info["options"].get("symbol")
        if symbol == "N":
            refused.append(
                (
                    site,
                    symbol,
                    info["options"].get("adds"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        ConjugationRule(),
        _AMINOPHENOL,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        symbol == "N" and added == "CCO" and count is None and not cleaved
        for _site, symbol, added, count, cleaved in refused
    )
    assert "*Oc1ccc(N)cc1" in products
    assert "*Nc1ccc(O)cc1" not in products


def test_phenol_acetyl_matches_old():
    """The target contains the acetyl conjugate. Canonical SMILES; order does not matter."""

    old = _old(OldAcetylation(as_star=False), _PHENOL)
    new = _new(Acetylation(as_star=False), _PHENOL)
    assert "CC(=O)Oc1ccccc1" in new
    assert all("C(=O)" in smiles for smiles in new)
    assert new == old


def test_filter_skips_nitrogen_and_keeps_the_phenol_acetyl():
    """The skip reads ``symbol``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        symbol = info["options"].get("symbol")
        if symbol == "N":
            refused.append(
                (
                    site,
                    symbol,
                    info["options"].get("adds"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        Acetylation(as_star=False),
        _AMINOPHENOL,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        symbol == "N" and added == "CCO" and count is None and not cleaved
        for _site, symbol, added, count, cleaved in refused
    )
    assert "CC(=O)Oc1ccc(N)cc1" in products
    assert "CC(=O)Nc1ccc(O)cc1" not in products


def _acetates_by_site(smiles):
    mol = Chem.MolFromSmiles(smiles)
    found = {}
    for product, info in Acetylation(as_star=False).metabolize(mol):
        site = info["site"]
        assert isinstance(site, frozenset) and len(site) == 1
        found[next(iter(site))] = next(iter(_fragments([product])))
    return mol, found


def test_pinned_acetate_is_the_stamped_atom():
    """The primary alcohol and the nitrogen each get their own acetate.

    The site names the atom. The first ``RunReactants`` product does not.
    """

    mol, found = _acetates_by_site("OCCN")
    nitrogen = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    alcohol = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    assert found[nitrogen] == _canon("CC(=O)NCCO")
    assert found[alcohol] == _canon("CC(=O)OCCN")

    mol, found = _acetates_by_site("CC(O)CO")
    primary = None
    secondary = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        carbon = next(iter(atom.GetNeighbors()))
        if carbon.GetTotalNumHs() == 2:
            primary = atom.GetIdx()
        elif carbon.GetTotalNumHs() == 1:
            secondary = atom.GetIdx()
    assert primary is not None and secondary is not None
    assert found[primary] == _canon("CC(=O)OCC(C)O")
    assert found[secondary] == _canon("CC(=O)OC(C)CO")


_HYDROXYBENZYL = "OCc1ccc(O)cc1"


def test_phenol_sulfate_matches_old():
    """The target contains the sulfate conjugate. Canonical SMILES; order does not matter."""

    old = _old(OldSulfation(as_star=False), _PHENOL)
    new = _new(Sulfation(as_star=False), _PHENOL)
    assert "O=S(=O)(O)Oc1ccccc1" in new
    assert all("S(=O)(O)" in smiles for smiles in new)
    assert new == old


def test_filter_skips_the_benzylic_alcohol_and_keeps_the_phenol_sulfate():
    """The skip reads ``partner_h``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        partner_h = info["options"].get("partner_h")
        if partner_h == 2:
            refused.append(
                (
                    site,
                    partner_h,
                    info["options"].get("symbol"),
                    info["options"].get("adds"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        Sulfation(as_star=False),
        _HYDROXYBENZYL,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        count == 2
        and symbol == "O"
        and added == "SOOO"
        and leave is None
        and not cleaved
        for _site, count, symbol, added, leave, cleaved in refused
    )
    assert "O=S(=O)(O)Oc1ccc(CO)cc1" in products
    assert "O=S(=O)(O)OCc1ccc(O)cc1" not in products


# Old `[#8:1][#6:2](=[O,N,P,S:3])` writes `=[#8:3]`. See DIVERGENCES.md.
_OLD_IMIDATE_ACID = {"O=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1"}


def test_phenol_glucuronide_matches_old():
    """The target contains the glucuronide. Canonical SMILES; order does not matter."""

    old = _old(OldGlucuronidation(as_star=False), _PHENOL)
    new = _new(Glucuronidation(as_star=False), _PHENOL)
    assert "O=C(O)C1OC(Oc2ccccc2)C(O)C(O)C1O" in new
    assert all("O=C(O)C1OC" in smiles for smiles in new)
    assert new == old


def test_filter_skips_the_benzylic_alcohol_and_keeps_the_phenol_glucuronide():
    """The skip reads ``partner_h``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        partner_h = info["options"].get("partner_h")
        if partner_h == 2:
            refused.append(
                (
                    site,
                    partner_h,
                    info["options"].get("symbol"),
                    info["options"].get("adds"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        Glucuronidation(as_star=False),
        _HYDROXYBENZYL,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        count == 2
        and symbol == "O"
        and added == "CCCCCCOOOOOO"
        and leave is None
        and not cleaved
        for _site, count, symbol, added, leave, cleaved in refused
    )
    assert "O=C(O)C1OC(Oc2ccc(CO)cc2)C(O)C(O)C1O" in products
    assert "O=C(O)C1OC(OCc2ccc(O)cc2)C(O)C(O)C1O" not in products


def test_imidate_glucuronide_keeps_the_nitrogen():
    """The old carbonyl OR rewrites nitrogen as oxygen. That string is not this glucuronide."""

    old = _old(OldGlucuronidation(as_star=False), "N=C(O)c1ccccc1")
    new = _new(Glucuronidation(as_star=False), "N=C(O)c1ccccc1")
    assert "N=C(OC1OC(C(=O)O)C(O)C(O)C1O)c1ccccc1" in new
    assert old - new == _OLD_IMIDATE_ACID
    assert new - old == set()


_STYRENE_OXIDE = "c1ccccc1C1OC1"
_EPICHLOROHYDRIN = "ClCC1CO1"
_GSH_ADDS = "CCCCCCCCCCNNNOOOOOOS"


def test_styrene_oxide_glutathione_matches_old():
    """The target contains the glutathione conjugate. Canonical SMILES; order does not matter."""

    old = _old(OldGlutathionation(as_star=False), _STYRENE_OXIDE)
    new = _new(Glutathionation(as_star=False), _STYRENE_OXIDE)
    assert "NC(CCC(=O)NC(CSC(CO)c1ccccc1)C(=O)NCC(=O)O)C(=O)O" in new
    assert all("C(=O)NCC(=O)O" in smiles for smiles in new)
    assert new == old


def test_filter_skips_the_alkyl_chloride_and_keeps_the_epoxide_glutathione():
    """The skip reads ``partner``. It does not ask which rule this is."""

    refused = []

    def filter_sites(site, info):
        partner = info["options"].get("partner")
        if partner == "Cl":
            refused.append(
                (
                    site,
                    partner,
                    info["options"].get("symbol"),
                    info["options"].get("adds"),
                    info["options"].get("removes"),
                    info["options"].get("leave_count"),
                    info["options"].get("cleaves"),
                )
            )
            return False
        return True

    products = _new(
        Glutathionation(as_star=False),
        _EPICHLOROHYDRIN,
        filter_sites=filter_sites,
    )
    assert refused
    assert all(
        partner == "Cl"
        and symbol == "C"
        and added == _GSH_ADDS
        and removed == "Cl"
        and leave is None
        and not cleaved
        for _site, partner, symbol, added, removed, leave, cleaved in refused
    )
    assert "NC(CCC(=O)NC(CSC(CO)CCl)C(=O)NCC(=O)O)C(=O)O" in products
    assert "NC(CCC(=O)NC(CSCC1CO1)C(=O)NCC(=O)O)C(=O)O" not in products


def test_cleaved_ring_bond_sets_breaks_ring():
    flags = set()
    for _product, info in NDealkylation().metabolize(Chem.MolFromSmiles("CN1CCCCC1")):
        options = info["options"]
        assert isinstance(options, dict)
        flags.add((options.get("leave_count"), options.get("breaks_ring")))
    assert (1, False) in flags
    assert (None, True) in flags


def test_benzene_and_phenol_quinone_match_old():
    """Quinone formation matches. The same product is also two steps, not a quinone step."""

    steps = RuleSet((Hydroxylation, Dehydrogenation), name="Poc")
    for smiles, depth in (("c1ccccc1", 3), ("Oc1ccccc1", 2)):
        old = _old(OldQuinone(), smiles)
        new = _new(QuinoneFormation(), smiles)
        assert new == old, smiles
        assert _canon(_QUINONE) in new
        enumerated = set()
        for _mol, info in bfs(smiles, steps, depth=depth):
            enumerated |= _fragments([info["csmi"]])
        assert _canon(_QUINONE) in enumerated
