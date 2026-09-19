"""The four easy pairs, both libraries, canonical fragment SMILES.

Order does not matter. A product only the old library emits is a failure
unless `src/xenosite/refactor_poc/DIVERGENCES.md` records both sets and the
reason. These four are not xfail.
"""

from rdkit import Chem

from xenosite.forest.rules import AzoSplitting as OldAzoSplitting
from xenosite.forest.rules import BenzodioxoleReduction as OldBenzodioxoleReduction
from xenosite.forest.rules import ConjugationRule as OldConjugationRule
from xenosite.forest.rules import Dealkylation as OldDealkylation
from xenosite.forest.rules import Hydroxylation as OldHydroxylation
from xenosite.forest.rules import NDealkylation as OldNDealkylation
from xenosite.forest.rules import NitroaromaticReduction as OldNitroaromaticReduction
from xenosite.forest.rules import QuinoneFormation as OldQuinone
from xenosite.forest.rules import ThiopheneSulfurOxidation as OldThiopheneSulfurOxidation
from xenosite.refactor_poc.find_path import bfs
from xenosite.refactor_poc.rules import (
    AzoSplitting,
    BenzodioxoleReduction,
    ConjugationRule,
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    NDealkylation,
    NitroaromaticReduction,
    QuinoneFormation,
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


def test_cleaved_ring_bond_sets_breaks_ring():
    flags = {
        (info["options"].get("leave_count"), info["options"].get("breaks_ring"))
        for _product, info in NDealkylation().metabolize(
            Chem.MolFromSmiles("CN1CCCCC1")
        )
    }
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
