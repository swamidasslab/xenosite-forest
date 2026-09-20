"""Remaining default-set rules, both libraries, canonical fragment SMILES.

Order does not matter. A product only the old library emits is a failure
unless `src/xenosite/forest/DIVERGENCES.md` records both sets and the
reason. These comparisons call each rule's metabolites method. They do not
call find_path.
"""

from rdkit import Chem

from xenosite._archive_forest.rules import Dehydration as OldDehydration
from xenosite._archive_forest.rules import Dehydrogenation as OldDehydrogenation
from xenosite._archive_forest.rules import Dephosphorylation as OldDephosphorylation
from xenosite._archive_forest.rules import Epoxidation as OldEpoxidation
from xenosite._archive_forest.rules import EpoxideOpening as OldEpoxideOpening
from xenosite._archive_forest.rules import Hydrogenation as OldHydrogenation
from xenosite._archive_forest.rules import Hydrolysis as OldHydrolysis
from xenosite._archive_forest.rules import NitrogenOxidation as OldNitrogenOxidation
from xenosite._archive_forest.rules import NitrogenReduction as OldNitrogenReduction
from xenosite._archive_forest.rules import (
    OxidativeDehalogenation as OldOxidativeDehalogenation,
)
from xenosite._archive_forest.rules import OxygenReduction as OldOxygenReduction
from xenosite._archive_forest.rules import (
    ReductiveDehalogenation as OldReductiveDehalogenation,
)
from xenosite._archive_forest.rules import SulfurOxidation as OldSulfurOxidation
from xenosite._archive_forest.rules import SulfurReduction as OldSulfurReduction
from xenosite.forest.rules import (
    Dehydration,
    Dehydrogenation,
    Dephosphorylation,
    Epoxidation,
    EpoxideOpening,
    Hydrogenation,
    Hydrolysis,
    NitrogenOxidation,
    NitrogenReduction,
    OxidativeDehalogenation,
    OxygenReduction,
    ReductiveDehalogenation,
    SulfurOxidation,
    SulfurReduction,
)

# Even path. Still C4H6. See DIVERGENCES.md.
_OLD_CUMULENE = {"C=C=CC"}

# Old `[#7D2:1]=[#8:2]>>([*:1].[*2])`. `[*2]` is a dummy, not oxygen.
_OLD_NITROSO_DUMMY = {"*"}

# Old RunReactants also cleaves a P-OH bond. See DIVERGENCES.md.
_OLD_METHYL_PHOSPHITE = {"CO[PH](=O)O", "O"}


def _fragments(mols):
    found = set()
    for mol in mols:
        if mol is None:
            continue
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
    for _site, metabolites in rule.metabolites(Chem.MolFromSmiles(smiles)):
        pieces.extend(metabolites)
    return _fragments(pieces)


def _new(rule, smiles):
    pieces = []
    for outcome in rule.metabolites(Chem.MolFromSmiles(smiles)):
        pieces.extend(outcome.products)
    return _fragments(pieces)


def test_ethene_and_ethyne_hydrogenation_match_old():
    """The one-bond SMARTS. These reactants have no pair path to swap."""

    assert _new(Hydrogenation(), "C=C") == _old(OldHydrogenation(), "C=C") == {"CC"}
    assert _new(Hydrogenation(), "C#C") == _old(OldHydrogenation(), "C#C") == {"C=C"}


def test_butadiene_hydrogenation_includes_2_butene():
    """1,4-hydrogenation is 2-butene. Cumulenes are false emits, not hydrogenation.

    Hard-fail regression: single-first / even bond-count ResonancePair walks
    must not emit ``C=C=C=C`` (C4H4) or ``C=C=CC`` (C4H6). Those are a
    Hydrogenation chemistry bug (unmasked when unique-edit stopped hiding
    them — not a dedup bug). Real products are C4H8 only. See DIVERGENCES.
    """

    old = _old(OldHydrogenation(), "C=CC=C")
    new = _new(Hydrogenation(), "C=CC=C")
    assert new == {"C=CCC", "CC=CC"}
    assert "C=C=CC" not in new
    assert "C=C=C=C" not in new
    assert old - new == _OLD_CUMULENE
    assert new - old == set()


def test_propene_epoxidation_matches_old():
    old = _old(OldEpoxidation(), "CC=C")
    new = _new(Epoxidation(), "CC=C")
    assert new == old == {"CC1CO1"}


def test_ethanol_dehydrogenation_matches_old():
    """One bond, not a conjugated pair."""

    old = _old(OldDehydrogenation(), "CCO")
    new = _new(Dehydrogenation(), "CCO")
    assert new == old == {"CC=O", "C=CO"}


def test_ethenediol_dehydrogenation_includes_glyoxal():
    """Both alcohols dehydrogenate to glyoxal. The ketene is the one-bond SMARTS."""

    old = _old(OldDehydrogenation(), "OC=CO")
    new = _new(Dehydrogenation(), "OC=CO")
    assert new == old == {"O=C=CO", "O=CC=O"}


def test_acetaldehyde_oxygen_reduction_matches_old():
    old = _old(OldOxygenReduction(), "CC=O")
    new = _new(OxygenReduction(), "CC=O")
    assert new == old == {"CCO"}


def test_ethylamine_nitrogen_oxidation_matches_old():
    old = _old(OldNitrogenOxidation(), "CCN")
    new = _new(NitrogenOxidation(), "CCN")
    assert new == old == {"CCNO", "CCN=O"}


def test_dimethyl_sulfide_sulfur_oxidation_matches_old():
    old = _old(OldSulfurOxidation(), "CSC")
    new = _new(SulfurOxidation(), "CSC")
    assert "CS(C)=O" in new
    assert new == old == {"CS(C)=O", "C[S+](C)[O-]", "C[SH](C)O"}


def test_chloroethane_oxidative_dehalogenation_matches_old():
    old = _old(OldOxidativeDehalogenation(), "CCCl")
    new = _new(OxidativeDehalogenation(), "CCCl")
    assert new == old == {"CCO", "CC(=O)O", "Cl"}


def test_nitromethane_reduction_matches_old():
    old = _old(OldNitrogenReduction(), "C[N+](=O)[O-]")
    new = _new(NitrogenReduction(), "C[N+](=O)[O-]")
    assert new == old == {"CN", "CN=O", "O"}


def test_nitrosomethane_reduction_replaces_the_dummy():
    """The old `[*2]` atom is not oxygen. The leaving group here is `O`."""

    old = _old(OldNitrogenReduction(), "CN=O")
    new = _new(NitrogenReduction(), "CN=O")
    assert "CN" in new
    assert old - new == _OLD_NITROSO_DUMMY
    assert new - old == {"O"}


def test_methyl_acetate_hydrolysis_matches_old():
    """The product set is both fragments of each pattern, not one endpoint."""

    old = _old(OldHydrolysis(), "CC(=O)OC")
    new = _new(Hydrolysis(), "CC(=O)OC")
    assert new == old == {"CC(=O)O", "CC=O", "CO"}


def test_ethanol_dehydration_matches_old():
    old = _old(OldDehydration(), "CCO")
    new = _new(Dehydration(), "CCO")
    assert new == old == {"C=C", "CC", "O"}


def test_dimethyl_disulfide_sulfur_reduction_matches_old():
    old = _old(OldSulfurReduction(), "CSSC")
    new = _new(SulfurReduction(), "CSSC")
    assert "CS" in new
    assert new == old == {"CS", "C", "CSS"}


def test_chloroethane_reductive_dehalogenation_matches_old():
    old = _old(OldReductiveDehalogenation(), "CCCl")
    new = _new(ReductiveDehalogenation(), "CCCl")
    assert new == old == {"CC", "C=C", "Cl"}


def test_ethylene_oxide_opening_matches_old():
    old = _old(OldEpoxideOpening(), "C1CO1")
    new = _new(EpoxideOpening(), "C1CO1")
    assert new == old == {"CCO", "OCCO"}


def test_methyl_phosphate_keeps_the_ester_cleavage():
    """Water and methyl phosphite are the old P-OH cleavage, not this ester."""

    old = _old(OldDephosphorylation(), "COP(=O)(O)O")
    new = _new(Dephosphorylation(), "COP(=O)(O)O")
    assert {"CO", "O=[PH](O)O"} <= new
    assert old - new == _OLD_METHYL_PHOSPHITE
    assert new - old == set()
