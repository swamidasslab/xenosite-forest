"""The four easy pairs, both libraries, canonical fragment SMILES.

Order does not matter. A product only the old library emits is a failure
unless `src/xenosite/refactor_poc/DIVERGENCES.md` records both sets and the
reason. These four are not xfail.
"""

from rdkit import Chem

from xenosite.forest.rules import Dealkylation as OldDealkylation
from xenosite.forest.rules import Hydroxylation as OldHydroxylation
from xenosite.forest.rules import QuinoneFormation as OldQuinone
from xenosite.refactor_poc.find_path import bfs
from xenosite.refactor_poc.rules import (
    Dealkylation,
    Dehydrogenation,
    Hydroxylation,
    QuinoneFormation,
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


def _new(rule, smiles):
    return _fragments(product for product, _info in rule.metabolize(Chem.MolFromSmiles(smiles)))


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
