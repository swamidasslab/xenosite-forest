"""RDKit names this proof of concept calls.

The ``TYPE_CHECKING`` branch is the type of those names. ``Mol._forest`` is
declared there and is not assigned, so molecules do not share one forest.
A constructor returns :class:`NoForestMol` when the object it builds has no
``_forest``. RDKit copies do not keep that attribute. A function that hands
back a molecule the caller already owned keeps that molecule's type.
At runtime the same names are the real RDKit objects.

Argument lists are the C++ signatures Boost printed when each used function
was called with the wrong arguments, narrowed to the overload the call site
uses. ``GetAtoms`` and ``GetBonds`` are Python wrappers on ``Mol``; Boost
does not wrap them.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import TYPE_CHECKING, Literal, overload

if TYPE_CHECKING:
    from xenosite.refactor_poc.records import Forest, TracingForest, UntracedForest

    class BondType:
        SINGLE: BondType
        DOUBLE: BondType
        TRIPLE: BondType
        AROMATIC: BondType
        UNSPECIFIED: BondType

    class Atom:
        def __init__(self, num: int) -> None: ...
        def GetAtomMapNum(self) -> int: ...
        def GetAtomicNum(self) -> int: ...
        def GetBonds(self) -> Iterator[Bond]: ...
        def GetFormalCharge(self) -> int: ...
        def GetIdx(self) -> int: ...
        def GetIsAromatic(self) -> bool: ...
        def GetIsotope(self) -> int: ...
        def GetNeighbors(self) -> Iterator[Atom]: ...
        def GetNumExplicitHs(self) -> int: ...
        def GetNumImplicitHs(self) -> int: ...
        def GetProp(self, key: str, autoConvert: bool = False) -> str: ...
        def GetSymbol(self) -> str: ...
        def GetTotalNumHs(self, includeNeighbors: bool = False) -> int: ...
        def HasProp(self, key: str) -> bool: ...
        def IsInRing(self) -> bool: ...
        def SetAtomMapNum(self, mapno: int, strict: bool = False) -> None: ...
        def SetAtomicNum(self, newNum: int) -> None: ...
        def SetFormalCharge(self, what: int) -> None: ...
        def SetIsAromatic(self, what: bool) -> None: ...
        def SetIsotope(self, what: int) -> None: ...
        def SetNoImplicit(self, what: bool) -> None: ...
        def SetNumExplicitHs(self, what: int) -> None: ...
        def SetProp(self, key: str, val: str) -> None: ...
        def UpdatePropertyCache(self, strict: bool = True) -> None: ...

    class Bond:
        def GetBeginAtomIdx(self) -> int: ...
        def GetBondType(self) -> BondType: ...
        def GetBondTypeAsDouble(self) -> float: ...
        def GetEndAtomIdx(self) -> int: ...
        def GetIsAromatic(self) -> bool: ...
        def GetOtherAtomIdx(self, idx: int) -> int: ...
        def IsInRing(self) -> bool: ...
        def SetBondType(self, bT: BondType) -> None: ...
        def SetIsAromatic(self, what: bool) -> None: ...

    class RingInfo:
        def AtomRings(self) -> tuple[tuple[int, ...], ...]: ...

    class Mol:
        _forest: None | Forest

        def __new__(
            cls, mol: Mol, quickCopy: bool = False, confId: int = -1
        ) -> NoForestMol: ...
        def __init__(
            self, mol: Mol, quickCopy: bool = False, confId: int = -1
        ) -> None: ...
        def GetAromaticAtoms(self) -> Iterator[Atom]: ...
        def GetAtomWithIdx(self, idx: int) -> Atom: ...
        def GetAtoms(self) -> Iterator[Atom]: ...
        def GetBondBetweenAtoms(self, idx1: int, idx2: int) -> Bond | None: ...
        def GetBonds(self) -> Iterator[Bond]: ...
        def GetNumAtoms(self, onlyHeavy: int = -1, onlyExplicit: bool = True) -> int: ...
        def GetNumHeavyAtoms(self) -> int: ...
        def GetProp(self, key: str, autoConvert: bool = False) -> str: ...
        def GetRingInfo(self) -> RingInfo: ...
        def GetSubstructMatches(
            self,
            query: Mol,
            uniquify: bool = True,
            useChirality: bool = False,
            useQueryQueryMatches: bool = False,
            maxMatches: int = 1000,
        ) -> tuple[tuple[int, ...], ...]: ...

    class NoTracingMol(Mol):
        """``atom_trace`` is absent. This does not say whether ``_forest`` exists."""

        # Instance attributes are invariant. These states are narrower on purpose.
        _forest: None | UntracedForest  # pyright: ignore[reportIncompatibleVariableOverride]

    class ForestMol(Mol):
        """``_forest`` is present. The trace may or may not be."""

        _forest: Forest  # pyright: ignore[reportIncompatibleVariableOverride]

    class NoForestMol(NoTracingMol):
        """No ``_forest``. A missing forest has no trace."""

        _forest: None  # pyright: ignore[reportIncompatibleVariableOverride]

    class ForestNoTracingMol(NoTracingMol, ForestMol):
        """``_forest`` is present and ``atom_trace`` is absent."""

        _forest: UntracedForest  # pyright: ignore[reportIncompatibleVariableOverride]

    class ForestTracingMol(ForestMol):
        """``_forest`` is present and ``atom_trace`` is initialized."""

        _forest: TracingForest  # pyright: ignore[reportIncompatibleVariableOverride]

    class RWMol(NoForestMol):
        def __new__(cls, m: Mol) -> RWMol: ...
        def __init__(self, m: Mol) -> None: ...
        def AddAtom(self, atom: Atom) -> int: ...
        def AddBond(
            self,
            beginAtomIdx: int,
            endAtomIdx: int,
            order: BondType = BondType.UNSPECIFIED,
        ) -> int: ...
        def GetMol(self) -> NoForestMol: ...
        def RemoveAtom(self, idx: int) -> None: ...
        def RemoveBond(self, idx1: int, idx2: int) -> None: ...

    class ChemicalReaction:
        def RunReactants(
            self, reactants: tuple[Mol, ...], maxProducts: int = 1000
        ) -> tuple[tuple[NoForestMol, ...], ...]: ...
        def _setImplicitPropertiesFlag(self, val: bool) -> None: ...

    class ResonanceMolSupplier:
        def __init__(
            self, mol: Mol, flags: int = 0, maxStructs: int = 1000
        ) -> None: ...
        def GetAtomConjGrpIdx(self, ai: int) -> int: ...
        def GetNumConjGrps(self) -> int: ...
        def __iter__(self) -> Iterator[NoForestMol | None]: ...

    class MCSResult:
        numAtoms: int
        canceled: bool
        smartsString: str

    class AtomCompare:
        CompareElements: AtomCompare

    class BondCompare:
        CompareAny: BondCompare

    class RingCompare:
        IgnoreRingFusion: RingCompare

    class SanitizeFlags:
        SANITIZE_ALL: int
        SANITIZE_SYMMRINGS: int

    KEKULE_ALL: int

    def CanonicalRankAtoms(
        mol: Mol,
        breakTies: bool = True,
        includeChirality: bool = True,
        includeIsotopes: bool = True,
        includeAtomMaps: bool = True,
        includeChiralPresence: bool = False,
    ) -> Sequence[int]: ...
    def DisableLog(spec: str) -> None: ...
    def FindMCS(
        mols: Sequence[Mol],
        maximizeBonds: bool = True,
        threshold: float = 1.0,
        timeout: int = 3600,
        verbose: bool = False,
        matchValences: bool = False,
        ringMatchesRingOnly: bool = False,
        completeRingsOnly: bool = False,
        matchChiralTag: bool = False,
        atomCompare: AtomCompare = AtomCompare.CompareElements,
        bondCompare: BondCompare = BondCompare.CompareAny,
        ringCompare: RingCompare = RingCompare.IgnoreRingFusion,
        seedSmarts: str = "",
    ) -> MCSResult: ...
    @overload
    def GetMolFrags(
        mol: Mol,
        asMols: Literal[True],
        sanitizeFrags: bool = True,
        frags: object = None,
        fragsMolAtomMapping: object = None,
    ) -> tuple[NoForestMol, ...]: ...
    @overload
    def GetMolFrags(
        mol: Mol,
        asMols: Literal[False] = False,
        sanitizeFrags: bool = True,
        frags: object = None,
        fragsMolAtomMapping: object = None,
    ) -> tuple[tuple[int, ...], ...]: ...
    def GetMolFrags(
        mol: Mol,
        asMols: bool = False,
        sanitizeFrags: bool = True,
        frags: object = None,
        fragsMolAtomMapping: object = None,
    ) -> tuple[NoForestMol, ...] | tuple[tuple[int, ...], ...]:
        raise AssertionError("rdkit_api.GetMolFrags is the real RDKit function at runtime")
    def MolFromSmarts(
        SMARTS: str,
        mergeHs: bool = False,
        replacements: dict[str, str] | None = None,
    ) -> NoForestMol | None: ...
    def MolFromSmiles(
        SMILES: str,
        sanitize: bool = True,
        replacements: dict[str, str] | None = None,
    ) -> NoForestMol | None: ...
    def MolToSmiles(
        mol: Mol,
        isomericSmiles: bool = True,
        kekuleSmiles: bool = False,
        rootedAtAtom: int = -1,
        canonical: bool = True,
        allBondsExplicit: bool = False,
        allHsExplicit: bool = False,
        doRandom: bool = False,
        ignoreAtomMapNumbers: bool = False,
    ) -> str: ...
    def ReactionFromSmarts(
        SMARTS: str,
        replacements: dict[str, str] | None = None,
        useSmiles: bool = False,
    ) -> ChemicalReaction: ...
    def RenumberAtoms(mol: Mol, newOrder: Sequence[int]) -> NoForestMol: ...
    def SanitizeMol(
        mol: Mol, sanitizeOps: int = SanitizeFlags.SANITIZE_ALL, catchErrors: bool = False
    ) -> int: ...

else:
    from rdkit import Chem, rdBase
    from rdkit.Chem import GetMolFrags as GetMolFrags
    from rdkit.Chem import SanitizeMol as SanitizeMol
    from rdkit.Chem import rdFMCS
    from rdkit.Chem.AllChem import CanonicalRankAtoms as CanonicalRankAtoms
    from rdkit.Chem.AllChem import MolToSmiles as MolToSmiles
    from rdkit.Chem.AllChem import ReactionFromSmarts as ReactionFromSmarts
    from rdkit.Chem.AllChem import RenumberAtoms as RenumberAtoms
    from rdkit.Chem.rdChemReactions import ChemicalReaction as ChemicalReaction
    from rdkit.Chem.rdchem import Atom as Atom
    from rdkit.Chem.rdchem import Bond as Bond
    from rdkit.Chem.rdchem import BondType as BondType
    from rdkit.Chem.rdchem import KEKULE_ALL as KEKULE_ALL
    from rdkit.Chem.rdchem import Mol as Mol
    from rdkit.Chem.rdchem import RWMol as RWMol
    from rdkit.Chem.rdchem import ResonanceMolSupplier as ResonanceMolSupplier

    AtomCompare = rdFMCS.AtomCompare
    BondCompare = rdFMCS.BondCompare
    RingCompare = rdFMCS.RingCompare
    FindMCS = rdFMCS.FindMCS
    DisableLog = rdBase.DisableLog
    MolFromSmiles = Chem.MolFromSmiles
    MolFromSmarts = Chem.MolFromSmarts
    SanitizeFlags = Chem.SanitizeFlags
    # Typing-only states. Runtime molecules stay RDKit's Mol; these names
    # exist so annotations can be imported. They are not a shared forest.
    NoForestMol = Mol
    NoTracingMol = Mol
    ForestMol = Mol
    ForestNoTracingMol = Mol
    ForestTracingMol = Mol
