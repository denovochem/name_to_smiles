from typing import Dict, Tuple
from dataclasses import dataclass, field
from enum import Enum, auto


class LigandType(Enum):
    """Classification of ligand charge types."""

    NEUTRAL = auto()
    ANIONIC = auto()
    CATIONIC = auto()


@dataclass
class LigandInfo:
    """
    Complete information about a ligand.

    Attributes:
        smiles: SMILES representation of the ligand
        denticity: Number of coordination sites (atoms that bind to metal)
        charge: Formal charge of the ligand
        aliases: Alternative names/abbreviations for this ligand
        description: Human-readable description
    """

    smiles: str
    denticity: int = 1
    charge: int = 0
    aliases: Tuple[str, ...] = field(default_factory=tuple)
    description: str = ""

    @property
    def ligand_type(self) -> LigandType:
        """Determine ligand type based on formal charge."""
        if self.charge < 0:
            return LigandType.ANIONIC
        elif self.charge > 0:
            return LigandType.CATIONIC
        return LigandType.NEUTRAL


@dataclass
class MetalInfo:
    """
    Information about a transition metal.

    Attributes:
        symbol: Element symbol (e.g., "Ir")
        name: Full element name (e.g., "Iridium")
        common_oxidation_states: List of typical oxidation states
        atomic_number: Atomic number of the element
    """

    symbol: str
    name: str
    common_oxidation_states: Tuple[int, ...]
    atomic_number: int


LIGAND_DATABASE: Dict[str, LigandInfo] = {
    # -------------------------------------------------------------------------
    # Monodentate Neutral Ligands
    # -------------------------------------------------------------------------
    "CO": LigandInfo(
        smiles="[C-]#[O+]",
        denticity=1,
        charge=0,
        aliases=("carbonyl",),
        description="Carbonyl ligand",
    ),
    "PPh3": LigandInfo(
        smiles="c1ccc(P(c2ccccc2)c3ccccc3)cc1",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine", "Ph3P"),
        description="Triphenylphosphine",
    ),
    "py": LigandInfo(
        smiles="c1ccncc1",
        denticity=1,
        charge=0,
        aliases=("pyridine", "Py"),
        description="Pyridine",
    ),
    "NH3": LigandInfo(
        smiles="N",
        denticity=1,
        charge=0,
        aliases=("ammonia", "ammine"),
        description="Ammonia/Ammine ligand",
    ),
    "H2O": LigandInfo(
        smiles="O",
        denticity=1,
        charge=0,
        aliases=("water", "aqua", "aquo"),
        description="Water/Aqua ligand",
    ),
    "MeCN": LigandInfo(
        smiles="CC#N",
        denticity=1,
        charge=0,
        aliases=("acetonitrile", "NCMe", "CH3CN"),
        description="Acetonitrile",
    ),
    "PMe3": LigandInfo(
        smiles="CP(C)C",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphine",),
        description="Trimethylphosphine",
    ),
    "PEt3": LigandInfo(
        smiles="CCP(CC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylphosphine",),
        description="Triethylphosphine",
    ),
    "PiPr3": LigandInfo(
        smiles="CC(C)P(C(C)C)C(C)C",
        denticity=1,
        charge=0,
        aliases=("triisopropylphosphine", "P(iPr)3"),
        description="Triisopropylphosphine",
    ),
    "PCy3": LigandInfo(
        smiles="C1(CCCCC1)P(C2CCCCC2)C3CCCCC3",
        denticity=1,
        charge=0,
        aliases=("tricyclohexylphosphine",),
        description="Tricyclohexylphosphine",
    ),
    "PtBu3": LigandInfo(
        smiles="CCCCP(CCCC)CCCC",
        denticity=1,
        charge=0,
        aliases=("tri-tert-butylphosphine", "P(tBu)3"),
        description="Tri-tert-butylphosphine",
    ),
    "P(OMe)3": LigandInfo(
        smiles="COP(C)(OC)OC",
        denticity=1,
        charge=0,
        aliases=("trimethylphosphite",),
        description="Trimethyl phosphite",
    ),
    "P(OEt)3": LigandInfo(
        smiles="O(P(OCC)OCC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylphosphite",),
        description="Triethyl phosphite",
    ),
    "P(OPh)3": LigandInfo(
        smiles="O(P(Oc1ccccc1)Oc2ccccc2)c3ccccc3",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphite",),
        description="Triphenyl phosphite",
    ),
    # N-Heterocyclic Carbenes (NHCs) - hugely important in modern catalysis
    "IMes": LigandInfo(
        smiles="Cc1cc(c(c(c1)C)N2C=CN([C]2)c3c(cc(cc3C)C)C)C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazol-2-ylidene",),
        description="IMes carbene",
    ),
    "IPr": LigandInfo(
        smiles="CC(C)c1cccc(C(C)C)c1N2[C]N(C=C2)c3c(cccc3C(C)C)C(C)C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene",),
        description="IPr carbene",
    ),
    "SIMes": LigandInfo(
        smiles="Cc1cc(C)c(c(C)c1)N2CCN([C]2)c3c(C)cc(C)cc3C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,4,6-trimethylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IMes carbene",
    ),
    "SIPr": LigandInfo(
        smiles="CC(C)C1=C(C(=CC=C1)C(C)C)N2CC[N+](=[C-]2)C3=C(C=CC=C3C(C)C)C(C)C",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(2,6-diisopropylphenyl)imidazolidin-2-ylidene",),
        description="Saturated IPr carbene",
    ),
    "ICy": LigandInfo(
        smiles="C1(CCCCC1)N1[C]N(C=C1)C1CCCCC1",
        denticity=1,
        charge=0,
        aliases=("1,3-dicyclohexylimidazol-2-ylidene",),
        description="ICy carbene",
    ),
    "ItBu": LigandInfo(
        smiles="C(C)(C)(C)N1[C]N(C=C1)C(C)(C)C",
        denticity=1,
        charge=0,
        aliases=("1,3-di-tert-butylimidazol-2-ylidene",),
        description="ItBu carbene",
    ),
    "IMe": LigandInfo(
        smiles="CN1[C]N(C=C1)C",
        denticity=1,
        charge=0,
        aliases=("1,3-dimethylimidazol-2-ylidene",),
        description="IMe carbene",
    ),
    "IAd": LigandInfo(
        smiles="C12(CC3CC(CC(C1)C3)C2)N2[C]N(C=C2)C23CC1CC(CC(C2)C1)C3",
        denticity=1,
        charge=0,
        aliases=("1,3-bis(adamantyl)imidazol-2-ylidene",),
        description="IAd carbene",
    ),
    # Other common neutral donors
    "THF": LigandInfo(
        smiles="C1CCOC1",
        denticity=1,
        charge=0,
        aliases=("tetrahydrofuran", "thf"),
        description="Tetrahydrofuran",
    ),
    "Et2O": LigandInfo(
        smiles="CCOCC",
        denticity=1,
        charge=0,
        aliases=("diethylether", "ether"),
        description="Diethyl ether",
    ),
    "DMF": LigandInfo(
        smiles="CN(C)C=O",
        denticity=1,
        charge=0,
        aliases=("dimethylformamide",),
        description="Dimethylformamide",
    ),
    "DMSO": LigandInfo(
        smiles="CS(=O)C",
        denticity=1,
        charge=0,
        aliases=("dimethylsulfoxide",),
        description="Dimethyl sulfoxide",
    ),
    "NMe3": LigandInfo(
        smiles="CN(C)C",
        denticity=1,
        charge=0,
        aliases=("trimethylamine",),
        description="Trimethylamine",
    ),
    "NEt3": LigandInfo(
        smiles="CCN(CC)CC",
        denticity=1,
        charge=0,
        aliases=("triethylamine", "TEA"),
        description="Triethylamine",
    ),
    "DMAP": LigandInfo(
        smiles="Cc1ccccc1N(C)C",
        denticity=1,
        charge=0,
        aliases=("4-dimethylaminopyridine",),
        description="4-Dimethylaminopyridine",
    ),
    # Isocyanides
    "CNtBu": LigandInfo(
        smiles="CC(C)(C)C#N",
        denticity=1,
        charge=0,
        aliases=("tert-butylisocyanide",),
        description="tert-Butyl isocyanide",
    ),
    "CNXyl": LigandInfo(
        smiles="CC1=C(C(=CC=C1)C)[N+]#[C-]",
        denticity=1,
        charge=0,
        aliases=("2,6-xylylisocyanide", "2,6-dimethylphenylisocyanide"),
        description="2,6-Xylyl isocyanide",
    ),
    "CNPh": LigandInfo(
        smiles="O=C=Nc1ccccc1",
        denticity=1,
        charge=0,
        aliases=("phenylisocyanide",),
        description="Phenyl isocyanide",
    ),
    "CNMe": LigandInfo(
        smiles="[C-]#[N+]C",
        denticity=1,
        charge=0,
        aliases=("methylisocyanide",),
        description="Methyl isocyanide",
    ),
    # Carbene ligands (non-NHC)
    "CHPh": LigandInfo(
        smiles="[CH]C1=CC=CC=C1",
        denticity=1,
        charge=0,
        aliases=("benzylidene", "phenylcarbene"),
        description="Benzylidene",
    ),
    # Olefins (η²)
    "ethylene": LigandInfo(
        smiles="C=C",
        denticity=1,
        charge=0,
        aliases=("C2H4", "eth"),
        description="Ethylene (η²)",
    ),
    # Dinitrogen and other small molecules
    "N2": LigandInfo(
        smiles="N#N",
        denticity=1,
        charge=0,
        aliases=("dinitrogen",),
        description="Dinitrogen",
    ),
    "NO": LigandInfo(
        smiles="N=O",
        denticity=1,
        charge=0,
        aliases=("nitrosyl",),
        description="Nitrosyl (neutral counting)",
    ),
    "CS": LigandInfo(
        smiles="C=S",
        denticity=1,
        charge=0,
        aliases=("thiocarbonyl",),
        description="Thiocarbonyl",
    ),
    "SO2": LigandInfo(
        smiles="O=S=O",
        denticity=1,
        charge=0,
        aliases=("sulfurdioxide",),
        description="Sulfur dioxide",
    ),
    "O2": LigandInfo(
        smiles="O=O",
        denticity=1,
        charge=0,
        aliases=("dioxygen",),
        description="Dioxygen",
    ),
    "H2": LigandInfo(
        smiles="H-H",
        denticity=1,
        charge=0,
        aliases=("dihydrogen",),
        description="Dihydrogen (η²-H2)",
    ),
    # -------------------------------------------------------------------------
    # Monodentate Anionic Ligands
    # -------------------------------------------------------------------------
    "Cl": LigandInfo(
        smiles="[Cl-]",
        denticity=1,
        charge=-1,
        aliases=("chloro", "chloride", "chlorido"),
        description="Chloride",
    ),
    "Br": LigandInfo(
        smiles="[Br-]",
        denticity=1,
        charge=-1,
        aliases=("bromo", "bromide", "bromido"),
        description="Bromide",
    ),
    "I": LigandInfo(
        smiles="[I-]",
        denticity=1,
        charge=-1,
        aliases=("iodo", "iodide", "iodido"),
        description="Iodide",
    ),
    "F": LigandInfo(
        smiles="[F-]",
        denticity=1,
        charge=-1,
        aliases=("fluoro", "fluoride", "fluorido"),
        description="Fluoride",
    ),
    "H": LigandInfo(
        smiles="[H-]",
        denticity=1,
        charge=-1,
        aliases=("hydrido", "hydride"),
        description="Hydride",
    ),
    "CN": LigandInfo(
        smiles="[C-]#N",
        denticity=1,
        charge=-1,
        aliases=("cyano", "cyanide", "cyanido"),
        description="Cyanide",
    ),
    "OAc": LigandInfo(
        smiles="CC([O-])=O",
        denticity=1,
        charge=-1,
        aliases=("acetato", "acetate", "OAc"),
        description="Acetate",
    ),
    "OMe": LigandInfo(
        smiles="[O-]C",
        denticity=1,
        charge=-1,
        aliases=("methoxo", "methoxide"),
        description="Methoxide",
    ),
    # Common anions
    "OH": LigandInfo(
        smiles="[OH-]",
        denticity=1,
        charge=-1,
        aliases=("hydroxo", "hydroxide", "hydroxido"),
        description="Hydroxide",
    ),
    "OEt": LigandInfo(
        smiles="[O-]C(C)C",
        denticity=1,
        charge=-1,
        aliases=("ethoxo", "ethoxide"),
        description="Ethoxide",
    ),
    "OiPr": LigandInfo(
        smiles="[O-]C(C)C",
        denticity=1,
        charge=-1,
        aliases=("isopropoxo", "isopropoxide"),
        description="Isopropoxide",
    ),
    "OtBu": LigandInfo(
        smiles="[O-]C(C)(C)C",
        denticity=1,
        charge=-1,
        aliases=("tert-butoxo", "tert-butoxide"),
        description="tert-Butoxide",
    ),
    "OPh": LigandInfo(
        smiles="[O-]c1ccccc1",
        denticity=1,
        charge=-1,
        aliases=("phenoxo", "phenoxide", "phenolate"),
        description="Phenoxide",
    ),
    # Alkyls
    "Me": LigandInfo(
        smiles="[CH3-]",
        denticity=1,
        charge=-1,
        aliases=("methyl", "CH3"),
        description="Methyl",
    ),
    "Et": LigandInfo(
        smiles="[CH2-]C",
        denticity=1,
        charge=-1,
        aliases=("ethyl", "C2H5"),
        description="Ethyl",
    ),
    "nBu": LigandInfo(
        smiles="[CH2-]CC",
        denticity=1,
        charge=-1,
        aliases=("n-butyl", "butyl"),
        description="n-Butyl",
    ),
    "Ph": LigandInfo(
        smiles="C1=CC=[C-]C=C1",
        denticity=1,
        charge=-1,
        aliases=("phenyl", "C6H5"),
        description="Phenyl",
    ),
    "Bn": LigandInfo(
        smiles="[CH2-]C1=CC=CC=C1",
        denticity=1,
        charge=-1,
        aliases=("benzyl",),
        description="Benzyl",
    ),
    "vinyl": LigandInfo(
        smiles="[C-]=C",
        denticity=1,
        charge=-1,
        aliases=("ethenyl",),
        description="Vinyl",
    ),
    "allyl": LigandInfo(
        smiles="[CH-]C=C",
        denticity=1,
        charge=-1,
        aliases=("η1-allyl", "propenyl", "C3H5"),
        description="Allyl (σ-bound)",
    ),
    "Np": LigandInfo(
        smiles="CC(C)(C)[CH2-]",
        denticity=1,
        charge=-1,
        aliases=("neopentyl",),
        description="Neopentyl",
    ),
    "Mes": LigandInfo(
        smiles="[c-]1c(C)cc(C)cc1C",
        denticity=1,
        charge=-1,
        aliases=("mesityl", "2,4,6-trimethylphenyl"),
        description="Mesityl",
    ),
    # Silyls
    "SiMe3": LigandInfo(
        smiles="C[Si-](C)C",
        denticity=1,
        charge=-1,
        aliases=("trimethylsilyl", "TMS"),
        description="Trimethylsilyl",
    ),
    "SiPh3": LigandInfo(
        smiles="C1(=CC=CC=C1)[Si-](C1=CC=CC=C1)C1=CC=CC=C1",
        denticity=1,
        charge=-1,
        aliases=("triphenylsilyl",),
        description="Triphenylsilyl",
    ),
    # Amides
    "NMe2": LigandInfo(
        smiles="C[N-]C",
        denticity=1,
        charge=-1,
        aliases=("dimethylamido",),
        description="Dimethylamide",
    ),
    "NEt2": LigandInfo(
        smiles="CC[N-]CC",
        denticity=1,
        charge=-1,
        aliases=("diethylamido",),
        description="Diethylamide",
    ),
    "NiPr2": LigandInfo(
        smiles="C(C)(C)[N-]C(C)C",
        denticity=1,
        charge=-1,
        aliases=("diisopropylamido",),
        description="Diisopropylamide",
    ),
    "NPh2": LigandInfo(
        smiles="C1(=CC=CC=C1)[N-]C1=CC=CC=C1",
        denticity=1,
        charge=-1,
        aliases=("diphenylamido",),
        description="Diphenylamide",
    ),
    "NTMS2": LigandInfo(
        smiles="C[Si](C)(C)[N-][Si](C)(C)C",
        denticity=1,
        charge=-1,
        aliases=("bis(trimethylsilyl)amido", "HMDS", "N(SiMe3)2"),
        description="Bis(trimethylsilyl)amide",
    ),
    "NHPh": LigandInfo(
        smiles="C1(=CC=CC=C1)[NH-]",
        denticity=1,
        charge=-1,
        aliases=("anilido", "phenylamido"),
        description="Anilide",
    ),
    # Other common anionic
    "SCN": LigandInfo(
        smiles="[S-]C#N",
        denticity=1,
        charge=-1,
        aliases=("thiocyanato", "thiocyanate"),
        description="Thiocyanate (S-bound)",
    ),
    "NCS": LigandInfo(
        smiles="SC#[N-]",
        denticity=1,
        charge=-1,
        aliases=("isothiocyanato", "isothiocyanate"),
        description="Isothiocyanate (N-bound)",
    ),
    "N3": LigandInfo(
        smiles="[N-]=[N+]=[N-]",
        denticity=1,
        charge=-1,
        aliases=("azido", "azide"),
        description="Azide",
    ),
    "NO2": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("nitrito", "nitrite"),
        description="Nitrite (N-bound nitro)",
    ),
    "ONO": LigandInfo(
        smiles="N(=O)[O-]",
        denticity=1,
        charge=-1,
        aliases=("nitrito-O",),
        description="Nitrite (O-bound nitrito)",
    ),
    "SH": LigandInfo(
        smiles="[SH-]",
        denticity=1,
        charge=-1,
        aliases=("mercapto", "sulfhydryl", "thiolato"),
        description="Hydrosulfide/Thiolate",
    ),
    "SPh": LigandInfo(
        smiles="C1(=CC=CC=C1)[S-]",
        denticity=1,
        charge=-1,
        aliases=("thiophenolate", "phenylthiolato"),
        description="Thiophenolate",
    ),
    "SMe": LigandInfo(
        smiles="[S-]C",
        denticity=1,
        charge=-1,
        aliases=("methylthiolate", "methanethiolato"),
        description="Methylthiolate",
    ),
    "StBu": LigandInfo(
        smiles="[S-]C(C)(C)C",
        denticity=1,
        charge=-1,
        aliases=("tert-butylthiolate",),
        description="tert-Butylthiolate",
    ),
    "OCN": LigandInfo(
        smiles="[O-]C#N",
        denticity=1,
        charge=-1,
        aliases=("cyanato", "cyanate"),
        description="Cyanate",
    ),
    "NCO": LigandInfo(
        smiles="[N-]=C=O",
        denticity=1,
        charge=-1,
        aliases=("isocyanato", "isocyanate"),
        description="Isocyanate",
    ),
    # Carboxylates
    "OBz": LigandInfo(
        smiles="C(C1=CC=CC=C1)(=O)[O-]",
        denticity=1,
        charge=-1,
        aliases=("benzoato", "benzoate"),
        description="Benzoate",
    ),
    "OPiv": LigandInfo(
        smiles="CC(C(=O)[O-])(C)C",
        denticity=1,
        charge=-1,
        aliases=("pivalato", "pivalate", "trimethylacetate"),
        description="Pivalate",
    ),
    "O2CCF3": LigandInfo(
        smiles="FC(C(=O)[O-])(F)F",
        denticity=1,
        charge=-1,
        aliases=("trifluoroacetato", "trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
    ),
    "formate": LigandInfo(
        smiles="C(=O)[O-]",
        denticity=1,
        charge=-1,
        aliases=("formato", "HCO2"),
        description="Formate",
    ),
    # Oxo and related
    "O": LigandInfo(
        smiles="[O-2]",
        denticity=1,
        charge=-2,
        aliases=("oxo", "oxide", "oxido"),
        description="Oxo",
    ),
    "S": LigandInfo(
        smiles="[S-2]",
        denticity=1,
        charge=-2,
        aliases=("sulfido", "sulfide"),
        description="Sulfido",
    ),
    "Se": LigandInfo(
        smiles="[Se-2]",
        denticity=1,
        charge=-2,
        aliases=("selenido", "selenide"),
        description="Selenido",
    ),
    "Te": LigandInfo(
        smiles="[Te-2]",
        denticity=1,
        charge=-2,
        aliases=("tellurido", "telluride"),
        description="Tellurido",
    ),
    "NR": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("imido",),
        description="Imido (generic)",
    ),
    "NAr": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("arylimido",),
        description="Aryl imido",
    ),
    "NtBu": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("tert-butylimido",),
        description="tert-Butyl imido",
    ),
    "NAd": LigandInfo(
        smiles="",
        denticity=1,
        charge=-2,
        aliases=("adamantylimido",),
        description="Adamantyl imido",
    ),
    # Borohydrides
    "BH4": LigandInfo(
        smiles="[BH4-]",
        denticity=1,
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride",
    ),
    # -------------------------------------------------------------------------
    # Bidentate Neutral Ligands
    # -------------------------------------------------------------------------
    "bpy": LigandInfo(
        smiles="c1ccc(-c2ccccn2)nc1",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyridine", "bipyridine", "bipy"),
        description="2,2'-Bipyridine",
    ),
    "dtbbpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(-c2cc(C(C)(C)C)ccn2)c1",
        denticity=2,
        charge=0,
        aliases=(
            "4,4'-di-tert-butyl-2,2'-bipyridine",
            "di-tert-butylbipyridine",
            "dtbpy",
        ),
        description="4,4'-Di-tert-butyl-2,2'-bipyridine",
    ),
    "phen": LigandInfo(
        smiles="c1cnc2c(c1)ccc1cccnc12",
        denticity=2,
        charge=0,
        aliases=("1,10-phenanthroline", "phenanthroline"),
        description="1,10-Phenanthroline",
    ),
    "en": LigandInfo(
        smiles="NCCN",
        denticity=2,
        charge=0,
        aliases=("ethylenediamine",),
        description="Ethylenediamine",
    ),
    "cod": LigandInfo(
        smiles="C1=CCCC=CCC1",
        denticity=2,
        charge=0,
        aliases=("1,5-cyclooctadiene", "cyclooctadiene", "COD"),
        description="1,5-Cyclooctadiene (η⁴)",
    ),
    "nbd": LigandInfo(
        smiles="C1C2C=CC1C=C2",
        denticity=2,
        charge=0,
        aliases=("norbornadiene", "2,5-norbornadiene"),
        description="Norbornadiene",
    ),
    "dppe": LigandInfo(
        smiles="c1ccc(P(CCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diphenylphosphino)ethane",),
        description="1,2-Bis(diphenylphosphino)ethane",
    ),
    "dppm": LigandInfo(
        smiles="c1ccc(P(CP(c2ccccc2)c2ccccc2)c2ccccc2)cc1",
        denticity=2,
        charge=0,
        aliases=("bis(diphenylphosphino)methane",),
        description="Bis(diphenylphosphino)methane",
    ),
    # Bipyridines and derivatives
    "4,4'-dmbpy": LigandInfo(
        smiles="CC1=CC(=NC=C1)C2=NC=CC(=C2)C",
        denticity=2,
        charge=0,
        aliases=("4,4'-dimethyl-2,2'-bipyridine", "dmb"),
        description="4,4'-Dimethyl-2,2'-bipyridine",
    ),
    "5,5'-dmbpy": LigandInfo(
        smiles="CC1=CN=C(C=C1)C2=NC=C(C=C2)C",
        denticity=2,
        charge=0,
        aliases=("5,5'-dimethyl-2,2'-bipyridine",),
        description="5,5'-Dimethyl-2,2'-bipyridine",
    ),
    "dCbpy": LigandInfo(
        smiles="C1=CN=C(C=C1C(=O)[O-])C2=NC=CC(=C2)C(=O)[O-]",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxy-2,2'-bipyridine",),
        description="4,4'-Dicarboxy-2,2'-bipyridine",
    ),
    "dCEbpy": LigandInfo(
        smiles="C1(C2=NC=CC(C(OCC)=O)=C2)=NC=CC(C(OCC)=O)=C1",
        denticity=2,
        charge=0,
        aliases=("4,4'-dicarboxyethyl-2,2'-bipyridine",),
        description="4,4'-Diethoxycarbonyl-2,2'-bipyridine",
    ),
    # Phenanthroline derivatives
    "dmp": LigandInfo(
        smiles="CC(=O)OI1(c2ccccc2C(=O)O1)(OC(=O)C)OC(=O)C",
        denticity=2,
        charge=0,
        aliases=("2,9-dimethyl-1,10-phenanthroline", "neocuproine"),
        description="2,9-Dimethyl-1,10-phenanthroline",
    ),
    "dpp": LigandInfo(
        smiles="C1=CC=C(C=C1)C2=NC3=C(C=CC4=C3N=C(C=C4)C5=CC=CC=C5)C=C2",
        denticity=2,
        charge=0,
        aliases=("2,9-diphenyl-1,10-phenanthroline", "bathocuproine"),
        description="2,9-Diphenyl-1,10-phenanthroline",
    ),
    "tmp": LigandInfo(
        smiles="CC1=CN=C2C(=C1C)C=CC3=C(C(=CN=C32)C)C",
        denticity=2,
        charge=0,
        aliases=("3,4,7,8-tetramethyl-1,10-phenanthroline",),
        description="3,4,7,8-Tetramethyl-1,10-phenanthroline",
    ),
    # Bipyrimidines and related
    "bpm": LigandInfo(
        smiles="C1=CN=C(N=C1)C2=NC=CC=N2",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrimidine", "bipyrimidine"),
        description="2,2'-Bipyrimidine",
    ),
    "bpz": LigandInfo(
        smiles="C1=CN=C(C=N1)C2=NC=CN=C2",
        denticity=2,
        charge=0,
        aliases=("2,2'-bipyrazine", "bipyrazine"),
        description="2,2'-Bipyrazine",
    ),
    # Diphosphines (very important in catalysis)
    "dppp": LigandInfo(
        smiles="P(c1ccccc1)(c2ccccc2)CCCP(c3ccccc3)c4ccccc4",
        denticity=2,
        charge=0,
        aliases=("1,3-bis(diphenylphosphino)propane",),
        description="1,3-Bis(diphenylphosphino)propane",
    ),
    "dppb": LigandInfo(
        smiles="C1=CC=C(C=C1)P(CCCCP(C2=CC=CC=C2)C3=CC=CC=C3)C4=CC=CC=C4",
        denticity=2,
        charge=0,
        aliases=("1,4-bis(diphenylphosphino)butane",),
        description="1,4-Bis(diphenylphosphino)butane",
    ),
    "dppf": LigandInfo(
        smiles="c1ccc(cc1)P(c2ccccc2)C34C5[Fe]3678912(C5C6C74)C3C8C9C1(C23)P(c1ccccc1)c1ccccc1",
        denticity=2,
        charge=0,
        aliases=("1,1'-bis(diphenylphosphino)ferrocene",),
        description="1,1'-Bis(diphenylphosphino)ferrocene",
    ),
    "dcpe": LigandInfo(
        smiles="C1CCC(CC1)P(CCP(C2CCCCC2)C3CCCCC3)C4CCCCC4",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dicyclohexylphosphino)ethane",),
        description="1,2-Bis(dicyclohexylphosphino)ethane",
    ),
    "dcpm": LigandInfo(
        smiles="C1CCC(CC1)P(CP(C2CCCCC2)C3CCCCC3)C4CCCCC4",
        denticity=2,
        charge=0,
        aliases=("bis(dicyclohexylphosphino)methane",),
        description="Bis(dicyclohexylphosphino)methane",
    ),
    "dmpe": LigandInfo(
        smiles="P(C)(C)CCP(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(dimethylphosphino)ethane",),
        description="1,2-Bis(dimethylphosphino)ethane",
    ),
    "depe": LigandInfo(
        smiles="CCP(CC)CCP(CC)CC",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diethylphosphino)ethane",),
        description="1,2-Bis(diethylphosphino)ethane",
    ),
    "dippe": LigandInfo(
        smiles="P(C(C)C)(CCP(C(C)C)C(C)C)C(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(diisopropylphosphino)ethane",),
        description="1,2-Bis(diisopropylphosphino)ethane",
    ),
    "dtbpe": LigandInfo(
        smiles="CC(C)(C)P(CCP(C(C)(C)C)C(C)(C)C)C(C)(C)C",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(di-tert-butylphosphino)ethane",),
        description="1,2-Bis(di-tert-butylphosphino)ethane",
    ),
    # BINAP and related chiral phosphines (critical in asymmetric catalysis)
    "BINAP": LigandInfo(
        smiles="c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3c5c6ccccc6ccc5P(c7ccccc7)c8ccccc8",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(diphenylphosphino)-1,1'-binaphthyl",),
        description="BINAP",
    ),
    "TolBINAP": LigandInfo(
        smiles="Cc1ccc(cc1)P(c2ccc(C)cc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccc(C)cc7)c8ccc(C)cc8",
        denticity=2,
        charge=0,
        aliases=("2,2'-bis(di-p-tolylphosphino)-1,1'-binaphthyl",),
        description="Tol-BINAP",
    ),
    "SEGPHOS": LigandInfo(
        smiles="C1OC2=C(O1)C(=C(C=C2)P(C3=CC=CC=C3)C4=CC=CC=C4)C5=C(C=CC6=C5OCO6)P(C7=CC=CC=C7)C8=CC=CC=C8",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(diphenylphosphino)-4,4'-bi-1,3-benzodioxole",),
        description="SEGPHOS",
    ),
    "DM-SEGPHOS": LigandInfo(
        smiles="Cc1cc(C)cc(c1)P(c2cc(C)cc(C)c2)c3ccc4OCOc4c3-c5c6OCOc6ccc5P(c7cc(C)cc(C)c7)c8cc(C)cc(C)c8",
        denticity=2,
        charge=0,
        aliases=("5,5'-bis(di(3,5-xylyl)phosphino)-4,4'-bi-1,3-benzodioxole",),
        description="DM-SEGPHOS",
    ),
    "DIFLUORPHOS": LigandInfo(
        smiles="FC1(F)Oc2ccc(P(c3ccccc3)c4ccccc4)c(c2O1)c5c6OC(F)(F)Oc6ccc5P(c7ccccc7)c8ccccc8",
        denticity=2,
        charge=0,
        aliases=(),
        description="DIFLUORPHOS",
    ),
    "MeO-BIPHEP": LigandInfo(
        smiles="COC1=CC=CC(=C1C1=C(C=CC=C1OC)P(C1=CC=CC=C1)C1=CC=CC=C1)P(C1=CC=CC=C1)C1=CC=CC=C1",
        denticity=2,
        charge=0,
        aliases=("6,6'-dimethoxy-2,2'-bis(diphenylphosphino)-1,1'-biphenyl",),
        description="MeO-BIPHEP",
    ),
    # Josiphos-type ligands
    "Josiphos": LigandInfo(
        smiles="[Fe].[CH]1[CH][CH][CH][CH]1.C[C@@H]([C]2[CH][CH][CH][C]2P(c3ccccc3)c4ccccc4)P(C5CCCCC5)C6CCCCC6",
        denticity=2,
        charge=0,
        aliases=(),
        description="Josiphos",
    ),
    # Other chiral ligands
    "CHIRAPHOS": LigandInfo(
        smiles="P(c1ccccc1)(c2ccccc2)[C@H]([C@@H](P(c3ccccc3)c4ccccc4)C)C",
        denticity=2,
        charge=0,
        aliases=("2,3-bis(diphenylphosphino)butane",),
        description="CHIRAPHOS",
    ),
    "DIOP": LigandInfo(
        smiles="CC1(O[C@H]([C@@H](O1)CP(c2ccccc2)c3ccccc3)CP(c4ccccc4)c5ccccc5)C",
        denticity=2,
        charge=0,
        aliases=(
            "2,3-O-isopropylidene-2,3-dihydroxy-1,4-bis(diphenylphosphino)butane",
        ),
        description="DIOP",
    ),
    "DuPhos": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)benzene",),
        description="DuPhos",
    ),
    "BPE": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("1,2-bis(phospholano)ethane",),
        description="BPE",
    ),
    # Diamines
    "tmeda": LigandInfo(
        smiles="CN(C)CCN(C)C",
        denticity=2,
        charge=0,
        aliases=("N,N,N',N'-tetramethylethylenediamine", "TMEDA"),
        description="TMEDA",
    ),
    "dach": LigandInfo(
        smiles="C1CCC(C(C1)N)N",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminocyclohexane", "chxn"),
        description="1,2-Diaminocyclohexane",
    ),
    "dpen": LigandInfo(
        smiles="N[C@H]([C@@H](N)c1ccccc1)c2ccccc2",
        denticity=2,
        charge=0,
        aliases=("1,2-diphenylethylenediamine",),
        description="1,2-Diphenylethylenediamine",
    ),
    "pn": LigandInfo(
        smiles="CC(N)CN",
        denticity=2,
        charge=0,
        aliases=("1,2-diaminopropane", "propylenediamine"),
        description="1,2-Diaminopropane",
    ),
    "bn": LigandInfo(
        smiles="NC(C)C(N)C",
        denticity=2,
        charge=0,
        aliases=("2,3-diaminobutane", "2,3-butanediamine"),
        description="2,3-Diaminobutane",
    ),
    # Diimine ligands
    "DAB": LigandInfo(
        smiles="N=CC=N",
        denticity=2,
        charge=0,
        aliases=("1,4-diazabutadiene",),
        description="1,4-Diazabutadiene",
    ),
    "Ar-DAB": LigandInfo(
        smiles="",
        denticity=2,
        charge=0,
        aliases=("N,N'-diaryl-1,4-diazabutadiene",),
        description="N,N'-Diaryl-1,4-diazabutadiene",
    ),
    "dpp-BIAN": LigandInfo(
        smiles="C(C)(C)C1=C(C(=CC=C1)C(C)C)C=1C(=C2C(C(C=3C=CC=C(C1)C32)=N)=N)C3=C(C=CC=C3C(C)C)C(C)C",
        denticity=2,
        charge=0,
        aliases=("bis(2,6-diisopropylphenyl)acenaphthenequinonediimine",),
        description="dpp-BIAN",
    ),
    # Mixed P,N donors
    "PHOX": LigandInfo(
        smiles="CC(C)(C)[C@@H]1N=C(C2=CC=CC=C2P(C3=CC=CC=C3)C4=CC=CC=C4)OC1",
        denticity=2,
        charge=0,
        aliases=("phosphinooxazoline",),
        description="Phosphinooxazoline",
    ),
    # Schiff bases (common in coordination chemistry)
    "salen-H2": LigandInfo(
        smiles="C(C=1C(O)=CC=CC1)=NCCN=CC=1C(O)=CC=CC1",
        denticity=2,
        charge=0,
        aliases=("N,N'-bis(salicylidene)ethylenediamine",),
        description="Salen (neutral form)",
    ),
    # Other bidentate neutrals
    "dbm": LigandInfo(
        smiles="C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=CC=C2",
        denticity=2,
        charge=0,
        aliases=("dibenzoylmethane",),
        description="Dibenzoylmethane (neutral)",
    ),
    "dme": LigandInfo(
        smiles="COCCOC",
        denticity=2,
        charge=0,
        aliases=("1,2-dimethoxyethane", "glyme", "DME"),
        description="1,2-Dimethoxyethane",
    ),
    "diglyme": LigandInfo(
        smiles="COCCOCCOC",
        denticity=2,
        charge=0,
        aliases=("diethylene glycol dimethyl ether",),
        description="Diglyme",
    ),
    "OPPh3": LigandInfo(
        smiles="O=P(c1ccccc1)(c2ccccc2)c3ccccc3",
        denticity=1,
        charge=0,
        aliases=("triphenylphosphine oxide",),
        description="Triphenylphosphine oxide",
    ),
    # -------------------------------------------------------------------------
    # Bidentate Anionic Ligands
    # -------------------------------------------------------------------------
    "acac": LigandInfo(
        smiles="CC(=O)[CH-]C(=O)C",
        denticity=2,
        charge=-1,
        aliases=("acetylacetonate", "acetylacetonato"),
        description="Acetylacetonate",
    ),
    "ppy": LigandInfo(
        smiles="[c-]1ccccc1-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("2-phenylpyridine", "phenylpyridine", "phenylpyridinato"),
        description="2-Phenylpyridinate (C^N cyclometalating)",
    ),
    "dfppy": LigandInfo(
        smiles="Fc1cc(F)c([c-]1)-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("2-(2,4-difluorophenyl)pyridine",),
        description="2-(2,4-Difluorophenyl)pyridinate",
    ),
    "F2ppy": LigandInfo(
        smiles="Fc1cc(F)c([c-]1)-c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("difluorophenylpyridine", "dFppy", "dF-ppy"),
        description="Difluorophenylpyridinate",
    ),
    "dF(CF3)ppy": LigandInfo(
        smiles="FC1=C(C=CC(=C1)F)C1=NC=C(C=C1)C(F)(F)F",
        denticity=2,
        charge=-1,
        aliases=("2-(2,4-Difluorophenyl)-5-(trifluoromethyl)pyridine", "dFCF3ppy"),
        description="2-(2,4-Difluorophenyl)-5-(trifluoromethyl)pyridine",
    ),
    "pic": LigandInfo(
        smiles="[O-]C(=O)c1ccccn1",
        denticity=2,
        charge=-1,
        aliases=("picolinate", "picolinato"),
        description="Picolinate",
    ),
    # Cyclometalating C^N ligands
    "bzq": LigandInfo(
        smiles="N1=CC=CC2=CC=C3C(=C12)C=CC=C3",
        denticity=2,
        charge=-1,
        aliases=("benzo[h]quinoline", "7,8-benzoquinoline"),
        description="Benzo[h]quinolinate",
    ),
    "thpy": LigandInfo(
        smiles="[S-]1C(=CC=C1)C1=NC=CC=C1",
        denticity=2,
        charge=-1,
        aliases=("2-(2-thienyl)pyridine", "thienylpyridine"),
        description="2-(2-Thienyl)pyridinate",
    ),
    "piq": LigandInfo(
        smiles="C1(=CC=CC=C1)C1=NC=CC2=CC=CC=C12",
        denticity=2,
        charge=-1,
        aliases=("1-phenylisoquinoline", "phenylisoquinoline"),
        description="1-Phenylisoquinolinate",
    ),
    "mppy": LigandInfo(
        smiles="CC1=CC=C(C=C1)C1=NC=CC=C1",
        denticity=2,
        charge=-1,
        aliases=("2-(4-methylphenyl)pyridine", "4-methyl-2-phenylpyridine"),
        description="2-(4-Methylphenyl)pyridinate",
    ),
    "btp": LigandInfo(
        smiles="S1C(=CC2=C1C=CC=C2)C2=NC=CC=C2",
        denticity=2,
        charge=-1,
        aliases=("2-benzothienylpyridine",),
        description="2-Benzothienylpyridinate",
    ),
    "pbpy": LigandInfo(
        smiles="C1(=CC=CC=C1)C1=CC=CC(=N1)C1=NC=CC=C1",
        denticity=2,
        charge=-1,
        aliases=("6-phenyl-2,2'-bipyridine",),
        description="6-Phenyl-2,2'-bipyridinate",
    ),
    # Beta-diketonates
    "hfac": LigandInfo(
        smiles="FC(F)(F)C(=O)CC(=O)C(F)(F)F",
        denticity=2,
        charge=-1,
        aliases=("hexafluoroacetylacetonate", "1,1,1,5,5,5-hexafluoroacetylacetonate"),
        description="Hexafluoroacetylacetonate",
    ),
    "tfac": LigandInfo(
        smiles="CC(=O)CC(=O)C(F)(F)F",
        denticity=2,
        charge=-1,
        aliases=("trifluoroacetylacetonate",),
        description="Trifluoroacetylacetonate",
    ),
    # "dbm": LigandInfo(
    #     smiles="C(C1=CC=CC=C1)(=O)CC(C1=CC=CC=C1)=O",
    #     denticity=2,
    #     charge=-1,
    #     aliases=("dibenzoylmethanate",),
    #     description="Dibenzoylmethanate",
    # ),
    "thd": LigandInfo(
        smiles="CC(C)(C(CC(C(C)(C)C)=O)=O)C",
        denticity=2,
        charge=-1,
        aliases=("2,2,6,6-tetramethyl-3,5-heptanedionate", "tmhd", "dpm"),
        description="2,2,6,6-Tetramethyl-3,5-heptanedionate",
    ),
    "fod": LigandInfo(
        smiles="FC(C(CC(C(C)(C)C)=O)=O)(C(C(F)(F)F)(F)F)F",
        denticity=2,
        charge=-1,
        aliases=("6,6,7,7,8,8,8-heptafluoro-2,2-dimethyl-3,5-octanedionate",),
        description="FOD",
    ),
    "trop": LigandInfo(
        smiles="C1=CC=C(C(=O)C=C1)O",
        denticity=2,
        charge=-1,
        aliases=("tropolonate",),
        description="Tropolonate",
    ),
    # Carboxylates (bridging/chelating)
    "OAc-bi": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("acetate-O,O'",),
        description="Acetate (bidentate)",
    ),
    "CO3": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("carbonato", "carbonate"),
        description="Carbonate",
    ),
    "SO4": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("sulfato", "sulfate"),
        description="Sulfate",
    ),
    # N,N chelates (anionic)
    "pz": LigandInfo(
        smiles="",
        denticity=1,
        charge=-1,
        aliases=("pyrazolato", "pyrazolate"),
        description="Pyrazolate",
    ),
    "pypz": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("3-(2-pyridyl)pyrazolate",),
        description="3-(2-Pyridyl)pyrazolate",
    ),
    "indazolato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("indazolate",),
        description="Indazolate",
    ),
    # N,O chelates
    "quinolinolate": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("8-hydroxyquinolinate", "oxinate", "Q"),
        description="8-Quinolinolate",
    ),
    "glycinato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("glycinate", "gly"),
        description="Glycinate",
    ),
    "alaninato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("alaninate", "ala"),
        description="Alaninate",
    ),
    "salicylate": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=("salicylato", "sal"),
        description="Salicylate",
    ),
    "oxalato-mono": LigandInfo(
        smiles="",
        denticity=2,
        charge=-1,
        aliases=(),
        description="Oxalate (monoanionic, monodentate)",
    ),
    # O,O chelates
    "catecholato": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("catecholate", "cat"),
        description="Catecholate",
    ),
    "semiquinone": LigandInfo(
        smiles="", denticity=2, charge=-1, aliases=("sq",), description="Semiquinone"
    ),
    # S,S chelates
    "S2CNMe2": LigandInfo(
        smiles="CN(C([S-])=S)C",
        denticity=2,
        charge=-1,
        aliases=("dimethyldithiocarbamate", "Me2dtc"),
        description="Dimethyldithiocarbamate",
    ),
    "S2CNEt2": LigandInfo(
        smiles="C(C)N(C([S-])=S)CC",
        denticity=2,
        charge=-1,
        aliases=("diethyldithiocarbamate", "Et2dtc", "dtc"),
        description="Diethyldithiocarbamate",
    ),
    "S2COEt": LigandInfo(
        smiles="C(C)OC(=S)[S-]",
        denticity=2,
        charge=-1,
        aliases=("ethylxanthate", "xanthate"),
        description="Ethyl xanthate",
    ),
    "S2PPh2": LigandInfo(
        smiles="C1(=CC=CC=C1)P(=S)([S-])C1=CC=CC=C1",
        denticity=2,
        charge=-1,
        aliases=("diphenylphosphinodithioate", "dtp"),
        description="Diphenylphosphinodithioate",
    ),
    "bdt": LigandInfo(
        smiles="C=1(C(=CC=CC1)[S-])[S-]",
        denticity=2,
        charge=-2,
        aliases=("1,2-benzenedithiolate", "benzene-1,2-dithiolate"),
        description="1,2-Benzenedithiolate",
    ),
    "mnt": LigandInfo(
        smiles="",
        denticity=2,
        charge=-2,
        aliases=("maleonitriledithiolate", "1,2-dicyanoethylene-1,2-dithiolate"),
        description="Maleonitriledithiolate",
    ),
    "dmit": LigandInfo(
        smiles="S=C1SC(=C(S1)[S-])[S-]",
        denticity=2,
        charge=-2,
        aliases=("2-thioxo-1,3-dithiole-4,5-dithiolate",),
        description="dmit",
    ),
    "edt": LigandInfo(
        smiles="C(C[S-])[S-]",
        denticity=2,
        charge=-2,
        aliases=("ethane-1,2-dithiolate", "1,2-ethanedithiolate"),
        description="Ethane-1,2-dithiolate",
    ),
    "tdt": LigandInfo(
        smiles="CC1=CC(=C(C=C1)[S-])[S-]",
        denticity=2,
        charge=-2,
        aliases=("toluene-3,4-dithiolate",),
        description="Toluene-3,4-dithiolate",
    ),
    # -------------------------------------------------------------------------
    # Tridentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Terpyridines
    "tpy": LigandInfo(
        smiles="c1ccnc(c1)c2cccc(n2)c3ccccn3",
        denticity=3,
        charge=0,
        aliases=("2,2':6',2''-terpyridine", "terpyridine", "terpy"),
        description="2,2':6',2''-Terpyridine",
    ),
    "ttpy": LigandInfo(
        smiles="CC1=CC=C(C=C1)C2=CC(=NC(=C2)C3=CC=CC=N3)C4=CC=CC=N4",
        denticity=3,
        charge=0,
        aliases=("4'-p-tolyl-2,2':6',2''-terpyridine",),
        description="4'-p-Tolyl-terpyridine",
    ),
    "tBu3tpy": LigandInfo(
        smiles="CC(C)(C)c1ccnc(c1)-c2cc(cc(n2)-c3cc(ccn3)C(C)(C)C)C(C)(C)C",
        denticity=3,
        charge=0,
        aliases=("4,4',4''-tri-tert-butyl-2,2':6',2''-terpyridine",),
        description="4,4',4''-Tri-tert-butylterpyridine",
    ),
    # Pincer ligands (hugely important class)
    "PNP": LigandInfo(
        smiles="C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=NC(=CC=C3)P(C4=CC=CC=C4)C5=CC=CC=C5",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)pyridine",),
        description="PNP pincer",
    ),
    "PCP": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(phosphino)aryl",),
        description="PCP pincer",
    ),
    "NCN": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(amino)aryl",),
        description="NCN pincer",
    ),
    "SCS": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(thio)aryl",),
        description="SCS pincer",
    ),
    "CNC": LigandInfo(
        smiles="",
        denticity=3,
        charge=0,
        aliases=("bis(NHC)pyridine",),
        description="CNC pincer (bis-NHC)",
    ),
    # PyBOX
    "PyBOX": LigandInfo(
        smiles="C1N=C(OC1)C1=NC(=CC=C1)C=1OCC(N1)",
        denticity=3,
        charge=0,
        aliases=("pyridinebisoxazoline", "pybox"),
        description="2,6-Bis(oxazolinyl)pyridine",
    ),
    "iPr-PyBOX": LigandInfo(
        smiles="C(C)(C)C1N=C(OC1)C1=NC(=CC=C1)C=1OCC(N1)C(C)C",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-isopropyl-2-oxazolinyl)pyridine",
    ),
    "Ph-PyBOX": LigandInfo(
        smiles="C1(=CC=CC=C1)C1N=C(OC1)C1=NC(=CC=C1)C=1OCC(N1)C1=CC=CC=C1",
        denticity=3,
        charge=0,
        aliases=(),
        description="2,6-Bis(4-phenyl-2-oxazolinyl)pyridine",
    ),
    # Triphosphines
    "triphos": LigandInfo(
        smiles="CC(CP(C1=CC=CC=C1)C2=CC=CC=C2)(CP(C3=CC=CC=C3)C4=CC=CC=C4)CP(C5=CC=CC=C5)C6=CC=CC=C6",
        denticity=3,
        charge=0,
        aliases=("1,1,1-tris(diphenylphosphinomethyl)ethane", "MeC(CH2PPh2)3"),
        description="Triphos",
    ),
    # Triamines
    "dien": LigandInfo(
        smiles="NCCNCCN",
        denticity=3,
        charge=0,
        aliases=("diethylenetriamine",),
        description="Diethylenetriamine",
    ),
    "tacn": LigandInfo(
        smiles="C1CNCCNCCN1",
        denticity=3,
        charge=0,
        aliases=("1,4,7-triazacyclononane",),
        description="1,4,7-Triazacyclononane",
    ),
    "Me3tacn": LigandInfo(
        smiles="CN1CCN(CCN(CC1)C)C",
        denticity=3,
        charge=0,
        aliases=("1,4,7-trimethyl-1,4,7-triazacyclononane",),
        description="1,4,7-Trimethyl-1,4,7-triazacyclononane",
    ),
    # Scorpionate-type (neutral)
    "Tp": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(pyrazolyl)borate", "trispyrazolylborate"),
        description="Hydrotris(pyrazolyl)borate",
    ),
    "Tp*": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-dimethylpyrazolyl)borate",),
        description="Hydrotris(3,5-dimethylpyrazolyl)borate",
    ),
    "TpiPr2": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("hydrotris(3,5-diisopropylpyrazolyl)borate",),
        description="Hydrotris(3,5-diisopropylpyrazolyl)borate",
    ),
    # Other
    "bpa": LigandInfo(
        smiles="N1=C(C=CC=C1)CNCC1=NC=CC=C1",
        denticity=3,
        charge=0,
        aliases=("bis(2-pyridylmethyl)amine",),
        description="Bis(2-pyridylmethyl)amine",
    ),
    "bpea": LigandInfo(
        smiles="N1=C(C=CC=C1)CN(CC1=NC=CC=C1)CC",
        denticity=3,
        charge=0,
        aliases=("N,N-bis(2-pyridylmethyl)ethylamine",),
        description="N,N-Bis(2-pyridylmethyl)ethylamine",
    ),
    "dap": LigandInfo(
        smiles="CC(=O)c1cccc(n1)C(=O)C",
        denticity=3,
        charge=0,
        aliases=("2,6-Diacetylpyridine",),
        description="2,6-Diacetylpyridine",
    ),
    # -------------------------------------------------------------------------
    # Tridentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Pincer anionic
    "PCP-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PCP pincer (cyclometalated)",
    ),
    "NCN-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="NCN pincer (cyclometalated)",
    ),
    "PNP-": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=(),
        description="PNP pincer (amido form)",
    ),
    # Bis(imino)pyridine
    "PDI": LigandInfo(
        smiles="",
        denticity=3,
        charge=-1,
        aliases=("bis(imino)pyridine", "pyridinediimine"),
        description="Pyridine-2,6-diimine (reduced form)",
    ),
    # Corroles
    "corrole": LigandInfo(
        smiles="",
        denticity=4,
        charge=-3,
        aliases=(),
        description="Corrole (tridentate in some counting)",
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Macrocycles
    "cyclam": LigandInfo(
        smiles="N1CCNCCCNCCNCCC1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetraazacyclotetradecane",),
        description="Cyclam",
    ),
    "cyclen": LigandInfo(
        smiles="N1CCNCCNCCNCC1",
        denticity=4,
        charge=0,
        aliases=("1,4,7,10-tetraazacyclododecane",),
        description="Cyclen",
    ),
    "Me4cyclam": LigandInfo(
        smiles="CN1CCCN(C)CCN(C)CCCN(C)CC1",
        denticity=4,
        charge=0,
        aliases=("1,4,8,11-tetramethyl-1,4,8,11-tetraazacyclotetradecane",),
        description="Tetramethylcyclam",
    ),
    # Linear tetradentate
    "trien": LigandInfo(
        smiles="NCCNCCNCCN",
        denticity=4,
        charge=0,
        aliases=("triethylenetetramine",),
        description="Triethylenetetramine",
    ),
    # Tetraphosphines
    "PP3": LigandInfo(
        smiles="",
        denticity=4,
        charge=0,
        aliases=("tris(2-(diphenylphosphino)ethyl)phosphine",),
        description="PP3",
    ),
    # -------------------------------------------------------------------------
    # Tetradentate Anionic Ligands
    # -------------------------------------------------------------------------
    # Porphyrins (crucial ligand class)
    "TPP": LigandInfo(
        smiles="C1(=CC=CC=C1)C=1C2=CC=C(N2)C(=C2C=CC(C(=C3C=CC(=C(C=4C=CC1N4)C4=CC=CC=C4)N3)C3=CC=CC=C3)=N2)C2=CC=CC=C2",
        denticity=4,
        charge=-2,
        aliases=("tetraphenylporphyrin", "5,10,15,20-tetraphenylporphyrin"),
        description="Tetraphenylporphyrin",
    ),
    "OEP": LigandInfo(
        smiles="C(C)C1=C2NC(=C1CC)C=C1C(=C(C(=N1)C=C1C(=C(C(N1)=CC=1C(=C(C(N1)=C2)CC)CC)CC)CC)CC)CC",
        denticity=4,
        charge=-2,
        aliases=("octaethylporphyrin", "2,3,7,8,12,13,17,18-octaethylporphyrin"),
        description="Octaethylporphyrin",
    ),
    "TMP": LigandInfo(
        smiles="C1(=C(C(=CC(=C1)C)C)C1=C2C=CC(C(=C3C=CC(=C(C=4C=CC(=C(C5=CC=C1N5)C5=C(C=C(C=C5C)C)C)N4)C4=C(C=C(C=C4C)C)C)N3)C3=C(C=C(C=C3C)C)C)=N2)C",
        denticity=4,
        charge=-2,
        aliases=("tetramesitylporphyrin",),
        description="Tetramesitylporphyrin",
    ),
    "por": LigandInfo(
        smiles="C12=CC=C(N1)C=C1C=CC(=N1)C=C1C=CC(N1)=CC=1C=CC(N1)=C2",
        denticity=4,
        charge=-2,
        aliases=("porphyrin", "porphyrinato"),
        description="Porphyrin (generic)",
    ),
    "TPFPP": LigandInfo(
        smiles="FC1=C(C(=C(C(=C1C1=C2C=CC(C(=C3C=CC(=C(C=4C=CC(=C(C5=CC=C1N5)C5=C(C(=C(C(=C5F)F)F)F)F)N4)C4=C(C(=C(C(=C4F)F)F)F)F)N3)C3=C(C(=C(C(=C3F)F)F)F)F)=N2)F)F)F)F",
        denticity=4,
        charge=-2,
        aliases=("tetrakis(pentafluorophenyl)porphyrin",),
        description="Tetrakis(pentafluorophenyl)porphyrin",
    ),
    # Phthalocyanines
    "Pc": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("phthalocyanine", "phthalocyaninato"),
        description="Phthalocyanine",
    ),
    # Salen-type
    "salen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-ethylenebis(salicylideneiminato)",),
        description="Salen",
    ),
    "salphen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-phenylenebis(salicylideneiminato)",),
        description="Salphen",
    ),
    "salophen": LigandInfo(
        smiles="", denticity=4, charge=-2, aliases=(), description="Salophen"
    ),
    "salcn": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("N,N'-cyclohexanebis(salicylideneiminato)",),
        description="Salen-cyclohexanediamine",
    ),
    "Jacobsen": LigandInfo(
        smiles="",
        denticity=4,
        charge=-2,
        aliases=("Jacobsen's salen",),
        description="Jacobsen's salen ligand",
    ),
    # -------------------------------------------------------------------------
    # Penta/Hexadentate Neutral Ligands
    # -------------------------------------------------------------------------
    # Pentadentate
    "Me6tren": LigandInfo(
        smiles="CN(C)CCN(CCN(C)C)CCN(C)C",
        denticity=4,
        charge=0,
        aliases=("tris(2-dimethylaminoethyl)amine",),
        description="Tris(2-dimethylaminoethyl)amine",
    ),
    "tpa": LigandInfo(
        smiles="c1ccnc(c1)CN(Cc2ccccn2)Cc3ccccn3",
        denticity=4,
        charge=0,
        aliases=("tris(2-pyridylmethyl)amine",),
        description="Tris(2-pyridylmethyl)amine",
    ),
    # Hexadentate
    "EDTA": LigandInfo(
        smiles="OC(=O)CN(CCN(CC(O)=O)CC(O)=O)CC(O)=O",
        denticity=6,
        charge=-4,
        aliases=("ethylenediaminetetraacetate", "edta"),
        description="Ethylenediaminetetraacetate",
    ),
    "DTPA": LigandInfo(
        smiles="C(CN(CC(=O)O)CC(=O)O)N(CCN(CC(=O)O)CC(=O)O)CC(=O)O",
        denticity=8,
        charge=-5,
        aliases=("diethylenetriaminepentaacetate",),
        description="Diethylenetriaminepentaacetate",
    ),
    "tpen": LigandInfo(
        smiles="C1=CC=NC(=C1)CN(CCN(CC2=CC=CC=N2)CC3=CC=CC=N3)CC4=CC=CC=N4",
        denticity=6,
        charge=0,
        aliases=("N,N,N',N'-tetrakis(2-pyridylmethyl)ethylenediamine",),
        description="TPEN",
    ),
    # -------------------------------------------------------------------------
    # η-Bonded Ligands
    # -------------------------------------------------------------------------
    "Cp": LigandInfo(
        smiles="[cH-]1cccc1",
        denticity=5,
        charge=-1,
        aliases=("cyclopentadienyl", "C5H5"),
        description="Cyclopentadienyl (η⁵)",
    ),
    "Cp*": LigandInfo(
        smiles="C[c-]1c(C)c(C)c(C)c1C",
        denticity=5,
        charge=-1,
        aliases=("pentamethylcyclopentadienyl", "C5Me5", "Cpstar"),
        description="Pentamethylcyclopentadienyl (η⁵)",
    ),
    # Arenes (η⁶)
    "benzene": LigandInfo(
        smiles="c1ccccc1",
        denticity=6,
        charge=0,
        aliases=("η6-benzene", "C6H6"),
        description="Benzene (η⁶)",
    ),
    "p-cymene": LigandInfo(
        smiles="c1cc(ccc1C(C)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-p-cymene", "4-isopropyltoluene", "p-cym", "p-Cym"),
        description="p-Cymene (η⁶)",
    ),
    "mesitylene": LigandInfo(
        smiles="Cc1cc(cc(c1)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-mesitylene", "1,3,5-trimethylbenzene"),
        description="Mesitylene (η⁶)",
    ),
    "hexamethylbenzene": LigandInfo(
        smiles="c1(c(c(c(c(c1C)C)C)C)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-C6Me6", "HMB"),
        description="Hexamethylbenzene (η⁶)",
    ),
    "toluene": LigandInfo(
        smiles="Cc1ccccc1",
        denticity=6,
        charge=0,
        aliases=("η6-toluene",),
        description="Toluene (η⁶)",
    ),
    "C6H3Me3": LigandInfo(
        smiles="Cc1cc(cc(c1)C)C",
        denticity=6,
        charge=0,
        aliases=("η6-trimethylbenzene",),
        description="Trimethylbenzene (η⁶)",
    ),
    # Allyl (η³)
    "η3-allyl": LigandInfo(
        smiles="[CH-]C=C",
        denticity=3,
        charge=-1,
        aliases=("η3-allyl", "η3-C3H5"),
        description="Allyl (η³)",
    ),
    "methallyl": LigandInfo(
        smiles="[CH-]C(C)=C",
        denticity=3,
        charge=-1,
        aliases=("η3-2-methylallyl", "η3-methallyl"),
        description="Methallyl (η³)",
    ),
    "crotyl": LigandInfo(
        smiles="[CH-]C=CC",
        denticity=3,
        charge=-1,
        aliases=("η3-crotyl", "η3-1-methylallyl"),
        description="Crotyl (η³)",
    ),
    "cinnamyl": LigandInfo(
        smiles="[CH-]C=CC1=CC=CC=C1",
        denticity=3,
        charge=-1,
        aliases=("η3-cinnamyl", "η3-3-phenylallyl"),
        description="Cinnamyl (η³)",
    ),
    # Tropylium and related
    "C7H7": LigandInfo(
        smiles="[C-]1=CC=CC=CC1",
        denticity=7,
        charge=1,
        aliases=("tropylium", "cycloheptatrienyl"),
        description="Tropylium (η⁷)",
    ),
    "cht": LigandInfo(
        smiles="C1=C\C/C=C\C=C1",
        denticity=6,
        charge=0,
        aliases=("cycloheptatriene", "η6-C7H8"),
        description="Cycloheptatriene (η⁶)",
    ),
    # Indenyl
    "Ind": LigandInfo(
        smiles="[C-]1C=CC2=CC=CC=C12",
        denticity=5,
        charge=-1,
        aliases=("indenyl", "η5-indenyl"),
        description="Indenyl (η⁵)",
    ),
    # Fluorenyl
    "Flu": LigandInfo(
        smiles="[C-]1=CC=CC=2C3=CC=CC=C3CC12",
        denticity=5,
        charge=-1,
        aliases=("fluorenyl", "η5-fluorenyl"),
        description="Fluorenyl (η⁵)",
    ),
    # Cyclooctatetraene
    "cot": LigandInfo(
        smiles="",
        denticity=8,
        charge=-2,
        aliases=("cyclooctatetraene", "η8-C8H8", "COT"),
        description="Cyclooctatetraene (η⁸)",
    ),
    # Butadiene
    "bd": LigandInfo(
        smiles="C=CC=C",
        denticity=4,
        charge=0,
        aliases=("butadiene", "η4-butadiene"),
        description="1,3-Butadiene (η⁴)",
    ),
    "isoprene": LigandInfo(
        smiles="CC(=C)C=C",
        denticity=4,
        charge=0,
        aliases=("η4-isoprene", "2-methylbutadiene"),
        description="Isoprene (η⁴)",
    ),
}


COUNTER_ION_DATABASE: Dict[str, LigandInfo] = {
    "PF6": LigandInfo(
        smiles="F[P-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluorophosphate",),
        description="Hexafluorophosphate",
    ),
    "BF4": LigandInfo(
        smiles="F[B-](F)(F)F",
        charge=-1,
        aliases=("tetrafluoroborate",),
        description="Tetrafluoroborate",
    ),
    "OTf": LigandInfo(
        smiles="[O-]S(=O)(=O)C(F)(F)F",
        charge=-1,
        aliases=("triflate", "trifluoromethanesulfonate", "CF3SO3"),
        description="Triflate",
    ),
    "ClO4": LigandInfo(
        smiles="[O-]Cl(=O)(=O)=O",
        charge=-1,
        aliases=("perchlorate",),
        description="Perchlorate",
    ),
    "SbF6": LigandInfo(
        smiles="F[Sb-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluoroantimonate",),
        description="Hexafluoroantimonate",
    ),
    "BArF": LigandInfo(
        smiles="FC(F)(F)c1cc([B-](c2cc(C(F)(F)F)cc(C(F)(F)F)c2)"
        "(c2cc(C(F)(F)F)cc(C(F)(F)F)c2)c2cc(C(F)(F)F)cc(C(F)(F)F)c2)"
        "cc(C(F)(F)F)c1",
        charge=-1,
        aliases=("BArF24", "tetrakis(3,5-bis(trifluoromethyl)phenyl)borate"),
        description="Tetrakis(3,5-bis(trifluoromethyl)phenyl)borate",
    ),
    "BAr4": LigandInfo(
        smiles="c1ccc([B-](c2ccccc2)(c2ccccc2)c2ccccc2)cc1",
        charge=-1,
        aliases=("tetraphenylborate", "BPh4"),
        description="Tetraphenylborate",
    ),
    "NO3": LigandInfo(
        smiles="[O-][N+](=O)[O-]",
        charge=-1,
        aliases=("nitrate",),
        description="Nitrate",
    ),
    "Cl": LigandInfo(
        smiles="[Cl-]", charge=-1, aliases=("chloride",), description="Chloride"
    ),
    "Br": LigandInfo(
        smiles="[Br-]", charge=-1, aliases=("bromide",), description="Bromide"
    ),
    "I": LigandInfo(
        smiles="[I-]", charge=-1, aliases=("iodide",), description="Iodide"
    ),
    "BPh4": LigandInfo(
        smiles="[B-](c1ccccc1)(c2ccccc2)(c3ccccc3)c4ccccc4",
        charge=-1,
        aliases=("tetraphenylborate",),
        description="Tetraphenylborate",
    ),
    "Al(OC(CF3)3)4": LigandInfo(
        smiles="C(C(F)(F)F)(C(F)(F)F)(C(F)(F)F)O[Al-](OC(C(F)(F)F)(C(F)(F)F)C(F)(F)F)(OC(C(F)(F)F)(C(F)(F)F)C(F)(F)F)OC(C(F)(F)F)(C(F)(F)F)C(F)(F)F",
        charge=-1,
        aliases=("perfluoro-tert-butoxide aluminate",),
        description="Perfluoro-tert-butoxide aluminate",
    ),
    "BAr4F": LigandInfo(
        smiles="FC(C=1C=C(C=C(C1)C(F)(F)F)[B-](C1=CC(=CC(=C1)C(F)(F)F)C(F)(F)F)(C1=CC(=CC(=C1)C(F)(F)F)C(F)(F)F)C1=CC(=CC(=C1)C(F)(F)F)C(F)(F)F)(F)F",
        charge=-1,
        aliases=("tetrakis(3,5-bis(trifluoromethyl)phenyl)borate", "BArF"),
        description="BArF",
    ),
    "NTf2": LigandInfo(
        smiles="C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",
        charge=-1,
        aliases=("bis(trifluoromethylsulfonyl)imide", "TFSI", "bistriflimide"),
        description="Bis(trifluoromethylsulfonyl)imide",
    ),
    "AsF6": LigandInfo(
        smiles="F[As-](F)(F)(F)(F)F",
        charge=-1,
        aliases=("hexafluoroarsenate",),
        description="Hexafluoroarsenate",
    ),
    "B(C6F5)4": LigandInfo(
        smiles="Fc1c(c(F)c(F)c(F)c1F)[B-](c2c(F)c(F)c(F)c(F)c2F)(c3c(F)c(F)c(F)c(F)c3F)c4c(F)c(F)c(F)c(F)c4F",
        charge=-1,
        aliases=("tetrakis(pentafluorophenyl)borate",),
        description="Tetrakis(pentafluorophenyl)borate",
    ),
    "HSO4": LigandInfo(
        smiles="OS(=O)(=O)[O-]",
        charge=-1,
        aliases=("hydrogensulfate", "bisulfate"),
        description="Hydrogen sulfate",
    ),
    "CF3CO2": LigandInfo(
        smiles="FC(F)(F)C(=O)[O-]",
        charge=-1,
        aliases=("trifluoroacetate", "TFA"),
        description="Trifluoroacetate",
    ),
    "MeSO3": LigandInfo(
        smiles="O=S(=O)([O-])C",
        charge=-1,
        aliases=("mesylate", "methanesulfonate"),
        description="Mesylate",
    ),
    "TsO": LigandInfo(
        smiles="CC1=CC=C(C=C1)S(=O)(=O)[O-]",
        charge=-1,
        aliases=("tosylate", "4-toluenesulfonate", "OTs"),
        description="Tosylate",
    ),
    "ReO4": LigandInfo(
        smiles="O=[Re](O)(O)(O)(O)",
        charge=-1,
        aliases=("perrhenate",),
        description="Perrhenate",
    ),
    "IO4": LigandInfo(
        smiles="O[I](O)(O)(O)O",
        charge=-1,
        aliases=("periodate",),
        description="Periodate",
    ),
    "BH4": LigandInfo(
        smiles="[BH4-]",
        charge=-1,
        aliases=("borohydride", "tetrahydroborate"),
        description="Borohydride (as counter ion)",
    ),
    "AlH4": LigandInfo(
        smiles="[AlH4-]",
        charge=-1,
        aliases=("aluminate", "tetrahydroaluminate"),
        description="Tetrahydroaluminate",
    ),
    "AlCl4": LigandInfo(
        smiles="Cl[Al-](Cl)(Cl)Cl",
        charge=-1,
        aliases=("tetrachloroaluminate",),
        description="Tetrachloroaluminate",
    ),
    "FeCl4": LigandInfo(
        smiles="Cl[Fe-](Cl)(Cl)Cl",
        charge=-1,
        aliases=("tetrachloroferrate",),
        description="Tetrachloroferrate",
    ),
    "CuCl2": LigandInfo(
        smiles="Cl[Cu-]Cl",
        charge=-1,
        aliases=("dichlorocuprate",),
        description="Dichlorocuprate",
    ),
    "ZnCl3": LigandInfo(
        smiles="Cl[Zn-](Cl)Cl",
        charge=-1,
        aliases=("trichlorozincate",),
        description="Trichlorozincate",
    ),
    "GaCl4": LigandInfo(
        smiles="Cl[Ga-](Cl)(Cl)Cl",
        charge=-1,
        aliases=("tetrachlorogallate",),
        description="Tetrachlorogallate",
    ),
    # Cationic counterions (for anionic complexes)
    "Na": LigandInfo(
        smiles="[Na+]", charge=1, aliases=("sodium",), description="Sodium"
    ),
    "K": LigandInfo(
        smiles="[K+]", charge=1, aliases=("potassium",), description="Potassium"
    ),
    "Li": LigandInfo(
        smiles="[Li+]", charge=1, aliases=("lithium",), description="Lithium"
    ),
    "Cs": LigandInfo(
        smiles="[Cs+]", charge=1, aliases=("cesium", "caesium"), description="Cesium"
    ),
    "NBu4": LigandInfo(
        smiles="CCCC[N+](CCCC)(CCCC)CCCC",
        charge=1,
        aliases=("tetrabutylammonium", "TBA", "nBu4N"),
        description="Tetrabutylammonium",
    ),
    "NEt4": LigandInfo(
        smiles="CC[N+](CC)(CC)CC",
        charge=1,
        aliases=("tetraethylammonium", "TEA", "Et4N"),
        description="Tetraethylammonium",
    ),
    "NMe4": LigandInfo(
        smiles="C[N+](C)(C)C",
        charge=1,
        aliases=("tetramethylammonium", "TMA", "Me4N"),
        description="Tetramethylammonium",
    ),
    "PPh4": LigandInfo(
        smiles="c1c(cccc1)[P+](c2ccccc2)(c3ccccc3)c4ccccc4",
        charge=1,
        aliases=("tetraphenylphosphonium", "Ph4P"),
        description="Tetraphenylphosphonium",
    ),
    "PPN": LigandInfo(
        smiles="C1=CC=C(C=C1)P(=N[P+](C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)(C5=CC=CC=C5)C6=CC=CC=C6",
        charge=1,
        aliases=("bis(triphenylphosphine)iminium", "Ph3P=N=PPh3"),
        description="Bis(triphenylphosphine)iminium",
    ),
    "Cp2Fe": LigandInfo(
        smiles="ferrocenium",
        charge=1,
        aliases=("ferrocenium", "Fc+"),
        description="Ferrocenium",
    ),
    "Cp2Co": LigandInfo(
        smiles="C1C=CC=C1.[CH-]1C=CC=C1.[Co+2]",
        charge=1,
        aliases=("cobaltocenium",),
        description="Cobaltocenium",
    ),
    "H3O": LigandInfo(
        smiles="", charge=1, aliases=("hydronium", "oxonium"), description="Hydronium"
    ),
    "NH4": LigandInfo(
        smiles="[NH4+]", charge=1, aliases=("ammonium",), description="Ammonium"
    ),
    "pyH": LigandInfo(
        smiles="[NH+]1=CC=CC=C1",
        charge=1,
        aliases=("pyridinium",),
        description="Pyridinium",
    ),
    "DMAH": LigandInfo(
        smiles="C[NH+](C1=CC=CC=C1)C",
        charge=1,
        aliases=("dimethylanilinium",),
        description="Dimethylanilinium",
    ),
}


METAL_DATABASE: Dict[str, MetalInfo] = {
    # Group 6
    "Cr": MetalInfo("Cr", "Chromium", (0, 2, 3, 6), 24),
    "Mo": MetalInfo("Mo", "Molybdenum", (0, 2, 4, 6), 42),
    "W": MetalInfo("W", "Tungsten", (0, 2, 4, 6), 74),
    # Group 7
    "Mn": MetalInfo("Mn", "Manganese", (0, 2, 3, 4, 7), 25),
    "Re": MetalInfo("Re", "Rhenium", (0, 1, 3, 5, 7), 75),
    # Group 8
    "Fe": MetalInfo("Fe", "Iron", (0, 2, 3), 26),
    "Ru": MetalInfo("Ru", "Ruthenium", (0, 2, 3, 4), 44),
    "Os": MetalInfo("Os", "Osmium", (0, 2, 3, 4, 6, 8), 76),
    # Group 9
    "Co": MetalInfo("Co", "Cobalt", (0, 2, 3), 27),
    "Rh": MetalInfo("Rh", "Rhodium", (0, 1, 2, 3), 45),
    "Ir": MetalInfo("Ir", "Iridium", (0, 1, 3, 4), 77),
    # Group 10
    "Ni": MetalInfo("Ni", "Nickel", (0, 2), 28),
    "Pd": MetalInfo("Pd", "Palladium", (0, 2, 4), 46),
    "Pt": MetalInfo("Pt", "Platinum", (0, 2, 4), 78),
    # Group 11
    "Cu": MetalInfo("Cu", "Copper", (1, 2), 29),
    "Ag": MetalInfo("Ag", "Silver", (1,), 47),
    "Au": MetalInfo("Au", "Gold", (1, 3), 79),
    # Group 12
    "Zn": MetalInfo("Zn", "Zinc", (2,), 30),
}
