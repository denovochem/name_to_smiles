import warnings
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

from rdkit import RDLogger

from placeholder_name.name_manipulation.split_names import (
    get_delimiter_split_dict,
    resolve_delimiter_split_dict,
)
from placeholder_name.resolvers.chemspipy_resolver import name_to_smiles_chemspipy
from placeholder_name.resolvers.cirpy_resolver import name_to_smiles_cirpy
from placeholder_name.resolvers.manual_resolver import name_to_smiles_manual
from placeholder_name.resolvers.opsin_resolver import name_to_smiles_opsin
from placeholder_name.resolvers.peptide_resolver import name_to_smiles_peptide
from placeholder_name.resolvers.pubchem_resolver import name_to_smiles_pubchem
from placeholder_name.resolvers.structural_formula_resolver import (
    name_to_smiles_structural_formula,
)
from placeholder_name.smiles_selector import SMILESSelector
from placeholder_name.utils.chem_utils import canonicalize_smiles
from placeholder_name.utils.logging_config import configure_logging, logger
from placeholder_name.utils.string_utils import clean_strings

# Configure loguru logging
configure_logging(level="WARNING")

# Ignore all RuntimeWarnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Disable rdkit warnings
RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]


class ChemicalNameResolver(ABC):
    """
    Abstract base class for chemical name-to-SMILES resolvers.

    Subclasses must implement the `name_to_smiles` method.
    """

    def __init__(self, resolver_type: str, resolver_name: str, resolver_weight: float):
        if not isinstance(resolver_type, str):
            raise TypeError("Invalid input: resolver_type must be a string.")
        self._resolver_type: str = resolver_type
        if not isinstance(resolver_name, str):
            raise TypeError("Invalid input: resolver_name must be a string.")
        self._resolver_name: str = resolver_name
        if not isinstance(resolver_weight, (int, float)):
            raise TypeError(
                "Invalid input: resolver_weight must be a number between 0-1000."
            )
        if resolver_weight < 0 or resolver_weight > 1000:
            raise ValueError(
                "Invalid input: resolver_weight must be a number between 0-1000."
            )
        self._resolver_weight: float = float(resolver_weight)

    @property
    def resolver_name(self) -> str:
        """Return resolver_name."""
        return self._resolver_name

    @property
    def resolver_weight(self) -> float:
        """Return resolver_weight."""
        return self._resolver_weight

    @abstractmethod
    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES strings.

        Args:
            compound_name_list: List of chemical names.

        Returns:
            Tuple of:
                - Dict mapping successful names to SMILES.
                - Dict mapping failed names to error messages.
        """
        pass


class OpsinNameResolver(ChemicalNameResolver):
    """
    Resolver using OPSIN via py2opsin.
    """

    def __init__(
        self,
        resolver_name: str,
        resolver_weight: float = 3,
        allow_acid: bool = False,
        allow_radicals: bool = True,
        allow_bad_stereo: bool = False,
        wildcard_radicals: bool = False,
        jar_fpath: str = "opsin-cli.jar",
    ):
        super().__init__("opsin", resolver_name, resolver_weight)
        self._allow_acid = allow_acid
        self._allow_radicals = allow_radicals
        self._allow_bad_stereo = allow_bad_stereo
        self._wildcard_radicals = wildcard_radicals

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        resolved_names, failure_message_dict = name_to_smiles_opsin(
            compound_name_list,
            allow_acid=self._allow_acid,
            allow_radicals=self._allow_radicals,
            allow_bad_stereo=self._allow_bad_stereo,
            wildcard_radicals=self._wildcard_radicals,
        )
        return resolved_names, failure_message_dict


class PubChemNameResolver(ChemicalNameResolver):
    """
    Resolver using PubChem via PubChemPy.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 2):
        super().__init__("pubchem", resolver_name, resolver_weight)

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using pubchem.
        """
        resolved_names = name_to_smiles_pubchem(compound_name_list)
        return resolved_names, {}


class CIRpyNameResolver(ChemicalNameResolver):
    """
    Resolver using Chemical Identity Resolver via CIRPy.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 1):
        super().__init__("cirpy", resolver_name, resolver_weight)

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using cirpy.
        """
        resolved_names = name_to_smiles_cirpy(compound_name_list)
        return resolved_names, {}


class ChemSpiPyResolver(ChemicalNameResolver):
    """
    Resolver using chemspipy.
    """

    def __init__(
        self, resolver_name: str, chemspider_api_key: str, resolver_weight: float = 3
    ):
        super().__init__("cirpy", resolver_name, resolver_weight)
        if chemspider_api_key:
            if not isinstance(chemspider_api_key, str):
                raise TypeError("Invalid input: chemspider_api_key must be a string.")
        self._chemspider_api_key = chemspider_api_key

    def name_to_smiles(
        self,
        compound_name_list: List[str],
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using ChemSpiPy.
        """
        resolved_names = name_to_smiles_chemspipy(
            compound_name_list, self._chemspider_api_key
        )
        return resolved_names, {}


class ManualNameResolver(ChemicalNameResolver):
    """
    Resolver using manually curated names and corresponding SMILES.
    """

    def __init__(
        self,
        resolver_name: str,
        provided_name_dict: dict | None = None,
        resolver_weight: float = 10,
    ):
        super().__init__("manual", resolver_name, resolver_weight)
        if provided_name_dict:
            if not isinstance(provided_name_dict, dict):
                raise TypeError(
                    "Invalid input: provided_name_dict must be a dictionary."
                )
            for k, v in provided_name_dict.items():
                if not isinstance(k, str) or not isinstance(v, str):
                    raise ValueError(
                        "Invalid input: keys and values in provided_name_dict must be strings."
                    )

        self._provided_name_dict = provided_name_dict

    def name_to_smiles(
        self,
        compound_name_list: List[str],
        provided_name_dict: Dict[str, str] | None = None,
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using manual name database.
        """
        if provided_name_dict is None:
            provided_name_dict = self._provided_name_dict
        resolved_names = name_to_smiles_manual(compound_name_list, provided_name_dict)
        return resolved_names, {}


class PeptideNameResolver(ChemicalNameResolver):
    """
    Resolver using peptide shorthand-to-IUPAC-like name, then resolved to SMILES
    with OPSIN via py2opsin.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 3):
        super().__init__("peptide", resolver_name, resolver_weight)

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using peptide name converter and OPSIN.
        """
        resolved_names, failure_message_dict = name_to_smiles_peptide(
            compound_name_list
        )
        return resolved_names, failure_message_dict


class StructuralFormulaNameResolver(ChemicalNameResolver):
    """
    Resolver using structural chemical formula (e.g. CH3CH2CH2COOH).
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 2):
        super().__init__("structural_formula", resolver_name, resolver_weight)

    def name_to_smiles(
        self, compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using structural formula converter.
        """
        resolved_names = name_to_smiles_structural_formula(compound_name_list)
        return resolved_names, {}


def clean_strings_and_return_mapping(
    compounds_list: List[str],
) -> Tuple[List[str], Dict[str, str]]:
    """
    Clean a list of strings by replacing certain characters and return a mapping from the original strings to their cleaned versions.

    Args:
        compounds_list (List[str]): A list of strings to clean.

    Returns:
        Tuple[List[str], Dict[str, str]]: A tuple containing the list of cleaned strings and a dictionary mapping the original strings to their cleaned versions.
    """
    cleaned_compounds_list = [clean_strings(ele) for ele in compounds_list]
    cleaned_compounds_dict = {
        k: v for k, v in zip(compounds_list, cleaned_compounds_list)
    }
    return cleaned_compounds_list, cleaned_compounds_dict


def get_resolvers_weight_dict(
    resolvers_list: List[ChemicalNameResolver],
) -> Dict[str, float]:
    """
    Get a dictionary mapping resolver names to their weights.

    Args:
        resolvers_list (List[ChemicalNameResolver]): A list of resolvers.

    Returns:
        Dict[str, float]: A dictionary mapping resolver names to their weights.
    """
    resolvers_weight_dict = {}
    for resolver in resolvers_list:
        resolvers_weight_dict[resolver.resolver_name] = resolver.resolver_weight
    return resolvers_weight_dict


def split_compounds_on_delimiters_and_return_mapping(
    compounds_list: List[str],
) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Split a list of compound names into individual parts based on the DELIMITERS list and store the split parts in a dictionary.

    Args:
        compounds_list (List[str]): The list of compound names to split.

    Returns:
        Tuple[List[str], Dict[str, List[str]]]: A tuple containing the list of split compound names and a dictionary mapping the original compound names to their split parts.
    """
    delimiter_split_dict: Dict[str, List[str]] = {}
    compounds_split_parts_list = []
    for compound in compounds_list:
        delimiter_split_dict, compound_split_parts = get_delimiter_split_dict(
            compound, delimiter_split_dict
        )
        compounds_split_parts_list.extend(compound_split_parts)
    compounds_and_split_parts_list = compounds_list + compounds_split_parts_list
    compounds_and_split_parts_list = list(set(compounds_and_split_parts_list))
    return compounds_and_split_parts_list, delimiter_split_dict


def resolve_compounds_using_resolvers(
    compounds_list: List[str],
    resolvers_list: List[ChemicalNameResolver],
    batch_size: int,
) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Resolve a list of compound names using a list of resolvers.

    Args:
        compounds_list (List[str]): A list of compound names to resolve.
        resolvers_list (List[ChemicalNameResolver]): A list of resolvers to use.
        batch_size (int): The number of compound names to process in each batch.

    Returns:
        Dict[str, Dict[str, Dict[str, str]]]: A dictionary mapping each resolver name to its output dictionary, which maps each compound name to its resolved SMILES string and error message.
    """
    resolvers_out_dict = {}
    for resolver in resolvers_list:
        for i in range(0, len(compounds_list), batch_size):
            chunk = compounds_list[i : i + batch_size]
            out, info_messages = resolver.name_to_smiles(chunk)
        resolvers_out_dict[resolver.resolver_name] = {
            "out": out,
            "info_messages": info_messages,
        }
    return resolvers_out_dict


def assemble_compounds_resolution_dict(
    compounds: List[str],
    resolvers_out_dict: Dict[str, Dict[str, Dict[str, str]]],
    cleaned_compounds_dict: Dict[str, str],
) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
    """
    Assemble a dictionary mapping each compound to its resolved SMILES string and resolvers.

    Args:
        compounds (List[str]): A list of compound names to resolve.
        resolvers_out_dict (Dict[str, Dict[str, Dict[str, str]]]): A dictionary mapping each resolver name to its output dictionary, which maps each compound name to its resolved SMILES string and error message.
        cleaned_compounds_dict (Dict[str, str]): A dictionary mapping each compound name to its cleaned version.

    Returns:
        Dict[str, Dict[str, Dict[str, List[str]]]]: A dictionary mapping each compound to its resolved SMILES string and resolvers.
    """

    def _assemble_compound_resolution_dict(
        compound, resolvers_out_dict, cleaned_compounds_dict
    ):
        """ """
        compound_cleaned = cleaned_compounds_dict.get(compound, "")
        compounds_out_dict[compound] = {
            "SMILES": "",
            "SMILES_source": [],
            "SMILES_dict": {},
            "info_messages": {},
        }

        for resolver, resolution_dict in resolvers_out_dict.items():
            compound_smiles = resolution_dict["out"].get(compound_cleaned, "")
            if resolution_dict["info_messages"].get(compound_cleaned, ""):
                compound_resolver_info_message = resolution_dict["info_messages"].get(
                    compound_cleaned, ""
                )
                compounds_out_dict[compound]["info_messages"][resolver] = (
                    compound_resolver_info_message
                )
            canonical_compound_smiles = canonicalize_smiles(compound_smiles)
            if not canonical_compound_smiles:
                continue
            if (
                canonical_compound_smiles
                not in compounds_out_dict[compound]["SMILES_dict"]
            ):
                compounds_out_dict[compound]["SMILES_dict"][
                    canonical_compound_smiles
                ] = [resolver]
            else:
                resolvers_list = compounds_out_dict[compound]["SMILES_dict"][
                    canonical_compound_smiles
                ]
                resolvers_list.append(resolver)
                compounds_out_dict[compound]["SMILES_dict"][
                    canonical_compound_smiles
                ] = resolvers_list

        return compounds_out_dict

    compounds_out_dict = {}
    for compound in compounds:
        compounds_out_dict = _assemble_compound_resolution_dict(
            compound, resolvers_out_dict, cleaned_compounds_dict
        )

    return compounds_out_dict


def assemble_split_compounds_resolution_dict(
    compounds_out_dict: Dict[str, Dict[str, Dict[str, List[str]]]],
    compounds: List[str],
    resolvers_out_dict: Dict[str, Dict[str, Dict[str, str]]],
    cleaned_compounds_dict: Dict[str, str],
    delimiter_split_dict: Dict[str, List[str]],
) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
    """
    Assemble a dictionary mapping each compound to its resolved SMILES string and resolvers.

    Args:
        compounds_out_dict (Dict[str, Dict[str, str]]): A dictionary mapping each compound to its resolved SMILES representations.
        compounds (List[str]): A list of compound names to resolve.
        resolvers_out_dict (Dict[str, Dict[str, str]]): A dictionary mapping each resolver name to its output dictionary, which maps each compound name to its resolved SMILES string and error message.
        cleaned_compounds_dict (Dict[str, str]): A dictionary mapping each compound name to its cleaned version.
        delimiter_split_dict (Dict[str, List[str]]): A dictionary mapping each compound name to its split parts.

    Returns:
        Dict[str, Dict[str, Dict[str, List[str]]]]: A dictionary mapping each compound to its resolved SMILES string and resolvers.
    """

    def _assemble_split_compound_resolution_dict_(
        compound,
        compounds_out_dict,
        resolvers_out_dict,
        cleaned_compounds_dict,
        delimiter_split_dict,
    ):
        """ """
        compound_cleaned = cleaned_compounds_dict.get(compound, "")
        if compound_cleaned not in delimiter_split_dict:
            return compounds_out_dict

        resolved_delimiter_split_dict = resolve_delimiter_split_dict(
            compound_cleaned, resolvers_out_dict, delimiter_split_dict
        )
        for (
            compound_smiles,
            resolvers_split_list,
        ) in resolved_delimiter_split_dict.items():
            canonical_compound_smiles = canonicalize_smiles(compound_smiles)
            if not canonical_compound_smiles:
                continue
            for resolver in resolvers_split_list:
                if (
                    canonical_compound_smiles
                    not in compounds_out_dict[compound]["SMILES_dict"]
                ):
                    compounds_out_dict[compound]["SMILES_dict"][
                        canonical_compound_smiles
                    ] = [resolver]
                else:
                    resolvers_list = compounds_out_dict[compound]["SMILES_dict"][
                        canonical_compound_smiles
                    ]
                    resolvers_list.append(resolver)
                    compounds_out_dict[compound]["SMILES_dict"][
                        canonical_compound_smiles
                    ] = resolvers_list

        return compounds_out_dict

    for compound in compounds:
        compounds_out_dict = _assemble_split_compound_resolution_dict_(
            compound,
            compounds_out_dict,
            resolvers_out_dict,
            cleaned_compounds_dict,
            delimiter_split_dict,
        )

    return compounds_out_dict


def select_smiles_with_criteria(
    compounds_out_dict,
    resolvers_weight_dict: Dict[str, float],
    resolvers_priority_order: List[str],
    smiles_selection_mode: str,
):
    """
    Select SMILES representation for each compound using the specified criteria.

    Args:
        compounds_out_dict: Dictionary of compound names to their resolved SMILES representations.
        resolvers_weight_dict (Dict[str, float]): Dictionary of resolver names to their weights.
        resolvers_priority_order (List[str]): List of priority order for resolvers
        smiles_selection_mode (str): The method to select the SMILES representation from multiple resolvers.

    Returns:
        Dictionary of compound names to their selected SMILES representations, SMILES source, and info messages.
    """
    selector = SMILESSelector(
        compounds_out_dict, resolvers_weight_dict, resolvers_priority_order
    )
    for k, v in compounds_out_dict.items():
        selected_smiles, selected_smiles_resolvers = selector.select_smiles(
            k, smiles_selection_mode
        )
        v["SMILES"] = selected_smiles
        v["SMILES_source"] = selected_smiles_resolvers

    return compounds_out_dict


def resolve_compounds_to_smiles(
    compounds_list: List[str],
    resolvers_list: List[ChemicalNameResolver] = [],
    smiles_selection_mode: str = "weighted",
    detailed_name_dict: bool = False,
    batch_size: int = 500,
    split_names_to_solve: bool = True,
) -> Dict[str, Dict[str, Dict[str, List[str]]]] | Dict[str, str]:
    """
    Resolve a list of compound names to their SMILES representations.

    Args:
        compounds_list (List[str]): A list of compound names.
        resolvers_list (List[ChemicalNameResolver], optional): A list of ChemicalNameResolver instances.
            Defaults to [].
        smiles_selection_mode (str, optional): The method to select the SMILES representation from multiple resolvers.
            Defaults to 'weighted'.
        detailed_name_dict (bool, optional): If True, returns a dictionary with detailed information about each compound.
            Defaults to False.
        batch_size (int, optional): The number of compounds to process in each batch. Defaults to 500.
        split_names_to_solve (bool, optional): Whether to split compound names on common delimiters to solve them as separate compounds.
            Can be used to solve otherwise unresolvable compound names such as BH3•THF. Defaults to True.

    Returns:
        Dict[str, Dict[str, Dict[str, List[str]]]] | Dict[str, str]: A dictionary mapping each compound to its SMILES representation and resolvers, or a simple dictionary mapping each compound to it's selected SMILES representation.
    """
    if not resolvers_list:
        resolvers_list = [
            PubChemNameResolver("pubchem_default"),
            OpsinNameResolver("opsin_default"),
            ManualNameResolver("manual_default"),
            PeptideNameResolver("peptide_default"),
            StructuralFormulaNameResolver("structural_formula_default"),
        ]

    if isinstance(compounds_list, str):
        compounds_list = [compounds_list]
    if not isinstance(compounds_list, list):
        raise ValueError(
            "Invalid input: compounds_list must be a string or a non-empty list of strings."
        )
    if isinstance(compounds_list, list):
        if len(compounds_list) == 0:
            raise ValueError(
                "Invalid input: compounds_list must be a string or a non-empty list of strings."
            )
        for compound in compounds_list:
            if not isinstance(compound, str):
                raise ValueError(
                    "Invalid input: compounds_list must be a string or a non-empty list of strings."
                )
    if len(compounds_list) != len(set(compounds_list)):
        logger.warning("Removing duplicate compound names from compounds_list.")
        compounds_list = list(set(compounds_list))

    if not isinstance(resolvers_list, list) or len(resolvers_list) == 0:
        raise ValueError(
            "Invalid input: resolvers_list must be a non-empty list of ChemicalNameResolver instances."
        )

    seen_resolvers = []
    for resolver in resolvers_list:
        if not isinstance(resolver, ChemicalNameResolver):
            raise ValueError(
                f"Invalid resolver: {resolver} is not an instance of ChemicalNameResolver."
            )
        if resolver.resolver_name in seen_resolvers:
            raise ValueError(f"Duplicate resolver name: {resolver.resolver_name}.")
        seen_resolvers.append(resolver.resolver_name)

    if not (isinstance(smiles_selection_mode, str) or callable(smiles_selection_mode)):
        raise ValueError(
            "Invalid input: smiles_selection_mode must be a string or function."
        )

    if not isinstance(detailed_name_dict, bool):
        raise ValueError("Invalid input: detailed_name_dict must be a bool.")

    if not isinstance(batch_size, int):
        raise TypeError("Invalid input: batch_size must be an integer.")
    if batch_size <= 0 or batch_size > 1000:
        raise ValueError("Invalid input: batch_size must be an integer between 1-1000.")

    if not isinstance(split_names_to_solve, bool):
        raise ValueError("Invalid input: split_names_to_solve must be a bool.")

    # Clean compound names (strip, remove/replace forbidden characters, etc.) and return a mapping dict
    cleaned_compounds_list, cleaned_compounds_dict = clean_strings_and_return_mapping(
        compounds_list
    )

    if split_names_to_solve:
        # Split compound names on delimiters, add split parts to compounds list
        # Return mapping between original compound names and split parts
        # Necessary to resolve names like H₂O•THF
        cleaned_compounds_list, delimiter_split_dict = (
            split_compounds_on_delimiters_and_return_mapping(cleaned_compounds_list)
        )

    # Resolve compounds and split compound names with resolvers
    resolvers_out_dict = resolve_compounds_using_resolvers(
        cleaned_compounds_list, resolvers_list, batch_size
    )

    # Assemble the resolution dictionary
    compounds_out_dict = assemble_compounds_resolution_dict(
        compounds_list, resolvers_out_dict, cleaned_compounds_dict
    )

    if split_names_to_solve:
        # Resolve compounds that were split with split_compounds_on_delimiters_and_return_mapping
        compounds_out_dict = assemble_split_compounds_resolution_dict(
            compounds_out_dict,
            compounds_list,
            resolvers_out_dict,
            cleaned_compounds_dict,
            delimiter_split_dict,
        )

    # Get the resolvers weight dict - needed for SMILESSelector
    resolvers_weight_dict = get_resolvers_weight_dict(resolvers_list)
    resolvers_priority_order = [resolver.resolver_name for resolver in resolvers_list]

    # Select "best" SMILES according to some criteria, add to resolution dict
    compounds_out_dict = select_smiles_with_criteria(
        compounds_out_dict,
        resolvers_weight_dict,
        resolvers_priority_order,
        smiles_selection_mode,
    )

    if not detailed_name_dict:
        return {k: v.get("SMILES", "") for k, v in compounds_out_dict.items()}

    return compounds_out_dict
