from abc import ABC, abstractmethod
from typing import List, Tuple, Dict
import warnings

from resolvers.split_name_resolver import get_delimiter_split_dict, resolve_delimiter_split_dict
from resolvers.manual_resolver import name_to_smiles_manual
from resolvers.opsin_resolver import name_to_smiles_opsin
from resolvers.pubchem_resolver import name_to_smiles_pubchem
from resolvers.cirpy_resolver import name_to_smiles_cirpy
from resolvers.peptide_resolver import name_to_smiles_peptide
from resolvers.structural_formula_resolver import name_to_smiles_structural_formula
from smiles_selector import SMILESSelector
from utils.chem_utils import canonicalize_smiles
from utils.logging_config import logger
from utils.string_utils import clean_strings

class ChemicalNameResolver(ABC):
    """
    Abstract base class for chemical name-to-SMILES resolvers.
    
    Subclasses must implement the `name_to_smiles` method.
    Optionally, they may expose a `config` property for runtime settings.
    """

    def __init__(self, resolver_type: str, resolver_name: str, resolver_weight: float):
        if not isinstance(resolver_type, str):
            raise TypeError("Invalid input: resolver_type must be a string.")
        self._resolver_type: str = resolver_type
        if not isinstance(resolver_name, str):
            raise TypeError("Invalid input: resolver_name must be a string.")
        self._resolver_name: str = resolver_name
        if not isinstance(resolver_weight, (int, float)):
            raise TypeError("Invalid input: resolver_weight must be a number between 0-1000.")
        if resolver_weight < 0 or resolver_weight > 1000:
            raise ValueError("Invalid input: resolver_weight must be a number between 0-1000.")
        self._resolver_weight: float = float(resolver_weight)

    @property
    @abstractmethod
    def config(self) -> Dict:
        """Return current configuration as a dictionary."""
        pass

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
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES strings.

        Args:
            chemical_name_list: List of chemical names.

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

    def __init__(self, resolver_name: str, allow_bad_stereo: bool = False, resolver_weight: float = 3):
        super().__init__("opsin", resolver_name, resolver_weight)
        self._allow_bad_stereo = allow_bad_stereo

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {
            "allow_bad_stereo": self._allow_bad_stereo,
        }

    def set_allow_bad_stereo(self, value: bool):
        """Update stereochemistry tolerance setting."""
        self._allow_bad_stereo = value

    def set_output_format(self, value: str):
        """Update output format (e.g., 'SMILES', 'InChI')."""

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        resolved_names, failure_message_dict = name_to_smiles_opsin(chemical_name_list)
        return resolved_names, failure_message_dict


class PubChemNameResolver(ChemicalNameResolver):
    """
    Resolver using PubChem via PubChemPy.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 2):
        super().__init__("pubchem", resolver_name, resolver_weight)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using pubchem.
        """
        resolved_names = name_to_smiles_pubchem(chemical_name_list)
        return resolved_names, None


class CIRPyNameResolver(ChemicalNameResolver):
    """
    Resolver using Chemical Identity Resolver via CIRPy.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 1):
        super().__init__("cirpy", resolver_name, resolver_weight)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using cirpy.
        """
        resolved_names = name_to_smiles_cirpy(chemical_name_list)
        return resolved_names, None


class ManualNameResolver(ChemicalNameResolver):
    """
    Resolver using manually curated names and corresponding SMILES. 
    """

    def __init__(self, resolver_name: str, provided_name_dict: dict = None, resolver_weight: float = 10):
        super().__init__("manual", resolver_name, resolver_weight)
        if provided_name_dict:
            if not isinstance(provided_name_dict, dict):
                raise TypeError("Invalid input: provided_name_dict must be a dictionary.")
            for k, v in provided_name_dict.items():
                if not isinstance(k, str) or not isinstance(v, str):
                    raise ValueError("Invalid input: keys and values in provided_name_dict must be strings.")

        self._provided_name_dict = provided_name_dict

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str],
        provided_name_dict: Dict[str, str] = None
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using manual name database.
        """
        if provided_name_dict is None:
            provided_name_dict = self._provided_name_dict
        resolved_names = name_to_smiles_manual(chemical_name_list, provided_name_dict)
        return resolved_names, None


class PeptideNameResolver(ChemicalNameResolver):
    """
    Resolver using peptide shorthand-to-IUPAC-like name, then resolved to SMILES
    with OPSIN via py2opsin.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 3):
        super().__init__("peptide", resolver_name, resolver_weight)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using peptide name converter and OPSIN.
        """
        resolved_names, failure_message_dict = name_to_smiles_peptide(chemical_name_list)
        return resolved_names, failure_message_dict


class StructuralFormulaNameResolver(ChemicalNameResolver):
    """
    Resolver using structural chemical formula (e.g. CH3CH2CH2COOH).
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 2):
        super().__init__("structural_formula", resolver_name, resolver_weight)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using structural formula converter.
        """
        resolved_names = name_to_smiles_structural_formula(chemical_name_list)
        return resolved_names, None


def resolve_compounds_to_smiles(
    compounds: List[str],
    resolvers: List[ChemicalNameResolver] = [],
    smiles_selection_mode: str = 'weighted',
    detailed_name_dict: bool = False,
    batch_size: int = 500
):

    if not resolvers:
        manual_resolver = ManualNameResolver('manual_default')
        opsin_resolver = OpsinNameResolver('opsin_default')
        pubchem_resolver =  PubChemNameResolver('pubchem_default')
        peptide_resolver = PeptideNameResolver('peptide_default')
        structural_formula_resolver = StructuralFormulaNameResolver('structural_formula_default')

        resolvers = [
            pubchem_resolver,
            opsin_resolver,
            manual_resolver,
            peptide_resolver,
            structural_formula_resolver
        ]

    if not isinstance(compounds, list):
        raise ValueError("Invalid input: compounds must be a string or a non-empty list of strings.")
    if isinstance(compounds, list):
        if len(compounds) == 0:
            raise ValueError("Invalid input: compounds must be a string or a non-empty list of strings.")
        for compound in compounds:
            if not isinstance(compound, str):
                raise ValueError("Invalid input: compounds must be a string or a non-empty list of strings.")
    if isinstance(compounds, str):
        compounds = [compounds]
    if len(compounds) != len(set(compounds)):
        warnings.warn("Removing duplicate compound names from input compounds list.")
        logger.info("Removing duplicate compound names from input compounds list.")
        compounds = list(set(compounds))

    if not isinstance(resolvers, list) or len(resolvers) == 0:
        raise ValueError("Invalid input: resolvers must be a non-empty list of ChemicalNameResolver instances.")
        
    seen_resolvers = []
    for resolver in resolvers:
        if not isinstance(resolver, ChemicalNameResolver):
            raise ValueError(f"Invalid resolver: {resolver} is not an instance of ChemicalNameResolver.")
        if resolver.resolver_name in seen_resolvers:
            raise ValueError(f"Duplicate resolver name: {resolver.resolver_name}.")
        seen_resolvers.append(resolver.resolver_name)

    if not isinstance(smiles_selection_mode, str):
        raise ValueError("Invalid input: smiles_selection_mode must be a string.")
    
    if not isinstance(detailed_name_dict, bool):
        raise ValueError("Invalid input: detailed_name_dict must be a bool.")
    
    if not isinstance(batch_size, int):
        raise TypeError("Invalid input: batch_size must be an integer.")
    if batch_size <= 0 or batch_size > 1000:
        raise ValueError("Invalid input: batch_size must be an integer between 1-1000.")

    cleaned_compounds = [clean_strings(ele) for ele in compounds]
    cleaned_compounds_dict = {k:v for k,v in zip(compounds, cleaned_compounds)}
    
    resolvers_weight_dict = {}
    for resolver in resolvers:
        resolvers_weight_dict[resolver.resolver_name] = resolver.resolver_weight

    delimiter_split_dict = {}
    compounds_split_parts = []
    for compound in cleaned_compounds:
        delimiter_split_dict, compound_split_parts = get_delimiter_split_dict(compound, delimiter_split_dict)
        compounds_split_parts.extend(compound_split_parts)
    compounds_and_split_parts = cleaned_compounds + compounds_split_parts
    compounds_and_split_parts = list(set(compounds_and_split_parts))
    
    resolvers_out_dict = {}
    for resolver in resolvers:
        full_resolver_dict = {}
        for i in range(0, len(compounds_and_split_parts), batch_size):
            chunk = compounds_and_split_parts[i:i + batch_size]
            out, _ = resolver.name_to_smiles(chunk)
            full_resolver_dict.update(out)

        resolvers_out_dict[resolver.resolver_name] = {
            "out": out
        }
    
    compounds_out_dict = {}
    for compound in compounds:

        compound_cleaned = cleaned_compounds_dict.get(compound, '')

        compounds_out_dict[compound] = {
            'SMILES': '',
            'SMILES_source': [],
            'SMILES_dict': {}
        }

        for resolver, resolution_dict in resolvers_out_dict.items():
            compound_smiles = resolution_dict['out'].get(compound_cleaned, '')
            canonical_compound_smiles = canonicalize_smiles(compound_smiles)
            if not canonical_compound_smiles:
                continue
            if canonical_compound_smiles not in compounds_out_dict[compound]['SMILES_dict']:
                compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles] = [resolver]
            else:
                resolvers_list = compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles]
                resolvers_list.append(resolver)
                compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles] = resolvers_list

        if compound_cleaned in delimiter_split_dict:
            resolved_delimiter_split_dict = resolve_delimiter_split_dict(compound_cleaned, resolvers_out_dict, delimiter_split_dict)
            for compound_smiles, resolvers in resolved_delimiter_split_dict.items():
                canonical_compound_smiles = canonicalize_smiles(compound_smiles)
                if not canonical_compound_smiles:
                    continue
                for resolver in resolvers:
                    if canonical_compound_smiles not in compounds_out_dict[compound]['SMILES_dict']:
                        compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles] = [resolver]
                    else:
                        resolvers_list = compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles]
                        resolvers_list.append(resolver)
                        compounds_out_dict[compound]['SMILES_dict'][canonical_compound_smiles] = resolvers_list
        
    
    assert len(compounds_out_dict) == len(compounds), "An error occurred, different input and output lengths" 

    selector = SMILESSelector(compounds_out_dict, resolvers_weight_dict)
    for k,v in compounds_out_dict.items():
        selected_smiles, selected_smiles_resolvers = selector.select_smiles(k, smiles_selection_mode)
        v['SMILES'] = selected_smiles
        v['SMILES_source'] = selected_smiles_resolvers
    
    if not detailed_name_dict:
        return {k:v.get('SMILES', '') for k,v in compounds_out_dict.items()}

    return compounds_out_dict


if __name__ == "__main__":    
    print(resolve_compounds_to_smiles(['benzene', 'aspirin', 'stab', 'Asp-Gly', 'CH3CH2OH', 'CH3CH(OH)COOH', 'H₂O', 'H₂O.THF'], detailed_name_dict=False))

    # print(resolve_compounds_to_smiles(['benzene', 'benzene'], detailed_name_dict=False))

    # print(resolve_compounds_to_smiles(['MeOH.benzene'], detailed_name_dict=False))