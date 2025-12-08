from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Any

from resolver.split_name_resolver import get_delimiter_split_dict, resolve_delimiter_split_dict
from resolvers.manual_resolver import name_to_smiles_manual
from resolvers.opsin_resolver import name_to_smiles_opsin
from resolvers.pubchem_resolver import name_to_smiles_pubchem
from resolvers.cirpy_resolver import name_to_smiles_cirpy
from resolvers.peptide_resolver import name_to_smiles_peptide
from resolvers.structural_formula_resolver import name_to_smiles_structural_formula
from utils.chem_utils import canonicalize_smiles

class ChemicalNameResolver(ABC):
    """
    Abstract base class for chemical name-to-SMILES resolvers.
    
    Subclasses must implement the `name_to_smiles` method.
    Optionally, they may expose a `config` property for runtime settings.
    """

    def __init__(self, resolver_type: str, resolver_name: str):
        self._resolver_type = resolver_type
        self._resolver_name = resolver_name

    @property
    @abstractmethod
    def config(self) -> Dict:
        """Return current configuration as a dictionary."""
        pass

    @property
    def resolver_name(self) -> str:
        """Return resolver_name."""
        return self._resolver_name

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
    Concrete resolver using OPSIN via py2opsin.
    """

    def __init__(self, resolver_name: str, allow_bad_stereo: bool = False, output_format: str = "SMILES"):
        super().__init__("opsin", resolver_name)
        self._allow_bad_stereo = allow_bad_stereo
        self._output_format = output_format

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {
            "allow_bad_stereo": self._allow_bad_stereo,
            "output_format": self._output_format,
        }

    def set_allow_bad_stereo(self, value: bool):
        """Update stereochemistry tolerance setting."""
        self._allow_bad_stereo = value

    def set_output_format(self, value: str):
        """Update output format (e.g., 'SMILES', 'InChI')."""
        self._output_format = value

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_opsin(chemical_name_list)


class PubChemNameResolver(ChemicalNameResolver):
    """
    """

    def __init__(self, resolver_name: str, output_format: str = "SMILES"):
        super().__init__("pubchem", resolver_name)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_pubchem(chemical_name_list)


class CIRPyNameResolver(ChemicalNameResolver):
    """
    """

    def __init__(self, resolver_name: str, output_format: str = "SMILES"):
        super().__init__("cirpy", resolver_name)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_cirpy(chemical_name_list)


class ManualNameResolver(ChemicalNameResolver):
    """
    """

    def __init__(self, resolver_name: str, output_format: str = "SMILES"):
        super().__init__("manual", resolver_name)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_manual(chemical_name_list)


class PeptideNameResolver(ChemicalNameResolver):
    """
    """

    def __init__(self, resolver_name: str, output_format: str = "SMILES"):
        super().__init__("peptide", resolver_name)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_peptide(chemical_name_list)


class StructuralFormulaNameResolver(ChemicalNameResolver):
    """
    """

    def __init__(self, resolver_name: str, output_format: str = "SMILES"):
        super().__init__("structural_formula", resolver_name)

    @property
    def config(self) -> Dict:
        """Return current configuration."""
        return {}

    def name_to_smiles(
        self,
        chemical_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Convert chemical names to SMILES using OPSIN.
        """
        return name_to_smiles_structural_formula(chemical_name_list)


def resolve_compounds_to_smiles(
    compounds: List[str],
    resolvers: List[ChemicalNameResolver] = [],
    return_type: str = 'full'
):

    if not resolvers:
        manual_resolver = ManualNameResolver('manual_default')
        opsin_resolver = OpsinNameResolver('opsin_default')
        pubchem_resolver =  PubChemNameResolver('pubchem_default')
        peptide_resolver = PeptideNameResolver('peptide_default')
        structural_formula_resolver = StructuralFormulaNameResolver('structural_formula_default')
        cirpy_resolver = CIRPyNameResolver('cirpy_default')

        resolvers = [
            pubchem_resolver,
            opsin_resolver,
            manual_resolver,
            peptide_resolver,
            structural_formula_resolver,
            cirpy_resolver
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

    if not isinstance(resolvers, list) or len(resolvers) == 0:
        raise ValueError("Invalid input: resolvers must be a non-empty list of ChemicalNameResolver instances.")
    for resolver in resolvers:
        if not isinstance(resolver, ChemicalNameResolver):
            raise ValueError(f"Invalid resolver: {resolver} is not an instance of ChemicalNameResolver.")

    delimiter_split_dict = {}
    compounds_split_parts = []
    for compound in compounds:
        delimiter_split_dict, compound_split_parts = get_delimiter_split_dict(compound, delimiter_split_dict)
        compounds_split_parts.extend(compound_split_parts)
    compounds_and_split_parts = compounds + compounds_split_parts
    
    # Apply each resolver
    resolvers_out_dict = {}
    for resolver in resolvers:
        # Call the resolver's method
        out = resolver.name_to_smiles(compounds_and_split_parts)

        resolvers_out_dict[resolver.resolver_name] = {
            "out": out
        }
    
    compounds_out_dict = {}
    for compound in compounds:
        compounds_out_dict[compound] = {
            'smiles': '',
            'smiles_dict': {}
        }
        for resolver, resolution_dict in resolvers_out_dict.items():
            compound_smiles = resolution_dict['out'].get(compound, '')
            canonical_compound_smiles = canonicalize_smiles(compound_smiles)
            if canonical_compound_smiles not in compounds_out_dict[compound]['smiles_dict']:
                compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles] = [resolver]
            else:
                resolvers_list = compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles]
                resolvers_list.append(resolver)
                compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles] = resolvers_list
        if compound in delimiter_split_dict:
            resolved_delimiter_split_dict = resolve_delimiter_split_dict(compound, resolvers_out_dict, delimiter_split_dict)
            for compound_smiles, resolvers in resolved_delimiter_split_dict.items():
                canonical_compound_smiles = canonicalize_smiles(compound_smiles)
                for resolver in resolvers:
                    if canonical_compound_smiles not in compounds_out_dict[compound]['smiles_dict']:
                        compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles] = [resolver]
                    else:
                        resolvers_list = compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles]
                        resolvers_list.append(resolver)
                        compounds_out_dict[compound]['smiles_dict'][canonical_compound_smiles] = resolvers_list
        
    
    assert len(compounds_out_dict) == len(compounds), "An error occurred, different input and output lengths" 

    
    if return_type == 'simple':
        return [v.get('smiles', '') for k,v in compounds_out_dict.items()]

    return compounds_out_dict


if __name__ == "__main__":    
    # print(resolve_compounds_to_smiles(['benzene', 'aspirin', 'stab', 'Asp-Gly', 'CH3CH2OH', 'CH3CH(OH)COOH'], return_type='simple'))

    # print(resolve_compounds_to_smiles(['benzene', 'benzene'], return_type='simple'))

    print(resolve_compounds_to_smiles(['MeOH.benzene']))