import cirpy
from placeholder_name.utils.logging_config import logger
from typing import List, Dict

def retrieve_cirpy_results(compound_name: str) -> str:
    """
    Retrieves the SMILES string for a given compound identifier using CIRpy.

    Args:
        compound_name (str): The compound name.

    Returns:
        str: The SMILES string of the compound, or an empty string if an exception occurs.
    """
    try:
        smiles = cirpy.resolve(compound_name, 'smiles')
        
    except Exception as e:
        logger.error(f"Exception in CIRpy query: {str(e)}")
        return ''
    
    return smiles

def name_to_smiles_cirpy(compound_name_list: List[str]) -> Dict[str, str]:
    """
    Converts a list of chemical names to their corresponding SMILES strings using CIRpy.

    Args:
        compound_name_list (List[str]): A list of chemical names to be converted.

    Returns:
        Dict[str, str]: A dictionary mapping each chemical name to its SMILES string.
    """
    cirpy_name_dict = {}
    for compound_name in compound_name_list:
        cirpy_name_dict[compound_name] = retrieve_cirpy_results(compound_name)
    return cirpy_name_dict