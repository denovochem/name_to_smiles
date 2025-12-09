import cirpy
from utils.logging_config import logger
from typing import List, Dict

def retrieve_cirpy_results(compound_identifier: str) -> str:
    """
    Retrieves the SMILES string for a given compound identifier using CIRpy.

    Args:
        compound_identifier (str): The identifier of the compound.

    Returns:
        str: The SMILES string of the compound, or an empty string if an exception occurs.
    """
    try:
        smiles = cirpy.resolve(compound_identifier, 'smiles')
        
    except Exception as e:
        logger.error(f"Exception in CIRpy query: {str(e)}")
        return ''
    
    return smiles

def name_to_smiles_cirpy(chemical_name_list: List[str]) -> Dict[str, str]:
    """
    """
    name_to_smiles_dict = {}
    for chemical_name in chemical_name_list:
        name_to_smiles_dict[chemical_name] = retrieve_cirpy_results(chemical_name)
    return name_to_smiles_dict