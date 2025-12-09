import cirpy
from utils.logging_config import logger
from typing import List

def retrieve_cirpy_results(compound_identifier: str):
    """
    """

    try:
        smiles = cirpy.resolve(compound_identifier, 'smiles')
        
    except Exception as e:
        print(f"Exception in CIRpy query: {str(e)}")
    
    return smiles

def name_to_smiles_cirpy(chemical_name_list: List[str]):
    """
    """
    name_to_smiles_dict = {}
    for chemical_name in chemical_name_list:
        name_to_smiles_dict[chemical_name] = retrieve_cirpy_results(chemical_name)
    return name_to_smiles_dict