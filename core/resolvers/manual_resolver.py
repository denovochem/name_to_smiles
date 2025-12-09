from pathlib import Path
from utils.logging_config import logger
from typing import Dict, List
import json

# Get the directory of the current file
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_MANUAL_NAME_DICT = BASE_DIR.parent / 'datafiles' / 'name_dicts' / 'manual_name_dict.json'

CUSTOM_NAME_DICT_DIR = BASE_DIR.parent / 'datafiles' / 'name_dicts' / 'custom_name_dicts'

def process_name_dict(
    compound_name_list: List[str], 
    name_dict: Dict[str, str], 
):
    """
    """
    processed_name_dict = {}
    compound_name_dict_lower = {ele.lower().strip():ele for ele in compound_name_list}
    name_dict_lower = {k.lower().strip():v for k,v in name_dict.items()}
    for compound, smiles in name_dict_lower.items():
        if smiles:
            if compound in compound_name_dict_lower:
                processed_name_dict[compound_name_dict_lower[compound]] = smiles

    return processed_name_dict

def name_to_smiles_manual(
    compound_name_list: List[str], 
    provided_name_dict: Dict[str, str] = None
):
    """
    """
    manual_name_dict = {}
    compound_name_dict_lower = {ele.lower():ele for ele in compound_name_list}
    if not provided_name_dict:
        with open(DEFAULT_MANUAL_NAME_DICT, 'rb') as handle:
            loaded_manual_name_dict = json.load(handle)
        logger.info(f"Loaded data from {DEFAULT_MANUAL_NAME_DICT}.")
    else:
        loaded_manual_name_dict = provided_name_dict
        logger.info(f"Using provided name dictionary.")

    manual_name_dict = process_name_dict(compound_name_list, loaded_manual_name_dict)

    return manual_name_dict