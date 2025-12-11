from importlib import resources
import json
from pathlib import Path
from typing import Dict, List

from placeholder_name.utils.logging_config import logger


# Get the directory of the current file
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_MANUAL_NAME_DICT = BASE_DIR.parent / 'datafiles' / 'name_dicts' / 'manual_name_dict.json'

def load_default_manual_name_dict() -> Dict[str, str]:
    """Load manual name dictionary from package data using importlib.resources."""
    # Open the file as text from within the package
    with resources.open_text('placeholder_name.datafiles.name_dicts', 'manual_name_dict.json') as f:
        return json.load(f)

def process_name_dict(
    compound_name_list: List[str], 
    name_dict: Dict[str, str], 
) -> Dict[str, str]:
    """
    Process a dictionary of compound names to SMILES by matching the lowercased and stripped
    compound names with their corresponding SMILES strings.

    Args:
        compound_name_list (List[str]): A list of compound names to process.
        name_dict (Dict[str, str]): A dictionary of compound names to SMILES strings.

    Returns:
        Dict[str, str]: A dictionary of processed compound names to SMILES strings.
    """
    processed_name_dict = {}
    compound_name_dict_lower = {ele.lower().strip():ele for ele in compound_name_list}
    name_dict_lower = {k.lower().strip():v for k,v in name_dict.items()}
    for compound_name, smiles in name_dict_lower.items():
        if not smiles:
            continue
        if compound_name in compound_name_dict_lower:
            processed_name_dict[compound_name_dict_lower[compound_name]] = smiles

    return processed_name_dict

def name_to_smiles_manual(
    compound_name_list: List[str], 
    provided_name_dict: Dict[str, str] = None
) -> Dict[str, str]:
    """
    Convert a list of compound names to their corresponding SMILES strings using a manual name dictionary.

    Args:
        compound_name_list (List[str]): A list of compound names to convert.
        provided_name_dict (Dict[str, str], optional): A manual name dictionary to use. Defaults to None.

    Returns:
        Dict[str, str]: A dictionary of converted compound names to SMILES strings.
    """
    manual_name_dict = {}
    compound_name_dict_lower = {ele.lower():ele for ele in compound_name_list}
    if not provided_name_dict:
        loaded_manual_name_dict = load_default_manual_name_dict()
        logger.info(f"Loaded data from {DEFAULT_MANUAL_NAME_DICT}.")
    else:
        loaded_manual_name_dict = provided_name_dict
        logger.info(f"Using provided name dictionary.")

    manual_name_dict = process_name_dict(compound_name_list, loaded_manual_name_dict)

    return manual_name_dict