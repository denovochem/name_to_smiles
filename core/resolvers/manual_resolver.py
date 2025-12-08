from pathlib import Path
from utils.logging_config import logger
import json

# Get the directory of the current file
BASE_DIR = Path(__file__).resolve().parent
MANUAL_NAME_DICT = BASE_DIR.parent / 'datafiles' / 'name_dicts' / 'manual_name_dict.json'

CUSTOM_NAME_DICT_DIR = BASE_DIR.parent / 'datafiles' / 'name_dicts' / 'custom_name_dicts'

def process_name_dict(name_dict, compound_name_list):
    processed_name_dict = {}
    compound_name_dict_lower = {ele.lower():ele for ele in compound_name_list}
    name_dict_lower = {k.lower():v for k,v in name_dict.items()}
    for compound, smiles in name_dict_lower.items():
        if smiles:
            if compound in compound_name_dict_lower:
                processed_name_dict[compound_name_dict_lower[compound]] = smiles

    return processed_name_dict

def name_to_smiles_manual(compound_name_list):
    
    manual_name_dict = {}
    compound_name_dict_lower = {ele.lower():ele for ele in compound_name_list}
    with open(MANUAL_NAME_DICT, 'rb') as handle:
        loaded_manual_name_dict = json.load(handle)
    logger.info(f"Loaded data from {MANUAL_NAME_DICT}")

    manual_name_dict = process_name_dict(loaded_manual_name_dict, compound_name_list)
    
    custom_name_dicts = {}
    for json_file in CUSTOM_NAME_DICT_DIR.glob("*.json"):
        if json_file.name == 'example_custom_dict.json':
            continue
        
        logger.info(f"Opening: {json_file.name}")
        with open(json_file, 'r', encoding='utf-8') as f:
            try:
                data = json.load(f)
                logger.info(f"Loaded {len(data)} items from {json_file.name}")
                processed_name_dict = process_name_dict(data['name_dict'], compound_name_list)
                custom_name_dicts[data['name_dict_name']] = processed_name_dict

            except json.JSONDecodeError as e:
                logger.info(f"Error decoding JSON in {json_file.name}: {e}")
            except Exception as e:
                logger.info(f"Error reading {json_file.name}: {e}")

                    
    return manual_name_dict