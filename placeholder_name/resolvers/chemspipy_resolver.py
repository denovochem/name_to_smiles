from typing import Dict, List

from chemspipy import ChemSpider

from placeholder_name.utils.logging_config import logger
from placeholder_name.utils.string_utils import filter_latin1_compatible


def name_to_smiles_chemspipy(
    compound_name_list: List[str],
    chemspider_api_key: str,
) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using ChemSpiPy.

    Args:
        compound_name_list (List[str]): List of compound names to convert to SMILES.
        chemspider_api_key (str): ChemSpider API key (https://developer.rsc.org/getting-started)

    Returns:
        Dict[str, str]: Dictionary of compound names to SMILES.
    """
    try:
        cs = ChemSpider(chemspider_api_key)
    except Exception as e:
        logger.warning(f"Error initializing ChemSpiPy: {e}")
        return {}

    chemspipy_name_dict = {}
    for compound_name in filter_latin1_compatible(compound_name_list):
        results = cs.search(compound_name)
        results.wait()
        results.ready()
        if len(results) == 0:
            continue
        c = results[0]
        if not c.smiles:
            continue
        chemspipy_name_dict[compound_name] = c.smiles

    return chemspipy_name_dict
