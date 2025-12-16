from typing import Dict, List

import pubchempy as pcp

from placeholder_name.utils.logging_config import logger
from placeholder_name.utils.string_utils import filter_latin1_compatible


def name_to_smiles_pubchem(compound_name_list: List[str]) -> Dict[str, str]:
    """
    Convert chemical names to SMILES using pubchem.

    Args:
        compound_name_list (List[str]): List of compound names to convert to SMILES.

    Returns:
        Dict[str, str]: Dictionary of compound names to SMILES.
    """
    pubchem_compounds = pcp.get_compounds(
        filter_latin1_compatible(compound_name_list), "name"
    )
    pubchem_name_dict = {
        k: v.smiles if v is not None else ""
        for k, v in zip(compound_name_list, pubchem_compounds)
    }

    if len(compound_name_list) != len(pubchem_name_dict):
        logger.warning(
            f"Mismatching lengths: "
            f"compound_name_list ({len(compound_name_list)}), "
            f"pubchem_name_dict ({len(pubchem_name_dict)})"
        )

    return pubchem_name_dict
