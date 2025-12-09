import pubchempy as pcp
from utils.string_utils import clean_strings, filter_latin1_compatible
from utils.logging_config import logger
from typing import List

def name_to_smiles_pubchem(compound_name_list: List[str]):
    """
    """
    pubchem_name_dict = {}
    pubchem_compounds = pcp.get_compounds(filter_latin1_compatible(compound_name_list), 'name')
    pubchem_name_dict = {k:v.smiles if v is not None else '' for k,v in zip(compound_name_list, pubchem_compounds)}

    return pubchem_name_dict