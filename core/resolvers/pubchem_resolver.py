import pubchempy as pcp
from utils.string_utils import clean_strings, filter_latin1_compatible



def name_to_smiles_pubchem(compound_name_list):
    pubchem_name_dict = {}

    cleaned_compound_name_list = [clean_strings(ele) for ele in compound_name_list]
    cleaned_strings_dict = {k:v for k,v in zip(cleaned_compound_name_list, compound_name_list)}
    if cleaned_strings_dict == {}:
        return {}
    
    print(filter_latin1_compatible(cleaned_compound_name_list))
    pubchem_compounds = pcp.get_compounds(filter_latin1_compatible(cleaned_compound_name_list), 'name')
    pubchem_name_dict = {k:v.smiles if v is not None else '' for k,v in zip(cleaned_compound_name_list, pubchem_compounds)}
                        
    return pubchem_name_dict