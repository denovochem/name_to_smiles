import cirpy

def retrieve_cirpy_results(compound_identifier):
    """
    Searches cirpy by compound identifier and identifier type, returning the top match for each property.
    Uses cirpy.resolve() which returns only the single best match for each requested property.

    Args:
        compound_identifier (str): The identifier for the compound being searched

    Returns:
        dict:  {}
    """

    try:
        # Get names (synonyms)
        smiles = cirpy.resolve(compound_identifier, 'smiles')
        
    except Exception as e:
        # Handle any exceptions during CIRpy interaction
        print(f"Exception in CIRpy query: {str(e)}")
    
    return smiles

def name_to_smiles_cirpy(chemical_name_list):
    name_to_smiles_dict = {}
    for chemical_name in chemical_name_list:
        name_to_smiles_dict[chemical_name] = retrieve_cirpy_results(chemical_name)
    return name_to_smiles_dict