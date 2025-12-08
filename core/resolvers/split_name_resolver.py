from itertools import product
from utils.constants import DELIMITERS

def get_delimiter_split_dict(compound_name, split_compound_name_dict):
    split_compound_name_list = []
    for delim in DELIMITERS:
        parts = compound_name.split(delim)
        if len(parts) > 1:
            split_compound_name_list.extend(parts)
            split_compound_name_dict[compound_name] = parts

    return split_compound_name_dict, split_compound_name_list

def resolve_delimiter_split_dict(compound, resolvers_out_dict, delimiter_split_dict):
    if compound not in delimiter_split_dict:
        return {}

    split_names = delimiter_split_dict[compound]

    all_parts_found = True
    smiles_parts_dict = {}
    
    for part in split_names:
        smiles_parts_dict[part] = {}
        part_found = False
        part_smiles = None

        for resolver, resolution_dict in resolvers_out_dict.items():
            part_smiles = resolution_dict['out'].get(part, '')
            canonicalized_part_smiles = canonicalize_smiles(part_smiles)

            if part_smiles:
                part_found = True
                if canonicalized_part_smiles in smiles_parts_dict[part]:
                    resolvers_list = smiles_parts_dict[part][canonicalized_part_smiles]
                    resolvers_list.append(resolver)
                    smiles_parts_dict[part][canonicalized_part_smiles] = resolvers_list
                else:
                    smiles_parts_dict[part][canonicalized_part_smiles] = [resolver]
                
        if not part_found:
            all_parts_found = False
            return {}

    result = {}

    compound_names = list(smiles_parts_dict.keys())
    compound_smiles_list = [list(smiles_dict.items()) for smiles_dict in smiles_parts_dict.values()]

    for smiles_combo in product(*compound_smiles_list):
        combined_smiles = '.'.join(smiles for smiles, _ in smiles_combo)
        
        resolver_lists = [resolvers for _, resolvers in smiles_combo]
        
        resolver_combos = product(*resolver_lists)
        
        resolver_strings = ['/'.join(resolver_combo) for resolver_combo in resolver_combos]
        
        if combined_smiles not in result:
            result[combined_smiles] = []
        result[combined_smiles].extend(resolver_strings)

    return result