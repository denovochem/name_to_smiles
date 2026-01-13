from itertools import product
from typing import Dict, List, Tuple

from placeholder_name.utils.chem_utils import canonicalize_smiles
from placeholder_name.utils.constants import DELIMITERS
# from placeholder_name.utils.logging_config import logger


def get_delimiter_split_dict(
    compound_name: str, split_compound_name_dict: Dict[str, List[str]]
) -> Tuple[Dict[str, List[str]], List[str]]:
    """
    Split a compound name into individual parts based on the DELIMITERS list and store the split parts in a dictionary.

    Args:
        compound_name (str): The compound name to split.
        split_compound_name_dict (Dict[str, List[str]]): The dictionary to store the split compound names.

    Returns:
        Dict[str, List[str]]: The dictionary containing the split compound names.
        List[str]: The list of split compound names.
    """
    split_compound_name_list = []
    for delim in DELIMITERS:
        parts = compound_name.split(delim)
        if len(parts) > 1:
            split_compound_name_list.extend(parts)
            split_compound_name_dict[compound_name] = parts

    return split_compound_name_dict, split_compound_name_list


def get_smiles_parts_from_name_parts(
    compound: str,
    resolvers_out_dict: Dict[str, Dict[str, Dict[str, str] | Dict[str, str]]],
    delimiter_split_dict: Dict[str, List[str]],
) -> Dict[str, Dict[str, List[str]]] | None:
    """
    Collect chemical name parts and corresponding SMILES and resolver names that generated the SMILES.

    Args:
        compound (str): The compound name to resolve.
        resolvers_out_dict (Dict[str, Dict[str, Dict[str, str] | Dict[str, str]]]): A dictionary containing the resolver names as keys and dictionaries as values.
            Each dictionary contains the input compound name as key and the resolved SMILES as value.
        delimiter_split_dict (Dict[str, List[str]]): A dictionary containing the compound name as key and a list of split compound names as value.

    Returns:
        Dict[str, Dict[str, List[str]]]: A dictionary mapping each compound to its SMILES representation.
    """
    split_names = delimiter_split_dict[compound]
    smiles_parts_dict: Dict[str, Dict[str, List[str]]] = {}
    for part in split_names:
        smiles_parts_dict[part] = {}
        part_found = False

        for resolver, resolution_dict in resolvers_out_dict.items():
            part_smiles = resolution_dict["out"].get(part)

            if not part_smiles:
                continue

            part_found = True

            canonicalized_part_smiles = canonicalize_smiles(part_smiles)

            if canonicalized_part_smiles in smiles_parts_dict[part]:
                resolvers_list = smiles_parts_dict[part][canonicalized_part_smiles]
                resolvers_list.append(resolver)
                smiles_parts_dict[part][canonicalized_part_smiles] = resolvers_list
            else:
                smiles_parts_dict[part][canonicalized_part_smiles] = [resolver]

        if not part_found:
            return None

    return smiles_parts_dict


def get_smiles_resolvers_combinations(
    compound_smiles_list: List[List[Tuple[str, List[str]]]],
) -> Dict[str, List[str]]:
    """
    Generate all possible combinations of SMILES strings and their corresponding resolvers.

    Args:
        compound_smiles_list (List[List[Tuple[str, List[str]]]]): A list of lists of tuples containing the SMILES strings and their corresponding resolvers.

    Returns:
        Dict[str, List[str]]: A dictionary mapping each SMILES string to a list of possible resolver combinations.
    """
    result: Dict[str, List[str]] = {}
    for smiles_combo in product(*compound_smiles_list):
        combined_smiles = ".".join(smiles for smiles, _ in smiles_combo)

        resolver_lists = [resolvers for _, resolvers in smiles_combo]

        resolver_combos = product(*resolver_lists)

        resolver_strings = [
            "/".join(resolver_combo) for resolver_combo in resolver_combos
        ]

        if combined_smiles not in result:
            result[combined_smiles] = []
        result[combined_smiles].extend(resolver_strings)

    return result


def resolve_delimiter_split_dict(
    compound: str,
    resolvers_out_dict: Dict[str, Dict[str, Dict[str, str] | Dict[str, str]]],
    delimiter_split_dict: Dict[str, List[str]],
) -> Dict[str, List[str]]:
    """
    Resolve a compound name into possible SMILES strings and their corresponding resolvers.

    Args:
        compound (str): The compound name to resolve.
        resolvers_out_dict (Dict[str, Dict[str, Dict[str, str] | Dict[str, str]]]): A dictionary containing the resolver names as keys and dictionaries as values.
            Each dictionary contains the input compound name as key and the resolved SMILES as value
        delimiter_split_dict (Dict[str, List[str]]): Dictionary of compound names to their split parts.

    Returns:
        Dict[str, List[str]]: A dictionary mapping each possible SMILES string to a list of possible resolver combinations.
    """
    if compound not in delimiter_split_dict:
        return {}

    smiles_parts_dict = get_smiles_parts_from_name_parts(
        compound, resolvers_out_dict, delimiter_split_dict
    )

    if smiles_parts_dict is None:
        return {}

    compound_smiles_list = [
        list(smiles_dict.items()) for smiles_dict in smiles_parts_dict.values()
    ]

    result: Dict[str, List[str]] = get_smiles_resolvers_combinations(
        compound_smiles_list
    )

    return result
