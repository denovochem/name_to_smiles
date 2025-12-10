from typing import List, Tuple, Dict
from py2opsin import py2opsin
from placeholder_name.utils.logging_config import logger


def name_to_smiles_opsin(
    compound_name_list: List[str],
    allow_acid: bool = False,
    allow_radicals: bool = True,
    allow_bad_stereo: bool = False,
    wildcard_radicals: bool = False
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Convert a list of chemical names to their corresponding SMILES representations using OPSIN.

    Args:
    compound_name_list (List[str]): A list of IUPAC or common chemical names to be converted.
    allow_acid (bool): If True, allow interpretation of acids.
    allow_radicals (bool): If True, enable radical interpretation.
    allow_bad_stereo (bool): If True, allow OPSIN to ignore uninterpretable stereochem.
    wildcard_radicals (bool): If True, output radicals as wildcards.

    Returns:
    Tuple[Dict[str, str], Dict[str, str]]:
        - First dict: Mapping from successfully converted chemical names to their SMILES strings.
        - Second dict: Mapping from chemical names that failed conversion to their error messages.

    Note:
        This function uses the `py2opsin` library to interface with OPSIN.
        Newline characters in input names are stripped to avoid CLI parsing issues.
    """
    opsin_name_dict: Dict[str, str] = {}
    failure_message_dict: Dict[str, str] = {}

    # Strip newlines to prevent CLI parsing issues
    sanitized_names = [compound_name.replace('\n', '') for compound_name in compound_name_list]

    smiles_strings, failure_messages = py2opsin(
        chemical_name=sanitized_names,
        output_format="SMILES",
        failure_analysis=True,
        allow_acid=allow_acid,
        allow_radicals=allow_radicals,
        allow_bad_stereo=allow_bad_stereo,
        wildcard_radicals=wildcard_radicals
    )

    if len(smiles_strings) != len(compound_name_list) or len(failure_messages) != len(compound_name_list):
        logger.warning(
            f"Mismatching lengths: "
            f"smiles_strings ({len(smiles_strings)}), "
            f"compound_name_list ({len(compound_name_list)}), "
            f"failure_messages ({len(failure_messages)})"
        )
        return {}, {}

    for compound_name, smiles, msg in zip(compound_name_list, smiles_strings, failure_messages):
        if smiles:
            opsin_name_dict[compound_name] = smiles
        if msg:
            failure_message_dict[compound_name] = msg

    return opsin_name_dict, failure_message_dict
