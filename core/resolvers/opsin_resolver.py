from typing import List, Tuple, Dict
from py2opsin import py2opsin
from utils.logging_config import logger


def name_to_smiles_opsin(
    chemical_name_list: List[str],
    allow_bad_stereo: bool = False
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Convert a list of chemical names to their corresponding SMILES representations using OPSIN.

    Args:
        chemical_name_list (List[str]): A list of IUPAC or common chemical names to be converted.
        allow_bad_stereo (bool): If True, allows OPSIN to generate SMILES even if stereochemistry
                                is invalid or ambiguous. Defaults to False.

    Returns:
        Tuple[Dict[str, str], Dict[str, str]]:
            - First dict: Mapping from successfully converted chemical names to their SMILES strings.
            - Second dict: Mapping from chemical names that failed conversion to their error messages.

    Note:
        This function uses the `py2opsin` library to interface with OPSIN.
        Newline characters in input names are stripped to avoid CLI parsing issues.
    """
    chemical_name_dict: Dict[str, str] = {}
    failure_message_dict: Dict[str, str] = {}

    # Strip newlines to prevent CLI parsing issues
    sanitized_names = [name.replace('\n', '') for name in chemical_name_list]

    smiles_strings, failure_messages = py2opsin(
        chemical_name=sanitized_names,
        output_format="SMILES",
        failure_analysis=True,
        allow_bad_stereo=allow_bad_stereo
    )

    if len(smiles_strings) != len(chemical_name_list) or len(failure_messages) != len(chemical_name_list):
        logger.warning(
            f"Mismatching lengths: "
            f"smiles_strings ({len(smiles_strings)}), "
            f"chemical_name_list ({len(chemical_name_list)}), "
            f"failure_messages ({len(failure_messages)})"
        )
        return {}, {}

    for name, smiles, msg in zip(chemical_name_list, smiles_strings, failure_messages):
        if smiles:
            chemical_name_dict[name] = smiles
        if msg:
            failure_message_dict[name] = msg

    return chemical_name_dict
