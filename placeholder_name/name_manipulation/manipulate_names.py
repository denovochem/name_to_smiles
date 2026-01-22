from typing import Any, Dict, List, Optional

from placeholder_name.name_manipulation.name_correction.dataclasses import (
    CorrectorConfig,
)
from placeholder_name.name_manipulation.name_correction.name_corrector import (
    ChemNameCorrector,
)
from placeholder_name.name_manipulation.peptide_shorthand_handler import (
    looks_like_peptide_shorthand,
    peptide_shorthand_to_iupac,
)
from placeholder_name.utils.logging_config import logger


corrector = ChemNameCorrector()


def correct_names(
    compounds_out_dict: Dict[str, Dict[str, Dict[str, List[str]]]],
    name_correction_config: Optional[CorrectorConfig] = None,
    resolve_peptide_shorthand: Optional[bool] = False,
) -> Dict[str, Dict[str, Any]]:
    """
    Correct chemical names using ChemNameCorrector and optionally resolve peptide shorthand.

    Args:
        compounds_out_dict: Dictionary mapping compound names to their resolution data
        name_correction_config: Configuration for name correction (e.g., allow_acid, allow_radicals)
        resolve_peptide_shorthand: Whether to attempt peptide shorthand expansion

    Returns:
        Dictionary mapping original names to their correction information
    """
    if not name_correction_config:
        corrector = ChemNameCorrector(config=name_correction_config)

    all_compound_correction_dict = {}
    names_to_correct = []
    for k, v in compounds_out_dict.items():
        v["name_correction_info"] = {}
        if not v.get("SMILES", ""):
            names_to_correct.append(k)

    corrected_names = corrector.correct_batch(names_to_correct)
    for name, correction_candidates in corrected_names.items():
        if name not in compounds_out_dict:
            continue
        compound_correction_dict: Dict[str, Any] = {}

        if resolve_peptide_shorthand and looks_like_peptide_shorthand(name):
            try:
                iupac_name = peptide_shorthand_to_iupac(name)
                compound_correction_dict["top_5"] = []
                compound_correction_dict["SMILES"] = ""
                compound_correction_dict["selected_name"] = iupac_name
                compound_correction_dict["name_manipulation_method"] = (
                    "peptide_shorthand_expansion"
                )

                all_compound_correction_dict[name] = compound_correction_dict

                continue

            except Exception as e:
                logger.info(f"Could not resolve peptide shorthand: {e}")
                pass

        compound_correction_dict["top_5"] = [
            [candidate.name, candidate.score] for candidate in correction_candidates[:5]
        ]
        for candidate in correction_candidates:
            if candidate.validation_result:
                resolved_smiles = candidate.validation_result
                compound_correction_dict["SMILES"] = resolved_smiles
                compound_correction_dict["selected_name"] = candidate.name
                compound_correction_dict["name_manipulation_method"] = "name_correction"
                break
            else:
                compound_correction_dict["SMILES"] = ""
                compound_correction_dict["selected_name"] = ""
                compound_correction_dict["name_manipulation_method"] = "name_correction"
        all_compound_correction_dict[name] = compound_correction_dict

    return all_compound_correction_dict
