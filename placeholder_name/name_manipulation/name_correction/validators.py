from typing import Dict, List, Optional, Tuple

from typing_extensions import Protocol

from placeholder_name.resolvers.opsin_resolver import name_to_smiles_opsin


class Validator(Protocol):
    """Protocol for external validation of chemical names."""

    def validate(self, name: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a chemical name.

        Args:
            name: Chemical name to validate

        Returns:
            Tuple of (is_valid, result) where result could be SMILES or None
        """
        ...

    def batch_validate(self, names: List[str]) -> Dict[str, Tuple[bool, Optional[str]]]:
        """
        Validate multiple chemical names.

        Args:
            names: List of chemical names to validate

        Returns:
            Dictionary mapping names to (is_valid, smiles_or_none) tuples
        """
        ...


class OPSINValidator(Validator):
    """
    Validator that uses OPSIN (via subprocess or API) to validate names.

    OPSIN must be installed and accessible via command line.
    """

    def __init__(self, opsin_path: str = "opsin"):
        """
        Initialize OPSIN validator.

        Args:
            opsin_path: Path to OPSIN executable
        """
        self.opsin_path = opsin_path

    def validate(self, name: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a chemical name using OPSIN.

        Args:
            name: Chemical name to validate

        Returns:
            Tuple of (is_valid, smiles_or_none)
        """
        smiles_dict, _ = name_to_smiles_opsin([name])
        if name in smiles_dict:
            return True, smiles_dict[name]
        return False, None

    def batch_validate(self, names: List[str]) -> Dict[str, Tuple[bool, Optional[str]]]:
        """
        Validate multiple chemical names using OPSIN.

        Args:
            names: List of chemical names to validate

        Returns:
            Dictionary mapping names to (is_valid, smiles_or_none) tuples
        """
        smiles_dict, _ = name_to_smiles_opsin(names)
        results: Dict[str, Tuple[bool, Optional[str]]] = {}
        for name in names:
            if name in smiles_dict:
                results[name] = (True, smiles_dict[name])
            else:
                results[name] = (False, None)

        return results
