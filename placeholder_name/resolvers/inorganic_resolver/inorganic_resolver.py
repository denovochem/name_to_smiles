from __future__ import annotations
from dataclasses import dataclass, field
import re
from typing import Dict, List, Optional, Tuple

from placeholder_name.resolvers.inorganic_resolver.inorganic_resolver_tokens import (
    COUNTER_ION_DATABASE,
    LIGAND_DATABASE,
    LigandInfo,
    METAL_DATABASE,
    MetalInfo,
)
from placeholder_name.utils.logging_config import logger


@dataclass
class ParsedLigand:
    """
    Represents a ligand as parsed from a complex name.

    Attributes:
        name: The ligand identifier as it appears in the name
        count: Number of this ligand coordinated to the metal
        modifiers: Any prefix/suffix modifiers (e.g., "dF(CF3)" in "dF(CF3)ppy")
    """

    name: str
    count: int = 1
    modifiers: List[str] = field(default_factory=list)

    def __repr__(self) -> str:
        return f"ParsedLigand(name='{self.name}', count={self.count})"


@dataclass
class ParsedComplex:
    """
    Complete parsed representation of a coordination complex.

    Attributes:
        metal: Metal symbol (e.g., "Ir", "Rh")
        ligands: List of ParsedLigand objects
        complex_charge: Overall charge of the complex ion
        multiplicity: Number of formula units (e.g., 2 for dimers like [IrCl(cod)]2)
        counter_ions: List of (ion_name, count) tuples
    """

    metal: str
    ligands: List[ParsedLigand]
    complex_charge: int = 0
    multiplicity: int = 1
    counter_ions: List[Tuple[str, int]] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"ParsedComplex(metal='{self.metal}', "
            f"ligands={self.ligands}, "
            f"charge={self.complex_charge}, "
            f"mult={self.multiplicity}, "
            f"counter_ions={self.counter_ions})"
        )


class ParserError(Exception):
    """Custom exception for parsing errors."""

    pass


class ComplexNameParser:
    """
    Parser for inorganic/organometallic complex names.

    Handles names in common formats such as:
        - [Metal(Ligand)n(Ligand2)m]charge
        - [Metal(Ligand)n]multiplicity
        - [Metal(Ligand)n]CounterIon

    Examples:
        >>> parser = ComplexNameParser()
        >>> result = parser.parse("[IrCl(cod)]2")
        >>> print(result.metal)
        'Ir'
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the parser with chemical databases.

        Args:
            ligand_db: Dictionary mapping ligand abbreviations to LigandInfo.
                      Uses LIGAND_DATABASE if None.
            metal_db: Dictionary mapping metal symbols to MetalInfo.
                     Uses METAL_DATABASE if None.
            counter_ion_db: Dictionary mapping counter ion names to LigandInfo.
                           Uses COUNTER_ION_DATABASE if None.
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

    def parse(self, name: str) -> ParsedComplex:
        """
        Parse an inorganic complex name into its components.

        Args:
            name: The complex name string (e.g., "[IrCl(cod)]2")

        Returns:
            ParsedComplex object containing all parsed components

        Raises:
            ParserError: If the name cannot be parsed
        """
        original_name = name.strip()
        working_name = original_name

        # Step 1: Extract trailing multiplicity (e.g., "]2" at end)
        multiplicity, working_name = self._extract_multiplicity(working_name)

        # Step 2: Extract counter ions (after the complex brackets)
        counter_ions, working_name = self._extract_counter_ions(working_name)

        # Step 3: Extract complex charge (e.g., "]+" or "]2-")
        complex_charge, working_name = self._extract_charge(working_name)

        # Step 4: Remove outer brackets
        working_name = self._strip_brackets(working_name)

        # Step 5: Extract metal symbol
        metal, ligand_string = self._extract_metal(working_name)

        # Step 6: Parse ligands
        ligands = self._parse_ligand_string(ligand_string)

        return ParsedComplex(
            metal=metal,
            ligands=ligands,
            complex_charge=complex_charge,
            multiplicity=multiplicity,
            counter_ions=counter_ions,
        )

    def _extract_multiplicity(self, name: str) -> Tuple[int, str]:
        """
        Extract multiplicity from end of name (e.g., "]2").

        Args:
            name: Complex name string

        Returns:
            Tuple of (multiplicity, remaining_string)
        """
        # Pattern: ends with ]<number> where number is the multiplicity
        match = re.search(r"\](\d+)$", name)
        if match:
            multiplicity = int(match.group(1))
            # Keep the closing bracket, remove the number
            return multiplicity, name[: match.start() + 1]
        return 1, name

    def _extract_counter_ions(self, name: str) -> Tuple[List[Tuple[str, int]], str]:
        """
        Extract counter ions from after the complex brackets.

        Args:
            name: Complex name string

        Returns:
            Tuple of (list of (ion_name, count) tuples, remaining_string)
        """
        counter_ions: List[Tuple[str, int]] = []

        # Find where the complex ends (after closing bracket)
        if not name.endswith("]"):
            # Look for content after the last ]
            match = re.search(r"\]([A-Za-z0-9()]+)$", name)
            if match:
                counter_string = match.group(1)
                remaining = name[: match.start() + 1]

                # Try to match known counter ions
                for ion_name in sorted(
                    self.counter_ion_db.keys(), key=len, reverse=True
                ):
                    # Look for the ion with optional count
                    ion_pattern = rf"({re.escape(ion_name)})(\d*)"
                    ion_match = re.search(ion_pattern, counter_string)
                    if ion_match:
                        count = int(ion_match.group(2)) if ion_match.group(2) else 1
                        counter_ions.append((ion_name, count))
                        counter_string = counter_string.replace(
                            ion_match.group(0), "", 1
                        )

                return counter_ions, remaining

        return counter_ions, name

    def _extract_charge(self, name: str) -> Tuple[int, str]:
        """
        Extract complex charge notation (e.g., "]+" or "]2-").

        Args:
            name: Complex name string

        Returns:
            Tuple of (charge, remaining_string)
        """
        # Pattern: ]<optional_number><+/->
        match = re.search(r"\](\d*)([+-])$", name)
        if match:
            charge_magnitude = int(match.group(1)) if match.group(1) else 1
            charge_sign = 1 if match.group(2) == "+" else -1
            charge = charge_magnitude * charge_sign
            return charge, name[: match.start() + 1]
        return 0, name

    def _strip_brackets(self, name: str) -> str:
        """
        Remove outer square brackets from complex name.

        Args:
            name: Complex name string

        Returns:
            String with outer brackets removed
        """
        name = name.strip()
        if name.startswith("[") and name.endswith("]"):
            return name[1:-1]
        return name

    def _extract_metal(self, name: str) -> Tuple[str, str]:
        """
        Extract metal symbol from beginning of name.

        Args:
            name: Complex name string (without outer brackets)

        Returns:
            Tuple of (metal_symbol, remaining_ligand_string)

        Raises:
            ParserError: If no known metal is found
        """
        # Try metals sorted by length (longest first) to handle cases like
        # "Ir" vs "I" (iodine)
        ## TODO: add check to see if metal symbol is in ligand tokens?
        potential_metal_symbols = {}
        for metal_symbol in sorted(self.metal_db.keys(), key=len, reverse=True):
            idx = name.find(metal_symbol)
            if idx != -1:
                potential_metal_symbols[metal_symbol] = idx

        if not potential_metal_symbols:
            raise ParserError(f"Could not identify metal in: {name}")

        if len(potential_metal_symbols) > 1:
            raise ParserError(f"Multiple potential metal symbols found in: {name}")

        best_idx = 1e10
        best_metal_symbol = ""
        for metal_symbol, idx in potential_metal_symbols.items():
            if idx < best_idx:
                best_metal_symbol = metal_symbol
                best_idx = idx
                continue
            if idx == best_idx:
                if len(metal_symbol) < best_metal_symbol:
                    continue
                best_metal_symbol = metal_symbol

        if not best_metal_symbol:
            raise ParserError(f"Could not identify metal in: {name}")

        return best_metal_symbol, name[:best_idx] + name[
            best_idx + len(best_metal_symbol) :
        ]

    def _parse_ligand_string(self, ligand_str: str) -> List[ParsedLigand]:
        """
        Parse the ligand portion of a complex name.

        Handles:
            - Parenthesized ligands: (cod), (ppy)2
            - Direct ligands: Cl, Cl2
            - Complex nested names: dF(CF3)ppy

        Args:
            ligand_str: String containing ligand specifications

        Returns:
            List of ParsedLigand objects
        """
        ligands: List[ParsedLigand] = []
        i = 0

        while i < len(ligand_str):
            # Skip whitespace
            if ligand_str[i].isspace():
                i += 1
                continue

            # Handle parenthesized ligands
            if ligand_str[i] == "(":
                ligand_name, count, new_i = self._parse_parenthesized_ligand(
                    ligand_str, i
                )
                ligands.append(ParsedLigand(name=ligand_name, count=count))
                i = new_i
            else:
                # Try to match a known ligand
                matched, new_i = self._try_match_known_ligand(ligand_str, i, ligands)
                if matched:
                    i = new_i
                else:
                    raise ParserError(f"could not match known ligand: {ligand_str}")
                    # print('WHAT')
                    # # Extract unknown token
                    # token, new_i = self._extract_unknown_token(ligand_str, i)
                    # if token:
                    #     name, count = self._split_name_and_count(token)
                    #     ligands.append(ParsedLigand(name=name, count=count))
                    # i = new_i

        return ligands

    def _parse_parenthesized_ligand(
        self, ligand_str: str, start: int
    ) -> Tuple[str, int, int]:
        """
        Parse a parenthesized ligand expression.

        Args:
            ligand_str: Full ligand string
            start: Index of opening parenthesis

        Returns:
            Tuple of (ligand_name, count, new_position)
        """
        # Find matching closing parenthesis
        depth = 1
        j = start + 1
        while j < len(ligand_str) and depth > 0:
            if ligand_str[j] == "(":
                depth += 1
            elif ligand_str[j] == ")":
                depth -= 1
            j += 1

        ligand_name = ligand_str[start + 1 : j - 1]

        # Check for count after closing parenthesis
        count = 1
        if j < len(ligand_str) and ligand_str[j].isdigit():
            count_str = ""
            while j < len(ligand_str) and ligand_str[j].isdigit():
                count_str += ligand_str[j]
                j += 1
            count = int(count_str)

        return ligand_name, count, j

    def _try_match_known_ligand(
        self, ligand_str: str, start: int, ligands: List[ParsedLigand]
    ) -> Tuple[bool, int]:
        """
        Try to match a known ligand at the current position.

        Args:
            ligand_str: Full ligand string
            start: Current position
            ligands: List to append matched ligand to

        Returns:
            Tuple of (matched: bool, new_position: int)
        """
        remaining = ligand_str[start:]

        # Sort by length (longest first) to match greedily
        for lig_name in sorted(self.ligand_db.keys(), key=len, reverse=True):
            if remaining.startswith(lig_name):
                j = start + len(lig_name)

                # Check for count after ligand
                count = 1
                if j < len(ligand_str) and ligand_str[j].isdigit():
                    count_str = ""
                    while j < len(ligand_str) and ligand_str[j].isdigit():
                        count_str += ligand_str[j]
                        j += 1
                    count = int(count_str)

                ligands.append(ParsedLigand(name=lig_name, count=count))
                return True, j

        return False, start

    def _extract_unknown_token(self, ligand_str: str, start: int) -> Tuple[str, int]:
        """
        Extract an unknown token until a delimiter is reached.

        Args:
            ligand_str: Full ligand string
            start: Current position

        Returns:
            Tuple of (token, new_position)
        """
        j = start
        delimiters = "()[]"

        while j < len(ligand_str) and ligand_str[j] not in delimiters:
            j += 1

        token = ligand_str[start:j] if j > start else ""
        return token, max(j, start + 1)

    def _split_name_and_count(self, token: str) -> Tuple[str, int]:
        """
        Split a token into name and trailing count.

        Args:
            token: Token string (e.g., "Cl2")

        Returns:
            Tuple of (name, count)
        """
        match = re.match(r"^(.+?)(\d+)$", token)
        if match:
            return match.group(1), int(match.group(2))
        return token, 1


class SMILESBuilderError(Exception):
    """Custom exception for SMILES building errors."""

    pass


class SMILESBuilder:
    """
    Builds SMILES strings from parsed complex data.

    Note:
        SMILES has inherent limitations for representing coordination
        compounds. This builder produces a valid SMILES that captures
        the molecular components but may not represent the true
        bonding topology or stereochemistry.

    Example:
        >>> builder = SMILESBuilder()
        >>> parsed = ParsedComplex(metal="Ir", ligands=[...], ...)
        >>> smiles = builder.build(parsed)
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the SMILES builder with chemical databases.

        Args:
            ligand_db: Ligand database. Uses LIGAND_DATABASE if None.
            metal_db: Metal database. Uses METAL_DATABASE if None.
            counter_ion_db: Counter ion database. Uses COUNTER_ION_DATABASE if None.
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

    def build(self, parsed: ParsedComplex) -> str:
        """
        Build a SMILES string from a parsed complex.

        Args:
            parsed: ParsedComplex object

        Returns:
            SMILES string representation

        Raises:
            SMILESBuilderError: If SMILES cannot be constructed
        """
        # Calculate metal oxidation state
        metal_charge = self._calculate_metal_charge(parsed)

        # Build metal center SMILES
        metal_smiles = self._format_metal_smiles(parsed.metal, metal_charge)

        # Build ligand SMILES
        ligand_smiles_list = self._build_ligand_smiles(parsed.ligands)

        # Combine into single complex unit
        complex_parts = [metal_smiles] + ligand_smiles_list
        single_unit_smiles = ".".join(complex_parts)

        # Handle multiplicity (dimers, etc.)
        if parsed.multiplicity > 1:
            full_smiles = ".".join([single_unit_smiles] * parsed.multiplicity)
        else:
            full_smiles = single_unit_smiles

        # Add counter ions
        for ion_name, ion_count in parsed.counter_ions:
            ion_smiles = self._get_counter_ion_smiles(ion_name)
            for _ in range(ion_count):
                full_smiles += "." + ion_smiles

        return full_smiles

    def _calculate_metal_charge(self, parsed: ParsedComplex) -> int:
        """
        Calculate the formal oxidation state of the metal center.

        Uses charge balance:
            metal_charge = complex_charge - sum(ligand_charges)

        Args:
            parsed: ParsedComplex object

        Returns:
            Calculated metal oxidation state
        """
        # Sum up ligand charges
        total_ligand_charge = 0
        for ligand in parsed.ligands:
            if ligand.name in self.ligand_db:
                lig_info = self.ligand_db[ligand.name]
                total_ligand_charge += lig_info.charge * ligand.count

        # Determine complex charge
        if parsed.counter_ions:
            # If counter ions present, calculate complex charge from them
            counter_ion_charge = 0
            for ion_name, ion_count in parsed.counter_ions:
                if ion_name in self.counter_ion_db:
                    counter_ion_charge += (
                        self.counter_ion_db[ion_name].charge * ion_count
                    )
            # Overall compound is neutral: complex_charge + counter_ion_charge = 0
            complex_charge = -counter_ion_charge
        else:
            complex_charge = parsed.complex_charge

        # Adjust for multiplicity
        complex_charge_per_unit = complex_charge // parsed.multiplicity

        # Calculate metal charge
        metal_charge = complex_charge_per_unit - total_ligand_charge

        return metal_charge

    def _format_metal_smiles(self, symbol: str, charge: int) -> str:
        """
        Format a metal symbol with charge as SMILES.

        Args:
            symbol: Metal element symbol
            charge: Formal charge

        Returns:
            SMILES representation (e.g., "[Ir+3]")
        """
        if charge == 0:
            return f"[{symbol}]"
        elif charge > 0:
            if charge == 1:
                return f"[{symbol}+]"
            else:
                return f"[{symbol}+{charge}]"
        else:  # charge < 0
            if charge == -1:
                return f"[{symbol}-]"
            else:
                return f"[{symbol}{charge}]"

    def _build_ligand_smiles(self, ligands: List[ParsedLigand]) -> List[str]:
        """
        Build SMILES strings for all ligands.

        Args:
            ligands: List of ParsedLigand objects

        Returns:
            List of SMILES strings
        """
        smiles_list: List[str] = []

        for ligand in ligands:
            lig_smiles = self._get_ligand_smiles(ligand.name)
            # Add SMILES for each copy of the ligand
            for _ in range(ligand.count):
                smiles_list.append(lig_smiles)

        return smiles_list

    def _get_ligand_smiles(self, name: str) -> str:
        """
        Get SMILES for a ligand by name.

        Attempts:
            1. Direct lookup in database
            2. Search aliases
            3. Match as suffix (for modified ligands)

        Args:
            name: Ligand name/abbreviation

        Returns:
            SMILES string

        Raises:
            SMILESBuilderError: If ligand not found
        """
        # Direct lookup
        if name in self.ligand_db:
            return self.ligand_db[name].smiles

        # Check aliases
        for lig_name, info in self.ligand_db.items():
            aliases_lower = [alias.lower() for alias in info.aliases]
            if name.lower() in aliases_lower:
                return info.smiles

        raise SMILESBuilderError(f"Unknown ligand: '{name}'")

    def _get_counter_ion_smiles(self, name: str) -> str:
        """
        Get SMILES for a counter ion by name.

        Args:
            name: Counter ion name

        Returns:
            SMILES string

        Raises:
            SMILESBuilderError: If counter ion not found
        """
        if name in self.counter_ion_db:
            return self.counter_ion_db[name].smiles

        # Check aliases
        for ion_name, info in self.counter_ion_db.items():
            if name in info.aliases:
                return info.smiles

        raise SMILESBuilderError(f"Unknown counter ion: '{name}'")


class InorganicNameToSMILES:
    """
    Main converter class for inorganic/organometallic names to SMILES.

    This class provides a high-level interface for converting coordination
    complex names to SMILES notation.

    Example:
        >>> converter = InorganicNameToSMILES()
        >>> smiles = converter.convert("[IrCl(cod)]2")
        >>> print(smiles)
        '[Ir+].[Cl-].C1=CCCC=CCC1.[Ir+].[Cl-].C1=CCCC=CCC1'
    """

    def __init__(
        self,
        ligand_db: Optional[Dict[str, LigandInfo]] = None,
        metal_db: Optional[Dict[str, MetalInfo]] = None,
        counter_ion_db: Optional[Dict[str, LigandInfo]] = None,
    ) -> None:
        """
        Initialize the converter with optional custom databases.

        Args:
            ligand_db: Custom ligand database (optional)
            metal_db: Custom metal database (optional)
            counter_ion_db: Custom counter ion database (optional)
        """
        self.ligand_db = ligand_db if ligand_db is not None else LIGAND_DATABASE
        self.metal_db = metal_db if metal_db is not None else METAL_DATABASE
        self.counter_ion_db = (
            counter_ion_db if counter_ion_db is not None else COUNTER_ION_DATABASE
        )

        self.parser = ComplexNameParser(
            self.ligand_db, self.metal_db, self.counter_ion_db
        )
        self.builder = SMILESBuilder(self.ligand_db, self.metal_db, self.counter_ion_db)

    def convert(self, name: str) -> str:
        """
        Convert an inorganic complex name to SMILES.

        Args:
            name: The complex name (e.g., "[IrCl(cod)]2")

        Returns:
            SMILES string representation

        Raises:
            ParserError: If the name cannot be parsed
            SMILESBuilderError: If SMILES cannot be constructed
        """
        parsed = self.parser.parse(name)
        return self.builder.build(parsed)

    def convert_with_details(self, name: str) -> Tuple[str, ParsedComplex]:
        """
        Convert and return both SMILES and parsed structure.

        Args:
            name: The complex name

        Returns:
            Tuple of (SMILES string, ParsedComplex object)
        """
        parsed = self.parser.parse(name)
        smiles = self.builder.build(parsed)
        return smiles, parsed

    def add_ligand(
        self,
        name: str,
        smiles: str,
        denticity: int = 1,
        charge: int = 0,
        aliases: Optional[Tuple[str, ...]] = None,
        description: str = "",
    ) -> None:
        """
        Add a new ligand to the database.

        Args:
            name: Ligand abbreviation (e.g., "dppe")
            smiles: SMILES representation
            denticity: Number of coordination sites
            charge: Formal charge of the ligand
            aliases: Alternative names
            description: Human-readable description
        """
        self.ligand_db[name] = LigandInfo(
            smiles=smiles,
            denticity=denticity,
            charge=charge,
            aliases=aliases if aliases else tuple(),
            description=description,
        )

    def add_counter_ion(
        self,
        name: str,
        smiles: str,
        charge: int = -1,
        aliases: Optional[Tuple[str, ...]] = None,
        description: str = "",
    ) -> None:
        """
        Add a new counter ion to the database.

        Args:
            name: Counter ion abbreviation
            smiles: SMILES representation
            charge: Formal charge (usually -1)
            aliases: Alternative names
            description: Human-readable description
        """
        self.counter_ion_db[name] = LigandInfo(
            smiles=smiles,
            charge=charge,
            aliases=aliases if aliases else tuple(),
            description=description,
        )

    def list_available_ligands(self) -> List[str]:
        """Return list of available ligand abbreviations."""
        return list(self.ligand_db.keys())

    def list_available_metals(self) -> List[str]:
        """Return list of available metal symbols."""
        return list(self.metal_db.keys())

    def list_available_counter_ions(self) -> List[str]:
        """Return list of available counter ion abbreviations."""
        return list(self.counter_ion_db.keys())


def name_to_smiles_inorganic_shorthand(
    names: List[str], strict: bool = True
) -> Dict[str, str]:
    """Convert multiple inorganic shorthand names to SMILES."""
    converter = InorganicNameToSMILES()
    name_smiles_dict = {}
    for name in names:
        try:
            smiles, _ = converter.convert_with_details(name)
        except Exception:
            continue
        if not smiles:
            continue
        name_smiles_dict[name] = smiles
    return name_smiles_dict
