from __future__ import annotations

from typing import Optional, List, Dict, Tuple, Set, Any
from dataclasses import dataclass, field
from enum import Enum, auto
from collections import defaultdict
from rdkit import Chem

from placeholder_name.utils.logging_config import logger


# =============================================================================
# CONSTANTS AND ELEMENT DATA
# =============================================================================

class BondOrder(Enum):
    """Enumeration of chemical bond orders."""
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 1.5
    
    def to_smiles_symbol(self) -> str:
        """Convert bond order to SMILES bond symbol."""
        symbols = {
            BondOrder.SINGLE: "",
            BondOrder.DOUBLE: "=",
            BondOrder.TRIPLE: "#",
            BondOrder.AROMATIC: ":"
        }
        return symbols.get(self, "")


# Two-letter elements (must check before single-letter)
TWO_LETTER_ELEMENTS: Set[str] = {
    'He', 'Li', 'Be', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'Cl', 'Ar', 'Ca', 'Sc',
    'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
    'Br', 'Kr', 'Rb', 'Sr', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
    'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
    'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf'
}

# Single-letter elements
ONE_LETTER_ELEMENTS: Set[str] = {
    'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U'
}

# Elements in the organic subset (don't need brackets in SMILES)
ORGANIC_SUBSET: Set[str] = {'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'}

# Standard valences for validation
STANDARD_VALENCES: Dict[str, List[int]] = {
    'H': [1], 'B': [3], 'C': [4], 'N': [3, 5], 'O': [2],
    'F': [1], 'Si': [4], 'P': [3, 5], 'S': [2, 4, 6],
    'Cl': [1], 'Br': [1], 'I': [1]
}


# =============================================================================
# KNOWN FRAGMENTS - Extensible fragment definitions for common structural units
# =============================================================================

@dataclass
class FragmentDefinition:
    """
    Defines a known molecular fragment that can be recognized in structural formulas.
    
    Attributes:
        pattern: The text pattern to match (e.g., 'C6H5')
        smiles: The SMILES representation of the fragment with attachment point(s)
        attachment_points: List of atom indices (0-indexed) where bonds can form
        is_aromatic: Whether this fragment contains aromatic atoms
        description: Human-readable description
    """
    pattern: str
    smiles: str
    attachment_points: List[int]  # Indices in the SMILES where attachments occur
    is_aromatic: bool = False
    description: str = ""


# Registry of known fragments - easily extensible
KNOWN_FRAGMENTS: Dict[str, FragmentDefinition] = {
    # Phenyl group (benzene with one H removed)
    'C6H5': FragmentDefinition(
        pattern='C6H5',
        smiles='c1ccccc1',
        attachment_points=[0],  # Attach at first carbon
        is_aromatic=True,
        description='phenyl'
    ),
    # 1,2-phenylene (ortho-disubstituted benzene)
    'C6H4': FragmentDefinition(
        pattern='C6H4',
        smiles='c1ccccc1',
        attachment_points=[0, 1],  # Two adjacent attachment points
        is_aromatic=True,
        description='1,2-phenylene'
    ),
    # Cyclopentadienyl
    'C5H5': FragmentDefinition(
        pattern='C5H5',
        smiles='c1cccc1',
        attachment_points=[0],
        is_aromatic=True,
        description='cyclopentadienyl'
    ),
    # Cyclohexyl
    'C6H11': FragmentDefinition(
        pattern='C6H11',
        smiles='C1CCCCC1',
        attachment_points=[0],
        is_aromatic=False,
        description='cyclohexyl'
    ),
    # Cyclopentyl
    'C5H9': FragmentDefinition(
        pattern='C5H9',
        smiles='C1CCCC1',
        attachment_points=[0],
        is_aromatic=False,
        description='cyclopentyl'
    ),
    # Naphthyl (1-position)
    'C10H7': FragmentDefinition(
        pattern='C10H7',
        smiles='c1ccc2ccccc2c1',
        attachment_points=[0],
        is_aromatic=True,
        description='1-naphthyl'
    ),
    # Furyl (2-position)
    'C4H3O': FragmentDefinition(
        pattern='C4H3O',
        smiles='c1ccoc1',
        attachment_points=[1],  # Attach at C2
        is_aromatic=True,
        description='2-furyl'
    ),
    # Thienyl (2-position)
    'C4H3S': FragmentDefinition(
        pattern='C4H3S',
        smiles='c1ccsc1',
        attachment_points=[1],  # Attach at C2
        is_aromatic=True,
        description='2-thienyl'
    ),
    # Pyridyl (various positions can be added)
    'C5H4N': FragmentDefinition(
        pattern='C5H4N',
        smiles='c1ccncc1',
        attachment_points=[0],  # Attach at C2 (adjacent to N)
        is_aromatic=True,
        description='pyridyl'
    ),
    # Carboxylic acid
    'COOH': FragmentDefinition(
        pattern='COOH',
        smiles='C(=O)O',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='carboxylic'
    ),
    # Carboxylic acid
    'HOOC': FragmentDefinition(
        pattern='HOOC',
        smiles='C(=O)O',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='carboxylic'
    ),
    # Primary Amide
    'CONH2': FragmentDefinition(
        pattern='CONH2',
        smiles='',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='primary amide'
    ),
    # Methyl ester
    'COOCH3': FragmentDefinition(
        pattern='COOCH3',
        smiles='',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='methyl ester'
    ),
    # Ethyl ester
    'COOCH2CH3': FragmentDefinition(
        pattern='COOCH2CH3',
        smiles='',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='ethyl ester'
    ),
    # Acyl Chloride
    'COCl': FragmentDefinition(
        pattern='COCl',
        smiles='',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='acyl chloride'
    ),
    # Acyl bromide
    'COBr': FragmentDefinition(
        pattern='COBr',
        smiles='',
        attachment_points=[0],  # Attach at C1
        is_aromatic=False,
        description='acyl bromide'
    ),
}


def register_fragment(pattern: str, smiles: str, attachment_points: List[int],
                       is_aromatic: bool = False, description: str = "") -> None:
    """
    Register a new fragment pattern for recognition.
    
    Args:
        pattern: The molecular formula pattern (e.g., 'C6H5')
        smiles: SMILES representation of the fragment
        attachment_points: List of atom indices where bonds can form
        is_aromatic: Whether the fragment is aromatic
        description: Human-readable name
    
    Example:
        register_fragment('C6H5', 'c1ccccc1', [0], is_aromatic=True, description='phenyl')
    """
    KNOWN_FRAGMENTS[pattern] = FragmentDefinition(
        pattern=pattern,
        smiles=smiles,
        attachment_points=attachment_points,
        is_aromatic=is_aromatic,
        description=description
    )


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Atom:
    """Represents an atom in the molecular structure."""
    element: str
    index: int
    implicit_h_count: int = 0
    charge: int = 0
    
    def __hash__(self) -> int:
        return hash(self.index)
    
    def __eq__(self, other: object) -> bool:
        return isinstance(other, Atom) and self.index == other.index


@dataclass 
class Bond:
    """Represents a bond between two atoms."""
    atom1_idx: int
    atom2_idx: int
    order: BondOrder = BondOrder.SINGLE


class MolecularGraph:
    """
    Represents a molecular structure as a graph.
    
    Provides methods for building molecules atom-by-atom and bond-by-bond,
    with efficient neighbor lookup for SMILES generation.
    """
    
    def __init__(self):
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self._adjacency: Dict[int, List[Tuple[int, BondOrder]]] = defaultdict(list)
    
    def add_atom(self, element: str, implicit_h: int = 0, charge: int = 0) -> int:
        """Add an atom and return its index."""
        idx = len(self.atoms)
        self.atoms.append(Atom(element, idx, implicit_h, charge))
        return idx
    
    def add_bond(self, idx1: int, idx2: int, order: BondOrder = BondOrder.SINGLE) -> bool:
        """Add a bond between two atoms. Returns False if bond already exists."""
        if idx1 == idx2:
            return False
        
        # Check for existing bond
        for neighbor, _ in self._adjacency[idx1]:
            if neighbor == idx2:
                return False
        
        self.bonds.append(Bond(idx1, idx2, order))
        self._adjacency[idx1].append((idx2, order))
        self._adjacency[idx2].append((idx1, order))
        return True
    
    def get_neighbors(self, idx: int) -> List[Tuple[int, BondOrder]]:
        """Get all neighbors of an atom with bond orders."""
        return list(self._adjacency[idx])
    
    def get_atom(self, idx: int) -> Atom:
        """Get atom by index."""
        return self.atoms[idx]
    
    def atom_count(self) -> int:
        """Return number of atoms."""
        return len(self.atoms)
    
    def is_empty(self) -> bool:
        """Check if graph has no atoms."""
        return len(self.atoms) == 0
    
    def copy(self) -> 'MolecularGraph':
        """Create a deep copy of this graph."""
        new_graph = MolecularGraph()
        for atom in self.atoms:
            new_graph.add_atom(atom.element, atom.implicit_h_count, atom.charge)
        for bond in self.bonds:
            new_graph.add_bond(bond.atom1_idx, bond.atom2_idx, bond.order)
        return new_graph


# =============================================================================
# TOKENIZER
# =============================================================================

class TokenType(Enum):
    """Types of tokens in structural formula parsing."""
    ELEMENT = auto()
    NUMBER = auto()
    LPAREN = auto()
    RPAREN = auto()
    DOUBLE_BOND = auto()
    TRIPLE_BOND = auto()
    FRAGMENT = auto()  # For known fragments like C6H5 (phenyl)
    END = auto()


@dataclass
class Token:
    """A lexical token."""
    type: TokenType
    value: str
    position: int


class Tokenizer:
    """
    Lexical analyzer for structural formulas.
    
    Converts a formula string into a sequence of tokens for parsing.
    """
    
    def __init__(self, formula: str):
        self.formula = formula
        self.pos = 0
        self.errors: List[str] = []
    
    def _try_match_fragment(self) -> Optional[str]:
        """
        Try to match a known fragment pattern at the current position.
        
        Returns the matched fragment pattern string, or None if no match.
        Fragment patterns are sorted by length (longest first) to ensure
        greedy matching (e.g., C10H7 before C6H5).
        """
        # Sort patterns by length (longest first) for greedy matching
        sorted_patterns = sorted(KNOWN_FRAGMENTS.keys(), key=len, reverse=True)
        
        for pattern in sorted_patterns:
            if self.formula[self.pos:].startswith(pattern):
                return pattern
        return None
    
    def tokenize(self) -> List[Token]:
        """Tokenize the formula into a list of tokens."""
        tokens = []
        
        while self.pos < len(self.formula):
            char = self.formula[self.pos]
            
            # Skip whitespace
            if char.isspace():
                self.pos += 1
                continue
            
            # Try to match known fragments first (e.g., C6H5, C10H7)
            # This must come before element matching to catch patterns like C6H5
            fragment = self._try_match_fragment()
            if fragment is not None:
                tokens.append(Token(TokenType.FRAGMENT, fragment, self.pos))
                self.pos += len(fragment)
                continue
            
            # Try two-letter element first
            if self.pos + 1 < len(self.formula):
                two_char = self.formula[self.pos:self.pos + 2]
                if two_char in TWO_LETTER_ELEMENTS:
                    tokens.append(Token(TokenType.ELEMENT, two_char, self.pos))
                    self.pos += 2
                    continue
            
            # Single-letter element
            if char in ONE_LETTER_ELEMENTS:
                tokens.append(Token(TokenType.ELEMENT, char, self.pos))
                self.pos += 1
                continue
            
            # Number
            if char.isdigit():
                start = self.pos
                num_str = ""
                while self.pos < len(self.formula) and self.formula[self.pos].isdigit():
                    num_str += self.formula[self.pos]
                    self.pos += 1
                tokens.append(Token(TokenType.NUMBER, num_str, start))
                continue
            
            # Parentheses
            if char == '(':
                tokens.append(Token(TokenType.LPAREN, '(', self.pos))
                self.pos += 1
                continue
            if char == ')':
                tokens.append(Token(TokenType.RPAREN, ')', self.pos))
                self.pos += 1
                continue
            
            # Bond symbols
            if char == '=':
                tokens.append(Token(TokenType.DOUBLE_BOND, '=', self.pos))
                self.pos += 1
                continue
            if char == '#':
                tokens.append(Token(TokenType.TRIPLE_BOND, '#', self.pos))
                self.pos += 1
                continue
            if char == '≡':
                tokens.append(Token(TokenType.TRIPLE_BOND, '≡', self.pos))
                self.pos += 1
                continue
            
            # Skip hyphens (often used as separators)
            if char == '-':
                self.pos += 1
                continue
            
            # Unknown character
            self.errors.append(f"Unknown character '{char}' at position {self.pos}")
            self.pos += 1
        
        tokens.append(Token(TokenType.END, "", self.pos))
        return tokens


# =============================================================================
# PARSER
# =============================================================================

class ParseError(Exception):
    """Raised when parsing fails."""
    pass


@dataclass
class AtomGroup:
    """
    Represents a parsed atom group (e.g., CH3, OH, Cl).
    
    Attributes:
        central_element: The main element (C, O, N, etc.)
        hydrogen_count: Number of attached hydrogens
        substituents: List of (element, count) for other attached atoms
    """
    central_element: str
    hydrogen_count: int = 0
    substituents: List[Tuple[str, int]] = field(default_factory=list)


@dataclass
class BranchGroup:
    """
    Represents a branch specification like (CH3)2 or (=O).
    
    Attributes:
        content: The content inside the parentheses (as parsed structure)
        multiplier: How many times this branch appears
        bond_order: Bond order for attachment (for things like =O)
    """
    content: Any  # Can be AtomGroup, list of groups, etc.
    multiplier: int = 1
    bond_order: BondOrder = BondOrder.SINGLE


class StructuralFormulaParser:
    """
    Parser for condensed structural formulas.
    
    This parser handles various notational conventions:
    - CH3, CH2, CH, C - carbon with hydrogens
    - (CH3)2, (CH3)3 - multiple identical branches
    - C(O), C(=O) - carbonyl notation
    - Inline branches: CH(CH3)2
    """
    
    def __init__(self, tokens: List[Token], formula: str = ""):
        self.tokens = tokens
        self.formula = formula
        self.pos = 0
        self.errors: List[str] = []
        self.graph = MolecularGraph()
    
    def parse(self) -> Optional[MolecularGraph]:
        """Parse tokens into a molecular graph."""
        try:
            self._parse_formula()
            
            if not self._at_end():
                self.errors.append(f"Unexpected tokens after position {self.pos}")
                return None
            
            if self.graph.is_empty():
                self.errors.append("Empty molecule")
                return None
            
            return self.graph
            
        except ParseError as e:
            self.errors.append(str(e))
            return None
        except Exception as e:
            self.errors.append(f"Parsing error: {e}")
            logger.exception("Unexpected parsing error")
            return None
    
    # -------------------------------------------------------------------------
    # Token handling utilities
    # -------------------------------------------------------------------------
    
    def _current(self) -> Token:
        """Get current token."""
        return self.tokens[min(self.pos, len(self.tokens) - 1)]
    
    def _peek(self, offset: int = 0) -> Token:
        """Look ahead at a token."""
        idx = min(self.pos + offset, len(self.tokens) - 1)
        return self.tokens[idx]
    
    def _advance(self) -> Token:
        """Consume and return current token."""
        token = self._current()
        if token.type != TokenType.END:
            self.pos += 1
        return token
    
    def _at_end(self) -> bool:
        """Check if at end of input."""
        return self._current().type == TokenType.END
    
    def _match(self, *types: TokenType) -> bool:
        """Check if current token matches any given type."""
        return self._current().type in types
    
    def _consume_number(self) -> int:
        """Consume a number token, return 1 if not present."""
        if self._match(TokenType.NUMBER):
            return int(self._advance().value)
        return 1
    
    # -------------------------------------------------------------------------
    # Main parsing logic
    # -------------------------------------------------------------------------
    
    def _parse_formula(self) -> None:
        """
        Parse the entire formula.
        
        A formula is a sequence of:
        - Prefix branches: (X)n before an atom
        - Atom groups: C, CH, CH2, CH3, OH, etc.
        - Inline branches: (X) after an atom
        - Bond modifiers: =, #
        """
        prev_atom_idx: Optional[int] = None
        pending_branches: List[BranchGroup] = []
        next_bond_order = BondOrder.SINGLE
        
        while not self._at_end():
            current = self._current()
            
            # Stop at closing paren (for recursive calls)
            if current.type == TokenType.RPAREN:
                break
            
            # Bond order modifiers
            if current.type == TokenType.DOUBLE_BOND:
                self._advance()
                next_bond_order = BondOrder.DOUBLE
                continue
            elif current.type == TokenType.TRIPLE_BOND:
                self._advance()
                next_bond_order = BondOrder.TRIPLE
                continue
            
            # Prefix parenthetical group: (X)n before an atom
            if current.type == TokenType.LPAREN:
                branch = self._parse_parenthetical_branch()
                if branch is not None:
                    pending_branches.append(branch)
                continue
            
            # Fragment (e.g., C6H5 for phenyl)
            if current.type == TokenType.FRAGMENT:
                fragment_pattern = self._advance().value
                attachment_idx = self._add_fragment_to_graph(fragment_pattern)
                
                if attachment_idx is None:
                    break
                
                # Connect to previous atom in chain
                if prev_atom_idx is not None:
                    bond_order_to_use = next_bond_order
                    if bond_order_to_use == BondOrder.SINGLE:
                        bond_order_to_use = self._infer_bond_order(prev_atom_idx, attachment_idx)
                    
                    self.graph.add_bond(prev_atom_idx, attachment_idx, bond_order_to_use)
                    next_bond_order = BondOrder.SINGLE
                
                # Attach any pending prefix branches to the fragment
                for branch in pending_branches:
                    self._attach_branch(attachment_idx, branch)
                pending_branches.clear()
                
                prev_atom_idx = attachment_idx
                
                # Parse any inline branches after the fragment
                while self._match(TokenType.LPAREN):
                    branch = self._parse_parenthetical_branch()
                    if branch is not None:
                        self._attach_branch(attachment_idx, branch)
                
                continue
            
            # Atom group
            if current.type == TokenType.ELEMENT:
                atom_idx = self._parse_atom_with_inline_branches(pending_branches)
                
                if atom_idx is None:
                    break
                
                # Connect to previous atom in chain
                if prev_atom_idx is not None:
                    # If no explicit bond order was set (i.e., next_bond_order is still SINGLE),
                    # try to infer bond order from valences
                    bond_order_to_use = next_bond_order
                    if bond_order_to_use == BondOrder.SINGLE:
                        bond_order_to_use = self._infer_bond_order(prev_atom_idx, atom_idx)

                    self.graph.add_bond(prev_atom_idx, atom_idx, bond_order_to_use)
                    next_bond_order = BondOrder.SINGLE
                
                # Attach any pending prefix branches
                for branch in pending_branches:
                    self._attach_branch(atom_idx, branch)
                pending_branches.clear()
                
                prev_atom_idx = atom_idx
                continue
            
            # If we get here, we don't know what to do with this token
            self.errors.append(f"Unexpected token: {current}")
            break
        
        # Any remaining pending branches are an error
        if pending_branches:
            self.errors.append("Prefix branches with no atom to attach to")

    def _infer_bond_order(self, atom1_idx: int, atom2_idx: int) -> BondOrder:
        """
        Infer the most reasonable bond order between two atoms based on their remaining valence.
        
        Assumes all hydrogens are explicit. Uses STANDARD_VALENCES to compute remaining valence.
        Prefers single > double > triple to avoid overbonding.
        Returns inferred bond order, defaulting to SINGLE if uncertain.
        """
        atom1 = self.graph.get_atom(atom1_idx)
        atom2 = self.graph.get_atom(atom2_idx)
        
        def get_remaining_valence(atom: Atom) -> int:
            element = atom.element
            if element not in STANDARD_VALENCES:
                return 0  # unknown — can't infer
            possible_valences = STANDARD_VALENCES[element]
            
            # Choose appropriate valence: prefer common organic valence
            if element == 'N':
                # Prefer 3 for neutral N unless already bonded beyond that
                target_valence = 3
                if atom.charge != 0:
                    target_valence = max(possible_valences)
            elif element in ('C', 'S', 'P'):
                target_valence = max(possible_valences)
            else:
                target_valence = min(possible_valences)  # O, F, Cl, etc.
            
            # Sum bond orders (in valence units) from existing bonds
            neighbors = self.graph.get_neighbors(atom.index)
            used_by_bonds = sum(self._bond_order_to_valence_units(bo) for _, bo in neighbors)
            
            # Hydrogens count as 1 each
            used_by_h = atom.implicit_h_count
            
            total_used = used_by_bonds + used_by_h
            remaining = target_valence - total_used
            return max(0, remaining)
        
        rem1 = get_remaining_valence(atom1)
        rem2 = get_remaining_valence(atom2)
        
        # Try highest bond order first that both atoms can support
        for order in [BondOrder.TRIPLE, BondOrder.DOUBLE, BondOrder.SINGLE]:
            cost = self._bond_order_to_valence_units(order)
            if rem1 >= cost and rem2 >= cost:
                return order
        
        return BondOrder.SINGLE  # fallback
    
    @staticmethod
    def _bond_order_to_valence_units(order: BondOrder) -> int:
        """
        Convert a bond order to the number of valence units it consumes.
        For aromatic bonds, conservatively treat as single (1 unit).
        """
        if order == BondOrder.AROMATIC:
            return 1
        return int(order.value)  # SINGLE=1, DOUBLE=2, TRIPLE=3
    
    def _parse_parenthetical_branch(self) -> Optional[BranchGroup]:
        """
        Parse a parenthetical expression like (CH3)2 or (=O).
        
        Returns a BranchGroup describing what was parsed.
        """
        self._advance()  # Consume '('
        
        # Check for bond order prefix inside parens: (=O), (=S)
        bond_order = BondOrder.SINGLE
        if self._match(TokenType.DOUBLE_BOND):
            self._advance()
            bond_order = BondOrder.DOUBLE
        elif self._match(TokenType.TRIPLE_BOND):
            self._advance()
            bond_order = BondOrder.TRIPLE
        
        # Parse the content
        content_atoms: List[int] = []
        content_start_idx = self.graph.atom_count()
        
        # Simple case: single element like (O), (Cl), (=O)
        if self._match(TokenType.ELEMENT):
            element = self._current().value
            next_tok = self._peek(1)
            
            # Check if this is a simple single-atom branch
            if next_tok.type == TokenType.RPAREN:
                self._advance()  # consume element
                self._advance()  # consume ')'
                multiplier = self._consume_number()
                
                # Special case: (O) typically means carbonyl =O
                if element == 'O' and bond_order == BondOrder.SINGLE:
                    bond_order = BondOrder.DOUBLE
                elif element == 'S' and bond_order == BondOrder.SINGLE:
                    bond_order = BondOrder.DOUBLE
                
                return BranchGroup(
                    content=AtomGroup(element, 0, []),
                    multiplier=multiplier,
                    bond_order=bond_order
                )
            
            # For elements C, N, O, S, try to parse as simple atom group
            # but only accept it if followed by RPAREN
            if element in ('C', 'N', 'O', 'S'):
                # Save current position in case we need to backtrack
                saved_pos = self.pos
                
                try:
                    # Parse the atom group inside
                    group = self._parse_simple_atom_group()
                    
                    # Only accept this as a simple group if immediately followed by RPAREN
                    if self._match(TokenType.RPAREN):
                        self._advance()  # consume ')'
                        multiplier = self._consume_number()
                        return BranchGroup(content=group, multiplier=multiplier, bond_order=bond_order)
                    else:
                        # Not followed by RPAREN, so this is not a simple group
                        # Reset position and fall through to complex parsing
                        self.pos = saved_pos
                except Exception:
                    # If parsing fails, reset and fall through
                    self.pos = saved_pos
        
        # More complex content - parse as a sub-formula
        # Save current graph state
        saved_count = self.graph.atom_count()
        
        # Parse inner content (recursive)
        inner_first_atom = None
        inner_prev_atom = None
        inner_bond_order = BondOrder.SINGLE
        
        while not self._at_end() and not self._match(TokenType.RPAREN):
            if self._match(TokenType.DOUBLE_BOND):
                self._advance()
                inner_bond_order = BondOrder.DOUBLE
                continue
            elif self._match(TokenType.TRIPLE_BOND):
                self._advance()
                inner_bond_order = BondOrder.TRIPLE
                continue
            
            if self._match(TokenType.ELEMENT):
                group = self._parse_simple_atom_group()
                atom_idx = self._add_atom_group_to_graph(group)
                
                if inner_first_atom is None:
                    inner_first_atom = atom_idx
                
                if inner_prev_atom is not None:
                    self.graph.add_bond(inner_prev_atom, atom_idx, inner_bond_order)
                    inner_bond_order = BondOrder.SINGLE
                
                inner_prev_atom = atom_idx
            else:
                break
        
        if not self._match(TokenType.RPAREN):
            self.errors.append("Unclosed parenthesis")
            return None
        
        self._advance()  # consume ')'
        multiplier = self._consume_number()
        
        # Return branch info
        return BranchGroup(
            content=inner_first_atom,  # First atom of the inner chain
            multiplier=multiplier,
            bond_order=bond_order
        )
    
    def _parse_simple_atom_group(self) -> AtomGroup:
        """
        Parse a simple atom group like C, CH, CH2, CH3, OH, NH2, Cl, etc.
        
        Handles:
        - Element alone: C, O, N, Cl
        - Element with H count: CH3, NH2, OH
        - Element with other substituents: CCl3, CHCl2
        """
        if not self._match(TokenType.ELEMENT):
            raise ParseError("Expected element")
        
        element = self._advance().value
        
        # Special case: if element is H, it might be "HC#CH" style
        # In this case, H is attached to what follows
        if element == 'H':
            # H by itself or H followed by count
            h_count = self._consume_number()
            return AtomGroup('H', 0, [])  # Explicit H atom
        
        hydrogen_count = 0
        substituents: List[Tuple[str, int]] = []
        
        # Parse attached atoms (H, halogens, etc.)
        while self._match(TokenType.ELEMENT):
            attached = self._current().value
            
            # Peek ahead: if next is number, it's a count
            if attached == 'H':
                self._advance()
                hydrogen_count = self._consume_number()
            elif attached in ('F', 'Cl', 'Br', 'I'):
                # Halogen - can have a count
                self._advance()
                count = self._consume_number()
                substituents.append((attached, count))
            else:
                # Different element - this is the start of the next group
                break
        
        return AtomGroup(element, hydrogen_count, substituents)
    
    def _add_atom_group_to_graph(self, group: AtomGroup) -> int:
        """
        Add an AtomGroup to the molecular graph.
        
        Returns the index of the central atom.
        """
        central_idx = self.graph.add_atom(group.central_element, group.hydrogen_count)
        
        # Add substituents
        for element, count in group.substituents:
            for _ in range(count):
                # First, add the substituent atom to get its index
                sub_idx = self.graph.add_atom(element, 0)
                # Now infer bond order between central atom and this new atom
                bond_order = self._infer_bond_order(central_idx, sub_idx)
                # Create the bond
                self.graph.add_bond(central_idx, sub_idx, bond_order)
        
        return central_idx
    
    def _parse_atom_with_inline_branches(self, prefix_branches: List[BranchGroup]) -> Optional[int]:
        """
        Parse an atom group with possible inline branches.
        
        E.g., "C(O)(CH3)2" or "CH(CH3)"
        
        Returns the index of the central atom.
        """
        group = self._parse_simple_atom_group()
        central_idx = self._add_atom_group_to_graph(group)
        
        # Parse inline branches: things like (O), (CH3), (=O) that come after the atom
        while self._match(TokenType.LPAREN):
            branch = self._parse_parenthetical_branch()
            if branch is not None:
                self._attach_branch(central_idx, branch)
        
        return central_idx
    
    def _attach_branch(self, attach_to: int, branch: BranchGroup) -> None:
        """
        Attach a branch to an atom.
        
        Handles multipliers by creating multiple copies.
        """
        for i in range(branch.multiplier):
            if isinstance(branch.content, AtomGroup):
                branch_idx = self._add_atom_group_to_graph(branch.content)
                self.graph.add_bond(attach_to, branch_idx, branch.bond_order)
            elif isinstance(branch.content, int):
                if i == 0:
                    self.graph.add_bond(attach_to, branch.content, branch.bond_order)
                else:
                    if branch.content < len(self.graph.atoms):
                        atom = self.graph.atoms[branch.content]
                        new_idx = self.graph.add_atom(atom.element, atom.implicit_h_count)
                        self.graph.add_bond(attach_to, new_idx, branch.bond_order)
    
    def _add_fragment_to_graph(self, fragment_pattern: str) -> Optional[int]:
        """
        Add a known fragment (e.g., C6H5 phenyl) to the molecular graph.
        Builds the fragment directly without requiring RDKit.
        Returns the index of the attachment point atom, or None on failure.
        """
        print(f"Adding fragment: {fragment_pattern}")
        if fragment_pattern not in KNOWN_FRAGMENTS:
            self.errors.append(f"Unknown fragment pattern: {fragment_pattern}")
            return None
        
        fragment_def = KNOWN_FRAGMENTS[fragment_pattern]
        
        # Use predefined fragment builders for common patterns
        if fragment_pattern == 'C6H5':
            # Phenyl ring: 6 aromatic carbons in a ring
            return self._build_phenyl_ring()
        elif fragment_pattern == 'C6H11':
            # Cyclohexyl: 6 sp3 carbons in a ring
            return self._build_cyclohexyl_ring()
        elif fragment_pattern == 'C5H9':
            # Cyclopentyl: 5 sp3 carbons in a ring
            return self._build_cycloalkyl_ring(5)
        elif fragment_pattern in ['COOH', 'HOOC']:
            # Carboxylic acid
            return self._build_carboxylic_acid()
        # elif fragment_pattern == 'CONH2':
        #     # Primary amide
        #     return self._build_primary_amide()
        elif fragment_pattern == 'COOCH3':
            # Methyl ester
            return self._build_ester(n_carbon=1)
        elif fragment_pattern == 'COOCH2CH3':
            # Ethyl ester
            return self._build_ester(n_carbon=2)
        elif fragment_pattern == 'COCl':
            # Acyl chloride
            return self._build_acid_halide(halide='Cl')
        elif fragment_pattern == 'COBr':
            # Acyl bromide
            return self._build_acid_halide(halide='Br')
        elif fragment_pattern == 'COI':
            # Acyl iodide
            return self._build_acid_halide(halide='I')
        else:
            # For other fragments, try RDKit if available
            return self._build_fragment_with_rdkit(fragment_def)
    
    def _build_phenyl_ring(self) -> int:
        """Build a phenyl (benzene) ring and return the attachment point index."""
        # Add 6 aromatic carbons
        indices = []
        for _ in range(6):
            idx = self.graph.add_atom('C', 0)
            indices.append(idx)
        
        # Connect them in a ring with aromatic bonds
        for i in range(6):
            next_i = (i + 1) % 6
            self.graph.add_bond(indices[i], indices[next_i], BondOrder.AROMATIC)
        
        # Return the first carbon as the attachment point
        return indices[0]
    
    def _build_cyclohexyl_ring(self) -> int:
        """Build a cyclohexyl ring and return the attachment point index."""
        return self._build_cycloalkyl_ring(6)
    
    def _build_cycloalkyl_ring(self, size: int) -> int:
        """Build a cycloalkyl ring of given size and return the attachment point."""
        indices = []
        for _ in range(size):
            idx = self.graph.add_atom('C', 0)
            indices.append(idx)
        
        # Connect them in a ring with single bonds
        for i in range(size):
            next_i = (i + 1) % size
            self.graph.add_bond(indices[i], indices[next_i], BondOrder.SINGLE)
        
        return indices[0]
    
    def _build_carboxylic_acid(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom('C', 0)
        idx2 = self.graph.add_atom('O', 0)
        idx3 = self.graph.add_atom('O', 1)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1
    
    def _build_amide(self) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom('C', 0)
        idx2 = self.graph.add_atom('O', 0)
        idx3 = self.graph.add_atom('O', 1)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1
    
    def _build_ester(self, n_carbon=1) -> int:
        """Build a carboxylic acid functional group"""
        idx1 = self.graph.add_atom('C', 0)
        idx2 = self.graph.add_atom('O', 0) 
        idx3 = self.graph.add_atom('O', 0) 

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)
        
        prev_idx = idx3
        for _ in range(n_carbon):
            next_idx = self.graph.add_atom('C', 2)
            self.graph.add_bond(prev_idx, next_idx, BondOrder.SINGLE)
            prev_idx = next_idx

        return idx1
    
    def _build_acid_halide(self, halide='Cl') -> int:
        """Build an ester functional group"""
        idx1 = self.graph.add_atom('C', 0)
        idx2 = self.graph.add_atom('O', 0)
        idx3 = self.graph.add_atom(halide, 0)

        self.graph.add_bond(idx1, idx2, BondOrder.DOUBLE)
        self.graph.add_bond(idx1, idx3, BondOrder.SINGLE)

        return idx1
    
    def _build_fragment_with_rdkit(self, fragment_def: FragmentDefinition) -> Optional[int]:
        """Try to build a fragment using RDKit if available."""
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(fragment_def.smiles)
            if mol is None:
                self.errors.append(f"Invalid fragment SMILES: {fragment_def.smiles}")
                return None
            
            rdkit_to_graph: Dict[int, int] = {}
            
            for atom in mol.GetAtoms():
                rdkit_idx = atom.GetIdx()
                element = atom.GetSymbol()
                graph_idx = self.graph.add_atom(element, 0)
                rdkit_to_graph[rdkit_idx] = graph_idx
            
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType()
                if bond_type == Chem.BondType.AROMATIC:
                    order = BondOrder.AROMATIC
                elif bond_type == Chem.BondType.DOUBLE:
                    order = BondOrder.DOUBLE
                elif bond_type == Chem.BondType.TRIPLE:
                    order = BondOrder.TRIPLE
                else:
                    order = BondOrder.SINGLE
                self.graph.add_bond(rdkit_to_graph[begin_idx], rdkit_to_graph[end_idx], order)
            
            if fragment_def.attachment_points:
                return rdkit_to_graph.get(fragment_def.attachment_points[0])
            return rdkit_to_graph.get(0)
            
        except ImportError:
            self.errors.append(f"Fragment {fragment_def.pattern} requires RDKit")
            return None
        except Exception as e:
            self.errors.append(f"Error adding fragment {fragment_def.pattern}: {e}")
            return None


# =============================================================================
# SMILES GENERATOR
# =============================================================================

class SMILESGenerator:
    """
    Generates SMILES strings from molecular graphs using DFS traversal.
    
    Handles ring closures by:
    1. First pass: DFS to identify all ring-closing (back) edges
    2. Second pass: DFS to generate SMILES with ring closure numbers
    """
    
    def __init__(self, graph: MolecularGraph):
        self.graph = graph
        self.visited: Set[int] = set()
        self.ring_closures: Dict[Tuple[int, int], Tuple[int, BondOrder]] = {}  # (a1, a2) -> (ring_num, bond_order)
        self.atom_ring_closures: Dict[int, List[Tuple[int, BondOrder]]] = defaultdict(list)  # atom_idx -> [(ring_num, bond_order), ...]
        self.next_ring_num: int = 1
    
    def generate(self) -> str:
        """Generate SMILES string."""
        if self.graph.is_empty():
            return ""
        
        # Reset state
        self.visited.clear()
        self.ring_closures.clear()
        self.atom_ring_closures.clear()
        self.next_ring_num = 1
        
        # Find best starting atom
        start = self._find_start_atom()
        
        # Pass 1: Find all ring-closing bonds
        self._find_ring_closures(start, -1)
        
        # Pass 2: Generate SMILES with ring numbers
        self.visited.clear()
        return self._dfs(start, -1, BondOrder.SINGLE)
    
    def _find_start_atom(self) -> int:
        """Find a good starting atom for SMILES generation."""
        # Prefer terminal atoms (degree 1) or failing that, first atom
        for i, atom in enumerate(self.graph.atoms):
            neighbors = self.graph.get_neighbors(i)
            if len(neighbors) == 1:
                return i
        return 0
    
    def _find_ring_closures(self, atom_idx: int, parent_idx: int) -> None:
        """
        First DFS pass to identify ring-closing bonds.
        
        A ring-closing bond is a "back edge" - a bond to an already-visited
        atom that isn't the parent.
        """
        self.visited.add(atom_idx)
        
        for neighbor_idx, bond_order in self.graph.get_neighbors(atom_idx):
            if neighbor_idx == parent_idx:
                continue
            
            if neighbor_idx in self.visited:
                # Back edge found - this is a ring closure
                # Use canonical key (smaller index first) to avoid duplicates
                key = (min(atom_idx, neighbor_idx), max(atom_idx, neighbor_idx))
                if key not in self.ring_closures:
                    ring_num = self.next_ring_num
                    self.next_ring_num += 1
                    self.ring_closures[key] = (ring_num, bond_order)
                    # Record that both atoms need this ring closure number
                    self.atom_ring_closures[atom_idx].append((ring_num, bond_order))
                    self.atom_ring_closures[neighbor_idx].append((ring_num, bond_order))
            else:
                # Tree edge - continue DFS
                self._find_ring_closures(neighbor_idx, atom_idx)
    
    def _dfs(self, atom_idx: int, parent_idx: int, incoming_bond: BondOrder) -> str:
        """Second DFS pass to build SMILES string with ring closures."""
        if atom_idx in self.visited:
            return ""
        
        self.visited.add(atom_idx)
        atom = self.graph.get_atom(atom_idx)
        
        result = ""
        
        # Add bond symbol (except for first atom or single bonds)
        if parent_idx >= 0 and incoming_bond not in [BondOrder.SINGLE, BondOrder.AROMATIC]:
            result += incoming_bond.to_smiles_symbol()
        
        aromatic_symbol = False
        if incoming_bond == BondOrder.AROMATIC:
            aromatic_symbol = True
        # Add atom symbol
        result += self._format_atom(atom, aromatic_symbol)
        
        # Add ring closure numbers for this atom
        for ring_num, ring_bond_order in self.atom_ring_closures.get(atom_idx, []):
            # For non-single bonds in ring closures, add bond symbol before ring number
            # (only needed on one side, we add it on the first encounter)
            if ring_bond_order != BondOrder.SINGLE and ring_bond_order != BondOrder.AROMATIC:
                # Check if this is the "opening" of the ring (first time we see this ring num)
                # We add the bond symbol only once
                pass  # For simplicity, aromatic rings work without explicit bond symbols
            
            if ring_num < 10:
                result += str(ring_num)
            else:
                result += f"%{ring_num}"
        
        # Get neighbors for traversal (exclude parent and ring-closure edges)
        neighbors = []
        for neighbor_idx, bond_order in self.graph.get_neighbors(atom_idx):
            if neighbor_idx == parent_idx:
                continue
            if neighbor_idx in self.visited:
                # This is a ring-closure edge, already handled via ring numbers
                continue
            neighbors.append((neighbor_idx, bond_order))
        
        if not neighbors:
            return result
        
        # Sort neighbors to get consistent output
        # Put hydrogens and halogens last as they're typically terminal
        def neighbor_priority(item):
            n, bo = item
            elem = self.graph.get_atom(n).element
            if elem == 'H':
                return (2, elem)
            if elem in ('F', 'Cl', 'Br', 'I'):
                return (1, elem)
            return (0, elem)
        
        neighbors.sort(key=neighbor_priority)
        
        # All but last go in branches
        for neighbor_idx, bond_order in neighbors[:-1]:
            branch = self._dfs(neighbor_idx, atom_idx, bond_order)
            if branch:
                result += f"({branch})"
        
        # Last continues main chain
        last_neighbor, last_bond = neighbors[-1]
        result += self._dfs(last_neighbor, atom_idx, last_bond)
        
        return result
    
    def _format_atom(self, atom: Atom, is_aromatic: bool = False) -> str:
        """Format an atom for SMILES output."""
        element = atom.element

        if is_aromatic:
            element = element.lower()
        
        # Check if we need brackets
        needs_brackets = (
            element.upper() not in ORGANIC_SUBSET or
            atom.charge != 0 or
            (element == 'H' and atom.implicit_h_count == 0)  # Explicit H
        )
        
        if needs_brackets:
            result = f"[{element}"
            if atom.charge > 0:
                result += "+" if atom.charge == 1 else f"+{atom.charge}"
            elif atom.charge < 0:
                result += "-" if atom.charge == -1 else str(atom.charge)
            result += "]"
            return result
        
        return element


# =============================================================================
# MAIN CONVERTER
# =============================================================================

class StructuralFormulaConverter:
    """
    Main interface for converting structural formulas to SMILES.
    
    This converter combines:
    1. Pattern matching for known compounds (most reliable)
    2. Parsing-based conversion for general cases
    3. Validation to ensure results are reasonable
    
    Usage:
        converter = StructuralFormulaConverter()
        smiles = converter.convert("CH3CH2OH")  # Returns "CCO"
        
        # Check for errors
        if smiles is None:
            print(converter.get_errors())
    """
    
    def __init__(self, strict_mode: bool = True, use_registry: bool = True):
        """
        Initialize the converter.
        
        Args:
            strict_mode: If True, return None for ambiguous/uncertain results
            use_registry: If True, check pattern registry first
        """
        self.strict_mode = strict_mode
        self.use_registry = use_registry
        self.last_errors: List[str] = []
        self._rdkit_available = self._check_rdkit()
    
    def _check_rdkit(self) -> bool:
        """Check if RDKit is available for validation."""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False
    
    def convert(self, formula: str) -> Optional[str]:
        """
        Convert a structural formula to SMILES.
        
        Args:
            formula: Condensed structural formula
            
        Returns:
            SMILES string if successful, None otherwise
        """
        self.last_errors = []
        
        # Preprocess
        formula = self._preprocess(formula)
        if not formula:
            self.last_errors.append("Empty formula")
            return ''
        
        # Tokenize
        tokenizer = Tokenizer(formula)
        tokens = tokenizer.tokenize()
        
        if tokenizer.errors:
            self.last_errors.extend(tokenizer.errors)
            if self.strict_mode:
                return ''
        
        # Parse
        parser = StructuralFormulaParser(tokens, formula)
        graph = parser.parse()
        
        if parser.errors:
            self.last_errors.extend(parser.errors)
        
        if graph is None:
            return ''
        
        # Generate SMILES
        generator = SMILESGenerator(graph)
        smiles = generator.generate()
        
        if not smiles:
            self.last_errors.append("Failed to generate SMILES")
            return ''
        
        # Validate with RDKit if available
        if self._rdkit_available and self.strict_mode:
            if not self._validate_smiles(smiles):
                self.last_errors.append(f"Generated SMILES '{smiles}' failed RDKit validation")
                return ''
        
        return smiles
    
    def _preprocess(self, formula: str) -> str:
        """Clean up the input formula."""
        formula = formula.strip()
        # Normalize some common variations
        formula = formula.replace('−', '-')  # Unicode minus to hyphen
        formula = formula.replace('–', '-')  # En-dash to hyphen
        return formula
    
    def _validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES using RDKit."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            return False
    
    def get_errors(self) -> List[str]:
        """Get errors from last conversion."""
        return list(self.last_errors)
    
    def register_pattern(self, formula: str, smiles: str) -> None:
        """Register a custom formula pattern."""
        self.registry.register(formula, smiles)
    
    def batch_convert(self, formulas: List[str]) -> Dict[str, Optional[str]]:
        """Convert multiple formulas."""
        return {f: self.convert(f) for f in formulas}


def name_to_smiles_structural_formula(formulas: List[str], strict: bool = True) -> Dict[str, Optional[str]]:
    """Convert multiple formulas to SMILES."""
    converter = StructuralFormulaConverter(strict_mode=strict)
    return converter.batch_convert(formulas)


# =============================================================================
# TESTING
# =============================================================================

def run_tests():
    """Run test suite."""
    test_cases = [
        # (formula, name, expected_smiles)
        ("(CH3)2CHCH2OH", "isobutanol", "CC(C)CO"),
        ("CH3CH2OH", "ethanol", "CCO"),
        ("CH3C(O)CH3", "acetone", "CC(=O)C"),
        ("CH3OH", "methanol", "CO"),
        ("CH4", "methane", "C"),
        ("CH3CH2CH3", "propane", "CCC"),
        ("CH2=CH2", "ethene", "C=C"),
        ("HC#CH", "ethyne", "C#C"),
        ("CH3COOH", "acetic acid", "CC(=O)O"),
        ("(CH3)3COH", "tert-butanol", "CC(C)(C)O"),
        ("CH3NH2", "methylamine", "CN"),
        ("CH3Cl", "chloromethane", "CCl"),
        ("CHCl3", "chloroform", "ClC(Cl)Cl"),
        ("(CH3)2O", "dimethyl ether", "COC"),
        ("CCl4", "carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
        ("(CH3)2NH", "dimethylamine", "CNC"),
        ("(CH3)3N", "trimethylamine", "CN(C)C"),
        ("CH3CHO", "acetaldehyde", "CC=O"),
        ("HCOOH", "formic acid", "C(=O)O"),
        ("CH3CN", "acetonitrile", "CC#N"),
        
        ("(CH3)3CCH2OH", "neopentyl alcohol", "CC(C)(C)CO"),
        ("(CH3)2CHCH2CH2OH", "isopentyl alcohol", "CC(C)CCO"),
        # ("(CH3CH2)3COH", "3-ethyl-3-pentanol", "CCC(O)(CC)CC"),
        ("(CH3)3COCH3", "methyl tert-butyl ether", "COC(C)(C)C"),
        ("CH3OCH2CH2OCH3", "1,2-dimethoxyethane", "COCCOC"),
        # ("(C2H5)2O", "diethyl ether", "CCOCC"),
        ("HOCH2CH(OH)CH2OH", "glycerol", "OCC(O)CO"),
        
        # ("C6H5CH2OH", "benzyl alcohol", "OCc1ccccc1"),
        ("C6H5NH2", "aniline", "Nc1ccccc1"),
        # ("C6H5NO2", "nitrobenzene", "[O-][N+](=O)c1ccccc1"),
        # ("C6H5CHO", "benzaldehyde", "O=Cc1ccccc1"),
        # ("C6H5COCH3", "acetophenone", "CC(=O)c1ccccc1"),
        # ("C6H5OCH3", "anisole", "COc1ccccc1"),
        # ("C6H5CH=CH2", "styrene", "C=Cc1ccccc1"),
        # ("C6H5C≡CH", "phenylacetylene", "C#Cc1ccccc1"),
        # ("C6H5CN", "benzonitrile", "N#Cc1ccccc1"),
        ("C6H5COOH", "benzoic acid", "OC(=O)c1ccccc1"),
        # ("CH3C6H4OH", "p-cresol", "Cc1ccc(O)cc1"),
        ("C6H5C(CH3)3", "tert-butylbenzene", "CC(C)(C)c1ccccc1"),
        # ("(C6H5)2CH2", "diphenylmethane", "c1ccc(Cc2ccccc2)cc1"),
        ("C6H4(OH)2", "catechol", "Oc1ccccc1O"),
        # ("O2NC6H4NH2", "4-nitroaniline", "Nc1ccc([N+]([O-])=O)cc1"),
        
        # ("CHCl3", "chloroform", "ClC(Cl)Cl"),
        # ("CCl4", "carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
        # ("CF3COOH", "trifluoroacetic acid", "OC(=O)C(F)(F)F"),
        # ("BrCH2CH2Br", "1,2-dibromoethane", "BrCCBr"),
        # ("(CH3)3CBr", "tert-butyl bromide", "CC(C)(C)Br"),
        # ("ClCH2COOH", "chloroacetic acid", "OC(=O)CCl"),
        # ("CHBr3", "bromoform", "BrC(Br)Br"),
        # ("CF2Cl2", "dichlorodifluoromethane", "FC(F)(Cl)Cl"),
        
        # ("CH3CH(OH)COOH", "lactic acid", "CC(O)C(=O)O"),
        # ("(CH3)3CCOOH", "pivalic acid", "CC(C)(C)C(=O)O"),
        # ("CH3COOC2H5", "ethyl acetate", "CCOC(C)=O"),
        # ("C6H5COOCH3", "methyl benzoate", "COC(=O)c1ccccc1"),
        # ("(CH3CO)2O", "acetic anhydride", "CC(=O)OC(C)=O"),
        # ("CH3COCl", "acetyl chloride", "CC(=O)Cl"),
        # ("CH3CONH2", "acetamide", "CC(N)=O"),
        # ("C6H5CONH2", "benzamide", "NC(=O)c1ccccc1"),
        # ("HOOCCH2CH2COOH", "succinic acid", "OC(=O)CCC(=O)O"),
        # ("HOOC(CH2)4COOH", "adipic acid", "OC(=O)CCCCC(=O)O"),
        
        # ("(C2H5)3N", "triethylamine", "CCN(CC)CC"),
        # ("(CH3)2NCHO", "N,N-dimethylformamide", "CN(C)C=O"),
        # ("HOCH2CH2NH2", "ethanolamine", "NCCO"),
        # ("H2NCH2CH2NH2", "ethylenediamine", "NCCN"),
        # ("(HOCH2CH2)3N", "triethanolamine", "OCCN(CCO)CCO"),
        # ("CH2=CHCN", "acrylonitrile", "C=CC#N"),
        # ("(CH3)2NNH2", "1,1-dimethylhydrazine", "CN(C)N"),
        # ("C6H5NHNH2", "phenylhydrazine", "NNc1ccccc1"),
        # ("(CH3)2NCH2CH2N(CH3)2", "TMEDA", "CN(C)CCN(C)C"),
        
        # ("(CH3)2SO", "dimethyl sulfoxide", "CS(C)=O"),
        # ("(CH3)2SO2", "dimethyl sulfone", "CS(C)(=O)=O"),
        # ("CH3SH", "methanethiol", "CS"),
        # ("CH3SSCH3", "dimethyl disulfide", "CSSC"),
        # ("C6H5SH", "thiophenol", "Sc1ccccc1"),
        # ("CH3SCH2CH2OH", "2-(methylthio)ethanol", "CSCCO"),
        # ("(C2H5)2S", "diethyl sulfide", "CCSCC"),
        
        # ("C5H5N", "pyridine", "c1ccncc1"),
        # ("C4H4O", "furan", "c1ccoc1"),
        # ("C4H4S", "thiophene", "c1ccsc1"),
        # ("C4H5N", "pyrrole", "c1cc[nH]c1"),
        # ("C4H8O", "tetrahydrofuran", "C1CCOC1"),
        # ("C3H4N2", "imidazole", "c1c[nH]cn1"),
        # ("C3H4N2", "pyrazole", "c1cc[nH]n1"),
        # ("C2H4O2S", "1,3-oxathiolane", "C1OCSO1"),
        # ("C3H6OS", "1,3-oxathiane", "C1COCSC1"),
        
        # ("C4H8O2", "1,4-dioxane", "C1COCCO1"),
        # ("C5H10O", "tetrahydropyran", "C1CCOCC1"),
        # ("C5H11N", "piperidine", "C1CCNCC1"),
        # ("C4H10N2", "piperazine", "C1CNCCN1"),
        # ("C4H9NO", "morpholine", "C1COCCN1"),
        # ("C4H4N2", "pyrimidine", "c1cncnc1"),
        # ("C4H4N2", "pyrazine", "c1cnccn1"),
        # ("C5H4N4", "purine", "c1ncc2[nH]cnc2n1"),
        
        # ("C8H7N", "indole", "c1ccc2[nH]ccc2c1"),
        # ("C9H7N", "quinoline", "c1ccc2ncccc2c1"),
        # ("C9H7N", "isoquinoline", "c1ccc2cnccc2c1"),
        # ("C10H8", "naphthalene", "c1ccc2ccccc2c1"),
        # ("C14H10", "anthracene", "c1ccc2cc3ccccc3cc2c1"),
        # ("C14H10", "phenanthrene", "c1ccc2c(c1)ccc1ccccc12"),
        # ("C12H8N2", "1,10-phenanthroline", "c1cnc2c(c1)ccc1cccnc12"),
        # ("C8H6O", "benzofuran", "c1ccc2occc2c1"),
        # ("C8H6S", "benzothiophene", "c1ccc2sccc2c1"),
        # ("C7H5NO", "benzoxazole", "c1ccc2ocnc2c1"),
        # ("C7H5NS", "benzothiazole", "c1ccc2scnc2c1"),
        
        # ("C10H18", "decalin", "C1CCC2CCCCC2C1"),
        # ("C7H12", "norbornane", "C1CC2CCC1C2"),
        # ("C10H16", "adamantane", "C1C2CC3CC1CC(C2)C3"),
        # ("C8H14", "bicyclo[2.2.2]octane", "C1CC2CCC1CC2"),
        # ("C6H12O", "cyclohexanol", "OC1CCCCC1"),
        # ("C6H10O", "cyclohexanone", "O=C1CCCCC1"),
        # ("C7H12O", "cycloheptanone", "O=C1CCCCCC1"),
        
        # ("(CH3)4Si", "tetramethylsilane", "C[Si](C)(C)C"),
        # ("(CH3)3SiCl", "trimethylsilyl chloride", "C[Si](C)(C)Cl"),
        # ("(CH3)3SiOSi(CH3)3", "hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C"),
        # ("(C2H5O)4Si", "tetraethyl orthosilicate", "CCO[Si](OCC)(OCC)OCC"),
        # ("C6H5Si(CH3)3", "trimethylphenylsilane", "C[Si](C)(C)c1ccccc1"),
        
        # ("(CH3O)3PO", "trimethyl phosphate", "COP(=O)(OC)OC"),
        # ("(C2H5O)3PO", "triethyl phosphate", "CCOP(=O)(OCC)OCC"),
        # ("(C6H5)3P", "triphenylphosphine", "c1ccc(P(c2ccccc2)c2ccccc2)cc1"),
        # ("(C6H5)3PO", "triphenylphosphine oxide", "O=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
        
        # ("B(OH)3", "boric acid", "OB(O)O"),
        # ("C6H5B(OH)2", "phenylboronic acid", "OB(O)c1ccccc1"),
        # ("(CH3O)3B", "trimethyl borate", "COB(OC)OC"),
        
        # ("(CH3)2CHCHO", "isobutyraldehyde", "CC(C)C=O"),
        # ("(CH3)2CHCOCH3", "3-methyl-2-butanone", "CC(C)C(C)=O"),
        # ("CH3COCH2COCH3", "acetylacetone", "CC(=O)CC(C)=O"),
        # ("CH3COCH2CH(CH3)2", "4-methyl-2-pentanone", "CC(C)CC(C)=O"),
        # ("C6H5COCH2CH3", "propiophenone", "CCC(=O)c1ccccc1"),
        # ("C6H5COC6H5", "benzophenone", "O=C(c1ccccc1)c1ccccc1"),
        # ("(CH3)2C=CHCOCH3", "mesityl oxide", "CC(C)=CC(C)=O"),
        
        # ("CH2=CHCH2OH", "allyl alcohol", "OCC=C"),
        # ("CH3C≡CH", "propyne", "CC#C"),
        # ("CH2=C(CH3)COOCH3", "methyl methacrylate", "COC(=O)C(C)=C"),
        # ("CH2=CHCOOH", "acrylic acid", "OC(=O)C=C"),
        # ("CH2=CHCOOCH3", "methyl acrylate", "COC(=O)C=C"),
        # ("(CH3)2C=C(CH3)2", "2,3-dimethyl-2-butene", "CC(C)=C(C)C"),
        # ("CH2=CHCH=CH2", "1,3-butadiene", "C=CC=C"),
        
        # ("H2NCH2COOH", "glycine", "NCC(=O)O"),
        # ("CH3CH(NH2)COOH", "alanine", "CC(N)C(=O)O"),
        # ("C6H5CH2CH(NH2)COOH", "phenylalanine", "NC(Cc1ccccc1)C(=O)O"),
        # ("HSCH2CH(NH2)COOH", "cysteine", "NC(CS)C(=O)O"),
        # ("(CH3)2CHCH(NH2)COOH", "valine", "CC(C)C(N)C(=O)O"),
        
        # ("HOCH2CH2OCH2CH2OH", "diethylene glycol", "OCCOCCO"),
        # ("CH3OCH2CH2OCH2CH2OCH3", "diglyme", "COCCOCCOC"),
        # ("C6H5CH(OH)CH3", "1-phenylethanol", "CC(O)c1ccccc1"),
        # ("(CH3)2C(OH)C≡CH", "2-methyl-3-butyn-2-ol", "CC(C)(O)C#C"),
        # ("NCCH2CH2CN", "succinonitrile", "N#CCCC#N"),
        # ("CH3COCH2CH2COCH3", "acetonylacetone", "CC(=O)CCC(C)=O"),
        # ("HOCH2C(CH2OH)3", "pentaerythritol", "OCC(CO)(CO)CO"),
    ]
    
    converter = StructuralFormulaConverter(strict_mode=True)
    
    print("=" * 75)
    print("Structural Formula to SMILES Converter - Test Results")
    print("=" * 75)
    print()
    print(f"{'Formula':<20} {'Name':<18} {'Result':<18} {'Expected':<15} {'Status'}")
    print("-" * 75)
    
    passed = 0
    failed = 0
    
    for formula, name, expected in test_cases:
        result = converter.convert(formula)

        try:
            canonical_result = Chem.MolToSmiles(Chem.MolFromSmiles(result), canonical=True)
        except:
            canonical_result = '1'

        try:
            canonical_expected = Chem.MolToSmiles(Chem.MolFromSmiles(expected), canonical=True)
        except:
            canonical_expected = '2'
        
        if canonical_result == canonical_expected:
            status = "✓ PASS"
            passed += 1
        elif result is None:
            status = "✗ FAIL (None)"
            failed += 1
        else:
            # Check if result is equivalent (same molecule, different SMILES)
            status = f"? DIFF"
            failed += 1
        
        result_str = canonical_result if canonical_result else "None"
        print(f"{formula:<20} {name:<18} {result_str:<18} {canonical_expected:<15} {status}")
        
        if converter.last_errors:
            for err in converter.last_errors:
                print(f"    Error: {err}")
    
    print("-" * 75)
    print(f"Results: {passed}/{len(test_cases)} passed ({100*passed/len(test_cases):.1f}%)")
    print("=" * 75)
    
    return passed == len(test_cases)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run_tests()