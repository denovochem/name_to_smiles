from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass, field
from enum import Enum, auto
from itertools import combinations, product
import re
from typing import ClassVar, Dict, FrozenSet, List, Optional, Protocol, Tuple, Union
import unicodedata


# =============================================================================
# CONSTANTS AND CONFIGURATION
# =============================================================================

class CorrectionType(Enum):
    """Enumeration of correction types for tracking provenance."""
    CHARACTER_SUBSTITUTION = auto()
    PUNCTUATION_RESTORATION = auto()
    GREEK_LETTER_NORMALIZATION = auto()
    UNICODE_NORMALIZATION = auto()
    SPACING_CORRECTION = auto()
    STEREOCHEMISTRY_FIX = auto()
    BRACKET_BALANCING = auto()
    CASE_CORRECTION = auto()
    LOCANT_FORMATTING = auto()
    CUSTOM = auto()


@dataclass(frozen=True)
class CorrectionRule:
    """
    Represents a single correction rule.
    
    Attributes:
        original: The original (erroneous) pattern
        replacement: The corrected replacement
        correction_type: Type of correction being applied
        priority: Higher priority rules are applied first (default: 0)
        context_pattern: Optional regex pattern that must match surrounding context
        description: Human-readable description of the rule
    """
    original: str
    replacement: str
    correction_type: CorrectionType
    priority: int = 0
    context_pattern: Optional[str] = None
    description: str = ""
    
    def __hash__(self) -> int:
        return hash((self.original, self.replacement, self.correction_type))


# Character substitution maps for OCR error correction
# Format: erroneous_char -> list of possible correct chars
OCR_SUBSTITUTIONS: Dict[str, List[str]] = {
    # Digit-letter confusions
    "0": ["O", "o"],
    "1": ["l", "I", "i"],
    "2": ["Z", "z"],
    "5": ["S", "s"],
    "6": ["G", "b"],
    "8": ["B"],
    "9": ["g", "q"],
    
    # Letter-digit confusions (reverse)
    "O": ["0"],
    "o": ["0"],
    "l": ["1", "I"],
    "I": ["1", "l"],
    "i": ["1"],
    "Z": ["2"],
    "z": ["2"],
    "S": ["5"],
    "s": ["5"],
    "G": ["6"],
    "B": ["8"],
    "g": ["9"],
    "q": ["9"],
    
    # Common OCR letter confusions
    "rn": ["m"],
    "m": ["rn"],
    "vv": ["w"],
    "w": ["vv"],
    "cl": ["d"],
    "d": ["cl"],
    "li": ["h"],
    "h": ["li"],
    "nn": ["m"],
    "ii": ["n"],
    "n": ["ii"],
    
    # Case confusions
    "C": ["c"],
    "c": ["C"],
}

# Greek letter mappings (various representations to standard form)
GREEK_LETTER_MAP: Dict[str, str] = {
    # Unicode Greek letters to text
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "ε": "epsilon",
    "ζ": "zeta",
    "η": "eta",
    "θ": "theta",
    "ι": "iota",
    "κ": "kappa",
    "λ": "lambda",
    "μ": "mu",
    "ν": "nu",
    "ξ": "xi",
    "ο": "omicron",
    "π": "pi",
    "ρ": "rho",
    "σ": "sigma",
    "τ": "tau",
    "υ": "upsilon",
    "φ": "phi",
    "χ": "chi",
    "ψ": "psi",
    "ω": "omega",
    
    # Uppercase Greek
    "Α": "Alpha",
    "Β": "Beta",
    "Γ": "Gamma",
    "Δ": "Delta",
    "Ε": "Epsilon",
    "Ζ": "Zeta",
    "Η": "Eta",
    "Θ": "Theta",
    "Ι": "Iota",
    "Κ": "Kappa",
    "Λ": "Lambda",
    "Μ": "Mu",
    "Ν": "Nu",
    "Ξ": "Xi",
    "Ο": "Omicron",
    "Π": "Pi",
    "Ρ": "Rho",
    "Σ": "Sigma",
    "Τ": "Tau",
    "Υ": "Upsilon",
    "Φ": "Phi",
    "Χ": "Chi",
    "Ψ": "Psi",
    "Ω": "Omega",
    
    # Common OCR errors for Greek
    "a": "alpha",  # Only in specific contexts
    "b": "beta",   # Only in specific contexts
    "g": "gamma",  # Only in specific contexts
}

# Unicode normalization map
UNICODE_NORMALIZATION_MAP: Dict[str, str] = {
    # Dashes and hyphens
    "–": "-",  # en dash
    "—": "-",  # em dash
    "−": "-",  # minus sign
    "‐": "-",  # hyphen
    "‑": "-",  # non-breaking hyphen
    "⁃": "-",  # hyphen bullet
    "﹣": "-",  # small hyphen-minus
    "－": "-",  # fullwidth hyphen-minus
    
    # Quotes
    "'": "'",
    "'": "'",
    """: '"',
    """: '"',
    "′": "'",  # prime
    "″": '"',  # double prime
    
    # Spaces
    "\u00a0": " ",  # non-breaking space
    "\u2002": " ",  # en space
    "\u2003": " ",  # em space
    "\u2009": " ",  # thin space
    "\u200a": " ",  # hair space
    "\u200b": "",   # zero-width space
    
    # Brackets
    "［": "[",
    "］": "]",
    "（": "(",
    "）": ")",
    "｛": "{",
    "｝": "}",
    
    # Other
    "×": "x",
    "·": ".",
    "•": ".",
}

# Common chemical prefixes for validation
CHEMICAL_PREFIXES: FrozenSet[str] = frozenset([
    # Multipliers
    "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca",
    "undeca", "dodeca", "mono", "bis", "tris", "tetrakis", "pentakis",
    
    # Alkyl groups
    "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl", "heptyl",
    "octyl", "nonyl", "decyl", "isopropyl", "isobutyl", "tert-butyl",
    "sec-butyl", "neopentyl", "cyclo", "cyclopropyl", "cyclobutyl",
    "cyclopentyl", "cyclohexyl",
    
    # Functional group prefixes
    "hydroxy", "amino", "nitro", "chloro", "bromo", "fluoro", "iodo",
    "cyano", "oxo", "thio", "mercapto", "carboxy", "formyl", "acetyl",
    "benzoyl", "phenyl", "benzyl", "vinyl", "allyl", "oxy", "alkoxy",
    "methoxy", "ethoxy", "propoxy",
    
    # Stereochemistry
    "cis", "trans", "endo", "exo", "syn", "anti", "meso", "ortho",
    "meta", "para", "sec", "tert", "neo", "iso", "n-", "o-", "m-", "p-",
])

# Common chemical suffixes for validation
CHEMICAL_SUFFIXES: FrozenSet[str] = frozenset([
    # Functional groups
    "ol", "al", "one", "oic", "amine", "amide", "nitrile", "ether",
    "ester", "acid", "aldehyde", "ketone", "alcohol",
    
    # Hydrocarbons
    "ane", "ene", "yne", "adiene", "atriene", "ylene",
    
    # Heteroatom
    "oxide", "sulfide", "phosphate", "sulfate", "nitrate", "carbonate",
    
    # Salts and derivatives
    "ate", "ite", "ide", "ium", "sodium", "potassium", "calcium",
    "hydrochloride", "hydrobromide", "sulfonate", "phosphonate",
])

# Stereochemistry patterns
STEREOCHEMISTRY_PATTERNS: List[Tuple[str, str]] = [
    # Patterns to fix (regex, replacement)
    (r'\(R\)', '(R)'),
    (r'\(S\)', '(S)'),
    (r'\(E\)', '(E)'),
    (r'\(Z\)', '(Z)'),
    (r'\bR-', '(R)-'),
    (r'\bS-', '(S)-'),
    (r'\bE-', '(E)-'),
    (r'\bZ-', '(Z)-'),
    (r'\bD-', 'D-'),
    (r'\bL-', 'L-'),
    (r'\bdl-', 'DL-'),
    (r'\brac-', 'rac-'),
    (r'\brel-', 'rel-'),
]


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class Correction:
    """
    Represents a single correction applied to a chemical name.
    
    Attributes:
        position: Character position where correction was applied
        original: Original text that was corrected
        replacement: Text that replaced the original
        correction_type: Type of correction applied
        description: Human-readable description of the correction
    """
    position: int
    original: str
    replacement: str
    correction_type: CorrectionType
    description: str = ""
    
    def __repr__(self) -> str:
        return (f"Correction(pos={self.position}, "
                f"'{self.original}'→'{self.replacement}', "
                f"type={self.correction_type.name})")


@dataclass
class CorrectionCandidate:
    """
    Represents a corrected chemical name candidate with scoring information.
    
    Attributes:
        name: The corrected chemical name
        original_name: The original (uncorrected) input name
        corrections: List of corrections that were applied
        score: Composite score indicating likelihood of correctness (0-1)
        score_components: Individual score components for debugging
        validated: Whether this candidate was validated by external tool
        validation_result: Result from external validation (e.g., SMILES)
    """
    name: str
    original_name: str
    corrections: List[Correction] = field(default_factory=list)
    score: float = 0.0
    score_components: Dict[str, float] = field(default_factory=dict)
    validated: bool = False
    validation_result: Optional[str] = None
    
    def __lt__(self, other: CorrectionCandidate) -> bool:
        """Enable sorting by score (higher is better)."""
        return self.score < other.score
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CorrectionCandidate):
            return NotImplemented
        return self.name == other.name
    
    def __hash__(self) -> int:
        return hash(self.name)
    
    @property
    def num_corrections(self) -> int:
        """Return the number of corrections applied."""
        return len(self.corrections)
    
    @property
    def correction_summary(self) -> str:
        """Return a human-readable summary of corrections."""
        if not self.corrections:
            return "No corrections applied"
        
        summaries = []
        for c in self.corrections:
            summaries.append(f"'{c.original}'→'{c.replacement}'")
        return "; ".join(summaries)
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return {
            "name": self.name,
            "original_name": self.original_name,
            "score": self.score,
            "score_components": self.score_components,
            "num_corrections": self.num_corrections,
            "corrections": [
                {
                    "position": c.position,
                    "original": c.original,
                    "replacement": c.replacement,
                    "type": c.correction_type.name,
                    "description": c.description,
                }
                for c in self.corrections
            ],
            "validated": self.validated,
            "validation_result": self.validation_result,
        }


@dataclass
class CorrectorConfig:
    """
    Configuration for the ChemNameCorrector.
    
    Attributes:
        max_candidates: Maximum number of candidates to generate
        max_corrections_per_candidate: Maximum corrections per candidate
        min_score_threshold: Minimum score to include candidate in results
        enable_character_substitution: Enable OCR character correction
        enable_punctuation_restoration: Enable missing punctuation detection
        enable_greek_normalization: Enable Greek letter conversion
        enable_unicode_normalization: Enable Unicode character normalization
        enable_spacing_correction: Enable space insertion/removal
        enable_stereochemistry_fixes: Enable stereochemistry prefix correction
        enable_bracket_balancing: Enable bracket matching correction
        enable_case_correction: Enable case sensitivity fixes
        custom_substitutions: Additional user-defined substitution rules
        custom_rules: Additional user-defined correction rules
    """
    max_candidates: int = 100
    max_corrections_per_candidate: int = 5
    min_score_threshold: float = 0.1
    enable_character_substitution: bool = True
    enable_punctuation_restoration: bool = True
    enable_greek_normalization: bool = True
    enable_unicode_normalization: bool = True
    enable_spacing_correction: bool = True
    enable_stereochemistry_fixes: bool = True
    enable_bracket_balancing: bool = True
    enable_case_correction: bool = True
    custom_substitutions: Dict[str, List[str]] = field(default_factory=dict)
    custom_rules: List[CorrectionRule] = field(default_factory=list)


# =============================================================================
# PROTOCOLS AND ABSTRACT CLASSES
# =============================================================================

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


class CorrectionStrategy(ABC):
    """
    Abstract base class for correction strategies.
    
    Subclass this to implement custom correction strategies.
    """
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of this strategy."""
        pass
    
    @property
    @abstractmethod
    def correction_type(self) -> CorrectionType:
        """Return the type of corrections this strategy produces."""
        pass
    
    @abstractmethod
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """
        Generate correction candidates from input text.
        
        Args:
            text: Input text to correct
            config: Corrector configuration
            
        Yields:
            Tuples of (corrected_text, list_of_corrections)
        """
        pass


# =============================================================================
# SCORING SYSTEM
# =============================================================================

class ChemicalNameScorer:
    """
    Scores chemical name candidates based on various heuristics.
    
    The scoring system uses multiple factors:
    - Bracket balance and validity
    - Presence of known chemical prefixes/suffixes
    - Valid locant patterns
    - Stereochemistry notation validity
    - Number of corrections (fewer is better)
    - Character distribution plausibility
    """
    
    # Weights for different scoring components
    WEIGHTS: ClassVar[Dict[str, float]] = {
        "bracket_balance": 0.20,
        "prefix_suffix_match": 0.15,
        "locant_validity": 0.15,
        "stereochemistry_validity": 0.10,
        "correction_penalty": 0.20,
        "character_plausibility": 0.10,
        "length_similarity": 0.10,
    }
    
    # Precompiled regex patterns
    LOCANT_PATTERN: ClassVar[re.Pattern] = re.compile(
        r'^\d+[,\d]*-|,\d+-|-\d+,'
    )
    
    STEREOCHEM_PATTERN: ClassVar[re.Pattern] = re.compile(
        r'^\(R\)-|^\(S\)-|^\(E\)-|^\(Z\)-|^D-|^L-|^DL-|^rac-|^rel-|'
        r'\(R\)|\(S\)|\(E\)|\(Z\)|'
        r'-\(R\)-|-\(S\)-|-\(E\)-|-\(Z\)-'
    )
    
    VALID_CHARS_PATTERN: ClassVar[re.Pattern] = re.compile(
        r'^[a-zA-Z0-9\-\(\)\[\]\{\},\.\'\"\s]+$'
    )
    
    def __init__(self, config: Optional[CorrectorConfig] = None):
        """
        Initialize the scorer.
        
        Args:
            config: Optional corrector configuration
        """
        self.config = config or CorrectorConfig()
    
    def score(self, candidate: CorrectionCandidate) -> CorrectionCandidate:
        """
        Calculate the composite score for a candidate.
        
        Args:
            candidate: The candidate to score
            
        Returns:
            The same candidate with updated score and score_components
        """
        components: Dict[str, float] = {}
        
        # Calculate individual components
        components["bracket_balance"] = self._score_bracket_balance(
            candidate.name
        )
        components["prefix_suffix_match"] = self._score_prefix_suffix(
            candidate.name
        )
        components["locant_validity"] = self._score_locants(candidate.name)
        components["stereochemistry_validity"] = self._score_stereochemistry(
            candidate.name
        )
        components["correction_penalty"] = self._score_correction_count(
            candidate.num_corrections
        )
        components["character_plausibility"] = self._score_characters(
            candidate.name
        )
        components["length_similarity"] = self._score_length_similarity(
            candidate.name, candidate.original_name
        )
        
        # Calculate weighted sum
        total_score = sum(
            components[key] * self.WEIGHTS[key]
            for key in self.WEIGHTS
        )
        
        candidate.score = total_score
        candidate.score_components = components
        
        return candidate
    
    def _score_bracket_balance(self, name: str) -> float:
        """
        Score based on bracket balance.
        
        Returns 1.0 if all brackets are balanced, 0.0 if severely unbalanced.
        """
        bracket_pairs = [("(", ")"), ("[", "]"), ("{", "}")]
        total_imbalance = 0
        
        for open_b, close_b in bracket_pairs:
            open_count = name.count(open_b)
            close_count = name.count(close_b)
            total_imbalance += abs(open_count - close_count)
        
        # Also check for proper nesting
        if not self._check_bracket_nesting(name):
            total_imbalance += 2
        
        # Convert imbalance to score (0-1)
        return max(0.0, 1.0 - (total_imbalance * 0.25))
    
    def _check_bracket_nesting(self, name: str) -> bool:
        """Check if brackets are properly nested."""
        stack = []
        bracket_map = {")": "(", "]": "[", "}": "{"}
        
        for char in name:
            if char in "([{":
                stack.append(char)
            elif char in ")]}":
                if not stack or stack[-1] != bracket_map[char]:
                    return False
                stack.pop()
        
        return len(stack) == 0
    
    def _score_prefix_suffix(self, name: str) -> float:
        """Score based on presence of known chemical prefixes/suffixes."""
        name_lower = name.lower()
        score = 0.5  # Base score
        
        # Check for prefixes
        for prefix in CHEMICAL_PREFIXES:
            if prefix in name_lower:
                score += 0.1
                break
        
        # Check for suffixes
        for suffix in CHEMICAL_SUFFIXES:
            if name_lower.endswith(suffix):
                score += 0.2
                break
        
        # Check if it looks like a chemical name pattern
        if re.search(r'\d+-', name):  # Has locant pattern
            score += 0.1
        
        if re.search(r'-\d+', name):  # Ends with locant
            score += 0.05
        
        return min(1.0, score)
    
    def _score_locants(self, name: str) -> float:
        """Score based on valid locant patterns."""
        # Check for proper locant formatting (e.g., "2,3-" or "1-")
        locant_matches = re.findall(r'\d+[,\d]*-', name)
        
        if not locant_matches:
            return 0.5  # Neutral score if no locants
        
        score = 0.7
        
        for match in locant_matches:
            # Check if locant numbers are reasonable (1-99 typically)
            numbers = re.findall(r'\d+', match)
            for num_str in numbers:
                num = int(num_str)
                if 1 <= num <= 50:
                    score += 0.05
                elif num > 100:
                    score -= 0.1
        
        return max(0.0, min(1.0, score))
    
    def _score_stereochemistry(self, name: str) -> float:
        """Score based on stereochemistry notation validity."""
        # Look for stereochemistry patterns
        stereo_matches = self.STEREOCHEM_PATTERN.findall(name)
        
        if not stereo_matches:
            return 0.5  # Neutral if no stereochemistry
        
        # Valid stereochemistry patterns boost score
        return min(1.0, 0.7 + len(stereo_matches) * 0.1)
    
    def _score_correction_count(self, num_corrections: int) -> float:
        """
        Score based on number of corrections (fewer is better).
        
        Follows the principle of minimal correction.
        """
        if num_corrections == 0:
            return 1.0
        elif num_corrections == 1:
            return 0.9
        elif num_corrections == 2:
            return 0.7
        elif num_corrections == 3:
            return 0.5
        else:
            return max(0.1, 1.0 - num_corrections * 0.15)
    
    def _score_characters(self, name: str) -> float:
        """Score based on character plausibility for chemical names."""
        if not name:
            return 0.0
        
        # Check for valid character set
        if not self.VALID_CHARS_PATTERN.match(name):
            return 0.3
        
        # Calculate letter/digit ratio (chemical names are mostly letters)
        letters = sum(1 for c in name if c.isalpha())
        digits = sum(1 for c in name if c.isdigit())
        total = len(name)
        
        if total == 0:
            return 0.0
        
        letter_ratio = letters / total
        
        # Ideal ratio is mostly letters with some digits
        if 0.6 <= letter_ratio <= 0.95:
            return 0.9
        elif 0.4 <= letter_ratio < 0.6:
            return 0.6
        elif letter_ratio > 0.95:
            return 0.8
        else:
            return 0.4
    
    def _score_length_similarity(
        self,
        corrected: str,
        original: str
    ) -> float:
        """Score based on length similarity to original."""
        if not original:
            return 0.5
        
        len_diff = abs(len(corrected) - len(original))
        max_len = max(len(corrected), len(original))
        
        if max_len == 0:
            return 1.0
        
        similarity = 1.0 - (len_diff / max_len)
        return max(0.0, similarity)


# =============================================================================
# CORRECTION STRATEGIES
# =============================================================================

class CharacterSubstitutionStrategy(CorrectionStrategy):
    """
    Strategy for correcting OCR character substitution errors.
    
    Handles common confusions like:
    - 0 ↔ O
    - 1 ↔ l ↔ I
    - rn ↔ m
    - etc.
    """
    
    @property
    def name(self) -> str:
        return "Character Substitution"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_SUBSTITUTION
    
    def __init__(
        self,
        substitutions: Optional[Dict[str, List[str]]] = None
    ):
        """
        Initialize with substitution map.
        
        Args:
            substitutions: Custom substitution map (default: OCR_SUBSTITUTIONS)
        """
        self.substitutions = substitutions or OCR_SUBSTITUTIONS.copy()
    
    def add_substitution(self, original: str, replacements: List[str]) -> None:
        """Add a custom substitution rule."""
        if original in self.substitutions:
            self.substitutions[original].extend(replacements)
        else:
            self.substitutions[original] = replacements
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by applying character substitutions."""
        # Find all positions where substitutions could apply
        substitution_points: List[Tuple[int, str, List[str]]] = []
        
        # Check single characters
        for i, char in enumerate(text):
            if char in self.substitutions:
                substitution_points.append(
                    (i, char, self.substitutions[char])
                )
        
        # Check multi-character patterns
        for pattern, replacements in self.substitutions.items():
            if len(pattern) > 1:
                start = 0
                while True:
                    pos = text.find(pattern, start)
                    if pos == -1:
                        break
                    substitution_points.append((pos, pattern, replacements))
                    start = pos + 1
        
        # Sort by position
        substitution_points.sort(key=lambda x: x[0])
        
        # Generate combinations of substitutions (limit to avoid explosion)
        max_subs = min(
            config.max_corrections_per_candidate,
            len(substitution_points)
        )
        
        # Yield original (no substitutions)
        yield text, []
        
        # Generate single substitutions
        for pos, original, replacements in substitution_points:
            for replacement in replacements:
                new_text = text[:pos] + replacement + text[pos + len(original):]
                correction = Correction(
                    position=pos,
                    original=original,
                    replacement=replacement,
                    correction_type=self.correction_type,
                    description=f"OCR correction: '{original}' → '{replacement}'"
                )
                yield new_text, [correction]
        
        # Generate combinations of 2 substitutions (if allowed)
        if max_subs >= 2 and len(substitution_points) >= 2:
            for combo in combinations(range(len(substitution_points)), 2):
                # Check for overlapping positions
                points = [substitution_points[i] for i in combo]
                if self._positions_overlap(points):
                    continue
                
                # Generate all combinations of replacements
                replacement_options = [
                    [(p[0], p[1], r) for r in p[2]] for p in points
                ]
                
                for replacement_combo in product(*replacement_options):
                    corrections = []
                    new_text = text
                    
                    # Apply in reverse order to maintain positions
                    for pos, original, replacement in sorted(
                        replacement_combo,
                        key=lambda x: x[0],
                        reverse=True
                    ):
                        new_text = (
                            new_text[:pos] +
                            replacement +
                            new_text[pos + len(original):]
                        )
                        corrections.append(Correction(
                            position=pos,
                            original=original,
                            replacement=replacement,
                            correction_type=self.correction_type,
                            description=f"OCR: '{original}' → '{replacement}'"
                        ))
                    
                    yield new_text, corrections
    
    @staticmethod
    def _positions_overlap(
        points: List[Tuple[int, str, List[str]]]
    ) -> bool:
        """Check if any substitution positions overlap."""
        ranges = [(p[0], p[0] + len(p[1])) for p in points]
        ranges.sort()
        
        for i in range(len(ranges) - 1):
            if ranges[i][1] > ranges[i + 1][0]:
                return True
        return False


class PunctuationRestorationStrategy(CorrectionStrategy):
    """
    Strategy for restoring missing punctuation in chemical names.
    
    Handles:
    - Missing hyphens between locants and substituent names
    - Missing commas between locants
    - Missing brackets around stereochemistry
    """
    
    @property
    def name(self) -> str:
        return "Punctuation Restoration"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.PUNCTUATION_RESTORATION
    
    # Patterns for missing punctuation
    MISSING_HYPHEN_PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        # digit followed directly by letter (missing hyphen after locant)
        (re.compile(r'(\d)([a-zA-Z])'), r'\1-\2', "digit-letter"),
        # letter followed by digit (missing hyphen before locant)
        (re.compile(r'([a-zA-Z])(\d)'), r'\1-\2', "letter-digit"),
    ]
    
    MISSING_COMMA_PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        # Adjacent digits that should be comma-separated locants
        (re.compile(r'(\d)\s+(\d)'), r'\1,\2', "spaced-digits"),
    ]
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by restoring missing punctuation."""
        yield text, []  # Original
        
        # Try hyphen restoration
        for pattern, replacement, desc in self.MISSING_HYPHEN_PATTERNS:
            matches = list(pattern.finditer(text))
            
            for match in matches:
                new_text = (
                    text[:match.start()] +
                    pattern.sub(replacement, match.group()) +
                    text[match.end():]
                )
                
                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added hyphen ({desc})"
                    )
                    yield new_text, [correction]
        
        # Try comma restoration
        for pattern, replacement, desc in self.MISSING_COMMA_PATTERNS:
            matches = list(pattern.finditer(text))
            
            for match in matches:
                new_text = (
                    text[:match.start()] +
                    pattern.sub(replacement, match.group()) +
                    text[match.end():]
                )
                
                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added comma ({desc})"
                    )
                    yield new_text, [correction]


class GreekLetterStrategy(CorrectionStrategy):
    """
    Strategy for normalizing Greek letter representations.
    
    Converts various representations to consistent format:
    - α → alpha
    - β → beta
    - etc.
    """
    
    @property
    def name(self) -> str:
        return "Greek Letter Normalization"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.GREEK_LETTER_NORMALIZATION
    
    def __init__(
        self,
        greek_map: Optional[Dict[str, str]] = None,
        to_unicode: bool = False
    ):
        """
        Initialize with Greek letter mapping.
        
        Args:
            greek_map: Custom Greek letter map (default: GREEK_LETTER_MAP)
            to_unicode: If True, convert text to Unicode Greek; if False,
                       convert Unicode to text
        """
        if to_unicode:
            # Reverse the map
            self.greek_map = {v: k for k, v in GREEK_LETTER_MAP.items()}
        else:
            self.greek_map = greek_map or GREEK_LETTER_MAP.copy()
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by normalizing Greek letters."""
        yield text, []  # Original
        
        corrections = []
        new_text = text
        
        # Apply all Greek letter substitutions
        for greek, replacement in self.greek_map.items():
            if greek in new_text:
                pos = new_text.find(greek)
                while pos != -1:
                    corrections.append(Correction(
                        position=pos,
                        original=greek,
                        replacement=replacement,
                        correction_type=self.correction_type,
                        description=f"Greek: '{greek}' → '{replacement}'"
                    ))
                    new_text = new_text[:pos] + replacement + new_text[pos + len(greek):]
                    pos = new_text.find(greek, pos + len(replacement))
        
        if corrections:
            yield new_text, corrections


class UnicodeNormalizationStrategy(CorrectionStrategy):
    """
    Strategy for normalizing Unicode characters to ASCII equivalents.
    
    Handles:
    - Various dash/hyphen characters → standard hyphen
    - Fancy quotes → standard quotes
    - Special spaces → standard space
    - Full-width brackets → standard brackets
    """
    
    @property
    def name(self) -> str:
        return "Unicode Normalization"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.UNICODE_NORMALIZATION
    
    def __init__(
        self,
        normalization_map: Optional[Dict[str, str]] = None
    ):
        """
        Initialize with normalization map.
        
        Args:
            normalization_map: Custom map (default: UNICODE_NORMALIZATION_MAP)
        """
        self.normalization_map = (
            normalization_map or UNICODE_NORMALIZATION_MAP.copy()
        )
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by normalizing Unicode characters."""
        corrections = []
        new_text = text
        
        # First, apply NFKC normalization
        normalized = unicodedata.normalize("NFKC", text)
        if normalized != text:
            corrections.append(Correction(
                position=0,
                original=text,
                replacement=normalized,
                correction_type=self.correction_type,
                description="NFKC normalization"
            ))
            new_text = normalized
        
        # Then apply custom normalizations
        for original, replacement in self.normalization_map.items():
            while original in new_text:
                pos = new_text.find(original)
                corrections.append(Correction(
                    position=pos,
                    original=original,
                    replacement=replacement,
                    correction_type=self.correction_type,
                    description=f"Unicode: '{repr(original)}' → '{replacement}'"
                ))
                new_text = new_text.replace(original, replacement, 1)
        
        if corrections:
            yield new_text, corrections
        
        yield text, []  # Also yield original


class SpacingCorrectionStrategy(CorrectionStrategy):
    """
    Strategy for correcting spacing issues in chemical names.
    
    Handles:
    - Missing spaces between words
    - Extra spaces within words
    - Spaces around hyphens
    """
    
    @property
    def name(self) -> str:
        return "Spacing Correction"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.SPACING_CORRECTION
    
    # Common word boundaries in chemical names
    WORD_BOUNDARIES: ClassVar[List[str]] = [
        "acid", "alcohol", "aldehyde", "amine", "amide", "ester",
        "ether", "ketone", "oxide", "sulfide", "chloride", "bromide",
        "iodide", "fluoride", "hydroxide", "carbonate", "sulfate",
        "nitrate", "phosphate",
    ]
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by correcting spacing."""
        yield text, []  # Original
        
        # Remove spaces around hyphens
        if " - " in text or "- " in text or " -" in text:
            new_text = re.sub(r'\s*-\s*', '-', text)
            if new_text != text:
                correction = Correction(
                    position=0,
                    original=text,
                    replacement=new_text,
                    correction_type=self.correction_type,
                    description="Removed spaces around hyphens"
                )
                yield new_text, [correction]
        
        # Check for missing spaces before common word endings
        for word in self.WORD_BOUNDARIES:
            # Pattern: letter + word without space
            pattern = re.compile(rf'([a-z])({word})\b', re.IGNORECASE)
            match = pattern.search(text)
            if match:
                new_text = pattern.sub(r'\1 \2', text)
                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=match.group(1) + " " + match.group(2),
                        correction_type=self.correction_type,
                        description=f"Added space before '{word}'"
                    )
                    yield new_text, [correction]
        
        # Collapse multiple spaces
        if "  " in text:
            new_text = re.sub(r'\s+', ' ', text)
            if new_text != text:
                correction = Correction(
                    position=0,
                    original=text,
                    replacement=new_text,
                    correction_type=self.correction_type,
                    description="Collapsed multiple spaces"
                )
                yield new_text, [correction]


class StereochemistryFixStrategy(CorrectionStrategy):
    """
    Strategy for fixing stereochemistry notation.
    
    Handles:
    - Missing parentheses around R/S/E/Z
    - Incorrect case
    - Missing hyphens
    """
    
    @property
    def name(self) -> str:
        return "Stereochemistry Fix"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.STEREOCHEMISTRY_FIX
    
    PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        (re.compile(r'\bR-(?!\))'), '(R)-', "Add parens to R"),
        (re.compile(r'\bS-(?!\))'), '(S)-', "Add parens to S"),
        (re.compile(r'\bE-(?!\))'), '(E)-', "Add parens to E"),
        (re.compile(r'\bZ-(?!\))'), '(Z)-', "Add parens to Z"),
        (re.compile(r'\(r\)'), '(R)', "Uppercase R"),
        (re.compile(r'\(s\)'), '(S)', "Uppercase S"),
        (re.compile(r'\(e\)'), '(E)', "Uppercase E"),
        (re.compile(r'\(z\)'), '(Z)', "Uppercase Z"),
        (re.compile(r'\bd-(?=[A-Za-z])'), 'D-', "Uppercase D"),
        (re.compile(r'\bl-(?=[A-Za-z])'), 'L-', "Uppercase L"),
    ]
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by fixing stereochemistry notation."""
        yield text, []  # Original
        
        for pattern, replacement, desc in self.PATTERNS:
            if pattern.search(text):
                new_text = pattern.sub(replacement, text)
                if new_text != text:
                    correction = Correction(
                        position=0,
                        original=text,
                        replacement=new_text,
                        correction_type=self.correction_type,
                        description=desc
                    )
                    yield new_text, [correction]


class BracketBalancingStrategy(CorrectionStrategy):
    """
    Strategy for balancing brackets in chemical names.
    
    Attempts to fix:
    - Missing closing brackets
    - Missing opening brackets
    - Mismatched bracket types
    """
    
    @property
    def name(self) -> str:
        return "Bracket Balancing"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.BRACKET_BALANCING
    
    BRACKET_PAIRS: ClassVar[List[Tuple[str, str]]] = [
        ("(", ")"),
        ("[", "]"),
        ("{", "}"),
    ]
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by balancing brackets."""
        yield text, []  # Original
        
        for open_b, close_b in self.BRACKET_PAIRS:
            open_count = text.count(open_b)
            close_count = text.count(close_b)
            
            if open_count > close_count:
                # Missing closing brackets - add at end
                diff = open_count - close_count
                new_text = text + close_b * diff
                correction = Correction(
                    position=len(text),
                    original="",
                    replacement=close_b * diff,
                    correction_type=self.correction_type,
                    description=f"Added {diff} closing '{close_b}'"
                )
                yield new_text, [correction]
            
            elif close_count > open_count:
                # Missing opening brackets - add at start
                diff = close_count - open_count
                new_text = open_b * diff + text
                correction = Correction(
                    position=0,
                    original="",
                    replacement=open_b * diff,
                    correction_type=self.correction_type,
                    description=f"Added {diff} opening '{open_b}'"
                )
                yield new_text, [correction]


class CaseNormalizationStrategy(CorrectionStrategy):
    """
    Strategy for normalizing case in chemical names.
    
    Chemical names typically:
    - Start with lowercase (unless it's a proper name)
    - Have specific case for stereochemistry (D-, L-, R, S, etc.)
    """
    
    @property
    def name(self) -> str:
        return "Case Normalization"
    
    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CASE_CORRECTION
    
    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates with case normalization."""
        yield text, []  # Original
        
        # Try lowercase version (common for chemical names)
        if text != text.lower():
            # Preserve stereochemistry markers
            new_text = self._smart_lowercase(text)
            if new_text != text:
                correction = Correction(
                    position=0,
                    original=text,
                    replacement=new_text,
                    correction_type=self.correction_type,
                    description="Smart lowercase conversion"
                )
                yield new_text, [correction]
    
    def _smart_lowercase(self, text: str) -> str:
        """
        Convert to lowercase while preserving stereochemistry markers.
        """
        # Patterns to preserve
        preserve_patterns = [
            r'\(R\)', r'\(S\)', r'\(E\)', r'\(Z\)',
            r'\bD-', r'\bL-', r'\bDL-',
            r'\bN-', r'\bO-', r'\bS-', r'\bC-',
        ]
        
        # Find all positions to preserve
        preserved: Dict[int, str] = {}
        for pattern in preserve_patterns:
            for match in re.finditer(pattern, text):
                for i in range(match.start(), match.end()):
                    preserved[i] = text[i]
        
        # Build new text
        result = []
        for i, char in enumerate(text):
            if i in preserved:
                result.append(preserved[i])
            else:
                result.append(char.lower())
        
        return "".join(result)


# =============================================================================
# VALIDATORS
# =============================================================================

class DummyValidator:
    """
    A dummy validator that always returns True.
    
    Use this for testing or when no external validation is needed.
    """
    
    def validate(self, name: str) -> Tuple[bool, Optional[str]]:
        """Always returns True with no result."""
        return True, None


class OPSINValidator:
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
        import subprocess
        
        try:
            result = subprocess.run(
                [self.opsin_path, "-osmi"],
                input=name,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            smiles = result.stdout.strip()
            if smiles and not smiles.startswith("Error"):
                return True, smiles
            return False, None
            
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            return False, None


class PubChemValidator:
    """
    Validator that uses PubChem's PUG REST API to validate names.
    
    Requires internet connection.
    """
    
    BASE_URL: ClassVar[str] = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
    )
    
    def __init__(self, timeout: int = 10):
        """
        Initialize PubChem validator.
        
        Args:
            timeout: Request timeout in seconds
        """
        self.timeout = timeout
    
    def validate(self, name: str) -> Tuple[bool, Optional[str]]:
        """
        Validate a chemical name using PubChem API.
        
        Args:
            name: Chemical name to validate
            
        Returns:
            Tuple of (is_valid, smiles_or_none)
        """
        import urllib.request
        import urllib.parse
        import json
        
        encoded_name = urllib.parse.quote(name)
        url = f"{self.BASE_URL}/{encoded_name}/property/CanonicalSMILES/JSON"
        
        try:
            request = urllib.request.Request(url)
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                data = json.loads(response.read().decode())
                
            smiles = (
                data.get("PropertyTable", {})
                .get("Properties", [{}])[0]
                .get("CanonicalSMILES")
            )
            
            if smiles:
                return True, smiles
            return False, None
            
        except Exception:
            return False, None


# =============================================================================
# MAIN CORRECTOR CLASS
# =============================================================================

class ChemNameCorrector:
    """
    Main class for correcting OCR errors in chemical names.
    
    This class orchestrates the correction process by:
    1. Applying configured correction strategies
    2. Generating candidate corrections
    3. Scoring candidates
    4. Optionally validating with external tools
    5. Returning ranked results
    
    Example:
        >>> corrector = ChemNameCorrector()
        >>> results = corrector.correct("2-ch1oropropanoic acid")
        >>> print(results[0].name)
        2-chloropropanoic acid
        
        >>> # With custom configuration
        >>> config = CorrectorConfig(max_candidates=50)
        >>> corrector = ChemNameCorrector(config)
        
        >>> # With external validation
        >>> validator = PubChemValidator()
        >>> results = corrector.correct("asprin", validator=validator)
    
    Attributes:
        config: Configuration for the corrector
        strategies: List of active correction strategies
        scorer: Scoring instance for ranking candidates
    """
    
    def __init__(
        self,
        config: Optional[CorrectorConfig] = None,
        strategies: Optional[List[CorrectionStrategy]] = None
    ):
        """
        Initialize the chemical name corrector.
        
        Args:
            config: Configuration object (uses defaults if None)
            strategies: Custom list of strategies (uses defaults if None)
        """
        self.config = config or CorrectorConfig()
        self.scorer = ChemicalNameScorer(self.config)
        
        if strategies is not None:
            self.strategies = strategies
        else:
            self.strategies = self._create_default_strategies()
    
    def _create_default_strategies(self) -> List[CorrectionStrategy]:
        """Create the default set of correction strategies based on config."""
        strategies: List[CorrectionStrategy] = []
        
        if self.config.enable_unicode_normalization:
            strategies.append(UnicodeNormalizationStrategy())
        
        if self.config.enable_greek_normalization:
            strategies.append(GreekLetterStrategy())
        
        if self.config.enable_character_substitution:
            char_strategy = CharacterSubstitutionStrategy()
            # Add any custom substitutions
            for orig, replacements in self.config.custom_substitutions.items():
                char_strategy.add_substitution(orig, replacements)
            strategies.append(char_strategy)
        
        if self.config.enable_punctuation_restoration:
            strategies.append(PunctuationRestorationStrategy())
        
        if self.config.enable_spacing_correction:
            strategies.append(SpacingCorrectionStrategy())
        
        if self.config.enable_stereochemistry_fixes:
            strategies.append(StereochemistryFixStrategy())
        
        if self.config.enable_bracket_balancing:
            strategies.append(BracketBalancingStrategy())
        
        if self.config.enable_case_correction:
            strategies.append(CaseNormalizationStrategy())
        
        return strategies
    
    def add_strategy(self, strategy: CorrectionStrategy) -> None:
        """
        Add a custom correction strategy.
        
        Args:
            strategy: The strategy to add
        """
        self.strategies.append(strategy)
    
    def remove_strategy(self, strategy_name: str) -> bool:
        """
        Remove a strategy by name.
        
        Args:
            strategy_name: Name of the strategy to remove
            
        Returns:
            True if strategy was found and removed, False otherwise
        """
        for i, strategy in enumerate(self.strategies):
            if strategy.name == strategy_name:
                self.strategies.pop(i)
                return True
        return False
    
    def correct(
        self,
        name: str,
        validator: Optional[Validator] = None,
        validate_all: bool = False
    ) -> List[CorrectionCandidate]:
        """
        Correct a chemical name and return ranked candidates.
        
        Args:
            name: The chemical name to correct
            validator: Optional external validator for candidates
            validate_all: If True, validate all candidates; if False,
                         stop after first valid candidate
        
        Returns:
            List of CorrectionCandidate objects, sorted by score (descending)
        """
        # Generate all candidates
        candidates = self._generate_all_candidates(name)
        
        # Remove duplicates while preserving best corrections
        unique_candidates = self._deduplicate_candidates(candidates)
        
        # Score all candidates
        scored_candidates = [
            self.scorer.score(candidate)
            for candidate in unique_candidates
        ]
        
        # Filter by minimum score threshold
        filtered_candidates = [
            c for c in scored_candidates
            if c.score >= self.config.min_score_threshold
        ]
        
        # Sort by score (descending)
        sorted_candidates = sorted(
            filtered_candidates,
            key=lambda c: c.score,
            reverse=True
        )
        
        # Limit to max candidates
        limited_candidates = sorted_candidates[:self.config.max_candidates]
        
        # Optionally validate with external tool
        if validator is not None:
            limited_candidates = self._validate_candidates(
                limited_candidates,
                validator,
                validate_all
            )
        
        return limited_candidates
    
    def correct_batch(
        self,
        names: List[str],
        validator: Optional[Validator] = None
    ) -> Dict[str, List[CorrectionCandidate]]:
        """
        Correct multiple chemical names.
        
        Args:
            names: List of chemical names to correct
            validator: Optional external validator
            
        Returns:
            Dictionary mapping original names to their candidates
        """
        results = {}
        for name in names:
            results[name] = self.correct(name, validator)
        return results
    
    def _generate_all_candidates(
        self,
        name: str
    ) -> List[CorrectionCandidate]:
        """Generate candidates from all strategies."""
        candidates: List[CorrectionCandidate] = []
        
        # Add original as a candidate
        candidates.append(CorrectionCandidate(
            name=name,
            original_name=name,
            corrections=[]
        ))
        
        # Apply each strategy
        current_texts = [(name, [])]  # (text, corrections)
        
        for strategy in self.strategies:
            new_texts = []
            
            for text, existing_corrections in current_texts:
                for new_text, new_corrections in strategy.generate_candidates(
                    text, self.config
                ):
                    combined_corrections = existing_corrections + new_corrections
                    
                    # Check max corrections limit
                    if len(combined_corrections) <= self.config.max_corrections_per_candidate:
                        new_texts.append((new_text, combined_corrections))
                        
                        candidates.append(CorrectionCandidate(
                            name=new_text,
                            original_name=name,
                            corrections=combined_corrections
                        ))
            
            # Update current texts for next strategy
            # Limit to prevent explosion
            current_texts = new_texts[:self.config.max_candidates * 2]
        
        return candidates
    
    def _deduplicate_candidates(
        self,
        candidates: List[CorrectionCandidate]
    ) -> List[CorrectionCandidate]:
        """Remove duplicate candidates, keeping the one with fewer corrections."""
        seen: Dict[str, CorrectionCandidate] = {}
        
        for candidate in candidates:
            if candidate.name not in seen:
                seen[candidate.name] = candidate
            else:
                # Keep the one with fewer corrections
                if candidate.num_corrections < seen[candidate.name].num_corrections:
                    seen[candidate.name] = candidate
        
        return list(seen.values())
    
    def _validate_candidates(
        self,
        candidates: List[CorrectionCandidate],
        validator: Validator,
        validate_all: bool
    ) -> List[CorrectionCandidate]:
        """Validate candidates using external validator."""
        validated_candidates = []
        
        for candidate in candidates:
            is_valid, result = validator.validate(candidate.name)
            candidate.validated = True
            
            if is_valid:
                candidate.validation_result = result
                # Boost score for valid candidates
                candidate.score = min(1.0, candidate.score + 0.3)
                validated_candidates.append(candidate)
                
                if not validate_all:
                    # Add remaining unvalidated candidates
                    remaining_idx = candidates.index(candidate) + 1
                    validated_candidates.extend(candidates[remaining_idx:])
                    break
            else:
                # Lower score for invalid candidates
                candidate.score = max(0.0, candidate.score - 0.2)
                validated_candidates.append(candidate)
        
        return sorted(validated_candidates, key=lambda c: c.score, reverse=True)
    
    def get_best_candidate(
        self,
        name: str,
        validator: Optional[Validator] = None
    ) -> Optional[CorrectionCandidate]:
        """
        Get the single best correction candidate.
        
        Args:
            name: Chemical name to correct
            validator: Optional external validator
            
        Returns:
            Best candidate, or None if no candidates found
        """
        candidates = self.correct(name, validator)
        return candidates[0] if candidates else None
    
    def explain_corrections(
        self,
        candidate: CorrectionCandidate
    ) -> str:
        """
        Generate a human-readable explanation of corrections.
        
        Args:
            candidate: The candidate to explain
            
        Returns:
            Multi-line string explaining all corrections
        """
        lines = [
            f"Original: {candidate.original_name}",
            f"Corrected: {candidate.name}",
            f"Score: {candidate.score:.3f}",
            f"Number of corrections: {candidate.num_corrections}",
            "",
            "Score components:",
        ]
        
        for component, value in candidate.score_components.items():
            lines.append(f"  - {component}: {value:.3f}")
        
        if candidate.corrections:
            lines.append("")
            lines.append("Corrections applied:")
            for i, correction in enumerate(candidate.corrections, 1):
                lines.append(
                    f"  {i}. [{correction.correction_type.name}] "
                    f"'{correction.original}' → '{correction.replacement}'"
                )
                if correction.description:
                    lines.append(f"     {correction.description}")
        
        if candidate.validated:
            lines.append("")
            lines.append(f"Validated: {candidate.validation_result or 'No result'}")
        
        return "\n".join(lines)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def create_custom_strategy(
    name: str,
    correction_type: CorrectionType,
    rules: List[Tuple[str, str, str]]
) -> CorrectionStrategy:
    """
    Factory function to create a custom correction strategy.
    
    Args:
        name: Name for the strategy
        correction_type: Type of corrections
        rules: List of (pattern, replacement, description) tuples
        
    Returns:
        A new CorrectionStrategy instance
        
    Example:
        >>> rules = [
        ...     ("phospherous", "phosphorus", "Common misspelling"),
        ...     ("flourine", "fluorine", "Common misspelling"),
        ... ]
        >>> strategy = create_custom_strategy(
        ...     "Spelling Corrections",
        ...     CorrectionType.CUSTOM,
        ...     rules
        ... )
        >>> corrector.add_strategy(strategy)
    """
    class CustomStrategy(CorrectionStrategy):
        def __init__(self):
            self._name = name
            self._correction_type = correction_type
            self._rules = rules
        
        @property
        def name(self) -> str:
            return self._name
        
        @property
        def correction_type(self) -> CorrectionType:
            return self._correction_type
        
        def generate_candidates(
            self,
            text: str,
            config: CorrectorConfig
        ) -> Iterator[Tuple[str, List[Correction]]]:
            yield text, []  # Original
            
            for pattern, replacement, description in self._rules:
                if pattern in text:
                    new_text = text.replace(pattern, replacement)
                    pos = text.find(pattern)
                    correction = Correction(
                        position=pos,
                        original=pattern,
                        replacement=replacement,
                        correction_type=self._correction_type,
                        description=description
                    )
                    yield new_text, [correction]
    
    return CustomStrategy()


def quick_correct(
    name: str,
    return_all: bool = False
) -> Union[str, List[CorrectionCandidate]]:
    """
    Quick correction function for simple use cases.
    
    Args:
        name: Chemical name to correct
        return_all: If True, return all candidates; if False, return best name
        
    Returns:
        Best corrected name (str) or list of all candidates
        
    Example:
        >>> quick_correct("2-ch1oropropanoic acid")
        '2-chloropropanoic acid'
    """
    corrector = ChemNameCorrector()
    candidates = corrector.correct(name)
    
    if return_all:
        return candidates
    
    return candidates[0].name if candidates else name


# =============================================================================
# MAIN / CLI
# =============================================================================

def main():
    """Command-line interface for the chemical name corrector."""
    import argparse
    import json
    import sys
    
    parser = argparse.ArgumentParser(
        description="Correct OCR errors in chemical names",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python chemname_corrector.py "2-ch1oropropanoic acid"
  python chemname_corrector.py -n 5 "2-ch1oropropanoic acid"
  python chemname_corrector.py --json "aspirin"
  echo "2-ch1oropropanoic acid" | python chemname_corrector.py -
        """
    )
    
    parser.add_argument(
        "name",
        nargs="?",
        help="Chemical name to correct (use '-' to read from stdin)"
    )
    
    parser.add_argument(
        "-n", "--num-candidates",
        type=int,
        default=5,
        help="Number of candidates to return (default: 5)"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output in JSON format"
    )
    
    parser.add_argument(
        "--explain",
        action="store_true",
        help="Show detailed explanation for best candidate"
    )
    
    parser.add_argument(
        "--validate-pubchem",
        action="store_true",
        help="Validate candidates against PubChem"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    # Get input
    if args.name == "-" or args.name is None:
        name = sys.stdin.read().strip()
    else:
        name = args.name
    
    if not name:
        parser.error("No chemical name provided")
    
    # Configure corrector
    config = CorrectorConfig(max_candidates=args.num_candidates * 2)
    corrector = ChemNameCorrector(config)
    
    # Set up validator if requested
    validator = None
    if args.validate_pubchem:
        validator = PubChemValidator()
    
    # Run correction
    candidates = corrector.correct(name, validator=validator)[:args.num_candidates]
    
    # Output results
    if args.json:
        output = {
            "original": name,
            "candidates": [c.to_dict() for c in candidates]
        }
        print(json.dumps(output, indent=2))
    
    elif args.explain and candidates:
        print(corrector.explain_corrections(candidates[0]))
    
    else:
        print(f"Original: {name}")
        print(f"Top {len(candidates)} candidates:")
        print("-" * 60)
        
        for i, candidate in enumerate(candidates, 1):
            validated_mark = " ✓" if candidate.validation_result else ""
            print(
                f"{i}. {candidate.name} "
                f"(score: {candidate.score:.3f}, "
                f"corrections: {candidate.num_corrections}){validated_mark}"
            )
            
            if args.verbose and candidate.corrections:
                for correction in candidate.corrections:
                    print(
                        f"   - {correction.correction_type.name}: "
                        f"'{correction.original}' → '{correction.replacement}'"
                    )


if __name__ == "__main__":
    main()