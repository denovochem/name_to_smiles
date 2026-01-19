from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional


@dataclass
class CorrectorConfig:
    """
    Configuration for the ChemNameCorrector.

    Attributes:
        max_candidates: Maximum number of candidates to generate
        max_corrections_per_candidate: Maximum corrections per candidate
        min_score_threshold: Minimum score to include candidate in results
        enable_character_substitution: Enable OCR character correction
        max_character_substitution_edits: Max number of substitution edits
        enable_punctuation_restoration: Enable missing punctuation detection
        enable_bracket_balancing: Enable bracket matching correction
        custom_substitutions: Additional user-defined substitution rules
        custom_rules: Additional user-defined correction rules
        enable_external_validation: Enable external validation of candidates
    """

    max_candidates: int = 100
    max_corrections_per_candidate: int = 5
    min_score_threshold: float = 0.1
    enable_locant_correction: bool = True
    enable_character_substitution: bool = True
    max_character_substitution_edits: int = 1
    enable_punctuation_restoration: bool = False
    enable_bracket_balancing: bool = False
    custom_substitutions: Dict[str, List[str]] = field(default_factory=dict)
    custom_rules: List[CorrectionRule] = field(default_factory=list)
    enable_external_validation: bool = True


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
        return (
            f"Correction(pos={self.position}, "
            f"'{self.original}'→'{self.replacement}', "
            f"type={self.correction_type.name})"
        )


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


class CorrectionType(Enum):
    """Enumeration of correction types for tracking provenance."""

    TOKEN_SEGMENTATION = auto()
    LOCANT_FIX = auto()
    CHARACTER_SUBSTITUTION = auto()
    CHARACTER_INSERTION = auto()
    CHARACTER_DELETION = auto()
    PUNCTUATION_RESTORATION = auto()
    STEREOCHEMISTRY_FIX = auto()
    BRACKET_BALANCING = auto()
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
        correction score: Higher correction scores indicate a less likely correction
        context_pattern: Optional regex pattern that must match surrounding context
        description: Human-readable description of the rule
    """

    original: str
    replacement: str
    correction_type: CorrectionType
    priority: int = 0
    correction_score: int = 0
    context_pattern: Optional[str] = None
    description: str = ""

    def __hash__(self) -> int:
        return hash((self.original, self.replacement, self.correction_type))
