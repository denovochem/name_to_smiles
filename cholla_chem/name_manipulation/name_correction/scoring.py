from __future__ import annotations
from importlib import resources
import json
import re
from typing import ClassVar, Dict, List, Optional

from Levenshtein import ratio as levenshtein_ratio
from cholla_chem.name_manipulation.name_correction.dataclasses import (
    CorrectionCandidate,
    CorrectorConfig,
)
from cholla_chem.name_manipulation.name_correction.regexes import PATTERNS


def load_chemical_name_tokens() -> List[str]:
    """Load chemical name tokens from package data using importlib.resources."""
    # Open the file as text from within the package
    with resources.open_text("cholla_chem.datafiles", "chemical_name_tokens.json") as f:
        return json.load(f)


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
        "correction_penalty": 0.20,
        "levenshtein_similarity": 0.25,
        "number_of_chemical_morphemes": 0.15,
        "number_of_regex_matches": 0.20,
    }

    def __init__(self, config: Optional[CorrectorConfig] = None):
        """
        Initialize the scorer.

        Args:
            config: Optional corrector configuration
        """
        self.config = config or CorrectorConfig()

        self.chemical_morpheme_list = load_chemical_name_tokens()

        self.chemical_morpheme_list = [
            ele for ele in self.chemical_morpheme_list if len(ele) > 1
        ]
        self.original_name_num_morphemes: Optional[int] = None
        self.original_name_num_regex_matches: Optional[int] = None

    def score(self, candidate: CorrectionCandidate) -> CorrectionCandidate:
        """
        Calculate the composite score for a candidate.

        Args:
            candidate: The candidate to score

        Returns:
            The same candidate with updated score and score_components
        """
        components: Dict[str, float] = {}

        if not self.original_name_num_morphemes:
            self.original_name_num_morphemes = self._get_number_of_chemical_morphemes(
                candidate.original_name
            )

        if not self.original_name_num_regex_matches:
            self.original_name_num_regex_matches = self._get_number_of_regex_matches(
                candidate.original_name
            )

        # Calculate individual components
        components["bracket_balance"] = self._score_bracket_balance(candidate.name)
        components["correction_penalty"] = self._score_correction_count(
            candidate.num_corrections
        )
        components["levenshtein_similarity"] = self._score_levenshtein_similarity(
            candidate.name, candidate.original_name
        )
        components["number_of_chemical_morphemes"] = self._score_chemical_morphemes(
            self._get_number_of_chemical_morphemes(candidate.name),
            self.original_name_num_morphemes,
        )
        components["number_of_regex_matches"] = self._score_regex_matches(
            self._get_number_of_regex_matches(candidate.name),
            self.original_name_num_regex_matches,
        )

        # Calculate weighted sum
        total_score = sum(components[key] * self.WEIGHTS[key] for key in self.WEIGHTS)

        candidate.score = total_score
        candidate.score_components = components

        self.original_name_num_regex_matches = None
        self.original_name_num_morphemes = None

        return candidate

    def _get_number_of_chemical_morphemes(self, name: str) -> int:
        """
        Get the number of chemical morphemes in the name.
        """
        num_morphemes = sum(
            1 for morpheme in self.chemical_morpheme_list if morpheme in name
        )

        if not num_morphemes:
            return 0

        return num_morphemes

    def _score_chemical_morphemes(
        self,
        corrected_name_num_morphemes: int,
        original_name_num_morphemes: Optional[int],
    ) -> float:
        """
        Score based on the number of chemical morphemes in the name.
        """
        if not original_name_num_morphemes:
            return 0

        return (
            corrected_name_num_morphemes / original_name_num_morphemes
            if original_name_num_morphemes > 0
            else 0
        )

    def _get_number_of_regex_matches(self, name: str) -> int:
        """
        Get the number of regex matches in the name.
        """
        num_regex_matches = sum(
            1 for pattern, _, _ in PATTERNS if re.search(pattern, name)
        )

        if not num_regex_matches:
            return 0

        return num_regex_matches

    def _score_regex_matches(
        self,
        corrected_name_num_regex_matches: int,
        original_name_num_regex_matches: Optional[int],
    ) -> float:
        """
        Score based on the number of regex matches in the name.
        """
        if not original_name_num_regex_matches:
            return 0

        return 1 - (
            corrected_name_num_regex_matches / original_name_num_regex_matches
            if original_name_num_regex_matches > 0
            else 0
        )

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

    def _score_levenshtein_similarity(self, corrected: str, original: str) -> float:
        """Score based on Levenshtein similarity to original."""
        return levenshtein_ratio(corrected, original)
