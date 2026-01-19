from __future__ import annotations
from abc import ABC, abstractmethod
from collections import deque
from collections.abc import Iterator
import itertools
import json
import os
from pathlib import Path
import re
from typing import ClassVar, Dict, List, Optional, Set, Tuple

from flashtext import KeywordProcessor
from placeholder_name.name_manipulation.name_correction.dataclasses import (
    Correction,
    CorrectionType,
    CorrectorConfig,
)
from placeholder_name.name_manipulation.name_correction.regexes import PATTERNS
from placeholder_name.utils.constants import (
    KEYBOARD_NEIGHBOR_SUBSTITUTIONS,
    OCR_SUBSTITUTIONS,
)
# from placeholder_name.utils.logging_config import logger

# Get the directory of the current file
BASE_DIR = Path(__file__).resolve().parent
CHEMICAL_NAME_TOKENS_PATH = os.path.abspath(
    BASE_DIR.parent.parent / "datafiles" / "chemical_name_tokens.json"
)


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
        self, text: str, config: CorrectorConfig
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


class CharacterSubstitutionStrategy(CorrectionStrategy):
    """
    Strategy for correcting OCR character substitution errors using Aho-Corasick (FlashText).
    """

    @property
    def name(self) -> str:
        return "Character Substitution"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.CHARACTER_SUBSTITUTION

    def __init__(
        self, max_edits: int = 1, substitutions: Optional[Dict[str, List[str]]] = None
    ):
        """
        Initialize with substitution map using FlashText for O(N) performance.
        """
        self.keyword_processor = KeywordProcessor()
        self.keyword_processor.non_word_boundaries = set()

        # We store metadata to handle 'OCR' vs 'Typo' descriptions if needed,
        # though standard FlashText maps keyword -> string.
        # Here we map: error_string -> correct_token

        def process_and_add(
            chemical_name_tokens: List[str],
            source_dict: Dict[str, List[str]],
            keyword_processor,
            max_edits: int,
            source_type_name: str,
        ):
            for chemical_name_token in chemical_name_tokens:
                generated_errors = self.generate_substitution_dict(
                    chemical_name_token, source_dict, max_edits=max_edits
                )

                for error in generated_errors:
                    if len(error) <= 2:
                        continue
                    if error in chemical_name_tokens:
                        continue

                    # OPTIMIZATION:
                    # We add the "squashed" (no space) version of the error to the processor.
                    # We store the correct token as the value.
                    clean_error_key = error.replace(" ", "")
                    keyword_processor.add_keyword(clean_error_key, chemical_name_token)

            return keyword_processor

        if not substitutions:
            with open(CHEMICAL_NAME_TOKENS_PATH, "rb") as f:
                chemical_name_tokens = json.load(f)

            # Process OCR Substitutions
            self.keyword_processor = process_and_add(
                chemical_name_tokens,
                OCR_SUBSTITUTIONS,
                self.keyword_processor,
                max_edits,
                "OCR",
            )

            # Process Typo Substitutions
            self.keyword_processor = process_and_add(
                chemical_name_tokens,
                KEYBOARD_NEIGHBOR_SUBSTITUTIONS,
                self.keyword_processor,
                max_edits,
                "Typo",
            )

        else:
            # Handle manual substitutions injection
            for error, replacement_list in substitutions.items():
                if replacement_list:
                    clean_error_key = error.replace(" ", "")
                    # Taking the first replacement as the primary one
                    self.keyword_processor.add_keyword(
                        clean_error_key, replacement_list[0]
                    )

    def generate_substitution_dict(
        self, word: str, substitutions: Dict[str, List[str]], max_edits: int = 1
    ) -> Set[str]:
        results = set()
        queue = deque([(word, 0)])
        visited = {word}

        while queue:
            current_s, depth = queue.popleft()
            if depth >= max_edits:
                continue

            for i in range(len(current_s)):
                for key, replacements in substitutions.items():
                    key_len = len(key)
                    if i + key_len > len(current_s):
                        continue

                    if current_s[i : i + key_len] == key:
                        for replacement in replacements:
                            prefix = current_s[:i]
                            suffix = current_s[i + key_len :]
                            new_candidate = prefix + replacement + suffix

                            if new_candidate not in visited:
                                visited.add(new_candidate)
                                results.add(new_candidate)
                                queue.append((new_candidate, depth + 1))
        return results

    def generate_candidates(
        self,
        text: str,
        config: CorrectorConfig,  # Assuming this config has a 'max_substitutions' or similar
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """
        Generate candidates by finding all matches and applying combinations of them.
        """

        # 1. SQUASH & MAP
        clean_chars = []
        idx_map = []

        for i, char in enumerate(text):
            if not char.isspace():
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)

        if not clean_text:
            return

        # 2. EXTRACT matches
        # matches = [('1,2-dichlorobenzene', 5, 12), ...]
        raw_matches = self.keyword_processor.extract_keywords(
            clean_text, span_info=True
        )

        if not raw_matches:
            return

        # 3. PREPARE Match Objects with Original Indices
        # We process these into a structured list so we can combine them easily
        mapped_matches = []

        for correct_token, start, end in raw_matches:
            # Map back to original text coordinates
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_match_string = text[orig_start_index:orig_end_index]

            match_data = {
                "replacement": correct_token,
                "start": orig_start_index,
                "end": orig_end_index,
                "original": original_match_string,
            }
            mapped_matches.append(match_data)

        # 4. GENERATE COMBINATIONS
        # By default, FlashText returns non-overlapping matches, so we don't need
        # complex overlap checks for the raw output.

        # Determine how many simultaneous corrections allowed.
        # Default to all found if not specified in config.
        max_simultaneous = 100

        # We iterate from 1 up to max_simultaneous
        for r in range(1, min(len(mapped_matches), max_simultaneous) + 1):
            # itertools.combinations produces tuples of matches: (match1, match2)
            for combo in itertools.combinations(mapped_matches, r):
                # Sort combo by start index DESCENDING.
                # This is crucial: we must replace from end to start so
                # earlier indices don't shift when we modify the string.
                sorted_combo = sorted(combo, key=lambda x: x["start"], reverse=True)

                # Check for overlaps within this specific combination
                # (Paranoia check: FlashText usually handles this, but strictly speaking
                # if you have custom logic, this ensures safety)
                if self._check_overlap(sorted_combo):
                    continue

                # Apply corrections
                new_text_chars = list(text)  # Mutable char list
                corrections_list = []

                for match in sorted_combo:
                    # Replace in the char list
                    # We replace the slice with the new string
                    new_text_chars[match["start"] : match["end"]] = list(
                        match["replacement"]
                    )

                    # Create correction object
                    correction = Correction(
                        position=match["start"],
                        original=match["original"],
                        replacement=match["replacement"],
                        correction_type=self.correction_type,
                        description=f"OCR: '{match['original']}' → '{match['replacement']}'",
                    )
                    # Because we process in reverse for string building,
                    # we might want to insert at 0 to keep the list in reading order,
                    # or just append and reverse later.
                    corrections_list.append(correction)

                # Reassemble string
                corrected_text = "".join(new_text_chars)

                # Sort corrections back to reading order (optional, good for UI)
                corrections_list.sort(key=lambda x: x.position)

                yield corrected_text, corrections_list

    def _check_overlap(self, sorted_reverse_matches: List[Dict]) -> bool:
        """
        Check if any matches in the combination overlap.
        Input is sorted by start index DESCENDING.
        """
        # range (start, end)
        # B is before A in the list (because reverse sort), so B has higher index
        for i in range(len(sorted_reverse_matches) - 1):
            later_match = sorted_reverse_matches[i]
            earlier_match = sorted_reverse_matches[i + 1]

            # If the earlier match ends after the later match starts -> Overlap
            if earlier_match["end"] > later_match["start"]:
                return True
        return False

    def add_substitution(self, original: str, replacements: List[str]) -> None:
        """Add a custom substitution rule dynamically."""
        if replacements:
            clean_key = original.replace(" ", "")
            self.keyword_processor.add_keyword(clean_key, replacements[0])

    @staticmethod
    def _positions_overlap(points: List[Tuple[int, str, List[str]]]) -> bool:
        """Check if any substitution positions overlap."""
        ranges = [(p[0], p[0] + len(p[1])) for p in points]
        ranges.sort()
        for i in range(len(ranges) - 1):
            if ranges[i][1] > ranges[i + 1][0]:
                return True
        return False


class TokenSegmentationStrategy(CorrectionStrategy):
    """
    Strategy for repairing segmented/broken chemical tokens using a list of valid sub-word tokens.

    It ignores 'noise' characters (hyphens, spaces, newlines) within a token
    to see if the sequence forms a valid known token.

    Example:
    - Tokens: ["dodec", "one"]
    - Input: "do-dec-2-one"
    - Logic:
      1. "do-dec" matches skeleton "dodec".
      2. "one" matches skeleton "one".
    - Result: "dodec-2-one" (The hyphens *between* tokens are preserved, hyphens *inside* tokens are removed).
    """

    @property
    def name(self) -> str:
        return "Token Segmentation"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.TOKEN_SEGMENTATION

    def __init__(self, min_token_length: int = 3):
        """
        min_token_length: Ignore tokens shorter than this to prevent aggressive
                          glueing of small words like "at" or "in".
        """
        self.keyword_processor = KeywordProcessor()
        self.keyword_processor.non_word_boundaries = set()
        self.min_token_length = min_token_length

        # Load tokens from your path
        try:
            with open(CHEMICAL_NAME_TOKENS_PATH, "rb") as f:
                chemical_tokens = json.load(f)
            chemical_tokens = [ele for ele in chemical_tokens if len(ele) > 1]
        except (NameError, FileNotFoundError):
            # Fallback for testing if constant not defined
            chemical_tokens = []

        for token in chemical_tokens:
            # We skip very short tokens to avoid false positives
            # (e.g. glueing "a n" -> "an" might not be desired in text)
            if len(token) < self.min_token_length:
                continue

            # 1. Create the skeleton (what the token looks like stripped)
            # Usually, the token list is already clean, so skeleton == token.
            # But we normalize just in case the source list has noise.
            skeleton = self._make_skeleton(token)

            # 2. Add to FlashText
            # Key: skeleton, Value: Clean Token
            self.keyword_processor.add_keyword(skeleton, token)

    def _make_skeleton(self, text: str) -> str:
        """
        Strip noise to create the search key.
        We remove hyphens, spaces, newlines to allow matching across them.
        """
        return (
            text.replace("-", "").replace("\n", "").replace(" ", "").replace("\xad", "")
        )

    def generate_candidates(
        self, text: str, config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        # 1. SQUASH & MAP
        # Create a "clean view" of the text while remembering original positions
        clean_chars = []
        idx_map = []

        # Characters we want to ignore when searching for tokens
        noise_chars = {"-", "\n", " ", "\xad", "\r", "\t"}

        for i, char in enumerate(text):
            if char not in noise_chars:
                clean_chars.append(char)
                idx_map.append(i)

        clean_text = "".join(clean_chars)

        if not clean_text:
            return

        # 2. EXTRACT matches
        # This finds valid tokens inside the squashed text
        matches = self.keyword_processor.extract_keywords(clean_text, span_info=True)

        if not matches:
            return

        mapped_matches = []

        # 3. CALCULATE original positions
        for correct_token, start, end in matches:
            orig_start_index = idx_map[start]
            last_char_idx = idx_map[end - 1]
            orig_end_index = last_char_idx + 1

            original_snippet = text[orig_start_index:orig_end_index]

            # CRITICAL: Only correct if the text is actually different.
            # If input is "dodec" and token is "dodec", do nothing.
            # If input is "do-dec", replace with "dodec".
            if original_snippet == correct_token:
                continue

            match_data = {
                "replacement": correct_token,
                "start": orig_start_index,
                "end": orig_end_index,
                "original": original_snippet,
            }
            mapped_matches.append(match_data)

        # 4. APPLY corrections (Batch)
        # Apply from end to start to maintain indices
        if mapped_matches:
            mapped_matches.sort(key=lambda x: x["start"], reverse=True)

            new_text_chars = list(text)
            corrections_list = []

            last_start = float("inf")

            for match in mapped_matches:
                # Basic overlap check
                if match["end"] > last_start:
                    continue

                last_start = match["start"]

                # Perform substitution in the character list
                new_text_chars[match["start"] : match["end"]] = list(
                    match["replacement"]
                )

                correction = Correction(
                    position=match["start"],
                    original=match["original"],
                    replacement=match["replacement"],
                    correction_type=self.correction_type,
                    description=f"Merged: '{match['original']}' → '{match['replacement']}'",
                )
                corrections_list.append(correction)

            if corrections_list:
                corrected_text = "".join(new_text_chars)
                corrections_list.sort(key=lambda x: x.position)
                yield corrected_text, corrections_list


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
        (re.compile(r"(\d)([a-zA-Z])"), r"\1-\2", "digit-letter"),
        # letter followed by digit (missing hyphen before locant)
        (re.compile(r"([a-zA-Z])(\d)"), r"\1-\2", "letter-digit"),
    ]

    MISSING_COMMA_PATTERNS: ClassVar[List[Tuple[re.Pattern, str, str]]] = [
        # Adjacent digits that should be comma-separated locants
        (re.compile(r"(\d)\s+(\d)"), r"\1,\2", "spaced-digits"),
    ]

    def generate_candidates(
        self, text: str, config: CorrectorConfig
    ) -> Iterator[Tuple[str, List[Correction]]]:
        """Generate candidates by restoring missing punctuation."""
        yield text, []  # Original

        # Try hyphen restoration
        for pattern, replacement, desc in self.MISSING_HYPHEN_PATTERNS:
            matches = list(pattern.finditer(text))

            for match in matches:
                new_text = (
                    text[: match.start()]
                    + pattern.sub(replacement, match.group())
                    + text[match.end() :]
                )

                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added hyphen ({desc})",
                    )
                    yield new_text, [correction]

        # Try comma restoration
        for pattern, replacement, desc in self.MISSING_COMMA_PATTERNS:
            matches = list(pattern.finditer(text))

            for match in matches:
                new_text = (
                    text[: match.start()]
                    + pattern.sub(replacement, match.group())
                    + text[match.end() :]
                )

                if new_text != text:
                    correction = Correction(
                        position=match.start(),
                        original=match.group(),
                        replacement=pattern.sub(replacement, match.group()),
                        correction_type=self.correction_type,
                        description=f"Added comma ({desc})",
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
        self, text: str, config: CorrectorConfig
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
                    description=f"Added {diff} closing '{close_b}'",
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
                    description=f"Added {diff} opening '{open_b}'",
                )
                yield new_text, [correction]


class LocantCorrectionStrategy(CorrectionStrategy):
    """
    Regex-based strategy for correcting chemical locants, stereochemistry,
    and punctuation errors (l->1, .->, missing commas).
    """

    @property
    def name(self) -> str:
        return "Locant Regex Correction"

    @property
    def correction_type(self) -> CorrectionType:
        return CorrectionType.LOCANT_FIX

    def generate_candidates(
        self, text: str, config: object
    ) -> Iterator[Tuple[str, List[Correction]]]:
        # We apply these patterns sequentially to the string.
        # Order matters: Fix characters -> Fix punctuation -> Fix formatting.

        # In this implementation, we apply the corrections sequentially
        # to build one final "Clean" string, yielding each step.

        current_text = text

        for pattern, replacement, desc in PATTERNS:
            # We use a loop to catch all instances of one error type (e.g. multiple 'l's)
            # before moving to the next pattern type.

            matches = list(re.finditer(pattern, current_text))
            if not matches:
                continue

            corrections_batch = []

            # We must apply replacements from right to left (reverse)
            # so indices don't shift for the remaining matches in this batch.
            for match in reversed(matches):
                start, end = match.span()
                original_str = match.group()

                # specific handling for regex group replacement (like \1,\2)
                if isinstance(replacement, str) and "\\" in replacement:
                    new_str = match.expand(replacement)
                else:
                    new_str = replacement

                # Create Correction Object
                correction = Correction(
                    position=start,
                    original=original_str,
                    replacement=new_str,
                    correction_type=self.correction_type,
                    description=desc,
                )
                corrections_batch.append(correction)

                # Mutate the string
                current_text = current_text[:start] + new_str + current_text[end:]

            # Since we processed reverse, reverse back for the UI/List
            corrections_batch.reverse()

            yield current_text, corrections_batch
