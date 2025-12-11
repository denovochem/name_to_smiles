import os
import sys

import pytest

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.name_manipulation.split_names import get_delimiter_split_dict  # noqa: E402
from placeholder_name.utils.constants import DELIMITERS  # noqa: E402


def test_no_delimiters_in_name_returns_empty_list_and_no_entry():
    compound_name = "SodiumChloride"
    split_dict = {}

    result_dict, split_list = get_delimiter_split_dict(compound_name, split_dict)

    assert result_dict == {}
    assert split_list == []


def test_single_delimiter_splits_once_and_records_entry():
    compound_name = "NaCl/H2O"
    split_dict = {}

    result_dict, split_list = get_delimiter_split_dict(compound_name, split_dict)

    # Expect split on '/' only, as defined in DELIMITERS
    expected_parts = compound_name.split("/")
    assert result_dict == {compound_name: expected_parts}
    assert split_list == expected_parts


def test_multiple_same_delimiter_keeps_all_parts_and_records_once():
    compound_name = "A/B/C/D"
    split_dict = {}

    result_dict, split_list = get_delimiter_split_dict(compound_name, split_dict)

    expected_parts = compound_name.split("/")
    assert result_dict == {compound_name: expected_parts}
    assert split_list == expected_parts


def test_multiple_different_delimiters_accumulate_parts_in_order():
    # Use several delimiters from DELIMITERS to ensure coverage
    # For example: '/', ':', '.', '*'
    name1 = "NaCl/H2O"     # uses '/'
    name2 = "CuSO4:5H2O"  # uses ':'
    name3 = "Fe2O3.H2O"   # uses '.'

    split_dict = {}

    # Each call should return only the parts for the name processed in that call,
    # while the shared dict accumulates entries across calls.
    dict_after_1, list_after_1 = get_delimiter_split_dict(name1, split_dict)
    dict_after_2, list_after_2 = get_delimiter_split_dict(name2, dict_after_1)
    dict_after_3, list_after_3 = get_delimiter_split_dict(name3, dict_after_2)

    # Dictionary should contain entries for all processed names with their respective parts
    assert dict_after_3[name1] == name1.split("/")
    assert dict_after_3[name2] == name2.split(":")
    assert dict_after_3[name3] == name3.split(".")

    # Each returned list should only contain parts from the name for that specific call
    assert list_after_1 == name1.split("/")
    assert list_after_2 == name2.split(":" )
    assert list_after_3 == name3.split(".")


def test_existing_dict_is_updated_not_replaced():
    existing_name = "Existing/Compound"
    existing_parts = existing_name.split("/")
    split_dict = {existing_name: existing_parts}

    new_name = "New•Compound"  # '•' is in DELIMITERS

    result_dict, split_list = get_delimiter_split_dict(new_name, split_dict)

    # Existing entry should still be present
    assert existing_name in result_dict
    assert result_dict[existing_name] == existing_parts

    # New entry should also be added
    assert new_name in result_dict
    assert result_dict[new_name] == new_name.split("•")

    # Returned list should contain only the parts from the processed name
    assert split_list == new_name.split("•")
