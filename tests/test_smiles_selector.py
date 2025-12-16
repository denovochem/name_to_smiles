import os
import sys

import pytest

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.smiles_selector import SMILESSelector  # noqa: E402


def _make_selector(
    smiles_dict, weight_dict=None, priority_order=None, custom_strategies=None
):
    data = {
        "compound1": {
            "SMILES_dict": smiles_dict,
        }
    }
    return SMILESSelector(
        data=data,
        weight_dict=weight_dict or {},
        priority_order=priority_order or [],
        custom_strategies=custom_strategies,
    )


def test_select_smiles_default_uses_consensus():
    smiles_dict = {
        "CCO": ["r1", "r2"],
        "O": ["r1"],
    }
    selector = _make_selector(smiles_dict)

    smiles, sources = selector.select_smiles("compound1")

    assert smiles == "CCO"
    assert sources == ["r1", "r2"]


def test_strategy_consensus_tie_breaker_lexicographic():
    smiles_dict = {
        "CC": ["r1"],
        "CN": ["r2"],
    }
    selector = _make_selector(smiles_dict)

    smiles, sources = selector._strategy_consensus(smiles_dict)

    assert smiles == "CC"
    assert sources == ["r1"]


def test_strategy_ordered_priority_picks_first_matching_resolver():
    smiles_dict = {
        "CCO": ["resolverA"],
        "O": ["resolverB"],
    }
    priority_order = ["resolverB", "resolverA"]
    selector = _make_selector(smiles_dict, priority_order=priority_order)

    smiles, sources = selector.select_smiles("compound1", strategy="ordered")

    assert smiles == "O"
    assert sources == ["resolverB"]


def test_strategy_ordered_priority_falls_back_to_any_smiles_if_no_match():
    smiles_dict = {
        "CCO": ["resolverA"],
        "O": ["resolverC"],
    }
    priority_order = ["resolverB"]
    selector = _make_selector(smiles_dict, priority_order=priority_order)

    smiles, sources = selector.select_smiles("compound1", strategy="ordered")

    assert smiles in smiles_dict
    assert sources == smiles_dict[smiles]


def test_strategy_weighted_consensus_uses_weights():
    smiles_dict = {
        "CCO": ["A"],
        "O": ["B", "B"],
    }
    weight_dict = {"A": 2.0, "B": 1.0}
    selector = _make_selector(smiles_dict, weight_dict=weight_dict)

    smiles, sources = selector.select_smiles("compound1", strategy="weighted")

    assert smiles == "CCO"
    assert sources == ["A"]


def test_strategy_weighted_consensus_tie_breaker_lexicographic():
    smiles_dict = {
        "CC": ["A"],
        "CN": ["B"],
    }
    weight_dict = {"A": 1.0, "B": 1.0}
    selector = _make_selector(smiles_dict, weight_dict=weight_dict)

    smiles, _ = selector.select_smiles("compound1", strategy="weighted")

    assert smiles == "CC"


def test_strategy_shortest_smiles_current_behavior_lexicographic_min():
    smiles_dict = {
        "CCO": ["r1"],
        "CO": ["r2"],
        "CN": ["r3"],
    }
    selector = _make_selector(smiles_dict)

    smiles, sources = selector.select_smiles("compound1", strategy="shortest_smiles")

    assert smiles == min(smiles_dict.keys())
    assert sources == smiles_dict[smiles]


def test_strategy_longest_smiles_current_behavior_lexicographic_max():
    smiles_dict = {
        "CCO": ["r1"],
        "CO": ["r2"],
        "CN": ["r3"],
    }
    selector = _make_selector(smiles_dict)

    smiles, sources = selector.select_smiles("compound1", strategy="longest_smiles")

    assert smiles == max(smiles_dict.keys())
    assert sources == smiles_dict[smiles]


def test_strategy_fewest_fragments():
    smiles_dict = {
        "C.O": ["r1"],
        "C.O.N": ["r2"],
        "CO": ["r3"],
    }
    selector = _make_selector(smiles_dict)

    smiles, sources = selector.select_smiles("compound1", strategy="fewest_fragments")

    assert smiles == "CO"
    assert sources == ["r3"]


def test_select_smiles_filters_out_invalid_and_empty_smiles(monkeypatch):
    # Construct entries that should be filtered out: None, empty, whitespace, invalid SMILES
    smiles_dict = {
        None: ["r1"],
        "": ["r2"],
        "   ": ["r3"],
        "not_a_smiles": ["r4"],
    }

    selector = _make_selector(smiles_dict)

    # Use the real RDKit Chem.MolFromSmiles via the implementation; invalid entries should be removed
    smiles, sources = selector.select_smiles("compound1")

    assert smiles == ""
    assert sources == []


def test_select_smiles_raises_for_unknown_strategy():
    smiles_dict = {"CCO": ["r1"]}
    selector = _make_selector(smiles_dict)

    with pytest.raises(ValueError):
        selector.select_smiles("compound1", strategy="nonexistent_strategy")


def test_select_smiles_with_callable_strategy():
    smiles_dict = {"CCO": ["r1"]}

    def custom_strategy(smiles_dict_arg, **kwargs):
        return "CUSTOM", ["custom_source"]

    selector = _make_selector(smiles_dict)

    smiles, sources = selector.select_smiles("compound1", strategy=custom_strategy)

    assert smiles == "CUSTOM"
    assert sources == ["custom_source"]


def test_custom_strategy_registered_by_name():
    smiles_dict = {"CCO": ["r1"]}

    def custom_strategy(smiles_dict_arg, **kwargs):
        return "CCO", ["custom_source"]

    selector = _make_selector(
        smiles_dict, custom_strategies={"custom": custom_strategy}
    )

    smiles, sources = selector.select_smiles("compound1", strategy="custom")

    assert smiles == "CCO"
    assert sources == ["custom_source"]
