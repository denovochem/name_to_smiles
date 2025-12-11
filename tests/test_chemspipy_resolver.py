import os
import sys

import pytest

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.resolvers.chemspipy_resolver import (  # noqa: E402
    name_to_smiles_chemspipy,
)


def test_name_to_smiles_chemspipy_basic_mapping(monkeypatch):
    """name_to_smiles_chemspipy should map each name to the SMILES returned by ChemSpiPy."""

    captured = {
        "api_key": None,
        "search_terms": [],
    }

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    class FakeResults:
        def __init__(self, compound_name):
            self._items = [FakeCompound(f"SMILES_{compound_name}")]

        def wait(self):
            return None

        def ready(self):
            return True

        def __len__(self):
            return len(self._items)

        def __getitem__(self, index):
            return self._items[index]

    class FakeChemSpider:
        def __init__(self, api_key):
            captured["api_key"] = api_key

        def search(self, compound_name):
            captured["search_terms"].append(compound_name)
            return FakeResults(compound_name)

    def fake_filter_latin1_compatible(names):
        # Pass-through so we can easily assert on which names were searched
        return names

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.ChemSpider",
        FakeChemSpider,
        raising=True,
    )

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    compound_names = ["ethanol", "water", "acetone"]
    api_key = "TEST_KEY"

    result = name_to_smiles_chemspipy(compound_names, api_key)

    # ChemSpider should have been constructed with our API key
    assert captured["api_key"] == api_key

    # Each name should have been searched exactly once
    assert captured["search_terms"] == compound_names

    # Result dict should contain all names with the expected SMILES
    assert set(result.keys()) == set(compound_names)
    for name in compound_names:
        assert result[name] == f"SMILES_{name}"


def test_name_to_smiles_chemspipy_respects_filter_latin1(monkeypatch):
    """Only names returned by filter_latin1_compatible should be queried and appear in the result."""

    captured_search_terms = []

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    class FakeResults:
        def __init__(self, compound_name):
            self._items = [FakeCompound(f"SMILES_{compound_name}")]

        def wait(self):
            return None

        def ready(self):
            return True

        def __len__(self):
            return len(self._items)

        def __getitem__(self, index):
            return self._items[index]

    class FakeChemSpider:
        def __init__(self, api_key):  # pragma: no cover - behavior validated via side effects
            self.api_key = api_key

        def search(self, compound_name):
            captured_search_terms.append(compound_name)
            return FakeResults(compound_name)

    original_names = ["ethanol", "naïve", "water"]
    filtered_names = ["ethanol", "water"]

    def fake_filter_latin1_compatible(names):
        # Simulate dropping non-Latin-1-compatible entries
        assert names == original_names
        return filtered_names

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.ChemSpider",
        FakeChemSpider,
        raising=True,
    )

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    result = name_to_smiles_chemspipy(original_names, "KEY")

    # Only filtered names should be searched and appear in the result
    assert set(captured_search_terms) == set(filtered_names)
    assert set(result.keys()) == set(filtered_names)


def test_name_to_smiles_chemspipy_skips_empty_results_and_missing_smiles(monkeypatch):
    """Entries with no results or missing SMILES should not appear in the output dict."""

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    class EmptyResults:
        def wait(self):
            return None

        def ready(self):
            return True

        def __len__(self):
            return 0

    class SingleResult:
        def __init__(self, smiles):
            self._items = [FakeCompound(smiles)]

        def wait(self):
            return None

        def ready(self):
            return True

        def __len__(self):
            return len(self._items)

        def __getitem__(self, index):
            return self._items[index]

    def fake_filter_latin1_compatible(names):
        return names

    class FakeChemSpider:
        def __init__(self, api_key):  # pragma: no cover - construction side effects not needed
            self.api_key = api_key

        def search(self, compound_name):
            if compound_name == "no_results":
                return EmptyResults()
            if compound_name == "no_smiles":
                # Result has a compound but with empty smiles
                return SingleResult("")
            # Normal successful compound
            return SingleResult(f"SMILES_{compound_name}")

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.ChemSpider",
        FakeChemSpider,
        raising=True,
    )

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    names = ["no_results", "no_smiles", "good"]
    result = name_to_smiles_chemspipy(names, "KEY")

    # Only the "good" compound should be present
    assert set(result.keys()) == {"good"}
    assert result["good"] == "SMILES_good"


def test_name_to_smiles_chemspipy_init_failure_returns_empty(monkeypatch):
    """If ChemSpider construction fails, the function should return an empty dict."""

    class FailingChemSpider:
        def __init__(self, api_key):  # pragma: no cover - behavior validated via return value
            raise RuntimeError("Failed to initialize ChemSpider")

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.ChemSpider",
        FailingChemSpider,
        raising=True,
    )

    # filter_latin1_compatible should not matter here, but keep behavior simple
    def fake_filter_latin1_compatible(names):  # pragma: no cover - not reached when init fails
        return names

    monkeypatch.setattr(
        "placeholder_name.resolvers.chemspipy_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    result = name_to_smiles_chemspipy(["ethanol"], "KEY")

    assert result == {}
