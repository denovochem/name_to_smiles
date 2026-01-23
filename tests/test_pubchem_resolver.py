import os
import sys

# Ensure project root is on sys.path so we can import cholla_chem modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cholla_chem.resolvers.pubchem_resolver import (  # noqa: E402
    name_to_smiles_pubchem,
)


def test_name_to_smiles_pubchem_basic_mapping(monkeypatch):
    """name_to_smiles_pubchem should map each name to the SMILES returned by PubChem."""

    # Make filter_latin1_compatible a pass-through so we control the argument into get_compounds
    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    captured_args = {}

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    def fake_get_compounds(names, identifier_type):
        captured_args["names"] = names
        captured_args["identifier_type"] = identifier_type
        # Return matching FakeCompound objects
        return [FakeCompound(f"SMILES_{name}") for name in names]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pcp.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    compound_names = ["ethanol", "water", "acetone"]
    result = name_to_smiles_pubchem(compound_names)

    # Ensure get_compounds was called with the expected arguments
    assert captured_args["names"] == compound_names
    assert captured_args["identifier_type"] == "name"

    # Ensure the result dictionary contains all names with the expected SMILES
    assert set(result.keys()) == set(compound_names)
    for name in compound_names:
        assert result[name] == f"SMILES_{name}"


def test_name_to_smiles_pubchem_handles_none_results(monkeypatch):
    """If PubChem returns None for an entry, the corresponding SMILES should be an empty string."""

    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    def fake_get_compounds(names, identifier_type):
        # Return a list with one valid compound and one None
        return [FakeCompound("SMILES_valid"), None]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pcp.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    compound_names = ["valid_name", "missing_name"]
    result = name_to_smiles_pubchem(compound_names)

    assert result["valid_name"] == "SMILES_valid"
    assert result["missing_name"] == ""


def test_name_to_smiles_pubchem_logs_length_mismatch(monkeypatch):
    """If lengths mismatch, a warning should be logged."""

    def fake_filter_latin1_compatible(names):
        return names

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.filter_latin1_compatible",
        fake_filter_latin1_compatible,
        raising=True,
    )

    class FakeCompound:
        def __init__(self, smiles):
            self.smiles = smiles

    # Return fewer results than input names to trigger the warning
    def fake_get_compounds(names, identifier_type):
        return [FakeCompound("SMILES_only_first")]

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.pcp.get_compounds",
        fake_get_compounds,
        raising=True,
    )

    logged_warnings = []

    class FakeLogger:
        @staticmethod
        def warning(msg):
            logged_warnings.append(msg)

    monkeypatch.setattr(
        "cholla_chem.resolvers.pubchem_resolver.logger",
        FakeLogger,
        raising=True,
    )

    compound_names = ["first", "second"]
    result = name_to_smiles_pubchem(compound_names)

    # Because of zip, only the first name will have an entry
    assert set(result.keys()) == {"first"}
    assert result["first"] == "SMILES_only_first"

    # A warning about mismatching lengths should have been logged
    assert any("Mismatching lengths" in msg for msg in logged_warnings)
