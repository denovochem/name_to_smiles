import os
import sys

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.resolvers.cirpy_resolver import (  # noqa: E402
    retrieve_cirpy_results,
    name_to_smiles_cirpy,
)


def test_retrieve_cirpy_results_success(monkeypatch):
    """retrieve_cirpy_results should return the SMILES string from cirpy.resolve."""

    def fake_resolve(compound_name, identifier_type):
        assert identifier_type == "smiles"
        return f"SMILES_{compound_name}"

    monkeypatch.setattr(
        "placeholder_name.resolvers.cirpy_resolver.cirpy.resolve",
        fake_resolve,
        raising=True,
    )

    compound_name = "ethanol"
    result = retrieve_cirpy_results(compound_name)

    assert result == "SMILES_ethanol"


def test_retrieve_cirpy_results_exception_returns_empty_string(monkeypatch):
    """If cirpy.resolve raises, retrieve_cirpy_results should return an empty string."""

    def fake_resolve(
        *args, **kwargs
    ):  # pragma: no cover - behavior validated via return value
        raise RuntimeError("Test error from cirpy")

    monkeypatch.setattr(
        "placeholder_name.resolvers.cirpy_resolver.cirpy.resolve",
        fake_resolve,
        raising=True,
    )

    result = retrieve_cirpy_results("any_name")

    assert result == ""


def test_name_to_smiles_cirpy_uses_retrieve_for_each_name(monkeypatch):
    """name_to_smiles_cirpy should map each input name to its SMILES using CIRpy."""

    def fake_resolve(compound_name, identifier_type):
        assert identifier_type == "smiles"
        return f"SMILES_{compound_name}"

    monkeypatch.setattr(
        "placeholder_name.resolvers.cirpy_resolver.cirpy.resolve",
        fake_resolve,
        raising=True,
    )

    names = ["ethanol", "water", "acetone"]
    result_dict = name_to_smiles_cirpy(names)

    # All names should appear as keys with the expected SMILES values
    assert set(result_dict.keys()) == set(names)
    for name in names:
        assert result_dict[name] == f"SMILES_{name}"
