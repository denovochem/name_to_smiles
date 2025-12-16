import os
import sys

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.resolvers.manual_resolver import (  # noqa: E402
    load_default_manual_name_dict,
    process_name_dict,
    name_to_smiles_manual,
)


def test_process_name_dict_basic_behavior():
    """process_name_dict should match names case-insensitively and strip whitespace."""
    compound_names = [" Ethanol ", "WATER", "Acetone"]
    name_dict = {
        "ethanol": "C2H6O",
        " water  ": "H2O",
        "acetone": "C3H6O",
        # Unused or empty entries should not appear in the result
        "methane": "",
    }

    result = process_name_dict(compound_names, name_dict)

    # Keys in the result should preserve the original input formatting
    assert set(result.keys()) == {" Ethanol ", "WATER", "Acetone"}
    assert result[" Ethanol "] == "C2H6O"
    assert result["WATER"] == "H2O"
    assert result["Acetone"] == "C3H6O"


def test_name_to_smiles_manual_with_provided_dict(monkeypatch):
    """name_to_smiles_manual should use the provided dict and not call the default loader."""

    # Ensure the default loader is not used when a dict is provided
    def fake_loader():  # pragma: no cover - validated by not being called
        raise AssertionError(
            "load_default_manual_name_dict should not be called when provided_name_dict is given"
        )

    monkeypatch.setattr(
        "placeholder_name.resolvers.manual_resolver.load_default_manual_name_dict",
        fake_loader,
        raising=True,
    )

    names = ["Ethanol", "Water"]
    provided_dict = {
        "ethanol": "C2H6O",
        "water": "H2O",
    }

    result = name_to_smiles_manual(names, provided_name_dict=provided_dict)

    assert result == {"Ethanol": "C2H6O", "Water": "H2O"}


def test_name_to_smiles_manual_uses_default_dict(monkeypatch):
    """When no dict is provided, name_to_smiles_manual should use load_default_manual_name_dict."""

    # Provide a small fake default dictionary via monkeypatch
    def fake_loader():
        return {
            "ethanol": "C2H6O",
            "water": "H2O",
        }

    monkeypatch.setattr(
        "placeholder_name.resolvers.manual_resolver.load_default_manual_name_dict",
        fake_loader,
        raising=True,
    )

    # Also monkeypatch logger.info so we don't depend on real logging output
    log_messages = []

    def fake_logger_info(msg):
        log_messages.append(msg)

    monkeypatch.setattr(
        "placeholder_name.resolvers.manual_resolver.logger",
        type("_L", (), {"info": staticmethod(fake_logger_info)}),
        raising=True,
    )

    names = ["Ethanol", "Water"]
    result = name_to_smiles_manual(names)

    assert result == {"Ethanol": "C2H6O", "Water": "H2O"}

    # We expect at least one info log call to have been made
    assert any("Loaded data from" in msg for msg in log_messages)
