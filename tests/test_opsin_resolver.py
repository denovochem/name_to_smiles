import os
import sys

# Ensure project root is on sys.path so we can import core modules
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from placeholder_name.resolvers.opsin_resolver import (  # noqa: E402
    name_to_smiles_opsin,
)


def test_name_to_smiles_opsin_success(monkeypatch):
    """name_to_smiles_opsin should map each input name to its SMILES using OPSIN."""

    captured_args = {}

    def fake_py2opsin(
        chemical_name,
        output_format,
        failure_analysis,
        allow_acid,
        allow_radicals,
        allow_bad_stereo,
        wildcard_radicals,
    ):
        # Verify that we pass the expected options through
        captured_args["chemical_name"] = chemical_name
        captured_args["output_format"] = output_format
        captured_args["failure_analysis"] = failure_analysis
        captured_args["allow_acid"] = allow_acid
        captured_args["allow_radicals"] = allow_radicals
        captured_args["allow_bad_stereo"] = allow_bad_stereo
        captured_args["wildcard_radicals"] = wildcard_radicals

        smiles = [f"SMILES_{name}" for name in chemical_name]
        failures = [""] * len(chemical_name)
        return smiles, failures

    monkeypatch.setattr(
        "placeholder_name.resolvers.opsin_resolver.py2opsin",
        fake_py2opsin,
        raising=True,
    )

    names = ["ethanol", "water", "acetone"]
    result_smiles, result_failures = name_to_smiles_opsin(
        names,
        allow_acid=True,
        allow_radicals=False,
        allow_bad_stereo=True,
        wildcard_radicals=True,
    )

    # All names should appear as keys with the expected SMILES values
    assert set(result_smiles.keys()) == set(names)
    for name in names:
        assert result_smiles[name] == f"SMILES_{name}"

    # No failures should be recorded
    assert result_failures == {}

    # Ensure options were forwarded correctly to py2opsin
    assert captured_args["chemical_name"] == names
    assert captured_args["output_format"] == "SMILES"
    assert captured_args["failure_analysis"] is True
    assert captured_args["allow_acid"] is True
    assert captured_args["allow_radicals"] is False
    assert captured_args["allow_bad_stereo"] is True
    assert captured_args["wildcard_radicals"] is True


def test_name_to_smiles_opsin_records_failures(monkeypatch):
    """name_to_smiles_opsin should populate the failure dict when OPSIN returns messages."""

    def fake_py2opsin(chemical_name, **kwargs):
        smiles = [""] * len(chemical_name)
        failures = [f"Error for {name}" for name in chemical_name]
        return smiles, failures

    monkeypatch.setattr(
        "placeholder_name.resolvers.opsin_resolver.py2opsin",
        fake_py2opsin,
        raising=True,
    )

    names = ["bad1", "bad2"]
    result_smiles, result_failures = name_to_smiles_opsin(names)

    # No successful SMILES should be present
    assert result_smiles == {}

    # All names should have corresponding failure messages
    assert set(result_failures.keys()) == set(names)
    for name in names:
        assert result_failures[name] == f"Error for {name}"


def test_name_to_smiles_opsin_strips_newlines(monkeypatch):
    """name_to_smiles_opsin should strip newline characters before passing to py2opsin."""

    seen_chemical_names = []

    def fake_py2opsin(chemical_name, **kwargs):
        seen_chemical_names.extend(chemical_name)
        smiles = [f"SMILES_{name}" for name in chemical_name]
        failures = [""] * len(chemical_name)
        return smiles, failures

    monkeypatch.setattr(
        "placeholder_name.resolvers.opsin_resolver.py2opsin",
        fake_py2opsin,
        raising=True,
    )

    names_with_newlines = ["ethanol\n", "water\n"]
    result_smiles, result_failures = name_to_smiles_opsin(names_with_newlines)

    # py2opsin should see names without newlines
    assert seen_chemical_names == ["ethanol", "water"]

    # Returned mapping should still use the original input keys
    assert set(result_smiles.keys()) == set(names_with_newlines)
    assert result_failures == {}


def test_name_to_smiles_opsin_mismatched_lengths_logs_warning_and_returns_empty(
    monkeypatch,
):
    """If OPSIN returns mismatched lengths, the function should log a warning and return empty dicts."""

    def fake_py2opsin(chemical_name, **kwargs):
        # Deliberately return mismatched lengths
        smiles = ["SMILES_only_one"]
        failures = ["err1", "err2"]
        return smiles, failures

    monkeypatch.setattr(
        "placeholder_name.resolvers.opsin_resolver.py2opsin",
        fake_py2opsin,
        raising=True,
    )

    warnings = []

    class FakeLogger:
        @staticmethod
        def warning(msg):  # pragma: no cover - behavior validated via side effects
            warnings.append(msg)

    monkeypatch.setattr(
        "placeholder_name.resolvers.opsin_resolver.logger",
        FakeLogger,
        raising=True,
    )

    names = ["n1", "n2"]
    result_smiles, result_failures = name_to_smiles_opsin(names)

    # Both dicts should be empty on mismatch
    assert result_smiles == {}
    assert result_failures == {}

    # A warning should have been logged mentioning mismatching lengths
    assert len(warnings) == 1
    assert "Mismatching lengths" in warnings[0]
