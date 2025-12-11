# import os
# import sys

# import pytest

# # Ensure project root is on sys.path so we can import core modules
# PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
# if PROJECT_ROOT not in sys.path:
#     sys.path.insert(0, PROJECT_ROOT)

# from placeholder_name.main import (  # noqa: E402
#     ChemicalNameResolver,
#     assemble_compounds_resolution_dict,
#     assemble_split_compounds_resolution_dict,
#     clean_strings_and_return_mapping,
#     get_resolvers_weight_dict,
#     resolve_compounds_to_smiles,
#     resolve_compounds_using_resolvers,
#     select_smiles_with_criteria,
#     split_compounds_on_delimiters_and_return_mapping,
# )


# class DummyResolver(ChemicalNameResolver):
#     """Simple in-memory resolver used for orchestration tests."""

#     def __init__(self, resolver_name, mapping, info_mapping=None, weight=1.0):
#         super().__init__("dummy", resolver_name, weight)
#         self._mapping = mapping
#         self._info_mapping = info_mapping or {}

#     def name_to_smiles(self, compound_name_list):
#         out = {name: self._mapping.get(name, "") for name in compound_name_list}
#         info = {name: self._info_mapping.get(name, "") for name in compound_name_list}
#         # Drop empty info messages, mimicking behavior of real resolvers
#         info = {k: v for k, v in info.items() if v}
#         return out, info


# def test_clean_strings_and_return_mapping_basic(monkeypatch):
#     """clean_strings_and_return_mapping should call clean_strings on each item and map originals to cleaned."""

#     from placeholder_name import main as main_module

#     calls = []

#     def fake_clean_strings(s):
#         calls.append(s)
#         return s.strip().lower()

#     monkeypatch.setattr(
#         "placeholder_name.utils.string_utils.clean_strings",
#         fake_clean_strings,
#         raising=True,
#     )

#     original = [" Ethanol ", "WATER"]
#     cleaned_list, mapping = clean_strings_and_return_mapping(original)

#     assert calls == original
#     assert cleaned_list == ["ethanol", "water"]
#     assert mapping == {" Ethanol ": "ethanol", "WATER": "water"}


# def test_get_resolvers_weight_dict_uses_resolver_name_and_weight():
#     """get_resolvers_weight_dict should pull weights from resolver instances."""

#     r1 = DummyResolver("r1", {}, weight=1.5)
#     r2 = DummyResolver("r2", {}, weight=3.0)

#     result = get_resolvers_weight_dict([r1, r2])

#     assert result == {"r1": 1.5, "r2": 3.0}


# def test_split_compounds_on_delimiters_and_return_mapping_uses_helper(monkeypatch):
#     """split_compounds_on_delimiters_and_return_mapping should call get_delimiter_split_dict for each compound."""

#     seen = []

#     def fake_get_delimiter_split_dict(compound, current_dict):
#         seen.append(compound)
#         # For testing, pretend every compound splits into two parts
#         split_parts = [f"{compound}_a", f"{compound}_b"]
#         current_dict[compound] = split_parts
#         return current_dict, split_parts

#     monkeypatch.setattr(
#         "placeholder_name.name_manipulation.split_names.get_delimiter_split_dict",
#         fake_get_delimiter_split_dict,
#         raising=True,
#     )

#     compounds = ["A", "B"]
#     all_names, split_dict = split_compounds_on_delimiters_and_return_mapping(compounds)

#     # Helper should be called once per original compound
#     assert seen == compounds

#     # All originals plus unique split parts should be present
#     assert set(all_names) == {"A", "B", "A_a", "A_b", "B_a", "B_b"}

#     # Mapping should contain the split parts we defined
#     assert split_dict == {
#         "A": ["A_a", "A_b"],
#         "B": ["B_a", "B_b"],
#     }


# def test_resolve_compounds_using_resolvers_batches_and_collects_info():
#     """resolve_compounds_using_resolvers should iterate over resolvers and batches, returning out and info dicts."""

#     # Two dummy resolvers with different names and mappings
#     r1 = DummyResolver("res1", {"a": "S_a1", "b": "S_b1"}, {"a": "i1"})
#     r2 = DummyResolver("res2", {"a": "S_a2"}, {"b": "i2"})

#     compounds = ["a", "b"]
#     # Use small batch size to exercise batching loop; behavior should be identical
#     out = resolve_compounds_using_resolvers(compounds, [r1, r2], batch_size=1)

#     assert set(out.keys()) == {"res1", "res2"}
#     assert out["res1"]["out"] == {"a": "S_a1", "b": "S_b1"}
#     assert out["res1"]["info_messages"] == {"a": "i1"}
#     assert out["res2"]["out"] == {"a": "S_a2", "b": ""}
#     assert out["res2"]["info_messages"] == {"b": "i2"}


# def test_assemble_compounds_resolution_dict_uses_canonical_smiles(monkeypatch):
#     """assemble_compounds_resolution_dict should canonicalize SMILES and group resolvers per canonical form."""

#     def fake_canonicalize(smiles):
#         # Normalize to uppercase as a simple canonicalization
#         return smiles.upper() if smiles else ""

#     monkeypatch.setattr(
#         "placeholder_name.utils.chem_utils.canonicalize_smiles",
#         fake_canonicalize,
#         raising=True,
#     )

#     compounds = ["ethanol"]
#     cleaned_mapping = {"ethanol": "ethanol_clean"}

#     resolvers_out = {
#         "r1": {"out": {"ethanol_clean": "c2h6o"}, "info_messages": {}},
#         "r2": {"out": {"ethanol_clean": "C2H6O"}, "info_messages": {"ethanol_clean": "ok"}},
#     }

#     result = assemble_compounds_resolution_dict(compounds, resolvers_out, cleaned_mapping)

#     assert set(result.keys()) == {"ethanol"}
#     entry = result["ethanol"]
#     assert entry["SMILES"] == ""  # not yet selected
#     assert entry["SMILES_dict"] == {"C2H6O": ["r1", "r2"]}
#     assert entry["info_messages"] == {"r2": "ok"}


# def test_assemble_split_compounds_resolution_dict_merges_split_smiles(monkeypatch):
#     """assemble_split_compounds_resolution_dict should merge SMILES from split parts via resolve_delimiter_split_dict."""

#     def fake_canonicalize(smiles):
#         return smiles

#     monkeypatch.setattr(
#         "placeholder_name.utils.chem_utils.canonicalize_smiles",
#         fake_canonicalize,
#         raising=True,
#     )

#     # resolvers_out_dict not used directly by our fake but required by API
#     resolvers_out_dict = {}
#     cleaned_mapping = {"mix": "mix_clean"}
#     delimiter_split_dict = {"mix_clean": ["a", "b"]}

#     # Start with an entry that has an existing SMILES_dict
#     compounds_out_dict = {
#         "mix": {
#             "SMILES": "",
#             "SMILES_source": [],
#             "SMILES_dict": {},
#             "info_messages": {},
#         }
#     }

#     def fake_resolve_delimiter_split_dict(compound_cleaned, ro_dict, del_dict):
#         assert compound_cleaned == "mix_clean"
#         assert del_dict is delimiter_split_dict
#         # Pretend the split parts combine into two possible SMILES
#         return {
#             "S_mix1": ["r_split1"],
#             "S_mix2": ["r_split2", "r_split3"],
#         }

#     monkeypatch.setattr(
#         "placeholder_name.name_manipulation.split_names.resolve_delimiter_split_dict",
#         fake_resolve_delimiter_split_dict,
#         raising=True,
#     )

#     result = assemble_split_compounds_resolution_dict(
#         compounds_out_dict,
#         ["mix"],
#         resolvers_out_dict,
#         cleaned_mapping,
#         delimiter_split_dict,
#     )

#     entry = result["mix"]
#     assert entry["SMILES_dict"] == {
#         "S_mix1": ["r_split1"],
#         "S_mix2": ["r_split2", "r_split3"],
#     }


# def test_select_smiles_with_criteria_uses_selector(monkeypatch):
#     """select_smiles_with_criteria should delegate to SMILESSelector and populate SMILES and SMILES_source."""

#     # Prepare a minimal compounds_out_dict with one compound and two candidate SMILES
#     compounds_out = {
#         "ethanol": {
#             "SMILES": "",
#             "SMILES_source": [],
#             "SMILES_dict": {
#                 "S1": ["r1"],
#                 "S2": ["r2"],
#             },
#             "info_messages": {},
#         }
#     }

#     resolvers_weight = {"r1": 1.0, "r2": 2.0}
#     priority_order = ["r2", "r1"]

#     selected_calls = []

#     class FakeSelector:
#         def __init__(self, compounds_out_dict, weight_dict, priority_list):  # pragma: no cover - wiring only
#             selected_calls.append((compounds_out_dict, weight_dict, priority_list))

#         def select_smiles(self, compound, mode):
#             # Always pick S2 for this test
#             assert compound == "ethanol"
#             assert mode == "weighted"
#             return "S2", ["r2"]

#     monkeypatch.setattr(
#         "placeholder_name.main.SMILESSelector",
#         FakeSelector,
#         raising=True,
#     )

#     result = select_smiles_with_criteria(
#         compounds_out,
#         resolvers_weight,
#         priority_order,
#         smiles_selection_mode="weighted",
#     )

#     assert selected_calls  # FakeSelector was constructed
#     assert result["ethanol"]["SMILES"] == "S2"
#     assert result["ethanol"]["SMILES_source"] == ["r2"]


# def test_resolve_compounds_to_smiles_validates_inputs_and_resolvers(monkeypatch):
#     """resolve_compounds_to_smiles should validate input types and resolver instances."""

#     # Invalid compounds_list type
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(123, resolvers_list=[DummyResolver("r", {})])

#     # Empty list
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles([], resolvers_list=[DummyResolver("r", {})])

#     # Non-string in list
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(["ok", 1], resolvers_list=[DummyResolver("r", {})])

#     # resolvers_list must be non-empty list of ChemicalNameResolver
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(["a"], resolvers_list=[])

#     class NotAResolver:
#         pass

#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(["a"], resolvers_list=[NotAResolver()])

#     # Duplicate resolver names should raise
#     r1 = DummyResolver("dup", {})
#     r2 = DummyResolver("dup", {})
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(["a"], resolvers_list=[r1, r2])

#     # smiles_selection_mode must be str or callable
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             smiles_selection_mode=123,
#         )

#     # detailed_name_dict must be bool
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             detailed_name_dict="yes",  # type: ignore[arg-type]
#         )

#     # batch_size must be int between 1 and 1000
#     with pytest.raises(TypeError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             batch_size=1.5,  # type: ignore[arg-type]
#         )

#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             batch_size=0,
#         )

#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             batch_size=1001,
#         )

#     # split_names_to_solve must be bool
#     with pytest.raises(ValueError):
#         resolve_compounds_to_smiles(
#             ["a"],
#             resolvers_list=[DummyResolver("r", {})],
#             split_names_to_solve="yes",  # type: ignore[arg-type]
#         )


# def test_resolve_compounds_to_smiles_happy_path_minimal(monkeypatch):
#     """End-to-end: with a simple resolver and selector, should return final SMILES mapping."""

#     # Use a DummyResolver that returns a direct mapping
#     resolver = DummyResolver("dummy_res", {"a": "S_a"})

#     # Avoid depending on real cleaning/splitting behavior: make them pass-through
#     monkeypatch.setattr(
#         "placeholder_name.main.clean_strings_and_return_mapping",
#         lambda names: (names, {n: n for n in names}),
#         raising=True,
#     )

#     monkeypatch.setattr(
#         "placeholder_name.main.split_compounds_on_delimiters_and_return_mapping",
#         lambda names: (names, {}),
#         raising=True,
#     )

#     # Canonicalization is also identity for this test
#     monkeypatch.setattr(
#         "placeholder_name.utils.chem_utils.canonicalize_smiles",
#         lambda s: s,
#         raising=True,
#     )

#     # Simplify SMILESSelector behavior: always select the only available SMILES if present
#     class FakeSelector:
#         def __init__(self, compounds_out_dict, weight_dict, priority_list):  # pragma: no cover - wiring only
#             self._d = compounds_out_dict

#         def select_smiles(self, compound, mode):
#             entry = self._d[compound]
#             # Pick the first key from SMILES_dict, or empty if none
#             smiles = next(iter(entry["SMILES_dict"].keys()), "")
#             return smiles, entry["SMILES_dict"].get(smiles, [])

#     monkeypatch.setattr(
#         "placeholder_name.main.SMILESSelector",
#         FakeSelector,
#         raising=True,
#     )

#     result = resolve_compounds_to_smiles(["a"], resolvers_list=[resolver])

#     assert result == {"a": "S_a"}
