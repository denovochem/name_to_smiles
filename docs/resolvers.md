placeholder_name uses a variety of resolvers to convert natural language chemical names to SMILES. Resolvers can be initialized and passed to the function resolve_compounds_to_smiles to customize how compounds are resolved to SMILES:

```
from placeholder_name import OpsinNameResolver, PubChemNameResolver, CIRpyNameResolver

opsin_resolver = OpsinNameResolver('opsin', resolver_weight=4)
pubchem_resolver =  PubChemNameResolver('pubchem', resolver_weight=3)
cirpy_resolver = CIRpyNameResolver('cirpy', resolver_weight=2)

resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'], [opsin_resolver, pubchem_resolver, cirpy_resolver])
```

## OpsinNameResolver
This resolver uses a fork of the [py2opsin](https://github.com/csnbritt/py2opsin) library that returns the error message from OPSIN if a name cannot be resolved. This resolver can be configured with the following arguments:

Arguments:
- allow_acid (bool, optional): Allow interpretation of acids. Defaults to False.
- allow_radicals (bool, optional): Enable radical interpretation. Defaults to False.
- allow_bad_stereo (bool, optional): Allow OPSIN to ignore uninterpreatable stereochem. Defaults to False.
- wildcard_radicals (bool, optional): Output radicals as wildcards. Defaults to False.

Default weight for 'weighted' SMILES selection: 3

```
from placeholder_name import OpsinNameResolver
opsin_resolver = OpsinNameResolver(
                                    'opsin',
                                    allow_acid=False,
                                    allow_radicals: True,
                                    allow_bad_stereo: False,
                                    wildcard_radicals: False
                                )

resolved_smiles = resolve_compounds_to_smiles(['2-acetyloxybenzoic acid'], [opsin_resolver])
```

## PubChemNameResolver
This resolver uses a fork of the [PubChemPy](https://github.com/csnbritt/PubChemPy) library which implements batching with the Power User Gateway XML schema to significantly speedup SMILES resolutions.

Default weight for 'weighted' SMILES selection: 2

```
from placeholder_name import PubChemNameResolver
pubchem_resolver = PubChemNameResolver('pubchem')

resolved_smiles = resolve_compounds_to_smiles(['acetone'], [pubchem_resolver])
```

## CIRpyNameResolver
This resolver uses the python library [CIRpy](https://github.com/mcs07/CIRpy), a Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.

Default weight for 'weighted' SMILES selection: 1

```
from placeholder_name import CIRpyNameResolver
cirpy_resolver = CIRpyNameResolver('cirpy')

resolved_smiles = resolve_compounds_to_smiles(['acetone'], [cirpy_resolver])
```

## ManualNameResolver 
This resolver uses a dataset of manually curated names and their corresponding SMILES, especially focused on common names that are incorrectly resolved by other resolvers (e.g. 'NaH').

Default weight for 'weighted' SMILES selection: 10

```
from placeholder_name import ManualNameResolver
manual_resolver = ManualNameResolver('manual')

resolved_smiles = resolve_compounds_to_smiles(['NAH'], [manual_resolver])
```

ManualNameResolver can be initialized with a custom dictionary mapping chemical names to SMILES:
```
from placeholder_name import ManualNameResolver

custom_name_dict = {'Foobar': 'c1ccccc1'}
manual_resolver = ManualNameResolver('manual', provided_name_dict=custom_name_dict)

resolved_smiles = resolve_compounds_to_smiles(['Foobar'], [manual_resolver])
```

## PeptideNameResolver
This resolver converts shorthand peptide names (e.g. 'cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)') to an IUPAC-like name, then attempts to resolve to SMILES with py2opsin.

Default weight for 'weighted' SMILES selection: 3

```
from placeholder_name import PeptideNameResolver
peptide_shorthand_resolver = PeptideNameResolver('peptide')

resolved_smiles = resolve_compounds_to_smiles(['cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)'], [peptide_shorthand_resolver])
```

## StructuralFormulaNameResolver
This resolver converts simple structural chemical formulas (e.g. 'CH3CH2CH2COOH') to SMILES.

Default weight for 'weighted' SMILES selection: 2

```
from placeholder_name import StructuralFormulaNameResolver
structural_formula_resolver = StructuralFormulaNameResolver('structural_formula')

resolved_smiles = resolve_compounds_to_smiles(['CH3CH2CH2COOH'], [structural_formula_resolver])
```

## Custom Resolvers
This library also supports using custom resolvers. To use a custom resolver import the base class ChemicalNameResolver, and create a subclass with the format shown below. The name_to_smiles method is used to resolve the names to SMILES. In this example, the method resolves names using a simple lookup dictionary, but it can be also used to call an API, use other name-to-SMILES libraries, run an algorithm, etc. This method must return a tuple of dictionaries, where the first dictionary maps chemical names (strings) to SMILES (strings). The second dictionary returns information to the detailed_name_dict by mapping chemical names (strings) to some info message (strings).

```
from placeholder_name import resolve_compounds_to_smiles
from placeholder_name import ChemicalNameResolver

class MyCustomResolver(ChemicalNameResolver):
    """
    My custom resolver.
    """

    def __init__(self, resolver_name: str, resolver_weight: float = 1):
        super().__init__("example", resolver_name, resolver_weight)

    def name_to_smiles(
        self,
        compound_name_list: List[str]
    ) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Lookup chemical names from a dict.
        """
        lookup_dict = {
            'benzene': 'c1ccccc1'
        }

        resolved_names_dict = {}
        info_messages_dict = {}
        for compound_name in compound_name_list:
            resolved_smiles = lookup_dict.get(compound_name, '')
            resolved_names_dict[compound_name] = resolved_smiles
            if not resolved_smiles:
                info_messages_dict[compound_name] = 'Some info message.'

        return resolved_names_dict, info_messages_dict

my_custom_resolver = MyCustomResolver('example', resolver_weight=1)

resolved_smiles = resolve_compounds_to_smiles(['benzene', 'aspirin'], [my_custom_resolver], detailed_name_dict=True)
'{'benzene': {'SMILES': 'c1ccccc1', 'SMILES_source': ['example'], 'SMILES_dict': {'c1ccccc1': ['example']}, 'info_messages': {}}, 'aspirin': {'SMILES': '', 'SMILES_source': [], 'SMILES_dict': {}, 'info_messages': {'example': 'Some info message.'}}}'
```