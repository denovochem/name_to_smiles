Resolve chemical names to SMILES by passing a string or a list of strings:
```pycon
from placeholder_name import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(['aspirin'])

"{'aspirin': 'CC(=O)Oc1ccccc1C(=O)O'}"
```

See detailed information including which resolver returned which SMILES with detailed_name_dict=True:
```pycon
from placeholder_name import resolve_compounds_to_smiles

resolved_smiles = resolve_compounds_to_smiles(
    ['2-acetyloxybenzoic acid'], 
    detailed_name_dict=True
)

"{'2-acetyloxybenzoic acid': {
    'SMILES': 'CC(=O)Oc1ccccc1C(=O)O',
    'SMILES_source': ['pubchem_default', 'opsin_default'],
    'SMILES_dict': {
        'CC(=O)Oc1ccccc1C(=O)O': ['pubchem_default', 'opsin_default']
    },
    'info_messages': {}
}}"
```

## Advanced usage
Many aspects of the name-to-SMILES resolution process can be customized, including the resolvers that are used, the configuration of those resolvers, and the strategy used to pick the best SMILES.

In this example, we resolve chemical names with OPSIN, PubChem, and CIRPy, and use a custom consensus weighting approach to pick the best SMILES:
```pycon
from placeholder_name import resolve_compounds_to_smiles
from placeholder_name import (
    OpsinNameResolver, 
    PubChemNameResolver, 
    CIRpyNameResolver
)

opsin_resolver = OpsinNameResolver('opsin', resolver_weight=4)
pubchem_resolver =  PubChemNameResolver('pubchem', resolver_weight=3)
cirpy_resolver = CIRpyNameResolver('cirpy', resolver_weight=2)

resolved_smiles = resolve_compounds_to_smiles(
    ['2-acetyloxybenzoic acid'],
    [opsin_resolver, pubchem_resolver, cirpy_resolver],
    smiles_selection_mode='weighted',
    detailed_name_dict=True
)

"{'2-acetyloxybenzoic acid': {
    'SMILES': 'CC(=O)Oc1ccccc1C(=O)O',
    'SMILES_source': ['opsin', 'pubchem', 'cirpy'],
    'SMILES_dict': {
        'CC(=O)Oc1ccccc1C(=O)O': ['opsin', 'pubchem', 'cirpy']
    },
    'info_messages': {}
}}"
```

More information about advanced usage can be found in [Resolvers](resolvers.md), [Smiles Selection](smiles_selection.md), and [Name Correction/Editing](name_manipulation.md)