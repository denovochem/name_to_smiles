# placeholder_name
[![PyPI Version](https://img.shields.io/pypi/v/PubChemPy?logo=python&logoColor=%23ffffff)](https://pypi.python.org/pypi/PubChemPy)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://gitHub.com/denovochem/name_to_smiles/graphs/commit-activity)
[![License](https://img.shields.io/pypi/l/PubChemPy)](https://github.com/denovochem/name_to_smiles/blob/main/LICENSE)
[![Tests](https://img.shields.io/github/actions/workflow/status/denovochem/name_to_smiles/tests.yml?logo=github&logoColor=%23ffffff&label=tests)](https://github.com/denovochem/name_to_smiles/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/readthedocs/pubchempy?logo=readthedocs&logoColor=%23ffffff)](https://denovochem.github.io/name_to_smiles/)

This library is used for performant, comprehensive, and customizable name-to-SMILES conversions. 

This library can use the following existing name-to-SMILES resolvers:
- [py2opsin](https://github.com/csnbritt/py2opsin)
- [PubChemPy](https://github.com/csnbritt/PubChemPy)
- [CIRpy](https://github.com/mcs07/CIRpy)
- [ChemSpiPy](https://github.com/mcs07/ChemSpiPy)


This library also implements the following new resolvers:
- Manually curated dataset of common names not correctly resolved by other resolvers (e.g. NaH)
- Structural formula resolver (e.g. 'CH3CH2CH2COOH')
- Peptide shorthand resolver (e.g. 'cyclo(Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)')


The following string editing/manipulation strategies may be applied to compounds to assist with name-to-SMILES resolution:
- Splitting compounds on common delimiters (useful for mixtures of compounds, e.g. 'BH3•THF')


When resolvers disagree on the SMILES for a given compound, a variety of SMILES selection methods can be employed to determine the "best" SMILES for a given compound name. See the documentation for more details.

## Installation

Install placeholder_name with pip directly from this repo:

```shell
pip install git+https://github.com/denovochem/name_to_smiles.git
```

## Basic usage
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
from placeholder_name import OpsinNameResolver, PubChemNameResolver, CIRpyNameResolver

opsin_resolver = OpsinNameResolver(
    resolver_name='opsin', 
    resolver_weight=4
)
pubchem_resolver =  PubChemNameResolver(
    resolver_name='pubchem', 
    resolver_weight=3
)
cirpy_resolver = CIRpyNameResolver(
    resolver_name='cirpy', 
    resolver_weight=2
)

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

See documentation for more details. 

## Documentation
Full documentation is availible [here](https://denovochem.github.io/name_to_smiles/)

## Contributing

- Feature ideas and bug reports are welcome on the Issue Tracker.
- Fork the [source code](https://github.com/denovochem/name_to_smiles) on GitHub, make changes and file a pull request.

## License

placeholder_name is licensed under the [MIT license](https://github.com/denovochem/name_to_smiles/blob/main/LICENSE).