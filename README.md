# placeholder_name
[![PyPI Version](https://img.shields.io/pypi/v/PubChemPy?logo=python&logoColor=%23ffffff)](https://pypi.python.org/pypi/PubChemPy)
[![License](https://img.shields.io/pypi/l/PubChemPy)](https://github.com/mcs07/PubChemPy/blob/main/LICENSE)
[![Tests](https://img.shields.io/github/actions/workflow/status/mcs07/pubchempy/test.yml?logo=github&logoColor=%23ffffff&label=tests)](https://github.com/mcs07/PubChemPy/actions/workflows/test.yml)
[![Docs](https://img.shields.io/readthedocs/pubchempy?logo=readthedocs&logoColor=%23ffffff)](https://docs.pubchempy.org)

This library is used for performant and customizable name-to-SMILES conversions using a variety of existing and custom libraries. 

This library uses the following existing resolvers:
- OPSIN
- PubChem
- CIRPy

This library implements the following new resolvers:
- Manual database of common names not resolved by other services
- Structural formula resolver (e.g. CH3CH2CH2COOH)
- Peptide shorthand resolver (e.g. Asp-Arg-Val-Tyr-Ile-His-Pro-Phe)
- Mixtures (e.g. H20.THF)

## Installation

Install placeholder_name with pip:

```shell
pip install pubchempy
```

## Basic usage
```pycon
from placeholder_name import resolve_compounds_to_smiles
resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'])
''
```

## Advanced usage
```pycon
from placeholder_name import resolve_compounds_to_smiles
from placeholder_name import OpsinNameResolver, PubChemNameResolver, CIRPyNameResolver

opsin_resolver = OpsinNameResolver('opsin')
pubchem_resolver =  PubChemNameResolver('pubchem')
cirpy_resolver = CIRPyNameResolver('cirpy')

resolved_smiles = resolve_compounds_to_smiles(['MeOH.benzene'], [opsin_resolver, pubchem_resolver, cirpy_resolver])
''
```

## Documentation
Full documentation is availible at 

## Contributing

- Feature ideas and bug reports are welcome on the Issue Tracker.
- Fork the [source code](https://github.com/denovochem/name_to_smiles) on GitHub, make changes and file a pull request.

## License

PubChemPy is licensed under the [MIT license](https://github.com/denovochem/name_to_smiles/blob/main/LICENSE).
