## Quick Start

Install placeholder_name with pip directly from the github repo:

```shell
pip install git+https://github.com/denovochem/name_to_smiles.git

```

## Basic usage
Resolve chemical names to SMILES by passing a string or a list of strings:
```pycon
from placeholder_name import resolve_compounds_to_smiles
resolved_smiles = resolve_compounds_to_smiles(['aspirin'])
'{'aspirin': 'CC(=O)Oc1ccccc1C(=O)O'}'
```
See detailed information including which resolver returned which SMILES with detailed_name_dict=True:
```pycon
from placeholder_name import resolve_compounds_to_smiles
resolved_smiles = resolve_compounds_to_smiles(['2-acetyloxybenzoic acid'], detailed_name_dict=True)
'{'2-acetyloxybenzoic acid': {'SMILES': 'CC(=O)Oc1ccccc1C(=O)O', 'SMILES_source': ['pubchem_default', 'opsin_default'], 'SMILES_dict': {'CC(=O)Oc1ccccc1C(=O)O': ['pubchem_default', 'opsin_default']}, 'info_messages': {}}}'
```