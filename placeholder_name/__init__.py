"""name_to_smiles initialization."""
from .main import (
    CIRpyNameResolver,
    ChemicalNameResolver,
    ManualNameResolver,
    OpsinNameResolver,
    PeptideNameResolver,
    PubChemNameResolver,
    StructuralFormulaNameResolver,
    resolve_compounds_to_smiles,
)

__all__ = [
  'resolve_compounds_to_smiles',
  'ChemicalNameResolver',
  'ManualNameResolver',
  'OpsinNameResolver',
  'PubChemNameResolver',
  'PeptideNameResolver',
  'StructuralFormulaNameResolver',
  'CIRpyNameResolver',
]

__version__ = "0.0.1"
