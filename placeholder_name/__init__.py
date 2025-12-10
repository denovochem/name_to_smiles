"""name_to_smiles initialization."""
from .main import resolve_compounds_to_smiles, ManualNameResolver, OpsinNameResolver, PubChemNameResolver, PeptideNameResolver, StructuralFormulaNameResolver, CIRPyNameResolver

__all__ = [
  'resolve_compounds_to_smiles',
  'ManualNameResolver',
  'OpsinNameResolver',
  'PubChemNameResolver',
  'PeptideNameResolver',
  'StructuralFormulaNameResolver',
  'CIRPyNameResolver',
]

__version__ = "0.0.1"
