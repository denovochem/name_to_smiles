from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

enumerator = rdMolStandardize.TautomerEnumerator()

def canonicalize_smiles(
    smiles: str,
    isomeric: bool = True,
    remove_mapping: bool = True,
    canonicalize_tautomer: bool =True,
) -> str:
    """
    Converts SMILES strings to their canonical form using RDKit.

    Takes a SMILES string (potentially containing multiple fragments separated by periods),
    splits it into fragments, sorts them, and converts each to its canonical form. Handles
    atom mapping and isomeric SMILES options.

    Args:
        smiles (str): The input SMILES string to canonicalize
        isomeric (bool): Whether to retain isomeric information. Defaults to True
        canonicalize_tautomer (bool): Whether to use the canonical tautomer. Defaults to True
        remove_mapping (bool): Whether to remove atom mapping numbers. Defaults to True

    Returns:
        str: The canonicalized SMILES string. If conversion fails, returns the input string
            unchanged.
    """
    try:
        x = smiles.split('.')
        x = sorted(x)
        frags = []
        for i in x:
            m = Chem.MolFromSmiles(i)
            if canonicalize_tautomer:
                m = enumerator.Canonicalize(m)
            if remove_mapping:
                [a.SetAtomMapNum(0) for a in m.GetAtoms()]
            canonical_smiles_string = str(Chem.MolToSmiles(m, canonical=True, isomericSmiles=isomeric))
            frags.append(canonical_smiles_string)
        canonical_smiles_string = '.'.join(i for i in sorted(frags))
        return(canonical_smiles_string)
    except:
        return smiles

