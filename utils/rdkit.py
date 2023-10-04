from rdkit import Chem
from rdkit.Chem import Descriptors


def molecular_weight(smiles: str) -> float:
    """
    Calculate exact molecular weight for the given SMILES string

    Args:
        smiles: SMILES string for which to calculate molecular weight

    Returns:
         float: exact molecular weight
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        molwt = Descriptors.ExactMolWt(mol)

        return molwt
    except Exception:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return 9999.0
        mol.UpdatePropertyCache(strict=False)
        molwt = Descriptors.ExactMolWt(mol)

        return molwt
