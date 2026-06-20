"""
Parse SMILES strings to extract atom-level chemical features (donor, acceptor, hydrophobe, aromatic).
"""
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
import os

_FDEF_PATH = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
_FACTORY = ChemicalFeatures.BuildFeatureFactory(_FDEF_PATH)

_RDKIT_TO_TASK = {
    "Donor":          "Donor",
    "Acceptor":       "Acceptor",
    "Hydrophobe":     "Hydrophobe",
    "LumpedHydrophobe": "Hydrophobe",
    "Aromatic":       "Aromatic",
}

def extract_atom_features(mol):
    """
    Extracts chemical features from a molecule to match pharmacophore sites.
    """
    result = {f: set() for f in ("Donor", "Acceptor", "Hydrophobe", "Aromatic")}
    for feat in _FACTORY.GetFeaturesForMol(mol):
        task_fam = _RDKIT_TO_TASK.get(feat.GetFamily())
        if task_fam is None:
            continue
        for idx in feat.GetAtomIds():
            if mol.GetAtomWithIdx(idx).GetAtomicNum() != 1:
                result[task_fam].add(idx)
    
    # Convert sets to sorted lists for consistency
    return {k: sorted(list(v)) for k, v in result.items()}

def get_feature_map(mol):
    """
    Returns a mapping of feature types to lists of atom indices.
    """
    return extract_atom_features(mol)
