"""
Parse SMILES strings to extract atom-level chemical features (donor, acceptor, hydrophobe, aromatic).
"""

def extract_atom_features(mol):
    """
    Extracts chemical features from a molecule to match pharmacophore sites.
    """
    features = []
    
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        
        if atom.GetIsAromatic():
            features.append((idx, 'Aromatic'))
            continue
            
        if atom.GetSymbol() in ('O', 'N'):
            if atom.GetExplicitValence() < atom.GetTotalValence() or atom.GetNumH() > 0:
                features.append((idx, 'Donor'))
            
            features.append((idx, 'Acceptor'))
            continue
            
        if atom.GetSymbol() == 'C' or atom.GetSymbol() in ('Cl', 'Br', 'I', 'F'):
            features.append((idx, 'Hydrophobe'))
            
    return features

def get_feature_map(mol):
    """
    Returns a mapping of feature types to lists of atom indices.
    """
    atom_features = extract_atom_features(mol)
    feature_map = {
        'Acceptor': [],
        'Donor': [],
        'Hydrophobe': [],
        'Aromatic': []
    }
    
    for idx, feat in atom_features:
        if feat in feature_map:
            feature_map[feat].append(idx)
            
    return feature_map
