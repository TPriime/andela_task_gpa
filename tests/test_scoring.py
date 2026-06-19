import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from src.scoring import align_to_pharmacophore, check_clashes, compute_site_score

def test_compute_site_score():
    # Test exact match
    assert compute_site_score(np.array([0,0,0]), np.array([0,0,0])) == 1.0
    # Test distance 1.25
    assert np.isclose(compute_site_score(np.array([1.25,0,0]), np.array([0,0,0])), np.exp(-1))

def test_alignment_score():
    mol = Chem.MolFromSmiles('CC')
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    
    feature_map = {'Hydrophobe': [0, 1]}
    interaction_sites = [{'type': 'Hydrophobe', 'coords': [0, 0, 0]}]
    
    # Move atom 0 to origin for perfect match
    conf = mol.GetConformer()
    conf.SetAtomPosition(0, (0, 0, 0))
    
    score = align_to_pharmacophore(mol, interaction_sites, feature_map)
    assert score == 1.0

def test_clash_detection():
    mol = Chem.MolFromSmiles('C')
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    
    # Atom at origin
    conf = mol.GetConformer()
    conf.SetAtomPosition(0, (0, 0, 0))
    
    # Volume at origin, radius 1.2. Distance 0 < 1.3 -> Clash
    excluded_volumes = [{'coords': [0, 0, 0], 'radius': 1.2}]
    assert check_clashes(mol, excluded_volumes) is False
    
    # Volume far away -> No clash
    excluded_volumes = [{'coords': [10, 0, 0], 'radius': 1.2}]
    assert check_clashes(mol, excluded_volumes) is True

if __name__ == "__main__":
    test_compute_site_score()
    test_alignment_score()
    test_clash_detection()
    print("All scoring tests passed!")
