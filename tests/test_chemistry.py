import unittest
from rdkit import Chem
from src.chemistry import extract_atom_features, get_feature_map

class TestChemistry(unittest.TestCase):
    def test_extract_features_simple(self):
        # Ethanol: CC(O)
        # C(0): Hydrophobe, C(1): Hydrophobe, O(2): Acceptor/Donor
        mol = Chem.MolFromSmiles("CCO")
        features = extract_atom_features(mol)
        
        # Check if we have the expected features
        feat_types = [f[1] for f in features]
        self.assertIn('Hydrophobe', feat_types)
        self.assertIn('Acceptor', feat_types)
        self.assertIn('Donor', feat_types)

    def test_extract_features_aromatic(self):
        # Benzene: c1ccccc1
        mol = Chem.MolFromSmiles("c1ccccc1")
        features = extract_atom_features(mol)
        
        for idx, feat in features:
            self.assertEqual(feat, 'Aromatic')

    def test_feature_map(self):
        mol = Chem.MolFromSmiles("CCO")
        fmap = get_feature_map(mol)
        
        self.assertTrue(len(fmap['Hydrophobe']) >= 2)
        self.assertTrue(len(fmap['Acceptor']) >= 1)
        self.assertTrue(len(fmap['Donor']) >= 1)

if __name__ == "__main__":
    unittest.main()
