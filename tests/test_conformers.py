import unittest
from rdkit import Chem
from src.conformers import generate_conformers, get_molecular_weight, count_atoms_and_bonds

class TestConformerGeneration(unittest.TestCase):
    def setUp(self):
        self.test_smiles = {
            "ethanol": "CCO", # Ethanol C2H6O: ~46.07
            "benzene": "C1=CC=CC=C1",
            "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "invalid": "NOT_A_SMILES"
        }

    def test_generate_conformers_success(self):
        """Test that valid SMILES generate the requested number of conformers."""
        smiles = self.test_smiles["ethanol"]
        max_confs = 5
        conformers = generate_conformers(smiles, max_confs=max_confs)
        
        self.assertIsInstance(conformers, list)
        self.assertLessEqual(len(conformers), max_confs)
        if conformers:
            self.assertIsInstance(conformers[0], Chem.Mol)

    def test_generate_conformers_invalid_smiles(self):
        """Test that invalid SMILES raise a ValueError."""
        with self.assertRaises(ValueError):
            generate_conformers(self.test_smiles["invalid"])

    def test_get_molecular_weight(self):
        mw = get_molecular_weight(self.test_smiles["ethanol"])
        self.assertAlmostEqual(mw, 46.07, places=1)

    def test_get_molecular_weight_invalid(self):        
        # Invalid SMILES should return None
        self.assertIsNone(get_molecular_weight(self.test_smiles["invalid"]))

    def test_count_atoms_and_bonds(self):
        """Test atom and bond counting."""
        smiles = self.test_smiles["benzene"]
        info = count_atoms_and_bonds(smiles)
        
        self.assertIsNotNone(info)
        self.assertEqual(info['atoms'].get('C'), 6)
        self.assertEqual(info['num_bonds'], 6)

    def test_conformer_geometry_have_3D_coordinates(self):
        smiles = self.test_smiles["benzene"]
        conformers = generate_conformers(smiles, max_confs=1)
        
        if conformers:
            mol = conformers[0]
            conf = mol.GetConformer()
            # Check if coordinates are not all zero
            pos = conf.GetAtomPosition(0)
            self.assertTrue(any([pos.x != 0, pos.y != 0, pos.z != 0]))

if __name__ == "__main__":
    unittest.main()
