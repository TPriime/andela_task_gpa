from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def generate_conformers(smiles: str, max_confs: int = 10) -> list[Chem.Mol]:
    """Generate 3D conformers from SMILES using ETKDGv3."""
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES string")

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    
    conformers = []
    for _ in range(min(max_confs * 2, 50)):
        conf_id = AllChem.EmbedMolecule(mol_h, params=params)
        
        if conf_id == -1:
            continue
        
        uff_result = AllChem.UFFOptimizeMolecule(mol_h, maxIters=200)
        
        if uff_result == 0 and len(conformers) < max_confs:
            conformers.append(Chem.RemoveHs(mol_h))
    
    return conformers


def get_molecular_weight(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    try:
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        return round(float(mw), 2)
    except Exception as e:
        print(f"Warning: Could not calculate molecular weight for '{smiles}': {e}")
        return None


def count_atoms_and_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    atom_counts = {}
    
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        
        if symbol != 'H':
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    return {
        'atoms': atom_counts,
        'total_atoms_including_h': mol.GetNumAtoms(),
        'num_bonds': mol.GetNumBonds()
    }
