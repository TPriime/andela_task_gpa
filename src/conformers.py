from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def generate_conformers(smiles: str, max_confs: int = 200) -> list[Chem.Mol]:
    """Generate 3D conformers from SMILES using ETKDGv3 and MMFF optimization."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(f"Invalid SMILES string")

    # Attach Hydrogen atoms
    mol_h = Chem.AddHs(mol)
    
    # Use ETKDGv3 (the best) for initial embeddings and fallback on failure 
    # to others: with random coordinates, ETKDGv2, basic
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0
    
    num_confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=50, params=params)
    
    if num_confs == 0:
        params.randomCoords = True
        num_confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=max_confs, params=params)

    if num_confs == 0:
        params_v2 = AllChem.ETKDGv2()
        params_v2.randomCoords = True
        num_confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=max_confs, params=params_v2)

    # If still no conformers, try a very basic embedding
    if num_confs == 0:
        AllChem.EmbedMolecule(mol_h, randomCoord=True)
        num_confs = mol_h.GetNumConformers()
    
    AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=0)

   # Extract optimized conformers and remove Hydrogen atoms for final output 
    conformers = []
    for i in range(mol_h.GetNumConformers()):
        conf_mol = Chem.Mol(mol_h)
        conf = mol_h.GetConformer(i)
        new_conf = Chem.Conformer(conf)
        conf_mol.AddConformer(new_conf)
        conformers.append(Chem.RemoveHs(conf_mol))
            
    return conformers[:max_confs]

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
