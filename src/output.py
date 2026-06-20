from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def write_sdf(mols, scores, output_path):
    """
    Writes a list of RDKit molecule objects to an SDF file.
    Includes the alignment score and basic molecular properties in the SDF tags.
    """
    writer = Chem.SDWriter(output_path)
    
    for mol, score in zip(mols, scores):
        if mol is None:
            continue
            
        # Calculate properties for SDF tags
        mw = Descriptors.MolWt(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        mol.SetProp('MolecularWeight', f"{mw:.2f}")
        mol.SetProp('H-Bond Donors', str(h_donors))
        mol.SetProp('H-Bond Acceptors', str(h_acceptors))
        mol.SetProp('Score', f"{score:.6f}")
        
        writer.write(mol)
        
    writer.close()
