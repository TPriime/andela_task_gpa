import numpy as np

def compute_site_score(ligand_feature_pos, target_site_pos):
    dist = np.linalg.norm(ligand_feature_pos - target_site_pos)
    return np.exp(-(dist / 1.25)**2)

def align_to_pharmacophore(mol, interaction_sites, feature_map):
    """Align molecule conformer to pharmacophore sites and compute total score."""
    conf = mol.GetConformer()
    total_score = 0.0
    
    for site in interaction_sites:
        site_type = site.get('family')
        site_coords = np.array([site.get('x'), site.get('y'), site.get('z')])
        
        matching_indices = feature_map.get(site_type, [])
        if not matching_indices:
            continue
            
        max_site_score = 0.0
        for idx in matching_indices:
            atom_pos = np.array(conf.GetAtomPosition(idx))
            current_score = compute_site_score(atom_pos, site_coords)
            if current_score > max_site_score:
                max_site_score = current_score
        
        total_score += max_site_score

    return total_score

def check_clashes(mol, excluded_volumes):
    """Check for steric clashes between the molecule and exclusion volumes."""
    conf = mol.GetConformer()
    tolerance = 0.1
    
    for volume in excluded_volumes:
        vol_coords = np.array([volume.get('x'), volume.get('y'), volume.get('z')])
        vol_radius = volume.get('radius', 1.2)
        
        for atom_idx in range(mol.GetNumAtoms()):
            atom_pos = np.array(conf.GetAtomPosition(atom_idx))
            dist = np.linalg.norm(atom_pos - vol_coords)
            
            if dist < (vol_radius + tolerance):
                return False
                
    return True
