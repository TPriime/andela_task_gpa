import numpy as np

def compute_site_score(ligand_feature_pos, target_site_pos, weight=1.0):
    dist = np.linalg.norm(ligand_feature_pos - target_site_pos)
    return weight * np.exp(-(dist / 1.25)**2)

def align_to_pharmacophore(coords_heavy, feature_map, interaction_sites):
    """
    score = Σ_i  w_i * exp(-(d_i / 1.25)^2)
    d_i = min distance from site i to nearest matching heavy atom.
    """
    total_score = 0.0
    for site in interaction_sites:
        site_type = site.get('family')
        site_coords = np.array([site.get('x'), site.get('y'), site.get('z')])
        weight = site.get('weight', 1.0)
        
        matching_indices = feature_map.get(site_type, [])
        if not matching_indices:
            continue
            
        # Find minimum distance to any atom of the correct feature type
        diffs = coords_heavy[matching_indices] - site_coords
        d_min = np.sqrt((diffs**2).sum(axis=1)).min()
        
        total_score += weight * np.exp(-(d_min / 1.25)**2)

    return total_score

def check_clashes(coords_all, excluded_volumes):
    """
    Check for steric clashes between the molecule and exclusion volumes.
    Spec: No ligand atom (including H) may come within 1.2 Å of any exclusion center (with 0.1 Å tolerance).
    Clash threshold = 1.2 - 0.1 = 1.1 Å.
    """
    clash_threshold = 1.1
    
    for volume in excluded_volumes:
        vol_coords = np.array([volume.get('x'), volume.get('y'), volume.get('z')])
        
        diffs = coords_all - vol_coords
        dists = np.sqrt((diffs**2).sum(axis=1))
        
        if dists.min() < clash_threshold:
            return False
                
    return True
