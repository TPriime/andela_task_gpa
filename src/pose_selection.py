import numpy as np
from scipy.spatial.transform import Rotation
from scoring import align_to_pharmacophore, check_clashes
from rdkit import Chem

def select_best_pose(poses, interaction_sites, feature_map, excluded_volumes):
    """
    Exhaustively aligns poses to pharmacophore sites, filters for steric clashes, 
    and selects the one with the highest pharmacophore score.
    """
    best_pose = None
    max_score = -1.0

    for pose in poses:
        mol_h = Chem.AddHs(pose)
        coords_heavy = np.array(pose.GetConformer().GetPositions())
        
        best_pose_for_conf = None
        best_score_for_conf = -1.0
        
        site_data = []
        for site in interaction_sites:
            site_data.append({
                'family': site.get('family'),
                'coords': np.array([site.get('x'), site.get('y'), site.get('z')]),
                'weight': site.get('weight', 1.0)
            })

        for site in site_data:
            family = site['family']
            site_coords = site['coords']
            matching_indices = feature_map.get(family, [])
            
            for atom_idx in matching_indices:
                atom_pos = coords_heavy[atom_idx]
                
                for _ in range(50):
                    R = Rotation.random().as_matrix()
                    centered_coords = coords_heavy - atom_pos
                    rotated_coords = (R @ centered_coords.T).T + atom_pos
                    t = site_coords - rotated_coords[atom_idx]
                    final_coords_heavy = rotated_coords + t
                    
                    all_pos = np.array(mol_h.GetConformer().GetPositions())
                    all_pos_centered = all_pos - all_pos[atom_idx]
                    all_pos_rotated = (R @ all_pos_centered.T).T
                    all_pos_final = all_pos_rotated + site_coords
                    
                    if check_clashes(all_pos_final, excluded_volumes):
                        heavy_mask = np.array([mol_h.GetAtomWithIdx(i).GetAtomicNum() != 1 for i in range(mol_h.GetNumAtoms())])
                        heavy_pos_final = all_pos_final[heavy_mask]
                        
                        score = align_to_pharmacophore(heavy_pos_final, feature_map, interaction_sites)
                        if score > best_score_for_conf:
                            best_score_for_conf = score
                            new_conf = pose.GetConformer()
                            new_conf.SetPositions(heavy_pos_final)
                            best_pose_for_conf = pose

        if best_pose_for_conf and best_score_for_conf > max_score:
            max_score = best_score_for_conf
            best_pose = best_pose_for_conf
            
    return best_pose, max_score
