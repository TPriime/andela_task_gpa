import numpy as np
from scoring import align_to_pharmacophore, check_clashes

def select_best_pose(poses, interaction_sites, feature_map, excluded_volumes):
    """
    Filters poses for steric clashes and selects the one with the highest pharmacophore score.
    """
    best_pose = None
    max_score = -1.0

    for pose in poses:
        # Steric clash check
        if not check_clashes(pose, excluded_volumes):
            continue
            
        # Alignment scoring
        score = align_to_pharmacophore(pose, interaction_sites, feature_map)
        
        if score > max_score:
            max_score = score
            best_pose = pose
            
    return best_pose, max_score
