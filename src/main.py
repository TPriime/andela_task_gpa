import json
import os
import argparse
from conformers import generate_conformers
from chemistry import get_feature_map
from pose_selection import select_best_pose
from output import write_sdf

def process_targets(targets_path, output_path):
    """
    Main pipeline to process targets, align to pharmacophores, and output SDF.
    """
    if not os.path.exists(targets_path):
        print(f"Error: Targets file not found at {targets_path}")
        return

    with open(targets_path, 'r') as f:
        targets_data = json.load(f)

    best_poses = []

    for target_id, target in targets_data.items():
        smiles = target.get('smiles')
        interaction_sites = target.get('interaction_sites', [])
        excluded_volumes = target.get('excluded_volumes', [])
        
        if not smiles:
            continue

        print(f"Processing {target_id}: {smiles}")
        
        try:
            # 1. Generate conformers
            poses = generate_conformers(smiles)
            
            if not poses:
                print(f"  No valid conformers generated for {smiles}")
                continue
                
            # 2. Extract chemical features (using the first pose as reference for topology)
            feature_map = get_feature_map(poses[0])
            
            # 3. Select best pose based on alignment and clashes
            best_pose, score = select_best_pose(poses, interaction_sites, feature_map, excluded_volumes)
            
            if best_pose:
                print(f"  Best pose found with score: {score:.4f}")
                best_poses.append(best_pose)
            else:
                print(f"  No valid pose surviving clash detection for {smiles}")
                
        except Exception as e:
            print(f"  Error processing {smiles}: {e}")

    # 4. Write results to SDF
    if best_poses:
        write_sdf(best_poses, output_path)
        print(f"\nSuccessfully wrote {len(best_poses)} poses to {output_path}")
    else:
        print("\nNo valid poses found for any target.")

def main():
    parser = argparse.ArgumentParser(description="Geometric Pharmacophore Alignment Tool")
    parser.add_argument("-t", "--targets", default="./data/targets.json", help="Path to targets.json")
    parser.add_argument("-o", "--output", default="./results/docked_poses.sdf", help="Path to output SDF file")
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    process_targets(args.targets, args.output)

if __name__ == "__main__":
    main()
