# Ge‎o‎met‎r‎ic P‎h‎ar‎m‎ac‎o‎pho‎r‎e A‎li‎gnm‎e‎nt (G‎P‎A‎)

This project implements a cross-docking solution for placing small molecules into protein pockets defined by pharmacophore interaction sites and exclusion spheres. The pocket is represented without an explicit protein structure, relying instead on chemical features and steric boundaries.

## Project Description

The goal is to take a set of target molecules (provided in SMILES format) and find the best 3D pose for each within a predefined pharmacophore environment. Each target environment consists of:

- **Interaction Sites**: Points with specific chemical families (Donor, Acceptor, Hydrophobe, Aromatic) and weights.
- **Exclusion Volumes**: Spheres that define steric boundaries where no ligand atom should reside.

The alignment is performed by generating 3D conformers for each molecule and scoring them based on how well their atoms match the pharmacophore interaction sites using a Gaussian-like scoring function:

$$score = \sum_{i} w_i \cdot \exp\left(-\left(\frac{d_i}{1.25}\right)^2\right)$$

where $d_i$ is the minimum distance from site $i$ to the nearest matching ligand atom. Poses that violate the exclusion volumes (within a 0.1 Å tolerance of the 1.2 Å radius) are rejected.

## Installation and Setup

This project is containerized using Docker for ease of use and environment consistency.

### Prerequisites

- Docker installed on your system.
- Bash shell (for running the provided script).

### Building the Docker Image

To build the project image, run the following command from the root directory:

```bash
docker build -t gpa-image .
```

## Usage

The project provides a bash script `gpa` to run the docking pipeline within the Docker container.

### Basic Execution

To run the pipeline with the default configuration (targets from `./data/targets.json` and output to `./results/docked_poses.sdf`):

```bash
./gpa
```

### Advanced Usage

You can specify custom paths for the input targets file and the output SDF file using flags:

```bash
./gpa -t "./path/to/your/targets.json" -o "./path/to/your/output.sdf"
```

**Options:**

- `-t, --targets <path>`: Path to the `targets.json` file (default: `./data/targets.json`).
- `-o, --output <path>`: Path where the resulting `.sdf` file will be saved (default: `./results/docked_poses.sdf`).
- `-h, --help`: Show this help message.

### Example

```bash
./gpa -t "./data/custom_targets.json" -o "./results/my_docking_results.sdf"
```

### Running directly with docker

For the most minimal setup, you can run the entire pipeline directly using the built image without relying on the local `./gpa` script. This command mounts the current directory (`$(pwd)`) into the container and executes the main Python script:

```bash
docker run --rm -v "$(pwd)":/app -w /app gpa-image
```

Full example:

```
docker run --rm \
  -v "$(pwd)":/app \
  -w /app gpa-image \
  -u $(id -u):$(id -g)
  -t "./path/to/your/targets.json" -o "./path/to/your/output.sdf"
```

---

### Development testing via docker

```
docker run --rm -v "$(pwd)":/app -w /app --entrypoint python3 gpa-image -m unittest tests/*.py;
```

## Implementation Details

This section provides a detailed explanation of the algorithm used to align small molecules to pharmacophore interaction sites and exclusion volumes.

### 1. Overview

The core objective is to find a 3D pose for each molecule that maximizes its alignment with a set of predefined pharmacophore features while strictly avoiding steric clashes with exclusion spheres.

The pipeline consists of several modules:

- `src/conformers.py`: Generates 3D conformers using RDKit.
- `src/chemistry.py`: Extracts chemical features (Donor, Acceptor, etc.) from molecules.
- `src/scoring.py`: Implements the pharmacophore scoring and clash detection logic.
- `src/pose_selection.py`: Filters valid poses and selects the best one per target.
- `src/output.py`: Generates the final SDF output file.
- `src/main.py`: Orchestrates the entire pipeline.

#### Project Structure

- `data/`: Contains input data (`targets.json`).
- `results/`: Output directory for generated `.sdf` files.
- `src/`: Python source code.
- `tests/`: Unit tests for the various modules.
- `Dockerfile`: Docker configuration for the environment.
- `./gpa`: Bash wrapper for running the containerized application.

### 2. Step-by-Step Process

#### A. Conformer Generation

For each target molecule (provided as a SMILES string):

1.  **SMILES Parsing**: The SMILES string is parsed into an RDKit molecule object.
2.  **3D Embedding**: Multiple 3D conformers are generated using RDKit's `AllChem.EmbedMultipleConfs`.
3.  **Geometry Optimization**: Each conformer undergoes energy minimization using the Universal Force Field (UFF) or Merck Molecular Force Field (MMFF) to ensure realistic bond lengths and angles.

#### B. Chemical Feature Extraction

For each generated conformer, we identify atoms that correspond to the pharmacophore feature families:

- **Donor**: Atoms capable of donating a hydrogen bond.
- **Acceptor**: Atoms capable of accepting a hydrogen bond.
  _Note: The implementation prioritizes extraction based on prevalence (Acceptor > Aromatic > Hydrophobe >> Donor) as identified during data analysis._

#### C. Pharmacophore Scoring

Each conformer is evaluated against the interaction sites defined in `targets.json`. For each site $i$:

1.  Identify all atoms in the molecule that belong to the same feature family as site $i$.
2.  Calculate the distance $d_{i,j}$ from site $i$ to each matching atom $j$.
3.  Find the minimum distance: $d_i = \min(d_{i,j})$.
4.  Compute the contribution of site $i$ using the Gaussian scoring function:
    $$score_i = w_i \cdot \exp\left(-\left(\frac{d_i}{1.25}\right)^2\right)$$
    where $w_i$ is the weight of the site.
5.  The total score for the conformer is the sum of all site contributions:
    $$Score_{total} = \sum_{i} score_i$$

#### D. Steric Clash Detection

A pose is considered invalid if any atom in the molecule enters an exclusion volume.

- **Exclusion Volume**: A sphere defined by a center $(x, y, z)$ and a radius of $1.2\text{ \AA}$.
- **Tolerance**: An additional tolerance of $0.1\text{ \AA}$ is applied.
- **Condition**: If the distance from any atom to an exclusion center is less than $1.3\text{ \AA}$ ($1.2 + 0.1$), the pose is rejected.

#### E. Pose Selection

After scoring and filtering:

1.  All conformers that passed the clash detection are considered "valid".
2.  From the valid set, the conformer with the highest total score is selected as the best pose for that target.
3.  If no conformers pass the clash detection, the target is skipped or an error is reported.

### 3. Complexity and Performance

- **Time Complexity**: Primarily driven by the number of conformers generated and the number of atoms/sites to check. The scoring is $O(N_{atoms} \cdot N_{sites})$.
- **Space Complexity**: $O(N_{conformers} \cdot N_{atoms})$ to store all generated poses before selection.
