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

## Implementation Details

The pipeline consists of several modules:

- `src/conformers.py`: Generates 3D conformers using RDKit.
- `src/chemistry.py`: Extracts chemical features (Donor, Acceptor, etc.) from molecules.
- `src/scoring.py`: Implements the pharmacophore scoring and clash detection logic.
- `src/pose_selection.py`: Filters valid poses and selects the best one per target.
- `src/output.py`: Generates the final SDF output file.
- `src/main.py`: Orchestrates the entire pipeline.

## Project Structure

- `data/`: Contains input data (`targets.json`).
- `results/`: Output directory for generated `.sdf` files.
- `src/`: Python source code.
- `tests/`: Unit tests for the various modules.
- `Dockerfile`: Docker configuration for the environment.
- `./gpa`: Bash wrapper for running the containerized application.
