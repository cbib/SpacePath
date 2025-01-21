#!/bin/bash
########################## Slurm options
#SBATCH --job-name=kpca-whole
#SBATCH --output=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/slurmout/kpca_%j.out
#SBATCH --workdir=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/user-pipeline/scripts/whole-ellis
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --ntasks-per-node=20
#SBATCH --time=30:00:00
#SBATCH --exclusive
##################################################
scontrol show job $SLURM_JOB_ID

eval "$(conda shell.bash hook)"  # this solves commandnotfounderror when conda activate

module load Metabo/stable

cd "/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/user-pipeline/scripts/whole-ellis"

python3 2-adata-pathways.py




