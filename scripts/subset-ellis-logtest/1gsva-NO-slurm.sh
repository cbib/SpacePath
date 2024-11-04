#!/bin/bash
########################## Slurm options
#SBATCH --job-name=brain-ellis_gsva
#SBATCH --output=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/slurmout/ellisgsva_%j.out
#SBATCH --workdir=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/
#SBATCH --mail-user=juana7@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --ntasks-per-node=20
#SBATCH --time=29:00:00
#SBATCH --exclusive
##################################################
scontrol show job $SLURM_JOB_ID

eval "$(conda shell.bash hook)"  # this solves commandnotfounderror when conda activate

conda activate r_gsva

cd /mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/

DATAPATH=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data/brain-ellis
MATRIXFILE="brain-ellis-mx-clean-transposed.tsv"
JSONFILE="gmt_mz.json"

Rscript gsva_NOinLog.R $DATAPATH $MATRIXFILE $JSONFILE




