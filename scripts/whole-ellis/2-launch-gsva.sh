#!/bin/bash
########################## Slurm options
#SBATCH --job-name=gsva-whole
#SBATCH --output=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/slurmout/gsvaw_%j.out
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

DATAPATH=/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/data/whole-ellis
MATRIXFILE="whole-ellis_mx_f2.tsv"
JSONFILE="gmt_mz.json"

### running with script that does log and then gsva !
Rscript zonist/src/r/run_gsva.R $DATAPATH $MATRIXFILE $JSONFILE




