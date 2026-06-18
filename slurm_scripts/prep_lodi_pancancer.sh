#!/bin/bash
#SBATCH --job-name=prep_lodi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/prep_lodi-%j.out

module load gnu14/14.3.0
module load openblas/0.3.30
module load lapack/3.11.0

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations

RAW=$PROJ_ROOT/raw_data/lodi_pancancer

Rscript preprocessing/prep_lodi_pancancer.R

echo "Done. Exit code: $?"

# Submission:
# sbatch slurm_scripts/prep_lodi_pancancer.sh
