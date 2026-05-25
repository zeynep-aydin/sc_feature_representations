#!/bin/bash
#SBATCH --job-name=prep_tabula
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/prep_tabula-%j.out

module load gnu14/14.3.0
module load openblas/0.3.30
module load lapack/3.11.0

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations

conda run -n scdata python preprocessing/prep_tabula_sapiens.py --proj_root $PROJ_ROOT && \
Rscript preprocessing/prep_tabula_sapiens.R

echo "Done. Exit code: $?"

# Submission:
# sbatch slurm_scripts/prep_tabula_sapiens.sh
