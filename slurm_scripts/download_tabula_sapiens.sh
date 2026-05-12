#!/bin/bash
#SBATCH --job-name=dl_tabula
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=12:00:00
#SBATCH --mem=32000MB
#SBATCH --output=logs/dl_tabula-%j.out

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations

source $HOME/.conda/etc/profile.d/conda.sh
conda activate scdata

echo "Downloading Tabula Sapiens (10x only) via CellxGene Census..."
python preprocessing/download_tabula_sapiens.py

echo "Done. Exit code: $?"

# Submission:
# sbatch slurm_scripts/download_tabula_sapiens.sh
