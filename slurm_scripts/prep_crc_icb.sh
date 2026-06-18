#!/bin/bash
#SBATCH --job-name=prep_crc_icb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/prep_crc_icb-%j.out

module load gnu14/14.3.0
module load openblas/0.3.30
module load lapack/3.11.0

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations

RAW=$PROJ_ROOT/raw_data/crc_icb
mkdir -p "$RAW"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236581/suppl"
for f in GSE236581_barcodes.tsv.gz GSE236581_features.tsv.gz GSE236581_counts.mtx.gz GSE236581_CRC-ICB_metadata.txt.gz; do
  if [ ! -f "$RAW/$f" ]; then
    echo "Downloading $f..."
    curl -sS --retry 5 --retry-delay 5 -C - -o "$RAW/$f" "$BASE/$f"
    if [ $? -ne 0 ]; then
      echo "Download failed for $f" >&2
      exit 1
    fi
  fi
done

Rscript preprocessing/prep_crc_icb.R

echo "Done. Exit code: $?"

# Submission:
# sbatch slurm_scripts/prep_crc_icb.sh
