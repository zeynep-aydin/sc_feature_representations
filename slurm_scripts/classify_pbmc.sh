#!/bin/bash
#SBATCH --job-name=pbmc_classify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=15-00
#SBATCH --mem=16000MB
#SBATCH --output=slurm_logs/%x-%A_%a.out

module load gnu14/14.3.0
module load openblas/0.3.30
module load lapack/3.11.0

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations

ulimit -s unlimited
ulimit -l unlimited

SCRIPT=${PROJ_ROOT}/src/classification.R
RUN_ID=${SLURM_ARRAY_TASK_ID}
METHOD=${1:-baseline}
N_HVG=2000
N_DIM_VALUES=(64 128 512 1024)
TRAIN_PCTS=(50) # 70 80
LABEL_LEVELS=(l1 l2 l3)

echo "Run ID: $RUN_ID  Method: $METHOD"
echo "==============================================================================="

run_classification() {
    local label_level=$1
    local train_pct=$2
    local extra=$3
    local desc=$4

    echo "Running: pbmc, label=$label_level, train=${train_pct}%, $desc"
    Rscript "$SCRIPT" \
        -r "$RUN_ID" \
        -d pbmc \
        -t celltype \
        -a glmnet \
        -s "$train_pct" \
        --label_level "$label_level" \
        $extra
}

for label_level in "${LABEL_LEVELS[@]}"; do
for train_pct in "${TRAIN_PCTS[@]}"; do
    if [ "$METHOD" = "baseline" ]; then
        run_classification "$label_level" "$train_pct" "" "baseline"
        run_classification "$label_level" "$train_pct" "--n_hvg $N_HVG" "baseline, n_hvg=$N_HVG"

    elif [ "$METHOD" = "rff_lapl" ] || [ "$METHOD" = "rff_gauss" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$label_level" "$train_pct" \
                "-m $METHOD -n $n_dim" \
                "method=$METHOD, n_dim=$n_dim"
            run_classification "$label_level" "$train_pct" \
                "-m $METHOD -n $n_dim --n_hvg $N_HVG" \
                "method=$METHOD, n_dim=$n_dim, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "pca" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$label_level" "$train_pct" \
                "-m pca -n $n_dim" \
                "method=pca, n_dim=$n_dim"
            run_classification "$label_level" "$train_pct" \
                "-m pca -n $n_dim --n_hvg $N_HVG" \
                "method=pca, n_dim=$n_dim, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "scimilarity" ]; then
        run_classification "$label_level" "$train_pct" "-m scimilarity" "method=scimilarity"

    elif [ "$METHOD" = "scvi" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$label_level" "$train_pct" \
                "-m scvi -n $n_dim --max_epochs 10" \
                "method=scvi, n_dim=$n_dim, epochs=10"
            run_classification "$label_level" "$train_pct" \
                "-m scvi -n $n_dim --max_epochs 100" \
                "method=scvi, n_dim=$n_dim, epochs=100"
        done
    fi
done
done

echo "Done. Exit code: $?"

# Submission (8 donors -> exhaustive split, array=1):
# for method in baseline pca rff_lapl rff_gauss; do
#     sbatch --array=1 --job-name=pbmc_$method classify_pbmc.sh $method
# done
# sbatch --array=1 --gres=gpu:1 --mem=20GB --job-name=pbmc_scimilarity classify_pbmc.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --job-name=pbmc_scvi classify_pbmc.sh scvi
