#!/bin/bash
#SBATCH --job-name=zil_classify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=15-00
#SBATCH --mem=8000MB
#SBATCH --output=slurm_logs/zil-%A_%a.out

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

echo "Run ID: $RUN_ID  Method: $METHOD"
echo "==============================================================================="

run_classification() {
    local train_pct=$1
    local extra=$2
    local desc=$3

    echo "Running: zilionis_lung, train=${train_pct}%, $desc"
    Rscript "$SCRIPT" \
        -r "$RUN_ID" \
        -d zilionis_lung \
        -t celltype \
        -a glmnet \
        -s "$train_pct" \
        $extra
}

for train_pct in "${TRAIN_PCTS[@]}"; do
    if [ "$METHOD" = "baseline" ]; then
        run_classification "$train_pct" "" "baseline"
        run_classification "$train_pct" "--n_hvg $N_HVG" "baseline, n_hvg=$N_HVG"

    elif [ "$METHOD" = "rff_lapl" ] || [ "$METHOD" = "rff_gauss" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$train_pct" \
                "-m $METHOD -n $n_dim" \
                "method=$METHOD, n_dim=$n_dim"
            run_classification "$train_pct" \
                "-m $METHOD -n $n_dim --n_hvg $N_HVG" \
                "method=$METHOD, n_dim=$n_dim, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "pca" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$train_pct" \
                "-m pca -n $n_dim" \
                "method=pca, n_dim=$n_dim"
            run_classification "$train_pct" \
                "-m pca -n $n_dim --n_hvg $N_HVG" \
                "method=pca, n_dim=$n_dim, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "scimilarity" ]; then
        run_classification "$train_pct" "-m scimilarity" "method=scimilarity"

    elif [ "$METHOD" = "scvi" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$train_pct" \
                "-m scvi -n $n_dim --max_epochs 10" \
                "method=scvi, n_dim=$n_dim"
            run_classification "$train_pct" \
                "-m scvi -n $n_dim --max_epochs 100" \
                "method=scvi, n_dim=$n_dim"
        done
    fi
done

echo "Done. Exit code: $?"

# Submission (single run):
# for method in baseline pca rff_lapl rff_gauss; do
#     sbatch --array=1 --job-name=zil_$method classify_zilionis.sh $method
# done
# sbatch --array=1 --gres=gpu:1 --mem=20GB --job-name=zil_scimilarity classify_zilionis.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --job-name=zil_scvi classify_zilionis.sh scvi
