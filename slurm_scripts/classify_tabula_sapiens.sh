#!/bin/bash
#SBATCH --job-name=tabula_classify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=30-00
#SBATCH --mem=64G
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

if [[ "$METHOD" == rff_* ]]; then
    module load singularity/4.3.2
    module load cuda/12.8.0
    SIF=/scratch/zeynepaydin21/containers/r-ver_4.6.0.sif
    export SINGULARITY_BIND="/scratch/zeynepaydin21/containers/r_libs:/r_libs,/opt/ohpc/pub/compiler/cuda/12.8.0/lib64:/cuda/lib64"
    export SINGULARITY_NV=1
    export SINGULARITYENV_LD_LIBRARY_PATH=/cuda/lib64
    export SINGULARITYENV_PYTHON3=$(which python3)
    export R_LIBS=/r_libs
    RUNNER="singularity exec $SIF Rscript"
else
    RUNNER="Rscript"
fi

TASK_ARG=${2:-both}
N_HVG=2000
N_DIM_VALUES=(64 128 512 1024 2048)
TRAIN_PCTS=(50 70 80)
if [ "$TASK_ARG" = "both" ]; then TASKS=(celltype tissue); else TASKS=("$TASK_ARG"); fi

echo "Run ID: $RUN_ID  Method: $METHOD"
echo "==============================================================================="

run_classification() {
    local task=$1
    local train_pct=$2
    local extra=$3
    local desc=$4

    echo "Running: tabula_sapiens, task=$task, train=${train_pct}%, $desc"
    $RUNNER "$SCRIPT" \
        -r "$RUN_ID" \
        -d tabula_sapiens \
        -t "$task" \
        -a glmnet \
        -s "$train_pct" \
        $extra
}

for task in "${TASKS[@]}"; do
for train_pct in "${TRAIN_PCTS[@]}"; do
    if [ "$METHOD" = "baseline" ]; then
        run_classification "$task" "$train_pct" "" "baseline"
        run_classification "$task" "$train_pct" "--n_hvg $N_HVG" "baseline, n_hvg=$N_HVG"

    elif [ "$METHOD" = "rff_lapl" ] || [ "$METHOD" = "rff_gauss" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$task" "$train_pct" \
                "-m $METHOD -n $n_dim" \
                "method=$METHOD, n_dim=$n_dim"
            run_classification "$task" "$train_pct" \
                "-m $METHOD -n $n_dim --n_hvg $N_HVG" \
                "method=$METHOD, n_dim=$n_dim, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "pca" ]; then
        for n_dim in "${N_DIM_VALUES[@]}"; do
            run_classification "$task" "$train_pct" \
                "-m pca -n $n_dim" \
                "method=pca, n_dim=$n_dim"
        done

    # scimilarity: needs ENSEMBL->symbol conversion + 53 duplicate symbol strategy TBD
    # elif [ "$METHOD" = "scimilarity" ]; then
    #     run_classification "$task" "$train_pct" "-m scimilarity" "method=scimilarity"

    # scvi: get_scvi_embeddings.py assumes symbol input; tabula_sapiens has ENSEMBL IDs — needs update
    # elif [ "$METHOD" = "scvi" ]; then
    #     run_classification "$task" "$train_pct" "-m scvi" "method=scvi"

    fi
done
done

echo "Done. Exit code: $?"

# Submission (17 donors -> exhaustive split, array=1):
# for method in baseline pca; do
#     sbatch --array=1 --job-name=tabula_$method slurm_scripts/classify_tabula_sapiens.sh $method
# done
# for method in rff_lapl rff_gauss; do
#     sbatch --array=1 --gres=gpu:1 --job-name=tabula_$method slurm_scripts/classify_tabula_sapiens.sh $method
# done
