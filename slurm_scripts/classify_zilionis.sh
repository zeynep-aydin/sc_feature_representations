#!/bin/bash
#SBATCH --job-name=zil_classify
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=ai[18-26]
#SBATCH --partition=ai
#SBATCH --time=15-00
#SBATCH --mem=8000MB
#SBATCH --output=slurm_logs/%x-%A_%a.out

module load gnu14/14.3.0
module load openblas/0.3.30
module load lapack/3.11.0

cd /scratch/zeynepaydin21/sc_feature_representations/
export PROJ_ROOT=/scratch/zeynepaydin21/sc_feature_representations
export SCVI_MODEL_DIR=/scratch/zeynepaydin21/scvi
export SCIMILARITY_MODEL_DIR=/scratch/zeynepaydin21/scimilarity/models/model_v1.1
export GENE_MAP_TSV=${PROJ_ROOT}/data/gene_map/gene_map_grch38_v114.tsv

ulimit -s unlimited
ulimit -l unlimited

SCRIPT=${PROJ_ROOT}/src/classification.R
RUN_ID=${SLURM_ARRAY_TASK_ID}
METHOD=${1:-reference}
ALGO=${ALGO:-glmnet}

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
HVG_ONLY=${HVG_ONLY:-0}  # pca: skip the no-HVG run (already done), submit only --n_hvg
PCA_DIMS=(32 64 128 512)
RFF_DIMS=(2048 1024 512 128)   # highest first
TRAIN_PCTS=(40 60 80)
if [ "$TASK_ARG" = "both" ]; then TASKS=(celltype tissue); else TASKS=("$TASK_ARG"); fi

echo "Run ID: $RUN_ID  Method: $METHOD  Algorithm: $ALGO"
echo "==============================================================================="

FAILS=0

run_classification() {
    local task=$1
    local train_pct=$2
    local extra=$3
    local desc=$4

    echo "Running: zilionis_lung, task=$task, train=${train_pct}%, $desc"
    $RUNNER "$SCRIPT" \
        -r "$RUN_ID" \
        -d zilionis_lung \
        -t "$task" \
        -a "$ALGO" \
        -s "$train_pct" \
        $extra
    local status=$?
    if [ "$status" -ne 0 ]; then
        echo "FAILED (exit $status): zilionis_lung task=$task train=${train_pct}% $desc"
        FAILS=$((FAILS + 1))
    fi
}

for task in "${TASKS[@]}"; do
    if [ "$METHOD" = "reference" ]; then
        for train_pct in "${TRAIN_PCTS[@]}"; do
            run_classification "$task" "$train_pct" "" "reference"
            run_classification "$task" "$train_pct" "--n_hvg $N_HVG" "reference, n_hvg=$N_HVG"
        done

    elif [ "$METHOD" = "rff_lapl" ] || [ "$METHOD" = "rff_gauss" ]; then
        for n_dim in "${RFF_DIMS[@]}"; do
        for train_pct in "${TRAIN_PCTS[@]}"; do
            run_classification "$task" "$train_pct" \
                "-m $METHOD -n $n_dim" \
                "method=$METHOD, n_dim=$n_dim"
            run_classification "$task" "$train_pct" \
                "-m $METHOD -n $n_dim --n_hvg $N_HVG" \
                "method=$METHOD, n_dim=$n_dim, n_hvg=$N_HVG"
        done
        done

    elif [ "$METHOD" = "pca" ]; then
        for n_dim in "${PCA_DIMS[@]}"; do
        for train_pct in "${TRAIN_PCTS[@]}"; do
            if [ "$HVG_ONLY" != "1" ]; then
                run_classification "$task" "$train_pct" \
                    "-m pca -n $n_dim" \
                    "method=pca, n_dim=$n_dim"
            fi
            run_classification "$task" "$train_pct" \
                "-m pca -n $n_dim --n_hvg $N_HVG" \
                "method=pca, n_dim=$n_dim, n_hvg=$N_HVG"
        done
        done

    elif [ "$METHOD" = "scimilarity" ]; then
        for train_pct in "${TRAIN_PCTS[@]}"; do
            run_classification "$task" "$train_pct" "-m scimilarity" "method=scimilarity"
        done

    elif [ "$METHOD" = "scvi" ]; then
        for train_pct in "${TRAIN_PCTS[@]}"; do
            run_classification "$task" "$train_pct" "-m scvi" "method=scvi"
        done

    fi
done

if [ "$FAILS" -eq 0 ]; then
    echo "Done. All runs succeeded."
else
    echo "Done. $FAILS run(s) FAILED (see FAILED lines above)."
fi
exit $(( FAILS > 0 ? 1 : 0 ))

# Submission (algorithm via ALGO env; default glmnet):
# sbatch --array=1 --job-name=zilionis_reference slurm_scripts/classify_zilionis.sh reference
# sbatch --array=1 --job-name=zilionis_pca slurm_scripts/classify_zilionis.sh pca
# sbatch --array=1 --export=ALL,HVG_ONLY=1 --job-name=zilionis_pca_hvg slurm_scripts/classify_zilionis.sh pca  # no-HVG pca already done
# sbatch --array=1 --gres=gpu:1 --job-name=zilionis_rff_lapl slurm_scripts/classify_zilionis.sh rff_lapl
# sbatch --array=1 --gres=gpu:1 --job-name=zilionis_rff_gauss slurm_scripts/classify_zilionis.sh rff_gauss
# sbatch --array=1 --gres=gpu:1 --mem=20G --job-name=zilionis_scimilarity slurm_scripts/classify_zilionis.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --mem=20G --job-name=zilionis_scvi slurm_scripts/classify_zilionis.sh scvi
