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
export SCVI_MODEL_DIR=/scratch/zeynepaydin21/scvi
export SCIMILARITY_MODEL_DIR=/scratch/zeynepaydin21/scimilarity/models/model_v1.1
export GENE_MAP_TSV=${PROJ_ROOT}/data/gene_map/gene_map_grch38_v114.tsv

ulimit -s unlimited
ulimit -l unlimited

SCRIPT=${PROJ_ROOT}/src/classification.R
RUN_ID=${SLURM_ARRAY_TASK_ID}
METHOD=${1:-reference}
ALGO=${ALGO:-glmnet}
LEVEL=${LEVEL:-broad}   # celltype label level: broad (default, 31) | compartment (6) | fine (121)

# broad is the default celltype label (no flag); other levels select labels_<level>.qs
if [ "$LEVEL" = "broad" ]; then LEVEL_ARG=""; else LEVEL_ARG="--label_level $LEVEL"; fi

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
TRAIN_PCTS=(40 60 80)
if [ "$TASK_ARG" = "both" ]; then TASKS=(celltype tissue); else TASKS=("$TASK_ARG"); fi

echo "Run ID: $RUN_ID  Method: $METHOD  Algorithm: $ALGO"
echo "==============================================================================="

run_classification() {
    local task=$1
    local train_pct=$2
    local extra=$3
    local desc=$4

    # label level applies to celltype only; tissue loads tissue.qs regardless
    local lvl=""
    if [ "$task" = "celltype" ]; then lvl="$LEVEL_ARG"; fi

    echo "Running: tabula_sapiens, task=$task, level=$LEVEL, train=${train_pct}%, $desc"
    $RUNNER "$SCRIPT" \
        -r "$RUN_ID" \
        -d tabula_sapiens \
        -t "$task" \
        -a "$ALGO" \
        -s "$train_pct" \
        $lvl \
        $extra
}

for task in "${TASKS[@]}"; do
# tissue is level-independent; only run it once (with the default broad job) so it
# isn't duplicated across the compartment/fine level jobs.
if [ "$task" = "tissue" ] && [ "$LEVEL" != "broad" ]; then continue; fi
for train_pct in "${TRAIN_PCTS[@]}"; do
    if [ "$METHOD" = "reference" ]; then
        run_classification "$task" "$train_pct" "" "reference"
        run_classification "$task" "$train_pct" "--n_hvg $N_HVG" "reference, n_hvg=$N_HVG"

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

    elif [ "$METHOD" = "scimilarity" ]; then
        run_classification "$task" "$train_pct" "-m scimilarity" "method=scimilarity"

    elif [ "$METHOD" = "scvi" ]; then
        run_classification "$task" "$train_pct" "-m scvi" "method=scvi"

    fi
done
done

echo "Done. Exit code: $?"

# Submission (17 donors -> nested split varies by run_id, use array=1-N).
# Separate jobs per celltype label level via LEVEL env (broad | compartment | fine); ALGO env default glmnet.
# tissue task runs only with LEVEL=broad (level-independent), so compartment/fine jobs are celltype-only.
# for level in broad compartment fine; do
#   for method in reference pca; do
#     sbatch --array=1-30 --export=ALL,LEVEL=$level --job-name=tabula_${level}_$method slurm_scripts/classify_tabula_sapiens.sh $method
#   done
#   for method in rff_lapl rff_gauss; do
#     sbatch --array=1-30 --gres=gpu:1 --export=ALL,LEVEL=$level --job-name=tabula_${level}_$method slurm_scripts/classify_tabula_sapiens.sh $method
#   done
#   sbatch --array=1-30 --gres=gpu:1 --mem=64G --export=ALL,LEVEL=$level --job-name=tabula_${level}_scvi slurm_scripts/classify_tabula_sapiens.sh scvi
# done
