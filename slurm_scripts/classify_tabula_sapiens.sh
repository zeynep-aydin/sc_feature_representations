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
DEFAULT_LEVEL=broad             # granularity stored in labels.qs (selected with no --label_level flag)
LEVEL=${LEVEL:-$DEFAULT_LEVEL}  # celltype level: broad (31) | compartment (6) | fine (121)

# default level lives in labels.qs (no flag); any other level selects labels_<level>.qs
if [ "$LEVEL" = "$DEFAULT_LEVEL" ]; then LEVEL_ARG=""; else LEVEL_ARG="--label_level $LEVEL"; fi

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

    # label_level applies to the celltype task only; tissue always loads tissue.qs regardless of LEVEL
    local label_level_arg=""
    local level_desc="-"
    if [ "$task" = "celltype" ]; then
        label_level_arg="$LEVEL_ARG"
        level_desc="$LEVEL"
    fi

    echo "Running: tabula_sapiens, task=$task, level=$level_desc, train=${train_pct}%, $desc"
    $RUNNER "$SCRIPT" \
        -r "$RUN_ID" \
        -d tabula_sapiens \
        -t "$task" \
        -a "$ALGO" \
        -s "$train_pct" \
        $label_level_arg \
        $extra
    local status=$?
    if [ "$status" -ne 0 ]; then
        echo "FAILED (exit $status): tabula_sapiens task=$task level=$level_desc train=${train_pct}% $desc"
        FAILS=$((FAILS + 1))
    fi
}

for task in "${TASKS[@]}"; do
    # tissue is level-independent; run it only under the default level so it isn't
    # duplicated across the compartment/fine jobs. Announce the skip so a level job
    # that quietly drops tissue is visible in the log, not silent.
    if [ "$task" = "tissue" ] && [ "$LEVEL" != "$DEFAULT_LEVEL" ]; then
        echo "Skipping tissue task (level-independent; runs only under LEVEL=$DEFAULT_LEVEL, current LEVEL=$LEVEL)"
        continue
    fi
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
            run_classification "$task" "$train_pct" \
                "-m pca -n $n_dim" \
                "method=pca, n_dim=$n_dim"
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

# Submission (17 donors -> nested split varies by run_id, use array=1-N).
# Separate jobs per celltype label level via LEVEL env (broad | compartment | fine); ALGO env default glmnet.
# tissue task runs only with LEVEL=broad (level-independent), so compartment/fine jobs are celltype-only.
#
# LEVEL=broad (default; includes tissue task):
# sbatch --array=1 --export=ALL,LEVEL=broad --job-name=tabula_broad_reference slurm_scripts/classify_tabula_sapiens.sh reference
# sbatch --array=1 --export=ALL,LEVEL=broad --job-name=tabula_broad_pca slurm_scripts/classify_tabula_sapiens.sh pca
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=broad --job-name=tabula_broad_rff_lapl slurm_scripts/classify_tabula_sapiens.sh rff_lapl
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=broad --job-name=tabula_broad_rff_gauss slurm_scripts/classify_tabula_sapiens.sh rff_gauss
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=broad --job-name=tabula_broad_scimilarity slurm_scripts/classify_tabula_sapiens.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=broad --job-name=tabula_broad_scvi slurm_scripts/classify_tabula_sapiens.sh scvi
#
# LEVEL=compartment (celltype task only):
# sbatch --array=1 --export=ALL,LEVEL=compartment --job-name=tabula_compartment_reference slurm_scripts/classify_tabula_sapiens.sh reference
# sbatch --array=1 --export=ALL,LEVEL=compartment --job-name=tabula_compartment_pca slurm_scripts/classify_tabula_sapiens.sh pca
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=compartment --job-name=tabula_compartment_rff_lapl slurm_scripts/classify_tabula_sapiens.sh rff_lapl
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=compartment --job-name=tabula_compartment_rff_gauss slurm_scripts/classify_tabula_sapiens.sh rff_gauss
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=compartment --job-name=tabula_compartment_scimilarity slurm_scripts/classify_tabula_sapiens.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=compartment --job-name=tabula_compartment_scvi slurm_scripts/classify_tabula_sapiens.sh scvi
#
# LEVEL=fine (celltype task only):
# sbatch --array=1 --export=ALL,LEVEL=fine --job-name=tabula_fine_reference slurm_scripts/classify_tabula_sapiens.sh reference
# sbatch --array=1 --export=ALL,LEVEL=fine --job-name=tabula_fine_pca slurm_scripts/classify_tabula_sapiens.sh pca
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=fine --job-name=tabula_fine_rff_lapl slurm_scripts/classify_tabula_sapiens.sh rff_lapl
# sbatch --array=1 --gres=gpu:1 --export=ALL,LEVEL=fine --job-name=tabula_fine_rff_gauss slurm_scripts/classify_tabula_sapiens.sh rff_gauss
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=fine --job-name=tabula_fine_scimilarity slurm_scripts/classify_tabula_sapiens.sh scimilarity
# sbatch --array=1 --gres=gpu:1 --mem=128G --export=ALL,LEVEL=fine --job-name=tabula_fine_scvi slurm_scripts/classify_tabula_sapiens.sh scvi
