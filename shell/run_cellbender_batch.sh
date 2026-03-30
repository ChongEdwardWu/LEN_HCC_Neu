#!/usr/bin/env bash
set -euo pipefail

###################### User-configurable parameters ######################
BASE_ROOT="path_to_cellranger_outputs"
OUT_DIR="path_to_cellbender_output"
THREADS=32
CHECKPOINT_MIN=120
CONDA_ENV="cellbender"
CONDA_INIT="path_to_conda_sh"
#######################################################################
mkdir -p "$OUT_DIR/"


# -------- 1. Batch log ---------
MASTER_LOG="$OUT_DIR/cellbender_batch.log"
echo "=== CellBender batch started $(date) ===" | tee  "$MASTER_LOG"

# -------- 2. Activate the conda environment --------
if [[ -f "$CONDA_INIT" ]]; then
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
else
    echo "ERROR: Conda init file not found at $CONDA_INIT" | tee -a "$MASTER_LOG"
    exit 1
fi
echo "Activated conda env: $(conda info --envs | grep '*' | awk '{print $1}')" | tee -a "$MASTER_LOG"

# -------- 3. Iterate through samples --------
find "$BASE_ROOT" -type d -path "*/outs/raw_feature_bc_matrix" \
| sort | while read -r INPUT; do
    SAMPLE_DIR="$(basename "$(dirname "$(dirname  "$INPUT")")")"
    SAMPLE="${SAMPLE_DIR%}"
    mkdir -p "$OUT_DIR/${SAMPLE}/"
    OUT="$OUT_DIR/${SAMPLE}/${SAMPLE}_filtered.h5"
    LOG="$OUT_DIR/${SAMPLE}/${SAMPLE}.log"

    echo ">>> [$SAMPLE] START $(date)" | tee  -a "$MASTER_LOG" "$LOG"
    cellbender remove-background \
        --cpu-threads "$THREADS" \
        --checkpoint-mins "$CHECKPOINT_MIN" \
        --input  "$INPUT" \
        --output "$OUT" \
        >> "$LOG" 2>&1
    CKPT_SRC="$(dirname "$OUT")/ckpt.tar.gz"
    if [[ -f "$CKPT_SRC" ]]; then
    mv "$CKPT_SRC" "$OUT_DIR/${SAMPLE}/${SAMPLE}_ckpt.tar.gz"
fi
    echo ">>> [$SAMPLE] DONE  $(date)"  | tee -a "$MASTER_LOG" "$LOG"
done

echo "=== All samples finished $(date) ===" | tee -a "$MASTER_LOG"
