#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_phase73_virtual_knockout_v1
BASE71=/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1
BASE64=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1
BASE67=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1
GATE=/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh
RSCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

mkdir -p "$ROOT/logs_v1" "$ROOT/dataset_runs_v1" "$ROOT/results_v1" "$ROOT/independent_validation_v1"

run_dataset() {
  local dataset="$1"
  local rds="$2"
  local metadata="$3"
  local output="$ROOT/dataset_runs_v1/$dataset"
  if [[ -f "$output/DATASET_COMPLETE" ]]; then
    echo "SKIP_COMPLETED_DATASET $dataset"
    return
  fi
  "$GATE" | tee "$ROOT/logs_v1/memory_gate_${dataset}_v1.log"
  "$RSCRIPT" --vanilla "$ROOT/01_run_phase73_dataset_v1.R" \
    "$dataset" "$rds" "$metadata" "$output" \
    2>&1 | tee "$ROOT/logs_v1/${dataset}_run_v1.log"
}

run_dataset \
  GSE207493 \
  "$BASE71/inputs_v1/seurat/GSE207493_annotated_local_v2.rds" \
  "$BASE64/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz"

run_dataset \
  GSE264586 \
  "$BASE71/inputs_v1/seurat/GSE264586_annotated_local_v2.rds" \
  "$BASE67/inputs_v1/metadata_v2/GSE264586_phase67_metadata_v1.tsv.gz"

run_dataset \
  GSE294109 \
  "$BASE71/inputs_v1/seurat/GSE294109_annotated_local_v2.rds" \
  "$BASE67/inputs_v1/metadata_v2/GSE294109_phase67_metadata_v1.tsv.gz"

"$GATE" | tee "$ROOT/logs_v1/memory_gate_finalize_v1.log"
"$PYTHON" "$ROOT/02_finalize_phase73_v1.py" 2>&1 | tee "$ROOT/logs_v1/finalize_v1.log"

"$GATE" | tee "$ROOT/logs_v1/memory_gate_validate_v1.log"
"$PYTHON" "$ROOT/03_validate_phase73_v1.py" 2>&1 | tee "$ROOT/logs_v1/validate_v1.log"

date --iso-8601=seconds > "$ROOT/RUN_COMPLETE"
echo PHASE73_BACKGROUND_PIPELINE_PASS
