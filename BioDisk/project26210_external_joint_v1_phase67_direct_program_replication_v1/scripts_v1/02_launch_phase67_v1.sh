#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
SCRIPT="$ROOT/scripts_v1/01_run_phase67_direct_replication_v1.py"
H5AD=/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad
METADATA="$ROOT/inputs_v1/metadata_v2"
SIGNATURES="$ROOT/inputs_v1/signatures_v1"
OUTPUT="$ROOT/results_v3"
LOG="$ROOT/logs_v1/run_phase67_v3.log"

mkdir -p "$ROOT/logs_v1"
exec > >(tee "$LOG") 2>&1

on_error() {
  code=$?
  printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "$code" "$(date -Is)" > "$ROOT/PHASE67_V3_FAIL_v1"
  exit "$code"
}
trap on_error ERR

printf 'status=RUNNING\nstarted_at=%s\npython=%s\n' "$(date -Is)" "$PY" > "$ROOT/PHASE67_V3_RUNNING_v1"
printf '%s\n' "$PY $SCRIPT --h5ad $H5AD --metadata-dir $METADATA --input-dir $SIGNATURES --cell-cycle-genes $SIGNATURES/seurat_cell_cycle_genes_v1.tsv --output $OUTPUT" > "$ROOT/logs_v1/command_v3.txt"

/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh
"$PY" "$SCRIPT" \
  --h5ad "$H5AD" \
  --metadata-dir "$METADATA" \
  --input-dir "$SIGNATURES" \
  --cell-cycle-genes "$SIGNATURES/seurat_cell_cycle_genes_v1.tsv" \
  --output "$OUTPUT"

printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "$ROOT/PHASE67_V3_PASS_v1"
