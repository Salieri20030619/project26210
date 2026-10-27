#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_phase60_v1
SOURCE=/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
F1=/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/cnv_clone_program_function_v1/rank7_program_top100_genes_v1.tsv
SCRIPTS="$ROOT/scripts"
LOGS="$ROOT/logs"
STATUS="$ROOT/status"
MANIFEST="$ROOT/00_manifest"
QC="$ROOT/01_qc"
VERIFY="$ROOT/02_verification"
INTEGRATION="$ROOT/03_integration_annotation"

mkdir -p "$LOGS" "$STATUS" "$MANIFEST" "$QC" "$VERIFY" "$INTEGRATION"
export OMP_NUM_THREADS=24
export OPENBLAS_NUM_THREADS=24
export MKL_NUM_THREADS=24
export NUMEXPR_NUM_THREADS=24
export PYTHONUNBUFFERED=1

on_error() {
  code=$?
  printf 'status\tFAILED\nexit_code\t%s\nfailed_at\t%s\n' "$code" "$(date -Iseconds)" > "$STATUS/FAILED_v1.tsv"
  exit "$code"
}
trap on_error ERR

printf 'status\tRUNNING\nstarted_at\t%s\nhost\t%s\npid\t%s\n' "$(date -Iseconds)" "$(hostname)" "$$" > "$STATUS/RUNNING_v1.tsv"
rm -f "$STATUS/FAILED_v1.tsv" "$STATUS/RUN_COMPLETE_v1.tsv"

echo "[$(date -Iseconds)] STAGE manifest"
"$PY" "$SCRIPTS/01_build_library_manifest_v1.py" \
  --download-root "$SOURCE/downloads" \
  --metadata-root "$SOURCE/metadata" \
  --output "$MANIFEST"
printf 'stage\tmanifest\nstatus\tPASS\ncompleted_at\t%s\n' "$(date -Iseconds)" > "$STATUS/STAGE_00_PASS.tsv"

echo "[$(date -Iseconds)] STAGE qc_build"
"$PY" "$SCRIPTS/02_qc_build_anndata_v1.py" \
  --manifest "$MANIFEST/library_manifest_v1.tsv" \
  --output "$QC"
printf 'stage\tqc_build\nstatus\tPASS\ncompleted_at\t%s\n' "$(date -Iseconds)" > "$STATUS/STAGE_01_PASS.tsv"

echo "[$(date -Iseconds)] STAGE qc_verify"
"$PY" "$SCRIPTS/03_verify_qc_v1.py" \
  --manifest "$MANIFEST/library_manifest_v1.tsv" \
  --qc-root "$QC" \
  --output "$VERIFY/qc_v1"
printf 'stage\tqc_verify\nstatus\tPASS\ncompleted_at\t%s\n' "$(date -Iseconds)" > "$STATUS/STAGE_02_PASS.tsv"

echo "[$(date -Iseconds)] STAGE integration_annotation_screen"
"$PY" "$SCRIPTS/04_integrate_annotate_screen_v1.py" \
  --qc-root "$QC" \
  --f1-signature "$F1" \
  --output "$INTEGRATION"
printf 'stage\tintegration_annotation_screen\nstatus\tPASS\ncompleted_at\t%s\n' "$(date -Iseconds)" > "$STATUS/STAGE_03_PASS.tsv"

echo "[$(date -Iseconds)] STAGE integration_verify"
"$PY" "$SCRIPTS/05_verify_integration_v1.py" \
  --integration-root "$INTEGRATION" \
  --qc-root "$QC" \
  --output "$VERIFY/integration_v1"
printf 'stage\tintegration_verify\nstatus\tPASS\ncompleted_at\t%s\n' "$(date -Iseconds)" > "$STATUS/STAGE_04_PASS.tsv"

printf 'status\tPASS\ncompleted_at\t%s\nhost\t%s\npid\t%s\n' "$(date -Iseconds)" "$(hostname)" "$$" > "$STATUS/RUN_COMPLETE_v1.tsv"
rm -f "$STATUS/RUNNING_v1.tsv"
echo "[$(date -Iseconds)] PHASE60_OVERNIGHT_PASS"
