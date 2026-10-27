#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1
SCENIC_ROOT="$PROJECT/results_v1/GSE207493/tf_regulon_de_novo_v1/scenic_de_novo_v1"
VIPER_ROOT="$PROJECT/results_v1/GSE207493/tf_regulon_de_novo_v1/collectri_viper_v1"
INTEGRATED_ROOT="$PROJECT/results_v1/GSE207493/tf_regulon_de_novo_v1/integrated_primary_v1"
SCRIPTS="$PROJECT/scripts_v1"
LOG="$PROJECT/logs_v1/gse207493_de_novo_pipeline_resume_v1.log"
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RSCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
WATCH_PID="${1:-}"

exec >> "$LOG" 2>&1
echo "$(date -Is) RESUME_WATCHER_STARTED watch_pid=$WATCH_PID"
if [[ -n "$WATCH_PID" ]]; then
  while kill -0 "$WATCH_PID" 2>/dev/null; do
    sleep 30
  done
fi
sleep 15

if [[ ! -f "$SCENIC_ROOT/RUN_COMPLETE" ]]; then
  echo "$(date -Is) SCENIC_RESUME_STARTED"
  "$SCRIPTS/34_run_gse207493_de_novo_scenic_v1.sh"
fi
if [[ ! -f "$SCENIC_ROOT/ANALYSIS_COMPLETE" ]]; then
  echo "$(date -Is) SCENIC_ANALYSIS_STARTED"
  "$PY" "$SCRIPTS/37_analyze_gse207493_de_novo_scenic_v1.py"
fi
if [[ ! -f "$VIPER_ROOT/ACTIVITY_COMPLETE" ]]; then
  echo "$(date -Is) COLLECTRI_VIPER_STARTED"
  "$RSCRIPT" --vanilla "$SCRIPTS/35_run_gse207493_collectri_viper_de_novo_v1.R"
fi
if [[ ! -f "$VIPER_ROOT/RUN_COMPLETE" ]]; then
  echo "$(date -Is) COLLECTRI_VIPER_FINALIZE_STARTED"
  "$PY" "$SCRIPTS/36_finalize_collectri_viper_de_novo_v1.py"
fi
if [[ ! -f "$INTEGRATED_ROOT/RUN_COMPLETE" ]]; then
  echo "$(date -Is) PRIMARY_INTEGRATION_STARTED"
  "$PY" "$SCRIPTS/38_integrate_gse207493_de_novo_tf_regulon_v1.py"
fi
echo "$(date -Is) DE_NOVO_TF_PIPELINE_V1_PASS"
