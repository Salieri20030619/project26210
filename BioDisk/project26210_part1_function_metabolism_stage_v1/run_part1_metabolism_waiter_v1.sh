#!/usr/bin/env bash
set -euo pipefail
PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
BASE="$ROOT/09_metabolism_trajectory"
OUT="$BASE/metabolism_v1"
SCENIC="$ROOT/07_pyscenic_v1"
STATUS="$OUT/background_status_v1.tsv"
R453=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript
PY_BATCH=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PY_SCFEA="$PROJECT/envs/part2_6_2_7_20260726_v1/python/bin/python"
printf 'time\tstage\tstatus\tdetail\n' >"$STATUS"
record(){ printf '%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" >>"$STATUS"; }
trap 'rc=$?; record controller FAIL "exit_code=$rc"; exit "$rc"' ERR
if [[ "${START_WITHOUT_TF_WAIT:-0}" == "1" ]]; then
  record wait_for_tf BYPASSED "manual concurrent launch after resource check: 256 cores, 419 GiB available, 0% iowait"
else
  record wait_for_tf STARTED "$SCENIC/OVERNIGHT_COMPLETE"
  while [[ ! -f "$SCENIC/OVERNIGHT_COMPLETE" ]]; do
    sleep 60
  done
  record wait_for_tf PASS scenic_and_metaviper_complete
fi
record auc STARTED KEGG_AUCell_11950_cells
"$R453" "$BASE/run_part1_metabolism_auc_v1.R" >"$OUT/logs/run_part1_metabolism_auc_v1.log" 2>&1
record auc PASS completed
record scfea_input STARTED current_raw_UMI
"$PY_BATCH" "$BASE/prepare_part1_scfea_input_v1.py" >"$OUT/logs/prepare_part1_scfea_input_v1.log" 2>&1
record scfea_input PASS completed
"$PROJECT/check_memory_gate.sh" part1_scfea_m168_v1 >"$OUT/qc/memory_gate_scfea_model_v1.txt"
record scfea_model STARTED M168_100_epochs_16_threads
"$PY_SCFEA" "$BASE/run_part1_scfea_m168_v1.py" --epochs 100 --threads 16 --out-dir "$OUT" >"$OUT/logs/run_part1_scfea_m168_v1.log" 2>&1
record scfea_model PASS completed
record scfea_analysis STARTED patient_paired
"$R453" "$BASE/analyze_part1_scfea_flux_v1.R" >"$OUT/logs/analyze_part1_scfea_flux_v1.log" 2>&1
record scfea_analysis PASS completed
record verification STARTED integrity_and_statistics
"$PY_BATCH" "$BASE/verify_part1_metabolism_v1.py" >"$OUT/logs/verify_part1_metabolism_v1.log" 2>&1
record verification PASS completed
record controller PASS all_stages_complete
echo PART1_METABOLISM_V1_PASS
