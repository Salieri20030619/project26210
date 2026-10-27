#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
BASE="$ROOT/09_metabolism_trajectory"
OUT="$BASE/trajectory_v1"
PY_BASE=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PY_TRAJ=/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1/bin/python
RS_TRAJ=/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1/bin/Rscript

on_exit() {
  status=$?
  if [[ $status -ne 0 && ! -e "$OUT/RUN_COMPLETE" ]]; then
    printf 'failed_at\texit_code\n%s\t%s\n' "$(date --iso-8601=seconds)" "$status" > "$OUT/failure_status_v1.tsv"
    touch "$OUT/RUN_FAILED"
  fi
}
trap on_exit EXIT

mkdir -p "$OUT/logs" "$OUT/qc"
{
  printf 'started_at\t%s\n' "$(date --iso-8601=seconds)"
  printf 'base_python\t%s\n' "$PY_BASE"
  printf 'trajectory_python\t%s\n' "$PY_TRAJ"
  printf 'trajectory_Rscript\t%s\n' "$RS_TRAJ"
  "$PY_BASE" --version 2>&1
  "$PY_TRAJ" --version 2>&1
  "$RS_TRAJ" -e 'cat(R.version.string, "\n")' 2>&1
} > "$OUT/qc/executable_provenance_v1.tsv"

"$PROJECT/check_memory_gate.sh" part1_trajectory_master_input_v1 | tee "$OUT/qc/memory_gate_master_input_v1.tsv"
if [[ ! -e "$OUT/input/INPUT_COMPLETE" ]]; then
  "$PY_BASE" "$BASE/prepare_part1_trajectory_input_v1.py"
fi
if [[ ! -e "$OUT/R_METHODS_COMPLETE" ]]; then
  "$RS_TRAJ" --no-init-file "$BASE/run_part1_trajectory_r_methods_v1.R"
fi
if [[ ! -e "$OUT/PYTHON_METHODS_COMPLETE" ]]; then
  "$PY_TRAJ" "$BASE/run_part1_trajectory_python_methods_v1.py"
fi
"$PROJECT/check_memory_gate.sh" part1_trajectory_final_verification_v1 | tee "$OUT/qc/memory_gate_final_verification_v1.tsv"
"$PY_TRAJ" "$BASE/verify_part1_trajectory_v1.py"
printf 'PART1_TRAJECTORY_V1_ALL_PASS\n'
