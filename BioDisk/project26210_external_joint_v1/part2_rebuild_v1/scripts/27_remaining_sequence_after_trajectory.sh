#!/usr/bin/env bash
set -euo pipefail
R=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1
Q=/mnt/BioDisk/ludingyu/project26210
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python
PYS="$Q/.venv-pyscenic/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
STAGE=dynamic_statistics
trap 'rc=$?; printf "stage=%s exit=%s\n" "$STAGE" "$rc" > "$R/logs/phase148_sequence_exit.txt"' EXIT
"$PY" "$R/scripts/20_trajectory_statistics.py" > "$R/logs/phase148_dynamics_v1.log" 2>&1
STAGE=regulon_input
"$PY" "$R/scripts/21_prepare_regulon_input.py" > "$R/logs/phase148_regulon_input_v1.log" 2>&1
"$PYS" "$R/scripts/22_build_new_regulon_looms.py" > "$R/logs/phase148_regulon_looms_v1.log" 2>&1
STAGE=new_pyscenic_network
bash "$R/scripts/23_run_new_pyscenic.sh" > "$R/logs/phase148_pyscenic_v1.log" 2>&1
STAGE=regulon_target_registry
"$PYS" "$R/scripts/24_extract_new_regulon_targets.py" > "$R/logs/phase148_regulon_registry_v1.log" 2>&1
STAGE=regulon_patient_statistics
"$PY" "$R/scripts/25_new_regulon_statistics.py" > "$R/logs/phase148_regulon_statistics_v1.log" 2>&1
STAGE=crosswalk_dependencies
"$PY" "$R/scripts/26_crosswalk_and_dependencies.py" > "$R/logs/phase148_crosswalk_v1.log" 2>&1
STAGE=complete_computation_await_model_review
printf '%s\n' "$STAGE"
