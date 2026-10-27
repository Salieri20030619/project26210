#!/usr/bin/env bash
set -euo pipefail
O=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase213_immune_restart_execution_20260919_v1
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
stage=I2
trap 'rc=$?; printf "stage=%s exit=%s\n" "$stage" "$rc" > "$O/I2_I3_EXIT.txt"' EXIT
"$PY" "$O/scripts/02_main_statistics.py" > "$O/I2.log" 2>&1
stage=I3
"$PY" "$O/scripts/03_robustness.py" > "$O/I3.log" 2>&1
stage=COMPLETE_I2_I3
