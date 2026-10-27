#!/usr/bin/env bash
set -Eeuo pipefail
P=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase257_TYMP_lysosome_cholesterol_validation_20260926_v1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 's=$?; printf "%s\n" "$s" > "$P/analysis_EXIT_STATUS"; if [ "$s" -eq 0 ]; then date --iso-8601=seconds > "$P/analysis_COMPLETE"; fi' EXIT
bash "$P/code/run_spatial.sh"
/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$P/code/paired.py" > "$P/logs/paired.log" 2>&1
