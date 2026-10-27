#!/usr/bin/env bash
set -Eeuo pipefail
P=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase256_spatial_region_function_20260926_v1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
G=/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh
trap 's=$?; printf "%s\n" "$s" > "$P/extract_EXIT_STATUS"; if [ "$s" -eq 0 ]; then date --iso-8601=seconds > "$P/extract_COMPLETE"; fi' EXIT
"$G" >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$P/code/export_libraries.R" "$P" > "$P/logs/export_libraries.log" 2>&1
"$G" >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$P/code/extract.py" > "$P/logs/extract.log" 2>&1
