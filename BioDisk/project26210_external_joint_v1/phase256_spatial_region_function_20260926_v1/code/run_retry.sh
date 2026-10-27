#!/usr/bin/env bash
set -Eeuo pipefail
P=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase256_spatial_region_function_20260926_v1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 's=$?; printf "%s\n" "$s" > "$P/retry_EXIT_STATUS"; if [ "$s" -eq 0 ]; then date --iso-8601=seconds > "$P/retry_COMPLETE"; fi' EXIT
/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$P/code/retry_gsea_NA.R" "$P" > "$P/logs/retry_gsea.log" 2>&1
bash "$P/code/run_validation.sh"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python "$P/code/summarize.py" > "$P/logs/summary.log" 2>&1
