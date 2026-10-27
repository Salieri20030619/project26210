#!/usr/bin/env bash
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u phase291_spatial_null_TYMP_matched_regions_20260930_v1/code/analyze.py > phase291_spatial_null_TYMP_matched_regions_20260930_v1/logs/analyze.log 2>&1
rc=$?
printf '%s\n' "$rc" > phase291_spatial_null_TYMP_matched_regions_20260930_v1/logs/analyze.exit
exit "$rc"
