#!/usr/bin/env bash
set -Eeuo pipefail
P=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase259_spatial_GSEA_reference_style_20260926_v4
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
trap 's=$?; printf "%s\n" "$s" > "$P/build_EXIT_STATUS"; if [ "$s" -eq 0 ]; then date --iso-8601=seconds > "$P/build_COMPLETE"; fi' EXIT
/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -u "$P/code/make_figure.py" > "$P/logs/build.log" 2>&1
/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >> "$P/logs/memory_gate.log"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/pdftoppm -r 150 -singlefile -png "$P/delivery/26210_Spatial_A4_GSEA_v4.pdf" "$P/qc/final_pdf_render" > "$P/logs/render.log" 2>&1
