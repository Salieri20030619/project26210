#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
SCRIPT="${ROOT}/run_glmf12_cell2location_reference_v3.py"
OUT="${ROOT}/glmf12_cell2location_reference_v3"
SAMPLE="${1:-}"

case "${SAMPLE}" in
  GLMF1|GLMF2) ;;
  *)
    printf 'Usage: %s GLMF1|GLMF2\n' "$0" >&2
    exit 2
    ;;
esac

mkdir -p "${OUT}/logs"
bash "/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh" \
  > "${OUT}/logs/${SAMPLE}_launcher_memory_gate.log" 2>&1

printf 'started_at\t%s\nsample\t%s\nlauncher_pid\t%s\nthreads\t32\n' \
  "$(date --iso-8601=seconds)" "${SAMPLE}" "$$" \
  > "${OUT}/${SAMPLE}_launcher_status.tsv"

set +e
env MPLCONFIGDIR="/tmp/project26210_glmf_${SAMPLE}_mpl_v3" \
    NUMBA_CACHE_DIR="/tmp/project26210_glmf_${SAMPLE}_numba_v3" \
    OMP_NUM_THREADS=32 \
    MKL_NUM_THREADS=32 \
    OPENBLAS_NUM_THREADS=32 \
    "${PYTHON}" "${SCRIPT}" --sample "${SAMPLE}" \
    > "${OUT}/logs/${SAMPLE}.log" 2>&1
rc=$?
set -e

printf 'finished_at\t%s\nexit_code\t%s\n' "$(date --iso-8601=seconds)" "${rc}" \
  >> "${OUT}/${SAMPLE}_launcher_status.tsv"
exit "${rc}"
