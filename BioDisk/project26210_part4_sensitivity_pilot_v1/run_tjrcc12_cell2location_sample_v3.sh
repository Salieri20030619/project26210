#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1"
OUT="${ROOT}/tjrcc12_cell2location_reference_v3"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
SCRIPT="${ROOT}/run_tjrcc12_cell2location_reference_v3.py"
SAMPLE="${1:-}"

case "${SAMPLE}" in
  H46_T|J38_T|R114_T|R29_T|R51_T|S15_T|X49_T|X98_T|Y12_T|Y27_T|Y7_T|Z43_T) ;;
  *) printf 'Invalid sample: %s\n' "${SAMPLE}" >&2; exit 2 ;;
esac

mkdir -p "${OUT}/logs"
bash "/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh" \
  > "${OUT}/logs/${SAMPLE}_launcher_memory_gate.log" 2>&1
printf 'started_at\t%s\nsample\t%s\nlauncher_pid\t%s\nthreads\t12\n' \
  "$(date --iso-8601=seconds)" "${SAMPLE}" "$$" \
  > "${OUT}/${SAMPLE}_launcher_status.tsv"

set +e
env MPLCONFIGDIR="/tmp/project26210_${SAMPLE}_mpl_refv3" \
    NUMBA_CACHE_DIR="/tmp/project26210_${SAMPLE}_numba_refv3" \
    OMP_NUM_THREADS=12 \
    MKL_NUM_THREADS=12 \
    OPENBLAS_NUM_THREADS=12 \
    "${PYTHON}" "${SCRIPT}" --sample "${SAMPLE}" \
    > "${OUT}/logs/${SAMPLE}.log" 2>&1
rc=$?
set -e
printf 'finished_at\t%s\nexit_code\t%s\n' "$(date --iso-8601=seconds)" "${rc}" \
  >> "${OUT}/${SAMPLE}_launcher_status.tsv"
exit "${rc}"
