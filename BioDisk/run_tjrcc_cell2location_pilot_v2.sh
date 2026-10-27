#!/usr/bin/env bash
set -u

PROJECT="/mnt/BioDisk/ludingyu/project26210"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
SCRIPT="${PROJECT}/external/ccRCC_spatial/run_tjrcc_cell2location_pilot_v2.py"
OUT="${PROJECT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/cell2location_pilot_v2"
LOGS="${OUT}/logs"

mkdir -p "${LOGS}"
printf '%s\n' "$$" > "${OUT}/launcher.pid"
printf 'started_at\t%s\npid\t%s\nexecution\tsequential_cpu\n' \
  "$(date --iso-8601=seconds)" "$$" > "${OUT}/launcher_status.tsv"

overall=0
for sample in R51_T X49_T X98_T; do
  printf '%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${sample}" "started" \
    >> "${OUT}/sample_transitions.tsv"
  env MPLCONFIGDIR=/tmp/project26210_mpl \
      NUMBA_CACHE_DIR=/tmp/project26210_numba \
      OMP_NUM_THREADS=32 \
      MKL_NUM_THREADS=32 \
      OPENBLAS_NUM_THREADS=32 \
      "${PYTHON}" "${SCRIPT}" --sample "${sample}" \
      > "${LOGS}/${sample}.log" 2>&1
  rc=$?
  if [ "${rc}" -eq 0 ]; then
    state="completed"
  else
    state="failed_rc_${rc}"
    overall=1
  fi
  printf '%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${sample}" "${state}" \
    >> "${OUT}/sample_transitions.tsv"
done

printf 'finished_at\t%s\nexit_code\t%s\n' "$(date --iso-8601=seconds)" "${overall}" \
  >> "${OUT}/launcher_status.tsv"
exit "${overall}"
