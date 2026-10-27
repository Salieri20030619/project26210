#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 PROJECT_ROOT WORK_ROOT CELL2LOC_PYTHON" >&2
  exit 2
fi
PROJECT_ROOT="$1"
WORK_ROOT="$2"
PYTHON="$3"
SCRIPT="${WORK_ROOT}/run_tjrcc12_cell2location_aux_resume_v4r1.py"
OUT="${WORK_ROOT}/tjrcc12_cell2location_reference_v4r1_aux_resume"
SAMPLES=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)
MAX_JOBS=4

if [[ ! -f "${OUT}/run_metadata_v4.json" ]]; then
  echo "Preflight metadata is missing" >&2
  exit 1
fi
mkdir -p "${OUT}/logs"
exec 9>"${OUT}/controlled_launcher_v4r1.lock"
if ! flock -n 9; then
  echo "Another TJRCC12 v4r1 launcher holds the lock" >&2
  exit 3
fi
PROCESS_TABLE="${OUT}/master_processes_v4r1.tsv"
COMPLETION_TABLE="${OUT}/master_completion_v4r1.tsv"
if [[ -e "${PROCESS_TABLE}" || -e "${COMPLETION_TABLE}" ]]; then
  echo "No-overwrite gate: launcher tables already exist" >&2
  exit 1
fi
printf 'sample\tpid\tstarted_at\n' >"${PROCESS_TABLE}.part"
mv "${PROCESS_TABLE}.part" "${PROCESS_TABLE}"
printf 'sample\tpid\texit_code\tfinished_at\n' >"${COMPLETION_TABLE}.part"
mv "${COMPLETION_TABLE}.part" "${COMPLETION_TABLE}"

declare -A PID_TO_SAMPLE=()
running=0
overall=0

reap_one() {
  local finished_pid rc sample
  set +e
  wait -n -p finished_pid
  rc=$?
  set -e
  sample="${PID_TO_SAMPLE[${finished_pid}]}"
  printf '%s\t%s\t%s\t%s\n' "${sample}" "${finished_pid}" "${rc}" "$(date --iso-8601=seconds)" \
    >>"${COMPLETION_TABLE}"
  unset 'PID_TO_SAMPLE['"${finished_pid}"']'
  running=$((running - 1))
  if [[ "${rc}" -ne 0 ]]; then overall=1; fi
}

for sample in "${SAMPLES[@]}"; do
  while [[ "${running}" -ge "${MAX_JOBS}" ]]; do reap_one; done
  nohup env \
    MPLCONFIGDIR="/tmp/project26210_${sample}_mpl_auxv4r1" \
    NUMBA_CACHE_DIR="/tmp/project26210_${sample}_numba_auxv4r1" \
    OMP_NUM_THREADS=12 MKL_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12 \
    "${PYTHON}" "${SCRIPT}" --project-root "${PROJECT_ROOT}" --work-root "${WORK_ROOT}" --sample "${sample}" \
    >"${OUT}/logs/${sample}.log" 2>&1 < /dev/null &
  pid=$!
  PID_TO_SAMPLE["${pid}"]="${sample}"
  running=$((running + 1))
  printf '%s\t%s\t%s\n' "${sample}" "${pid}" "$(date --iso-8601=seconds)" >>"${PROCESS_TABLE}"
done
while [[ "${running}" -gt 0 ]]; do reap_one; done

if [[ "${overall}" -eq 0 ]]; then
  date --iso-8601=seconds >"${OUT}/ALL_SAMPLES_RUN_COMPLETE.part"
  mv "${OUT}/ALL_SAMPLES_RUN_COMPLETE.part" "${OUT}/ALL_SAMPLES_RUN_COMPLETE"
else
  date --iso-8601=seconds >"${OUT}/ALL_SAMPLES_RUN_FAILED.part"
  mv "${OUT}/ALL_SAMPLES_RUN_FAILED.part" "${OUT}/ALL_SAMPLES_RUN_FAILED"
fi
exit "${overall}"
