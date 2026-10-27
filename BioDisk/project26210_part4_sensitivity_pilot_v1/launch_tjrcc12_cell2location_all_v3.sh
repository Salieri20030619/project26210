#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1"
OUT="${ROOT}/tjrcc12_cell2location_reference_v3"
LAUNCHER="${ROOT}/run_tjrcc12_cell2location_sample_v3.sh"
SAMPLES=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)

mkdir -p "${OUT}/logs"
printf 'sample\tpid\tstarted_at\n' > "${OUT}/master_processes_v3.tsv"
pids=()
samples_started=()

for sample in "${SAMPLES[@]}"; do
  "${LAUNCHER}" "${sample}" &
  pid=$!
  pids+=("${pid}")
  samples_started+=("${sample}")
  printf '%s\t%s\t%s\n' "${sample}" "${pid}" "$(date --iso-8601=seconds)" \
    >> "${OUT}/master_processes_v3.tsv"
  sleep 10
done

printf 'sample\tpid\texit_code\tfinished_at\n' > "${OUT}/master_completion_v3.tsv"
overall=0
for index in "${!pids[@]}"; do
  sample="${samples_started[$index]}"
  pid="${pids[$index]}"
  set +e
  wait "${pid}"
  rc=$?
  set -e
  printf '%s\t%s\t%s\t%s\n' "${sample}" "${pid}" "${rc}" "$(date --iso-8601=seconds)" \
    >> "${OUT}/master_completion_v3.tsv"
  if [[ "${rc}" -ne 0 ]]; then overall=1; fi
done
exit "${overall}"
