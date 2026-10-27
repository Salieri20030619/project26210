#!/usr/bin/env bash
set -euo pipefail

WORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT26210_ROOT:?Set PROJECT26210_ROOT to the server project root}"
PYTHON="${CELL2LOC_PYTHON:?Set CELL2LOC_PYTHON to the cell2location Python executable}"
SCRIPT="${WORK_ROOT}/run_glmf12_cell2location_resume_v4.py"
OUT="${WORK_ROOT}/glmf12_cell2location_reference_v4_resume"

mkdir -p "${OUT}/logs"
exec 9>"${OUT}/serial_launcher_v4.lock"
if ! flock -n 9; then
  printf 'Another v4 serial launcher already holds %s\n' "${OUT}/serial_launcher_v4.lock" >&2
  exit 3
fi

"${PYTHON}" "${SCRIPT}" \
  --project-root "${PROJECT_ROOT}" \
  --work-root "${WORK_ROOT}" \
  --preflight-only \
  >"${OUT}/logs/preflight_v4.log" 2>&1

for sample in GLMF1 GLMF2; do
  printf '%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${sample}" "start" \
    >>"${OUT}/serial_launcher_status_v4.tsv"
  set +e
  "${PYTHON}" "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --work-root "${WORK_ROOT}" \
    --sample "${sample}" \
    >"${OUT}/logs/${sample}.log" 2>&1
  rc=$?
  set -e
  printf '%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${sample}" "finish" "${rc}" \
    >>"${OUT}/serial_launcher_status_v4.tsv"
  if [[ "${rc}" -ne 0 ]]; then
    exit "${rc}"
  fi
done
