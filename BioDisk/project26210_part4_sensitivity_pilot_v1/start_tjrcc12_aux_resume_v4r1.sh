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
LAUNCHER="${WORK_ROOT}/launch_tjrcc12_aux_resume_v4r1.sh"
OUT="${WORK_ROOT}/tjrcc12_cell2location_reference_v4r1_aux_resume"

if [[ -e "${OUT}" ]]; then
  echo "No-overwrite gate: ${OUT} already exists" >&2
  exit 1
fi
"${PYTHON}" "${SCRIPT}" --project-root "${PROJECT_ROOT}" --work-root "${WORK_ROOT}" --preflight-only
PID_FILE="${OUT}/controlled_launcher_v4r1.pid"
nohup "${LAUNCHER}" "${PROJECT_ROOT}" "${WORK_ROOT}" "${PYTHON}" \
  >"${OUT}/logs/controlled_launcher_v4r1.log" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "${pid}" >"${PID_FILE}.part"
mv "${PID_FILE}.part" "${PID_FILE}"
echo "TJRCC12_AUX_V4R1_LAUNCHED pid=${pid} max_concurrent=4"
