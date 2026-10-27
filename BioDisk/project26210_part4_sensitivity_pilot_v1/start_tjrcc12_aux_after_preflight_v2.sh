#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 PROJECT_ROOT WORK_ROOT CELL2LOC_PYTHON" >&2
  exit 2
fi
PROJECT_ROOT="$1"
WORK_ROOT="$2"
PYTHON="$3"
LAUNCHER="${WORK_ROOT}/launch_tjrcc12_aux_resume_v4r1.sh"
OUT="${WORK_ROOT}/tjrcc12_cell2location_reference_v4r1_aux_resume"
PID_FILE="${OUT}/controlled_launcher_recovery_v2.pid"
LOG_FILE="${OUT}/logs/controlled_launcher_recovery_v2.log"

if [[ ! -f "${OUT}/run_metadata_v4.json" || ! -f "${OUT}/preflight_tjrcc12_v4r1.tsv" ]]; then
  echo "Validated preflight outputs are missing" >&2
  exit 1
fi
if [[ -e "${PID_FILE}" || -e "${LOG_FILE}" || -e "${OUT}/master_processes_v4r1.tsv" ]]; then
  echo "No-overwrite/duplicate-launch gate blocked recovery start" >&2
  exit 1
fi
if find "${OUT}" -mindepth 1 -maxdepth 1 -type d ! -name logs | grep -q .; then
  echo "Unexpected sample/output directory exists" >&2
  exit 1
fi
mkdir -p "${OUT}/logs"
nohup "${LAUNCHER}" "${PROJECT_ROOT}" "${WORK_ROOT}" "${PYTHON}" \
  >"${LOG_FILE}" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "${pid}" >"${PID_FILE}.part"
mv "${PID_FILE}.part" "${PID_FILE}"
echo "TJRCC12_AUX_RECOVERY_V2_LAUNCHED pid=${pid} max_concurrent=4"
