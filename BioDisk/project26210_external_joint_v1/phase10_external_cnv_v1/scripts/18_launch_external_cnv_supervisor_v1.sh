#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
  echo "Usage: $0 STAGE_ROOT INPUT_ROOT PREFLIGHT_DIR BUILDER_PID PYTHON RSCRIPT MEMORY_GATE" >&2
  exit 2
fi
STAGE_ROOT="$1"
PID_FILE="${STAGE_ROOT}/logs/full_cohort_supervisor_v1.pid"
if [[ -e "${PID_FILE}" ]]; then echo "Duplicate supervisor launch blocked" >&2; exit 1; fi
nohup "${STAGE_ROOT}/scripts/17_supervise_and_launch_external_cnv_cohorts_v1.sh" \
  "$@" >"${STAGE_ROOT}/logs/full_cohort_supervisor_launcher_v1.log" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "${pid}" >"${PID_FILE}.part"
mv "${PID_FILE}.part" "${PID_FILE}"
echo "EXTERNAL_CNV_SUPERVISOR_LAUNCHED pid=${pid}"
