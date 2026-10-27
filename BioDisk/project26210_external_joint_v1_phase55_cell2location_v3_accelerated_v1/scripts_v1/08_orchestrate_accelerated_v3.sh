#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 7 ]]; then
  echo "usage: $0 V2_ROOT V3_ROOT REFERENCE_H5AD SPATIAL_SECTIONS PYTHON REQUESTED_WORKERS MIN_AVAILABLE_GIB" >&2
  exit 64
fi

V2_ROOT="$1"
V3_ROOT="$2"
REFERENCE_H5AD="$3"
SPATIAL_SECTIONS="$4"
PYTHON_BIN="$5"
REQUESTED_WORKERS="$6"
MIN_AVAILABLE_GIB="$7"
V2_REFERENCE="${V2_ROOT}/run_v1/reference_model_v1"
V3_REFERENCE="${V3_ROOT}/run_v1/reference_model_v1"
SCRIPT_DIR="${V3_ROOT}/scripts_v1"
LOG_PATH="${V3_ROOT}/logs_v1/acceleration_orchestrator_v3.log"

mkdir -p "${V3_ROOT}/logs_v1"
exec >> "${LOG_PATH}" 2>&1

fail() {
  local message="$1"
  printf '[%s] FAIL: %s\n' "$(date -Is)" "${message}"
  printf 'status=FAIL\nfailed_at=%s\nmessage=%s\n' "$(date -Is)" "${message}" > "${V3_ROOT}/PHASE55_ACCELERATION_FAIL_v3"
  exit 2
}

for value in "${REQUESTED_WORKERS}" "${MIN_AVAILABLE_GIB}"; do
  [[ "${value}" =~ ^[0-9]+$ ]] || fail "Worker and memory arguments must be positive integers"
done
(( REQUESTED_WORKERS >= 1 )) || fail "Requested workers must be positive"
(( MIN_AVAILABLE_GIB >= 1 )) || fail "Minimum available GiB must be positive"

for required in \
  "${REFERENCE_H5AD}" \
  "${SCRIPT_DIR}/02_run_cell2location_v1.py" \
  "${SCRIPT_DIR}/03_launch_cell2location_background_v1.sh" \
  "${SCRIPT_DIR}/05_validate_cell2location_results_v2.py" \
  "${SCRIPT_DIR}/06_wait_and_validate_cell2location_v2.sh"; do
  [[ -f "${required}" ]] || fail "Missing required file: ${required}"
done

section_count="$(find "${SPATIAL_SECTIONS}" -mindepth 1 -maxdepth 1 -type d -name 'GSM*' | wc -l)"
[[ "${section_count}" -eq 24 ]] || fail "Expected 24 spatial sections, found ${section_count}"

printf 'status=WAITING_FOR_REFERENCE\npid=%s\nstarted_at=%s\nrequested_workers=%s\n' \
  "$$" "$(date -Is)" "${REQUESTED_WORKERS}" > "${V3_ROOT}/PHASE55_ACCELERATION_WAITING_v3"
printf '[%s] waiting for v2 reference checkpoint\n' "$(date -Is)"

while [[ ! -f "${V2_REFERENCE}/REFERENCE_MODEL_PASS_v1" ]]; do
  if [[ -f "${V2_ROOT}/PHASE55_BACKGROUND_FAIL_v1" ]]; then
    fail "V2 analysis failed before the reference checkpoint completed"
  fi
  sleep 10
done

printf '[%s] v2 reference checkpoint found\n' "$(date -Is)"
for required in \
  "${V2_REFERENCE}/reference_cell_state_signatures_v1.tsv.gz" \
  "${V2_REFERENCE}/reference_metadata_v1.json" \
  "${V2_REFERENCE}/reference_posterior_v1.h5ad"; do
  [[ -f "${required}" ]] || fail "Incomplete v2 reference checkpoint: ${required}"
done

available_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
minimum_kib="$(( MIN_AVAILABLE_GIB * 1024 * 1024 ))"
effective_workers="${REQUESTED_WORKERS}"
if (( available_kib < minimum_kib )); then
  if (( available_kib >= 150 * 1024 * 1024 )); then
    effective_workers=12
  elif (( available_kib >= 100 * 1024 * 1024 )); then
    effective_workers=8
  else
    fail "Available memory below 100 GiB; accelerated spatial launch withheld"
  fi
fi
printf '[%s] available_kib=%s effective_workers=%s\n' "$(date -Is)" "${available_kib}" "${effective_workers}"

mkdir -p "${V3_ROOT}/run_v1"
if [[ -e "${V3_REFERENCE}" || -L "${V3_REFERENCE}" ]]; then
  [[ "$(readlink -f "${V3_REFERENCE}")" == "$(readlink -f "${V2_REFERENCE}")" ]] \
    || fail "Existing v3 reference path does not resolve to the v2 checkpoint"
else
  ln -s "${V2_REFERENCE}" "${V3_REFERENCE}"
fi
printf 'status=PASS\nsource=%s\nlinked_at=%s\n' \
  "${V2_REFERENCE}" "$(date -Is)" > "${V3_ROOT}/REFERENCE_REUSE_PASS_v3"

setsid -f bash "${SCRIPT_DIR}/03_launch_cell2location_background_v1.sh" \
  "${REFERENCE_H5AD}" \
  "${SPATIAL_SECTIONS}" \
  "${V3_ROOT}" \
  "${PYTHON_BIN}" \
  "${effective_workers}" \
  > "${V3_ROOT}/bootstrap_v3.log" 2>&1 < /dev/null

for _ in $(seq 1 30); do
  [[ -f "${V3_ROOT}/PHASE55_BACKGROUND_RUNNING_v1" ]] && break
  [[ -f "${V3_ROOT}/PHASE55_BACKGROUND_FAIL_v1" ]] && fail "V3 analysis failed during launch"
  sleep 1
done
[[ -f "${V3_ROOT}/PHASE55_BACKGROUND_RUNNING_v1" ]] || fail "V3 analysis did not create a RUNNING marker"

setsid -f bash "${SCRIPT_DIR}/06_wait_and_validate_cell2location_v2.sh" \
  "${V3_ROOT}" \
  "${PYTHON_BIN}" \
  "${SCRIPT_DIR}/05_validate_cell2location_results_v2.py" \
  > "${V3_ROOT}/logs_v1/watcher_bootstrap_v3.log" 2>&1 < /dev/null

for _ in $(seq 1 30); do
  [[ -f "${V3_ROOT}/PHASE55_VALIDATION_WATCHER_RUNNING_v2" ]] && break
  sleep 1
done
[[ -f "${V3_ROOT}/PHASE55_VALIDATION_WATCHER_RUNNING_v2" ]] || fail "V3 validation watcher did not start"

printf 'status=PASS\nlaunched_at=%s\nrequested_workers=%s\neffective_workers=%s\nreference_source=%s\n' \
  "$(date -Is)" "${REQUESTED_WORKERS}" "${effective_workers}" "${V2_REFERENCE}" \
  > "${V3_ROOT}/PHASE55_ACCELERATION_LAUNCHED_v3"

v2_launcher_pid="$(awk -F= '$1 == "pid" {print $2}' "${V2_ROOT}/PHASE55_BACKGROUND_RUNNING_v1")"
[[ "${v2_launcher_pid}" =~ ^[0-9]+$ ]] || fail "Could not resolve the v2 launcher PID"
v2_command="$(ps -p "${v2_launcher_pid}" -o args= || true)"
[[ "${v2_command}" == *"${V2_ROOT}"* ]] || fail "V2 launcher PID command does not match the expected root"
v2_pgid="$(ps -p "${v2_launcher_pid}" -o pgid= | tr -d ' ')"
[[ "${v2_pgid}" =~ ^[0-9]+$ ]] || fail "Could not resolve the v2 process group"
[[ "${v2_pgid}" -ne "$$" ]] || fail "Refusing to terminate the orchestrator process group"

printf 'status=SUPERSEDED_BY_ACCELERATED_V3\nsuperseded_at=%s\nv3_root=%s\nreference_preserved=%s\n' \
  "$(date -Is)" "${V3_ROOT}" "${V2_REFERENCE}" > "${V2_ROOT}/PHASE55_ACCELERATION_SUPERSEDED_v2"
printf 'status=SUPERSEDED_BY_ACCELERATED_V3\nexit_code=0\nfailed_at=%s\n' \
  "$(date -Is)" > "${V2_ROOT}/PHASE55_BACKGROUND_FAIL_v1"

printf '[%s] terminating superseded v2 process group %s\n' "$(date -Is)" "${v2_pgid}"
kill -TERM -- "-${v2_pgid}" || true
sleep 10
if ps -g "${v2_pgid}" -o pid= | grep -q '[0-9]'; then
  printf '[%s] v2 process group still alive after TERM; sending KILL\n' "$(date -Is)"
  kill -KILL -- "-${v2_pgid}" || true
fi

printf '[%s] accelerated v3 launch complete\n' "$(date -Is)"
printf 'status=PASS\ncompleted_at=%s\nv3_root=%s\neffective_workers=%s\n' \
  "$(date -Is)" "${V3_ROOT}" "${effective_workers}" > "${V3_ROOT}/PHASE55_ACCELERATION_ORCHESTRATOR_PASS_v3"

