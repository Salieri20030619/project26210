#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 SERVER_ROOT PYTHON VALIDATOR" >&2
  exit 64
fi

SERVER_ROOT="$1"
PYTHON_BIN="$2"
VALIDATOR="$3"
LOG_PATH="${SERVER_ROOT}/logs_v1/postrun_validation_watcher_v2.log"
mkdir -p "$(dirname "${LOG_PATH}")"
exec >> "${LOG_PATH}" 2>&1

printf '[%s] watcher started pid=%s\n' "$(date -Is)" "$$"
printf 'status=WAITING\npid=%s\nstarted_at=%s\n' "$$" "$(date -Is)" > "${SERVER_ROOT}/PHASE55_VALIDATION_WATCHER_RUNNING_v2"

while true; do
  if [[ -f "${SERVER_ROOT}/PHASE55_BACKGROUND_FAIL_v1" ]]; then
    printf '[%s] analysis failed; validation not run\n' "$(date -Is)"
    printf 'status=ANALYSIS_FAIL\ncompleted_at=%s\n' "$(date -Is)" > "${SERVER_ROOT}/PHASE55_VALIDATION_WATCHER_FAIL_v2"
    exit 2
  fi
  if [[ -f "${SERVER_ROOT}/PHASE55_BACKGROUND_PASS_v1" ]]; then
    printf '[%s] analysis pass marker found; starting validation\n' "$(date -Is)"
    if "${PYTHON_BIN}" "${VALIDATOR}" --server-root "${SERVER_ROOT}"; then
      printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${SERVER_ROOT}/PHASE55_VALIDATION_WATCHER_PASS_v2"
      printf '[%s] validation passed\n' "$(date -Is)"
      exit 0
    else
      code=$?
      printf 'status=VALIDATION_FAIL\nexit_code=%s\ncompleted_at=%s\n' "${code}" "$(date -Is)" > "${SERVER_ROOT}/PHASE55_VALIDATION_WATCHER_FAIL_v2"
      printf '[%s] validation failed with exit code %s\n' "$(date -Is)" "${code}"
      exit "${code}"
    fi
  fi
  sleep 300
done
