#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_phase75_mta_spatial_transfer_v1"
PHASE54="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase54_gse175540_spatial_v1"
PHASE55="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase55_cell2location_v5_spawn_resume_20260821_v1"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
RESULTS="${ROOT}/results_v1"
VALIDATION="${ROOT}/independent_validation_v1"
LOGS="${ROOT}/logs_v1"

mkdir -p "${LOGS}"
exec > >(tee -a "${LOGS}/background_master_v1.log") 2>&1

on_error() {
  code=$?
  printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "${code}" "$(date -Is)" > "${ROOT}/RUN_FAILED"
  exit "${code}"
}
trap on_error ERR

available_kib=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
minimum_kib=$((50 * 1024 * 1024))
printf 'time=%s\nMemAvailable_kib=%s\nminimum_kib=%s\n' "$(date -Is)" "${available_kib}" "${minimum_kib}" > "${LOGS}/memory_gate_v1.log"
if [[ "${available_kib}" -lt "${minimum_kib}" ]]; then
  echo "Memory gate failed" >&2
  exit 70
fi

if [[ -e "${RESULTS}" || -e "${RESULTS}.part" || -e "${VALIDATION}" || -e "${VALIDATION}.part" ]]; then
  echo "Refusing existing Phase75 outputs" >&2
  exit 73
fi

printf 'status=RUNNING\npid=%s\nstarted_at=%s\n' "$$" "$(date -Is)" > "${ROOT}/RUNNING"
"${PYTHON}" "${ROOT}/01_run_phase75_spatial_transfer_v1.py" --phase54 "${PHASE54}" --phase55 "${PHASE55}" --out "${RESULTS}"
"${PYTHON}" "${ROOT}/02_validate_phase75_v1.py" --results "${RESULTS}" --out "${VALIDATION}"
printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${ROOT}/RUN_COMPLETE"
echo "PHASE75 BACKGROUND RUN PASS"
