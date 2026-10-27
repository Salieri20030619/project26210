#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_phase78_fine_state_spatial_v1"
MAPPING_PID="3893967"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"

printf 'status=WAITING\nmapping_pid=%s\nstarted_at=%s\n' "${MAPPING_PID}" "$(date -Is)" > "${ROOT}/PHASE78_POSTMAP_WATCHER_RUNNING_v1"
while kill -0 "${MAPPING_PID}" 2>/dev/null; do
  sleep 900
done

if [[ ! -f "${ROOT}/run_v1/PHASE78_MAPPING_PASS_v1" ]]; then
  printf 'status=FAIL\nreason=mapping_did_not_pass\nfailed_at=%s\n' "$(date -Is)" > "${ROOT}/PHASE78_POSTMAP_WATCHER_FAIL_v1"
  exit 1
fi

"${PYTHON}" "${ROOT}/scripts_v1/analyze_phase78_spatial_cross_compartment_v1.py" \
  --mapping-root "${ROOT}/run_v1" \
  --section-status "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase54_gse175540_spatial_v1/run_v1/section_status_v1.tsv" \
  --preregistration "${ROOT}/PHASE78_PREREGISTRATION_v1.md" \
  --output "${ROOT}/analysis_v1"

printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${ROOT}/PHASE78_POSTMAP_WATCHER_PASS_v1"
