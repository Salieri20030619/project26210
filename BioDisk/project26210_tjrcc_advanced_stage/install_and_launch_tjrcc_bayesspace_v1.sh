#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
STAGE_ROOT="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
RUN_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
PREFLIGHT_DIR="${RUN_ROOT}/preflight"
LOG_DIR="${RUN_ROOT}/logs"

mkdir -p "${SCRIPT_DIR}" "${PREFLIGHT_DIR}/preflight_v3" "${LOG_DIR}"
cp "${STAGE_ROOT}/run_tjrcc_bayesspace_sample_v1.R" "${SCRIPT_DIR}/"
cp "${STAGE_ROOT}/postprocess_tjrcc_bayesspace_boundary_v1.py" "${SCRIPT_DIR}/"
cp "${STAGE_ROOT}/run_tjrcc_bayesspace_sample_v1.sh" "${SCRIPT_DIR}/"
cp "${STAGE_ROOT}/run_tjrcc_bayesspace_cohort_v1.sh" "${SCRIPT_DIR}/"
cp "${STAGE_ROOT}/preflight_tjrcc_advanced_v1.py" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/initial_time_and_scope_estimate_v1.json" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/preflight_v1_failure_audit.md" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/bayesspace_v1_failure_audit.md" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/preflight_v3/"* "${PREFLIGHT_DIR}/preflight_v3/"
chmod 0755 "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.sh" "${SCRIPT_DIR}/run_tjrcc_bayesspace_cohort_v1.sh"
chmod 0644 "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.R" "${SCRIPT_DIR}/postprocess_tjrcc_bayesspace_boundary_v1.py"

sha256sum \
  "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.R" \
  "${SCRIPT_DIR}/postprocess_tjrcc_bayesspace_boundary_v1.py" \
  "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.sh" \
  "${SCRIPT_DIR}/run_tjrcc_bayesspace_cohort_v1.sh" \
  >"${PREFLIGHT_DIR}/bayesspace_script_sha256_v1.tsv"

"${PROJECT_ROOT}/check_memory_gate.sh" | tee "${LOG_DIR}/launch_memory_gate_v1.log"
echo "TJRCC_BAYESSPACE_V2_INSTALLED run_root=${RUN_ROOT}"
