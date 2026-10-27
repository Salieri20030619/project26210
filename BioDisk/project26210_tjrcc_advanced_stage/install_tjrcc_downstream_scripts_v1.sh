#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
STAGE_ROOT="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
RUN_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
PREFLIGHT_DIR="${RUN_ROOT}/preflight"

mkdir -p "${SCRIPT_DIR}" "${PREFLIGHT_DIR}/preflight_v3"
for file in \
  run_tjrcc_sparkx_spatial_programs_v1.R \
  run_tjrcc_hotspot_spatial_programs_v1.py \
  run_tjrcc_spatial_programs_sample_v1.sh \
  run_tjrcc_spatial_programs_cohort_v1.sh \
  prepare_tjrcc_misty_inputs_v1.py \
  run_tjrcc_misty_sample_v1.R \
  run_tjrcc_misty_cohort_v1.sh \
  run_tjrcc_commot_candidate_v1.py \
  run_tjrcc_commot_cohort_v1.sh \
  run_tjrcc_advanced_downstream_v1.sh \
  summarize_tjrcc_advanced_cohort_v1.py \
  verify_tjrcc_advanced_full_v1.py \
  run_tjrcc_cohort_synthesis_and_qc_v1.sh \
  run_tjrcc_recovery_after_misty_deldir_v1.sh \
  run_tjrcc_recovery_after_commot_missing_gene_v1.sh
do
  cp "${STAGE_ROOT}/${file}" "${SCRIPT_DIR}/${file}"
done
cp "${STAGE_ROOT}/preflight_tjrcc_advanced_v1.py" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/initial_time_and_scope_estimate_v1.json" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/preflight_v1_failure_audit.md" "${PREFLIGHT_DIR}/"
cp "${STAGE_ROOT}/preflight_v3/"* "${PREFLIGHT_DIR}/preflight_v3/"
chmod 0755 "${SCRIPT_DIR}"/*.sh
chmod 0644 "${SCRIPT_DIR}"/*.py "${SCRIPT_DIR}"/*.R
sha256sum "${SCRIPT_DIR}"/* >"${PREFLIGHT_DIR}/all_installed_script_sha256_v1.tsv"
echo "TJRCC_DOWNSTREAM_SCRIPTS_INSTALLED run_root=${RUN_ROOT}"
