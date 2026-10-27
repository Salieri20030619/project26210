#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage"
DEST="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
LOG_DEST="${PROJECT}/results/part1_rebuild_v1/logs"

mkdir -p "${DEST}" "${LOG_DEST}"

for script in \
  prepare_part1_cnv_unit_inputs_v1.py \
  prepare_part1_cnv_unit_inputs_v2.py \
  run_part1_infercnv_unit_v1.R \
  run_part1_copykat_unit_v1.R \
  run_part1_cnv_pilot_v1.sh \
  run_part1_cnv_pilot_v2.sh \
  verify_part1_cnv_pilot_v2.py \
  verify_part1_cnv_pilot_v2_verification_v2.py \
  run_part1_cnv_pilot_v2_verification_v1.sh \
  run_part1_cnv_pilot_v2_verification_v2.sh \
  run_part1_infercnv_cohort_v1.sh \
  run_part1_copykat_cohort_v1.sh \
  run_part1_malignancy_cnv_full_v1.sh \
  launch_part1_malignancy_cnv_full_v1.sh
do
  mode=0644
  case "${script}" in
    *.sh) mode=0755 ;;
  esac
  install -m "${mode}" "${STAGE}/${script}" "${DEST}/${script}"
done

cp -a "${STAGE}/cnv_pilot_v1" "${DEST}/"
cp -a "${STAGE}/cnv_pilot_v2" "${DEST}/"

test -f "${DEST}/cnv_pilot_v2/infercnv/GSE156632__RCC1/RUN_COMPLETE"
test -f "${DEST}/cnv_pilot_v2/copykat/RCC1t/RUN_COMPLETE"
test -f "${DEST}/cnv_pilot_v2/verification_v1/RUN_FAILED"
test -f "${DEST}/cnv_pilot_v2/verification_v2/RUN_COMPLETE"
grep -q 'PART1_CNV_PILOT_V2_ALL_PASS' "${DEST}/cnv_pilot_v2/cnv_pilot_v2.log"
grep -q 'PART1_CNV_PILOT_V2_VERIFICATION_V2_PASS' "${DEST}/cnv_pilot_v2/verification_v2.log"

echo "PART1_MALIGNANCY_CNV_PIPELINE_INSTALL_PASS"
