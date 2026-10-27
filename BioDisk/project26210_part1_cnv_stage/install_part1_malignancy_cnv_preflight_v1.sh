#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage"
DEST="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
LOG_DEST="${PROJECT}/results/part1_rebuild_v1/logs"

mkdir -p "${DEST}" "${LOG_DEST}"

for script in \
  audit_part1_malignancy_cnv_inputs_v1.py \
  audit_part1_malignancy_cnv_inputs_v2.py \
  run_part1_malignancy_cnv_input_audit_v1.sh \
  run_part1_malignancy_cnv_input_audit_v2.sh
do
  install -m 0644 "${STAGE}/${script}" "${DEST}/${script}"
done

cp -a "${STAGE}/malignancy_cnv_input_audit_v1" "${DEST}/"
cp -a "${STAGE}/malignancy_cnv_input_audit_v2" "${DEST}/"
install -m 0644 \
  "${STAGE}/malignancy_cnv_input_audit_v1.log" \
  "${LOG_DEST}/malignancy_cnv_input_audit_v1.log"
install -m 0644 \
  "${STAGE}/malignancy_cnv_input_audit_v2.log" \
  "${LOG_DEST}/malignancy_cnv_input_audit_v2.log"

test -f "${DEST}/malignancy_cnv_input_audit_v1/RUN_FAILED"
test -f "${DEST}/malignancy_cnv_input_audit_v2/RUN_COMPLETE"
grep -q 'PART1_MALIGNANCY_CNV_INPUT_AUDIT_V1_FAIL' \
  "${LOG_DEST}/malignancy_cnv_input_audit_v1.log"
grep -q 'PART1_MALIGNANCY_CNV_INPUT_AUDIT_V2_PASS' \
  "${LOG_DEST}/malignancy_cnv_input_audit_v2.log"

echo "PART1_MALIGNANCY_CNV_PREFLIGHT_INSTALL_PASS"
