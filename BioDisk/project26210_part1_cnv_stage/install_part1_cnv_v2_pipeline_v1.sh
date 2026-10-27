#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
FILES=(
  run_part1_infercnv_unit_v2.R
  run_part1_copykat_unit_v2.R
  run_part1_infercnv_cohort_v2.sh
  run_part1_copykat_cohort_v2.sh
  run_part1_cnv_v2_pilot.sh
  verify_part1_cnv_v2_pilot.py
  run_part1_malignancy_cnv_full_v2.sh
  run_part1_cnv_v2_pilot_then_full_v1.sh
  launch_part1_cnv_v2_pilot_then_full_v1.sh
  recover_part1_cnv_v2_terminal_v1.py
  verify_part1_cnv_v2_terminal_v1.py
  run_part1_cnv_v2_terminal_recovery_v1.sh
  summarize_part1_malignancy_cnv_v1.py
  run_part1_malignancy_cnv_summary_v1.sh
  verify_part1_malignancy_cnv_summary_v1.py
  run_part1_malignancy_cnv_summary_verification_v1.sh
  prepare_part1_mixed_cnv_sensitivity_inputs_v1.py
  run_part1_mixed_cnv_sensitivity_v1.sh
  summarize_part1_mixed_cnv_sensitivity_v1.py
  run_part1_mixed_cnv_sensitivity_summary_v1.sh
  build_part1_cnv_annotated_h5ad_v1.py
  run_part1_cnv_annotated_h5ad_v1.sh
  verify_part1_cnv_annotated_h5ad_v1.py
  run_part1_cnv_annotated_h5ad_verification_v1.sh
  verify_part1_mixed_cnv_sensitivity_v1.py
  run_part1_mixed_cnv_sensitivity_verification_v1.sh
  run_part1_mixed_cnv_infercnv_control_v1.sh
  summarize_part1_mixed_cnv_control_v1.py
  run_part1_mixed_cnv_control_summary_v2.sh
)

"${PROJECT}/check_memory_gate.sh"
for name in "${FILES[@]}"
do
  test -s "${STAGE}/${name}"
  cp "${STAGE}/${name}" "${DEST}/${name}"
  cmp -s "${STAGE}/${name}" "${DEST}/${name}"
done

/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript -e \
  "parse(file='${DEST}/run_part1_infercnv_unit_v2.R'); parse(file='${DEST}/run_part1_copykat_unit_v2.R'); cat('R_PARSE_PASS\\n')"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile \
  "${DEST}/verify_part1_cnv_v2_pilot.py" \
  "${DEST}/recover_part1_cnv_v2_terminal_v1.py" \
  "${DEST}/verify_part1_cnv_v2_terminal_v1.py" \
  "${DEST}/summarize_part1_malignancy_cnv_v1.py" \
  "${DEST}/verify_part1_malignancy_cnv_summary_v1.py" \
  "${DEST}/prepare_part1_mixed_cnv_sensitivity_inputs_v1.py" \
  "${DEST}/summarize_part1_mixed_cnv_sensitivity_v1.py" \
  "${DEST}/build_part1_cnv_annotated_h5ad_v1.py" \
  "${DEST}/verify_part1_cnv_annotated_h5ad_v1.py" \
  "${DEST}/verify_part1_mixed_cnv_sensitivity_v1.py" \
  "${DEST}/summarize_part1_mixed_cnv_control_v1.py"
bash -n \
  "${DEST}/run_part1_infercnv_cohort_v2.sh" \
  "${DEST}/run_part1_copykat_cohort_v2.sh" \
  "${DEST}/run_part1_cnv_v2_pilot.sh" \
  "${DEST}/run_part1_malignancy_cnv_full_v2.sh" \
  "${DEST}/run_part1_cnv_v2_pilot_then_full_v1.sh" \
  "${DEST}/launch_part1_cnv_v2_pilot_then_full_v1.sh" \
  "${DEST}/run_part1_cnv_v2_terminal_recovery_v1.sh" \
  "${DEST}/run_part1_malignancy_cnv_summary_v1.sh" \
  "${DEST}/run_part1_malignancy_cnv_summary_verification_v1.sh" \
  "${DEST}/run_part1_mixed_cnv_sensitivity_v1.sh" \
  "${DEST}/run_part1_mixed_cnv_sensitivity_summary_v1.sh" \
  "${DEST}/run_part1_cnv_annotated_h5ad_v1.sh" \
  "${DEST}/run_part1_cnv_annotated_h5ad_verification_v1.sh" \
  "${DEST}/run_part1_mixed_cnv_sensitivity_verification_v1.sh" \
  "${DEST}/run_part1_mixed_cnv_infercnv_control_v1.sh" \
  "${DEST}/run_part1_mixed_cnv_control_summary_v2.sh"

echo PART1_CNV_V2_PIPELINE_INSTALL_PASS
