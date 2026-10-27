#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/04_subclones_tymp"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage

"${PROJECT}/check_memory_gate.sh"
cp "${STAGE}/audit_part1_subclones_cnmf_inputs_v1.py" "${DEST}/audit_part1_subclones_cnmf_inputs_v1.py"
cp "${STAGE}/run_part1_subclones_cnmf_input_audit_v1.sh" "${DEST}/run_part1_subclones_cnmf_input_audit_v1.sh"
cp "${STAGE}/verify_part1_subclones_cnmf_input_audit_v1.py" "${DEST}/verify_part1_subclones_cnmf_input_audit_v1.py"
cp "${STAGE}/run_part1_subclones_cnmf_input_audit_verification_v1.sh" "${DEST}/run_part1_subclones_cnmf_input_audit_verification_v1.sh"
cmp -s "${STAGE}/audit_part1_subclones_cnmf_inputs_v1.py" "${DEST}/audit_part1_subclones_cnmf_inputs_v1.py"
cmp -s "${STAGE}/run_part1_subclones_cnmf_input_audit_v1.sh" "${DEST}/run_part1_subclones_cnmf_input_audit_v1.sh"
cmp -s "${STAGE}/verify_part1_subclones_cnmf_input_audit_v1.py" "${DEST}/verify_part1_subclones_cnmf_input_audit_v1.py"
cmp -s "${STAGE}/run_part1_subclones_cnmf_input_audit_verification_v1.sh" "${DEST}/run_part1_subclones_cnmf_input_audit_verification_v1.sh"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile "${DEST}/audit_part1_subclones_cnmf_inputs_v1.py"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python -m py_compile "${DEST}/verify_part1_subclones_cnmf_input_audit_v1.py"
bash -n "${DEST}/run_part1_subclones_cnmf_input_audit_v1.sh"
bash -n "${DEST}/run_part1_subclones_cnmf_input_audit_verification_v1.sh"
echo PART1_SUBCLONE_CNMF_AUDIT_INSTALL_PASS
