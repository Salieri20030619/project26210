#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_completion_stage_v1
PYSCENIC_PY="${PROJECT}/.venv-pyscenic/bin/python"
R_SCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript

"${PROJECT}/check_memory_gate.sh" part1_tf_terminal_verification_install_v1
test -x "${PYSCENIC_PY}"
test -x "${R_SCRIPT}"
"${PYSCENIC_PY}" --version
"${R_SCRIPT}" --version
cp "${STAGE}/verify_part1_pyscenic_terminal_v1.py" "${ROOT}/07_pyscenic_v1/verify_part1_pyscenic_terminal_v1.py"
cp "${STAGE}/verify_part1_metaviper_terminal_v1.R" "${ROOT}/08_metaviper_v1/verify_part1_metaviper_terminal_v1.R"
cmp -s "${STAGE}/verify_part1_pyscenic_terminal_v1.py" "${ROOT}/07_pyscenic_v1/verify_part1_pyscenic_terminal_v1.py"
cmp -s "${STAGE}/verify_part1_metaviper_terminal_v1.R" "${ROOT}/08_metaviper_v1/verify_part1_metaviper_terminal_v1.R"
"${PYSCENIC_PY}" -m py_compile "${ROOT}/07_pyscenic_v1/verify_part1_pyscenic_terminal_v1.py"
"${PROJECT}/check_memory_gate.sh" part1_pyscenic_terminal_verification_run_v1
"${PYSCENIC_PY}" "${ROOT}/07_pyscenic_v1/verify_part1_pyscenic_terminal_v1.py" \
  > "${ROOT}/07_pyscenic_v1/terminal_verification_v1.log" 2>&1
grep -F 'PART1_PYSCENIC_TERMINAL_VERIFICATION_V1_PASS' "${ROOT}/07_pyscenic_v1/terminal_verification_v1.log"
"${PROJECT}/check_memory_gate.sh" part1_metaviper_terminal_verification_run_v1
"${R_SCRIPT}" "${ROOT}/08_metaviper_v1/verify_part1_metaviper_terminal_v1.R" \
  > "${ROOT}/08_metaviper_v1/terminal_verification_v1.log" 2>&1
grep -F 'PART1_METAVIPER_TERMINAL_VERIFICATION_V1_PASS' "${ROOT}/08_metaviper_v1/terminal_verification_v1.log"
echo PART1_TF_TERMINAL_VERIFICATION_V1_ALL_PASS
