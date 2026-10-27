#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/04_subclones_tymp"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh"
cp "${STAGE}/audit_part1_cnv_subclone_profiles_v1.py" "${DEST}/audit_part1_cnv_subclone_profiles_v1.py"
cp "${STAGE}/run_part1_cnv_subclone_preflight_v1.sh" "${DEST}/run_part1_cnv_subclone_preflight_v1.sh"
cmp -s "${STAGE}/audit_part1_cnv_subclone_profiles_v1.py" "${DEST}/audit_part1_cnv_subclone_profiles_v1.py"
cmp -s "${STAGE}/run_part1_cnv_subclone_preflight_v1.sh" "${DEST}/run_part1_cnv_subclone_preflight_v1.sh"
"${PYTHON}" -m py_compile "${DEST}/audit_part1_cnv_subclone_profiles_v1.py"
bash -n "${DEST}/run_part1_cnv_subclone_preflight_v1.sh"
echo PART1_CNV_SUBCLONE_PREFLIGHT_INSTALL_PASS
