#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/04_subclones_tymp"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
FILES=(
    verify_part1_formal_cnv_subclones_v1.py
    run_part1_formal_cnv_subclones_verification_v1.sh
    run_part1_cnv_clone_tymp_cnmf_association_v1.py
    run_part1_cnv_clone_tymp_cnmf_association_v1.sh
    verify_part1_cnv_clone_tymp_cnmf_association_v1.py
    run_part1_cnv_clone_tymp_cnmf_association_verification_v1.sh
    run_part1_cnv_clone_overnight_v1.sh
    launch_part1_cnv_clone_overnight_v1.sh
)

"${PROJECT}/check_memory_gate.sh"
for file in "${FILES[@]}"; do
    cp "${STAGE}/${file}" "${DEST}/${file}"
    cmp -s "${STAGE}/${file}" "${DEST}/${file}"
done
"${PYTHON}" -m py_compile \
    "${DEST}/verify_part1_formal_cnv_subclones_v1.py" \
    "${DEST}/run_part1_cnv_clone_tymp_cnmf_association_v1.py" \
    "${DEST}/verify_part1_cnv_clone_tymp_cnmf_association_v1.py"
bash -n "${DEST}/run_part1_formal_cnv_subclones_verification_v1.sh"
bash -n "${DEST}/run_part1_cnv_clone_tymp_cnmf_association_v1.sh"
bash -n "${DEST}/run_part1_cnv_clone_tymp_cnmf_association_verification_v1.sh"
bash -n "${DEST}/run_part1_cnv_clone_overnight_v1.sh"
bash -n "${DEST}/launch_part1_cnv_clone_overnight_v1.sh"
echo PART1_CNV_CLONE_OVERNIGHT_INSTALL_PASS
