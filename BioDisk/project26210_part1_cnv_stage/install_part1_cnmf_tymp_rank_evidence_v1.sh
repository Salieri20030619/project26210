#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/05_cnmf"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
RSCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript

"${PROJECT}/check_memory_gate.sh"
cp "${STAGE}/analyze_part1_cnmf_tymp_rank_evidence_v1.R" "${DEST}/analyze_part1_cnmf_tymp_rank_evidence_v1.R"
cp "${STAGE}/run_part1_cnmf_tymp_rank_evidence_v1.sh" "${DEST}/run_part1_cnmf_tymp_rank_evidence_v1.sh"
cp "${STAGE}/verify_part1_cnmf_tymp_rank_evidence_v1.R" "${DEST}/verify_part1_cnmf_tymp_rank_evidence_v1.R"
cp "${STAGE}/run_part1_cnmf_tymp_rank_evidence_verification_v1.sh" "${DEST}/run_part1_cnmf_tymp_rank_evidence_verification_v1.sh"
cmp -s "${STAGE}/analyze_part1_cnmf_tymp_rank_evidence_v1.R" "${DEST}/analyze_part1_cnmf_tymp_rank_evidence_v1.R"
cmp -s "${STAGE}/run_part1_cnmf_tymp_rank_evidence_v1.sh" "${DEST}/run_part1_cnmf_tymp_rank_evidence_v1.sh"
cmp -s "${STAGE}/verify_part1_cnmf_tymp_rank_evidence_v1.R" "${DEST}/verify_part1_cnmf_tymp_rank_evidence_v1.R"
cmp -s "${STAGE}/run_part1_cnmf_tymp_rank_evidence_verification_v1.sh" "${DEST}/run_part1_cnmf_tymp_rank_evidence_verification_v1.sh"
bash -n "${DEST}/run_part1_cnmf_tymp_rank_evidence_v1.sh"
bash -n "${DEST}/run_part1_cnmf_tymp_rank_evidence_verification_v1.sh"
"${RSCRIPT}" -e "parse(file='${DEST}/analyze_part1_cnmf_tymp_rank_evidence_v1.R'); cat('R_PARSE_PASS\\n')" >/dev/null
"${RSCRIPT}" -e "parse(file='${DEST}/verify_part1_cnmf_tymp_rank_evidence_v1.R'); cat('R_VERIFY_PARSE_PASS\\n')" >/dev/null
echo PART1_CNMF_TYMP_RANK_EVIDENCE_INSTALL_PASS
