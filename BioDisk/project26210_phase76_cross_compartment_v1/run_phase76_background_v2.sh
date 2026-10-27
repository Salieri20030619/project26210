#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_phase76_cross_compartment_v1"
R="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript"
PY="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
P71="/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1"
META67="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2"
META207="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz"
RUNS="${ROOT}/dataset_runs_v1"
RESULTS="${ROOT}/results_v1"
VALIDATION="${ROOT}/independent_validation_v1"
LOGS="${ROOT}/logs_v1"

mkdir -p "${LOGS}" "${RUNS}"
exec > >(tee -a "${LOGS}/background_master_v2.log") 2>&1
on_error(){ code=$?; printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "${code}" "$(date -Is)" > "${ROOT}/RUN_FAILED"; exit "${code}"; }
trap on_error ERR

available_kib=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
minimum_kib=$((50*1024*1024))
printf 'time=%s\nMemAvailable_kib=%s\nminimum_kib=%s\n' "$(date -Is)" "${available_kib}" "${minimum_kib}" > "${LOGS}/memory_gate_v2.log"
[[ "${available_kib}" -ge "${minimum_kib}" ]] || exit 70

for target in "${RUNS}/GSE207493" "${RUNS}/GSE264586" "${RUNS}/GSE294109" "${RESULTS}" "${VALIDATION}"; do
  [[ ! -e "${target}" && ! -e "${target}.part" ]] || { echo "Refusing existing output ${target}" >&2; exit 73; }
done

printf 'status=RUNNING\npid=%s\nstarted_at=%s\n' "$$" "$(date -Is)" > "${ROOT}/RUNNING"
"${R}" "${ROOT}/01_run_phase76_dataset_v2.R" GSE207493 "${P71}/inputs_v1/seurat/GSE207493_annotated_local_v2.rds" "${META207}" "${RUNS}/GSE207493"
"${R}" "${ROOT}/01_run_phase76_dataset_v2.R" GSE264586 "${P71}/inputs_v1/seurat/GSE264586_annotated_local_v2.rds" "${META67}/GSE264586_phase67_metadata_v1.tsv.gz" "${RUNS}/GSE264586"
"${R}" "${ROOT}/01_run_phase76_dataset_v2.R" GSE294109 "${P71}/inputs_v1/seurat/GSE294109_annotated_local_v2.rds" "${META67}/GSE294109_phase67_metadata_v1.tsv.gz" "${RUNS}/GSE294109"
"${PY}" "${ROOT}/02_finalize_phase76_v1.py" --runs "${RUNS}" --out "${RESULTS}"
"${PY}" "${ROOT}/03_validate_phase76_v1.py" --results "${RESULTS}" --out "${VALIDATION}"
printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${ROOT}/RUN_COMPLETE"
echo "PHASE76 BACKGROUND RUN PASS"
