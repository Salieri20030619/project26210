#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 STAGE_ROOT INPUT_ROOT RSCRIPT MEMORY_GATE" >&2
  exit 2
fi
STAGE_ROOT="$1"
INPUT_ROOT="$2/infercnv"
RSCRIPT="$3"
MEMORY_GATE="$4"
UNIT_SCRIPT="${STAGE_ROOT}/scripts/run_external_infercnv_unit_v1.R"
OUTPUT_ROOT="${STAGE_ROOT}/results/full_v2/infercnv"
LOG_ROOT="${STAGE_ROOT}/logs/full_v2/infercnv"
STATUS_ROOT="${STAGE_ROOT}/results/full_v2/status"
mkdir -p "${OUTPUT_ROOT}" "${LOG_ROOT}" "${STATUS_ROOT}"
FAIL_FILE="${STATUS_ROOT}/infercnv_failed_units_v2.tsv"
if [[ -e "${FAIL_FILE}" ]]; then echo "No-overwrite gate: ${FAIL_FILE}" >&2; exit 1; fi
printf 'unit_id\tstatus\ttimestamp\n' >"${FAIL_FILE}"
completed=0
failed=0
for input_dir in "${INPUT_ROOT}"/*; do
  [[ -d "${input_dir}" ]] || continue
  unit_id="$(basename "${input_dir}")"
  output_dir="${OUTPUT_ROOT}/${unit_id}"
  log_file="${LOG_ROOT}/${unit_id}.log"
  if [[ -e "${output_dir}" ]]; then
    printf '%s\tduplicate_output_dir\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${FAIL_FILE}"
    failed=$((failed + 1))
    continue
  fi
  mkdir -p "${output_dir}"
  if "${MEMORY_GATE}" && "${RSCRIPT}" --vanilla "${UNIT_SCRIPT}" \
      --input-dir "${input_dir}" --output-dir "${output_dir}" --threads 8 >"${log_file}" 2>&1; then
    if [[ -f "${output_dir}/RUN_COMPLETE" ]]; then
      completed=$((completed + 1))
    else
      printf '%s\tmissing_RUN_COMPLETE\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${FAIL_FILE}"
      failed=$((failed + 1))
    fi
  else
    printf '%s\tfailed\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${FAIL_FILE}"
    failed=$((failed + 1))
  fi
done
printf 'completed\tfailed\n%s\t%s\n' "${completed}" "${failed}" >"${STATUS_ROOT}/infercnv_cohort_counts_v2.tsv"
if [[ "${failed}" -eq 0 && "${completed}" -eq 49 ]]; then
  date --iso-8601=seconds >"${STATUS_ROOT}/INFERCNV_RUN_COMPLETE"
  echo "EXTERNAL_INFERCNV_COHORT_V2_PASS completed=${completed}"
else
  echo "EXTERNAL_INFERCNV_COHORT_V2_FAILED completed=${completed} failed=${failed}" >&2
  exit 1
fi
