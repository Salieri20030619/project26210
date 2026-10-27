#!/usr/bin/env bash
set -uo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 STAGE_ROOT INPUT_ROOT RSCRIPT MEMORY_GATE" >&2
  exit 2
fi
STAGE_ROOT="$1"
INPUT_ROOT="$2/copykat"
RSCRIPT="$3"
MEMORY_GATE="$4"
UNIT_SCRIPT="${STAGE_ROOT}/scripts/run_external_copykat_unit_v1.R"
OUTPUT_ROOT="${STAGE_ROOT}/results/full_v1/copykat"
LOG_ROOT="${STAGE_ROOT}/logs/full_v1/copykat"
STATUS_ROOT="${STAGE_ROOT}/results/full_v1/status"
mkdir -p "${OUTPUT_ROOT}" "${LOG_ROOT}" "${STATUS_ROOT}"
FAIL_FILE="${STATUS_ROOT}/copykat_failed_units_v1.tsv"
NOT_ESTIMABLE_FILE="${STATUS_ROOT}/copykat_not_estimable_units_v1.tsv"
if [[ -e "${FAIL_FILE}" || -e "${NOT_ESTIMABLE_FILE}" ]]; then echo "No-overwrite gate: cohort status exists" >&2; exit 1; fi
printf 'unit_id\tstatus\ttimestamp\n' >"${FAIL_FILE}"
printf 'unit_id\tstatus\ttimestamp\n' >"${NOT_ESTIMABLE_FILE}"
completed=0
not_estimable=0
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
      --input-dir "${input_dir}" --output-dir "${output_dir}" --n-cores 8 >"${log_file}" 2>&1; then
    if [[ -f "${output_dir}/RUN_COMPLETE" ]]; then
      completed=$((completed + 1))
    elif [[ -f "${output_dir}/RUN_NOT_ESTIMABLE" ]]; then
      printf '%s\tnot_estimable\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${NOT_ESTIMABLE_FILE}"
      not_estimable=$((not_estimable + 1))
    else
      printf '%s\tmissing_terminal_marker\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${FAIL_FILE}"
      failed=$((failed + 1))
    fi
  else
    printf '%s\tfailed\t%s\n' "${unit_id}" "$(date --iso-8601=seconds)" >>"${FAIL_FILE}"
    failed=$((failed + 1))
  fi
done
printf 'completed\tnot_estimable\tfailed\n%s\t%s\t%s\n' "${completed}" "${not_estimable}" "${failed}" >"${STATUS_ROOT}/copykat_cohort_counts_v1.tsv"
if [[ "${failed}" -eq 0 && $((completed + not_estimable)) -eq 62 ]]; then
  date --iso-8601=seconds >"${STATUS_ROOT}/COPYKAT_RUN_COMPLETE"
  echo "EXTERNAL_COPYKAT_COHORT_PASS completed=${completed} not_estimable=${not_estimable}"
else
  echo "EXTERNAL_COPYKAT_COHORT_FAILED completed=${completed} not_estimable=${not_estimable} failed=${failed}" >&2
  exit 1
fi
