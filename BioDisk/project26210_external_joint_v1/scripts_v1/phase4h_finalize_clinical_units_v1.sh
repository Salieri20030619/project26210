#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
OUT="${ROOT}/phase4h_finalize_clinical_units_v1"

export PHASE4H_ROOT="${ROOT}"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase4h_finalize_clinical_units_v1.R"

printf '\n===== STATUS =====\n'
cat "${OUT}/PHASE4H_COMPLETE"

printf '\n===== CORE QC =====\n'
column -s $'\t' -t \
    "${OUT}/qc/phase4h_qc.tsv"

printf '\n===== GSE294109 TREATMENT COUNTS =====\n'
column -s $'\t' -t \
    "${OUT}/qc/gse294109_treatment_counts.tsv"

printf '\n===== GSE294109 PATIENT-LEVEL TABLE =====\n'
awk -F'\t' '
NR == 1 || $1 == "GSE294109"
' "${OUT}/tables/patient_level_clinical_phase4h.tsv" |
column -s $'\t' -t

printf '\n===== RESIDUAL TREATMENT REVIEW =====\n'
residual_file="${OUT}/tables/residual_treatment_review.tsv"

printf 'rows=%s\n' \
    "$(( $(wc -l < "${residual_file}") - 1 ))"

column -s $'\t' -t "${residual_file}"

printf '\n===== UPDATE FIELD COUNTS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") dataset_col = i
        if ($i == "field") field_col = i
        if ($i == "update_decision") decision_col = i
    }
    next
}
{
        key = $dataset_col SUBSEP $field_col SUBSEP $decision_col
        count[key]++
}
END {
    for (x in count) {
        split(x, y, SUBSEP)
        print y[1] "\t" y[2] "\t" y[3] "\t" count[x]
    }
}
' "${OUT}/tables/phase4h_reviewed_updates.tsv" |
sort |
column -s $'\t' -t

printf '\n===== OUTPUT FILES =====\n'
find "${OUT}" \
    -maxdepth 2 \
    -type f \
    -printf '%s\t%p\n' |
sort -k2,2

printf '\n===== REPORT =====\n'
cat "${OUT}/PHASE4H_REPORT.md"
