#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
OUT="${ROOT}/phase4i_finalize_residual_treatment_v1"

export PHASE4I_ROOT="${ROOT}"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase4i_finalize_residual_treatment_v1.R"

printf '\n===== STATUS =====\n'
cat "${OUT}/PHASE4I_COMPLETE"

printf '\n===== CORE QC =====\n'
column -s $'\t' -t \
    "${OUT}/qc/phase4i_qc.tsv"

printf '\n===== RESIDUAL FINAL ADJUDICATION =====\n'
column -s $'\t' -t \
    "${OUT}/tables/residual_treatment_final_adjudication.tsv"

printf '\n===== GSE202374 PATIENT TABLE =====\n'
awk -F'\t' '
NR == 1 || $1 == "GSE202374"
' "${OUT}/tables/patient_level_clinical_phase4i.tsv" |
column -s $'\t' -t

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
    for (key in count) {
        split(key, parts, SUBSEP)
        print \
            parts[1] "\t" \
            parts[2] "\t" \
            parts[3] "\t" \
            count[key]
    }
}
' "${OUT}/tables/phase4i_reviewed_updates.tsv" |
sort |
column -s $'\t' -t

printf '\n===== OUTPUT FILES =====\n'
find "${OUT}" \
    -maxdepth 2 \
    -type f \
    -printf '%s\t%p\n' |
sort -k2,2

printf '\n===== REPORT =====\n'
cat "${OUT}/PHASE4I_REPORT.md"
