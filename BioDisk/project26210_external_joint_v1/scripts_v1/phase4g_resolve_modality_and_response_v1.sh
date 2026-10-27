#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
OUT="${ROOT}/phase4g_resolve_modality_and_response_v1"

export PHASE4G_ROOT="${ROOT}"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase4g_resolve_modality_and_response_v1.R"

printf '\n===== STATUS =====\n'
cat "${OUT}/PHASE4G_COMPLETE"

printf '\n===== QC =====\n'
column -s $'\t' -t \
    "${OUT}/qc/phase4g_qc.tsv"

printf '\n===== OVERRIDE SUMMARY =====\n'
column -s $'\t' -t \
    "${OUT}/tables/phase4g_override_summary.tsv"

printf '\n===== RESPONSE AUDIT =====\n'
column -s $'\t' -t \
    "${OUT}/tables/response_candidate_audit.tsv"

printf '\n===== TREATMENT AUDIT COUNTS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "audit_decision") {
            decision_col = i
        }
    }
    next
}
{
    count[$decision_col]++
}
END {
    for (x in count) {
        print x "\t" count[x]
    }
}
' "${OUT}/tables/treatment_candidate_audit.tsv" |
sort |
column -s $'\t' -t

printf '\n===== REJECTED TREATMENT VALUES =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") dataset_col = i
        if ($i == "object_name") object_col = i
        if ($i == "raw_key") key_col = i
        if ($i == "raw_value") value_col = i
        if ($i == "audit_decision") decision_col = i
    }

    print "dataset\tobject_name\traw_key\traw_value\taudit_decision"
    next
}

$decision_col ~ /^REJECT/ {
    print \
        $dataset_col "\t" \
        $object_col "\t" \
        $key_col "\t" \
        $value_col "\t" \
        $decision_col
}
' "${OUT}/tables/treatment_candidate_audit.tsv" |
column -s $'\t' -t

printf '\n===== PROPOSED PHASE 4E OVERRIDES =====\n'
override_file="${OUT}/tables/phase4e_manual_clinical_overrides_PHASE4G_ACCEPT_CANDIDATES.tsv"

printf 'rows=%s\n' \
    "$(( $(wc -l < "${override_file}") - 1 ))"

column -s $'\t' -t "${override_file}"

printf '\n===== REPORT =====\n'
cat "${OUT}/PHASE4G_REPORT.md"
