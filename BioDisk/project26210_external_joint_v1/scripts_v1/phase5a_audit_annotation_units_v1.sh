#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
OUT="${ROOT}/phase5a_audit_annotation_units_v1"

export PHASE5A_ROOT="${ROOT}"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase5a_audit_annotation_units_v1.R"

printf '\n===== STATUS =====\n'
cat "${OUT}/PHASE5A_COMPLETE"

printf '\n===== STRUCTURAL QC =====\n'
column -s $'\t' -t \
    "${OUT}/qc/phase5a_structural_qc.tsv"

printf '\n===== DATASET MAINLINE READINESS =====\n'
column -s $'\t' -t \
    "${OUT}/tables/dataset_mainline_readiness.tsv"

printf '\n===== SELECTED ANNOTATION FIELDS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") d = i
        if ($i == "annotation_field") a = i
        if ($i == "read_status") r = i
    }
    next
}
{
    key = $d SUBSEP $a SUBSEP $r
    count[key]++
}
END {
    print "dataset\tannotation_field\tread_status\tobjects"
    for (key in count) {
        split(key, part, SUBSEP)
        print part[1] "\t" part[2] "\t" \
              part[3] "\t" count[key]
    }
}
' "${OUT}/tables/object_annotation_audit.tsv" |
sort -t $'\t' -k1,1 -k2,2 |
column -s $'\t' -t

printf '\n===== LABEL MAPPING STATUS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") d = i
        if ($i == "mapping_status") s = i
        if ($i == "cells") c = i
    }
    next
}
{
    key = $d SUBSEP $s
    labels[key]++
    cells[key] += $c
}
END {
    print "dataset\tmapping_status\tlabels\tcells"
    for (key in labels) {
        split(key, part, SUBSEP)
        print part[1] "\t" part[2] "\t" \
              labels[key] "\t" cells[key]
    }
}
' "${OUT}/tables/annotation_crosswalk_candidate.tsv" |
sort -t $'\t' -k1,1 -k2,2 |
column -s $'\t' -t

printf '\n===== UNMAPPED LABELS =====\n'
awk -F'\t' '
NR == 1 {
    print
    next
}
$5 == "UNMAPPED_REVIEW_REQUIRED"
' "${OUT}/tables/annotation_crosswalk_candidate.tsv" |
column -s $'\t' -t

printf '\n===== PATIENT UNIT STATUS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") d = i
        if ($i == "patient_status") s = i
    }
    next
}
{
    key = $d SUBSEP $s
    count[key]++
}
END {
    print "dataset\tpatient_status\tobject_rows"
    for (key in count) {
        split(key, part, SUBSEP)
        print part[1] "\t" part[2] "\t" count[key]
    }
}
' "${OUT}/tables/patient_analysis_units_candidate.tsv" |
sort -t $'\t' -k1,1 -k2,2 |
column -s $'\t' -t

printf '\n===== QC FAILURE FILE COUNTS =====\n'
for qc_file in \
    object_read_failures.tsv \
    cell_count_mismatches.tsv \
    objects_without_annotation.tsv
do
    file_path="${OUT}/qc/${qc_file}"

    row_count=$(
        awk '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
        ' "${file_path}"
    )

    printf '%s\t%s\n' \
        "${qc_file}" \
        "${row_count}"
done

printf '\n===== OUTPUT FILES =====\n'
find "${OUT}" \
    -maxdepth 2 \
    -type f \
    -printf '%s\t%p\n' |
sort -k2,2

printf '\n===== REPORT =====\n'
cat "${OUT}/PHASE5A_REPORT.md"
