#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
MANIFEST="${ROOT}/phase4e_clinical_metadata_harmonization_v1/tables/object_clinical_manifest.tsv"
OUT="${ROOT}/phase4f_geo_clinical_evidence_v1"
RAW="${OUT}/raw_geo"
DATASET_FILE="${OUT}/datasets.txt"

mkdir -p "${RAW}" "${OUT}/tables" "${OUT}/qc"

if [[ ! -f "${MANIFEST}" ]]; then
    printf 'ERROR: missing manifest: %s\n' "${MANIFEST}" >&2
    exit 1
fi

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "dataset") {
            dataset_col = i
        }
    }

    if (!dataset_col) {
        exit 2
    }

    next
}

$dataset_col != "" {
    print $dataset_col
}
' "${MANIFEST}" |
tr -d '\r' |
sort -u > "${DATASET_FILE}"

dataset_count=$(wc -l < "${DATASET_FILE}")

if [[ "${dataset_count}" -ne 5 ]]; then
    printf 'ERROR: expected 5 datasets, observed %s\n' \
        "${dataset_count}" >&2
    cat "${DATASET_FILE}" >&2
    exit 1
fi

download_soft() {
    local gse="$1"
    local prefix
    local url
    local destination
    local attempt

    prefix="${gse%???}nnn"
    url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/soft/${gse}_family.soft.gz"
    destination="${RAW}/${gse}_family.soft.gz"

    if [[ -s "${destination}" ]] &&
       gzip -t "${destination}" 2>/dev/null
    then
        printf 'REUSE_VALID=%s\n' "${destination}"
        return 0
    fi

    rm -f "${destination}" "${destination}.partial"

    for attempt in 1 2 3 4 5; do
        printf 'DOWNLOAD dataset=%s attempt=%s url=%s\n' \
            "${gse}" "${attempt}" "${url}"

        if curl \
            --fail \
            --location \
            --retry 20 \
            --retry-delay 5 \
            --retry-all-errors \
            --connect-timeout 30 \
            --max-time 7200 \
            --continue-at - \
            --output "${destination}.partial" \
            "${url}"
        then
            if gzip -t "${destination}.partial"; then
                mv \
                    "${destination}.partial" \
                    "${destination}"

                printf 'DOWNLOAD_OK=%s\n' "${destination}"
                return 0
            fi
        fi

        rm -f "${destination}.partial"
        sleep $((attempt * 10))
    done

    printf 'ERROR: failed to download valid GEO SOFT for %s\n' \
        "${gse}" >&2
    return 1
}

while IFS= read -r gse
do
    [[ -n "${gse}" ]] || continue
    download_soft "${gse}"
done < "${DATASET_FILE}"

export PHASE4F_ROOT="${ROOT}"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase4f_geo_clinical_evidence_v1.R"

printf '\n===== PHASE 4F STATUS =====\n'
find "${OUT}" -maxdepth 1 \
    -type f \
    \( \
        -name 'PHASE4F_GEO_EVIDENCE_COMPLETE' \
        -o -name 'PHASE4F_INCOMPLETE' \
    \) \
    -print \
    -exec cat {} \;

printf '\n===== QC =====\n'
column -s $'\t' -t \
    "${OUT}/qc/phase4f_qc.tsv"

printf '\n===== DATASET SUMMARY =====\n'
column -s $'\t' -t \
    "${OUT}/tables/phase4f_dataset_summary.tsv"

printf '\n===== OBJECT MATCH COUNTS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "match_status") status_col = i
        if ($i == "match_method") method_col = i
    }
    next
}

{
    status[$status_col]++
    method[$method_col]++
}

END {
    for (x in status) {
        print "status", x, status[x]
    }

    for (x in method) {
        if (x != "") {
            print "method", x, method[x]
        }
    }
}
' "${OUT}/tables/object_geo_exact_match.tsv" |
sort

printf '\n===== CLINICAL CANDIDATE COUNTS =====\n'
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "field") field_col = i
        if ($i == "candidate_status") status_col = i
    }
    next
}

{
    count[$field_col SUBSEP $status_col]++
}

END {
    for (x in count) {
        split(x, y, SUBSEP)
        print y[1], y[2], count[x]
    }
}
' "${OUT}/tables/object_clinical_candidates.tsv" |
sort

printf '\n===== REVIEW OVERRIDE PREVIEW =====\n'
review_file="${OUT}/tables/phase4e_manual_clinical_overrides_REVIEW_REQUIRED.tsv"

printf 'rows=%s\n' \
    "$(( $(wc -l < "${review_file}") - 1 ))"

head -n 31 "${review_file}"
