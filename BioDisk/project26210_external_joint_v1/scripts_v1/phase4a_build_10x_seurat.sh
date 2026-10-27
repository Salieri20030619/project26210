#!/usr/bin/env bash

set -Eeuo pipefail

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export VECLIB_MAXIMUM_THREADS=2

WORK_ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
INPUT_ROOT="${WORK_ROOT}/external_inputs_v1"
PHASE3_ROOT="${WORK_ROOT}/phase3_input_audit_20260808_215808"
METADATA="${PHASE3_ROOT}/geo_sample_metadata_summary.tsv"

SCRIPT_ROOT="${WORK_ROOT}/scripts_v1"
R_SCRIPT="${SCRIPT_ROOT}/phase4a_build_one_seurat.R"

OUT_ROOT="${WORK_ROOT}/phase4a_10x_seurat_v1"
OBJECT_DIR="${OUT_ROOT}/objects"
SUMMARY_DIR="${OUT_ROOT}/sample_summaries"
DONE_DIR="${OUT_ROOT}/done"
FAILED_DIR="${OUT_ROOT}/failed"
LOG_DIR="${OUT_ROOT}/sample_logs"
MANIFEST="${OUT_ROOT}/sample_manifest.tsv"

LOCK_DIR="${OUT_ROOT}.lock"
COMPLETE_MARKER="${OUT_ROOT}/PHASE4A_COMPLETE"

if [[ -s "${COMPLETE_MARKER}" ]]; then
    echo "Phase 4A is already complete:"
    echo "${COMPLETE_MARKER}"
    exit 0
fi

if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
    echo "ERROR: another Phase 4A process may already be running."
    echo "LOCK=${LOCK_DIR}"
    exit 1
fi

cleanup_lock() {
    rmdir "${LOCK_DIR}" 2>/dev/null || true
}

trap cleanup_lock EXIT

mkdir -p \
    "${OUT_ROOT}" \
    "${OBJECT_DIR}" \
    "${SUMMARY_DIR}" \
    "${DONE_DIR}" \
    "${FAILED_DIR}" \
    "${LOG_DIR}"

exec > >(
    tee -a "${OUT_ROOT}/phase4a_driver.log"
) 2>&1

echo "START=$(date --iso-8601=seconds)"
echo "HOST=$(hostname)"
echo "USER=$(id -un)"
echo "CONDA_ENV=${CONDA_DEFAULT_ENV:-NONE}"
echo "WORK_ROOT=${WORK_ROOT}"
echo "OUT_ROOT=${OUT_ROOT}"

if [[ "${CONDA_DEFAULT_ENV:-}" != "R_4.5.2" ]]; then
    echo "ERROR: Phase 4A requires conda environment R_4.5.2"
    exit 1
fi

for required_path in \
    "${INPUT_ROOT}" \
    "${METADATA}" \
    "${R_SCRIPT}"
do
    if [[ ! -e "${required_path}" ]]; then
        echo "ERROR: required path does not exist: ${required_path}"
        exit 1
    fi
done

available_kb="$(
    df -k --output=avail "${WORK_ROOT}" |
    tail -n 1 |
    tr -d ' '
)"

available_gb=$((available_kb / 1024 / 1024))

echo "AVAILABLE_DISK_GB=${available_gb}"

if (( available_gb < 50 )); then
    echo "ERROR: fewer than 50 GB are available."
    echo "Phase 4A was not started."
    exit 1
fi

Rscript --vanilla -e '
required <- c("Matrix", "Seurat")
missing <- required[
    !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "))
}
cat("R_PACKAGE_PREFLIGHT=PASS\n")
cat("R_VERSION=", R.version.string, "\n", sep = "")
cat(
    "SEURAT_VERSION=",
    as.character(packageVersion("Seurat")),
    "\n",
    sep = ""
)
'

Rscript --vanilla -e 'sessionInfo()' \
    > "${OUT_ROOT}/sessionInfo.txt" 2>&1

manifest_unsorted="${OUT_ROOT}/sample_manifest.unsorted.tsv"
manifest_sorted="${OUT_ROOT}/sample_manifest.sorted.tsv"

printf \
    'dataset\tsample_key\tgsm\tmatrix\tfeatures\tbarcodes\n' \
    > "${manifest_unsorted}"

declare -A observed_counts=()

add_dataset() {
    local dataset="$1"
    local search_root="$2"

    if [[ ! -d "${search_root}" ]]; then
        echo "ERROR: dataset directory not found: ${search_root}"
        exit 1
    fi

    while IFS= read -r matrix_path; do
        local matrix_dir
        local matrix_name
        local prefix
        local sample_key
        local feature_path
        local barcode_path
        local gsm

        matrix_dir="$(dirname "${matrix_path}")"
        matrix_name="$(basename "${matrix_path}")"
        prefix="${matrix_name%matrix.mtx.gz}"

        if [[ -n "${prefix}" ]]; then
            sample_key="${prefix%_}"
        else
            sample_key="$(basename "${matrix_dir}")"
        fi

        feature_path="${matrix_dir}/${prefix}features.tsv.gz"

        if [[ ! -s "${feature_path}" ]]; then
            feature_path="${matrix_dir}/${prefix}genes.tsv.gz"
        fi

        barcode_path="${matrix_dir}/${prefix}barcodes.tsv.gz"

        if [[ ! -s "${feature_path}" ]]; then
            echo "ERROR: missing features for ${matrix_path}"
            exit 1
        fi

        if [[ ! -s "${barcode_path}" ]]; then
            echo "ERROR: missing barcodes for ${matrix_path}"
            exit 1
        fi

        gsm="$(
            printf '%s\n' "${matrix_path}" |
            grep -oE 'GSM[0-9]+' |
            head -n 1 ||
            true
        )"

        if [[ -z "${gsm}" && "${dataset}" == "GSE207493" ]]; then
            gsm="$(
                awk -F $'\t' \
                    -v key="${sample_key}" '
                    $1 == "GSE207493" &&
                    $3 == ("scRNA-seq_" key) {
                        print $2
                        exit
                    }
                    ' "${METADATA}"
            )"
        fi

        if [[ -z "${gsm}" ]]; then
            gsm="NA"
        fi

        printf \
            '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${dataset}" \
            "${sample_key}" \
            "${gsm}" \
            "${matrix_path}" \
            "${feature_path}" \
            "${barcode_path}" \
            >> "${manifest_unsorted}"

            (( observed_counts["${dataset}"] += 1 ))
    done < <(
        find "${search_root}" \
            -type f \
            -name '*matrix.mtx.gz' |
        sort
    )
}

add_dataset \
    "GSE207493" \
    "${INPUT_ROOT}/GSE207493/scRNA-seq_extracted_v1"

add_dataset \
    "GSE294109" \
    "${INPUT_ROOT}/GSE294109"

add_dataset \
    "GSE264586" \
    "${INPUT_ROOT}/GSE264586"

add_dataset \
    "GSE202374" \
    "${INPUT_ROOT}/GSE202374"

{
    head -n 1 "${manifest_unsorted}"

    tail -n +2 "${manifest_unsorted}" |
    sort -t $'\t' -k1,1 -k2,2
} > "${manifest_sorted}"

mv "${manifest_sorted}" "${MANIFEST}"
rm -f "${manifest_unsorted}"

declare -A expected_counts=(
    ["GSE207493"]=19
    ["GSE294109"]=25
    ["GSE264586"]=19
    ["GSE202374"]=6
)

for dataset in \
    GSE207493 \
    GSE294109 \
    GSE264586 \
    GSE202374
do
    observed="${observed_counts[${dataset}]:-0}"
    expected="${expected_counts[${dataset}]}"

    echo "${dataset}_MATRICES=${observed}"

    if [[ "${observed}" -ne "${expected}" ]]; then
        echo \
            "ERROR: ${dataset} expected ${expected} matrices, found ${observed}"
        exit 1
    fi
done

manifest_rows="$(
    awk 'END {print NR - 1}' "${MANIFEST}"
)"

echo "TOTAL_10X_MATRICES=${manifest_rows}"

if [[ "${manifest_rows}" -ne 69 ]]; then
    echo "ERROR: expected exactly 69 standard 10x matrices."
    exit 1
fi

printf \
    'dataset\tsample_key\tgsm\treason\n' \
    > "${OUT_ROOT}/failed_samples.tsv"

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r \
    dataset sample_key gsm matrix_path feature_path barcode_path
do
    safe_sample="$(
        printf '%s' "${sample_key}" |
        sed 's/[^A-Za-z0-9._-]/_/g'
    )"

    object_id="${dataset}__${safe_sample}"
    object_path="${OBJECT_DIR}/${object_id}.rds"
    summary_path="${SUMMARY_DIR}/${object_id}.summary.tsv"
    done_path="${DONE_DIR}/${object_id}.done"
    failed_path="${FAILED_DIR}/${object_id}.failed"
    sample_log="${LOG_DIR}/${object_id}.log"

    if [[ -s "${done_path}" &&
          -s "${object_path}" &&
          -s "${summary_path}" ]]
    then
        echo "SKIP_COMPLETED=${object_id}"
        continue
    fi

    echo
    echo "PROCESSING=${object_id}"
    echo "START_SAMPLE_TIME=$(date --iso-8601=seconds)"

    rm -f "${failed_path}"

    if Rscript --vanilla "${R_SCRIPT}" \
        "${dataset}" \
        "${sample_key}" \
        "${gsm}" \
        "${matrix_path}" \
        "${feature_path}" \
        "${barcode_path}" \
        "${object_path}" \
        "${summary_path}" \
        "${done_path}" \
        "${object_id}" \
        > "${sample_log}" 2>&1
    then
        echo "PASS=${object_id}"
    else
        exit_code=$?

        printf \
            'status=FAIL\nexit_code=%s\nsample=%s\nlog=%s\n' \
            "${exit_code}" \
            "${object_id}" \
            "${sample_log}" \
            > "${failed_path}"

        printf \
            '%s\t%s\t%s\tRscript_exit_%s\n' \
            "${dataset}" \
            "${sample_key}" \
            "${gsm}" \
            "${exit_code}" \
            >> "${OUT_ROOT}/failed_samples.tsv"

        echo "FAIL=${object_id}"
        echo "LOG=${sample_log}"
    fi

    echo "END_SAMPLE_TIME=$(date --iso-8601=seconds)"
done

summary_output="${OUT_ROOT}/sample_qc_summary.tsv"

printf \
    'dataset\tsample_key\tgsm\tgenes\tcells\tnnz\tmedian_nCount\tmedian_nFeature\tmedian_percent_mt\tmedian_percent_ribo\tbasic_qc_pass_cells\tbasic_qc_pass_rate\tfeature_type_filter_applied\n' \
    > "${summary_output}"

while IFS= read -r summary_file; do
    tail -n +2 "${summary_file}" >> "${summary_output}"
done < <(
    find "${SUMMARY_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.summary.tsv' |
    sort
)

done_count="$(
    find "${DONE_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.done' |
    wc -l |
    awk '{print $1}'
)"

failed_count="$(
    find "${FAILED_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.failed' |
    wc -l |
    awk '{print $1}'
)"

object_count="$(
    find "${OBJECT_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.rds' |
    wc -l |
    awk '{print $1}'
)"

echo
echo "[SUMMARY]"
echo "EXPECTED_OBJECTS=69"
echo "DONE_OBJECTS=${done_count}"
echo "RDS_OBJECTS=${object_count}"
echo "FAILED_OBJECTS=${failed_count}"
echo "QC_SUMMARY=${summary_output}"
echo "END=$(date --iso-8601=seconds)"

sha256sum \
    "${R_SCRIPT}" \
    "$0" \
    > "${OUT_ROOT}/script_sha256.txt"

if [[ "${done_count}" -eq 69 &&
      "${object_count}" -eq 69 &&
      "${failed_count}" -eq 0 ]]
then
    echo "Calculating object SHA256 checksums ..."

    sha256sum "${OBJECT_DIR}"/*.rds |
    sort -k2,2 \
        > "${OUT_ROOT}/object_sha256.txt"

    {
        echo "status=PASS"
        echo "objects=69"
        echo "completed_at=$(date --iso-8601=seconds)"
        echo "qc_summary=${summary_output}"
    } > "${COMPLETE_MARKER}"

    echo "PHASE4A_COMPLETE"
else
    echo "PHASE4A_COMPLETE_WITH_FAILURES"
    exit 2
fi
