#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
MANIFEST="${ROOT}/phase4e_clinical_metadata_harmonization_v1/tables/object_clinical_manifest.tsv"
COMPLETE="${ROOT}/phase4e_clinical_metadata_harmonization_v1/PHASE4E_CLINICAL_METADATA_COMPLETE"
OUT="${ROOT}/local_transfer_v1"
PRE="${OUT}/RDS_FILES_PRE.tsv"
FINAL="${OUT}/RDS_FILES.tsv"
SUPPORT_LIST="${OUT}/SUPPORT_FILES.txt"
SUPPORT_ARCHIVE="${OUT}/project_support_v1.tar.gz"

mkdir -p "${OUT}"

if [[ ! -f "${COMPLETE}" ]]; then
    echo "ERROR: Phase 4E completion marker is missing:"
    echo "${COMPLETE}"
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: clinical manifest is missing:"
    echo "${MANIFEST}"
    exit 1
fi

export TRANSFER_ROOT="${ROOT}"
export TRANSFER_MANIFEST="${MANIFEST}"
export TRANSFER_PRE="${PRE}"

conda run --no-capture-output -n R_4.5.2 \
Rscript --vanilla - <<'RS'
options(stringsAsFactors = FALSE, scipen = 999)

root <- Sys.getenv("TRANSFER_ROOT")
manifest_file <- Sys.getenv("TRANSFER_MANIFEST")
out_file <- Sys.getenv("TRANSFER_PRE")

m <- read.delim(
  manifest_file,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  fill = TRUE
)

needed <- c(
  "dataset",
  "object_name",
  "source_file",
  "rds_read_status"
)

missing_cols <- setdiff(needed, names(m))

if (length(missing_cols)) {
  stop(
    "Missing manifest columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (nrow(m) != 95L) {
  stop("Expected 95 objects, observed ", nrow(m))
}

if (any(m$rds_read_status != "PASS")) {
  stop("Manifest contains non-PASS RDS records")
}

if (any(!file.exists(m$source_file))) {
  stop(
    "Missing source files:\n",
    paste(
      m$source_file[!file.exists(m$source_file)],
      collapse = "\n"
    )
  )
}

source_file <- normalizePath(
  m$source_file,
  winslash = "/",
  mustWork = TRUE
)

if (anyDuplicated(source_file)) {
  stop("Duplicate source RDS paths detected")
}

safe_name <- m$object_name
safe_name <- sub(
  "\\.rds$",
  "",
  safe_name,
  ignore.case = TRUE
)
safe_name <- gsub(
  "[<>:\"/\\\\|?*]",
  "_",
  safe_name
)
safe_name <- sub(
  "[. ]+$",
  "",
  safe_name
)

local_relpath <- file.path(
  "data",
  "rds",
  m$dataset,
  paste0(safe_name, ".rds")
)

local_relpath <- gsub(
  "\\\\",
  "/",
  local_relpath
)

if (anyDuplicated(tolower(local_relpath))) {
  stop(
    paste0(
      "Windows local filename collision detected. ",
      "Review object names before transfer."
    )
  )
}

size_bytes <- file.info(source_file)$size

out <- data.frame(
  dataset = m$dataset,
  object_name = m$object_name,
  source_file = source_file,
  size_bytes = sprintf("%.0f", size_bytes),
  local_relpath = local_relpath,
  stringsAsFactors = FALSE
)

write.table(
  out,
  out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat(
  "OBJECTS=", nrow(out), "\n",
  "TOTAL_CELLS=", sum(m$cells, na.rm = TRUE), "\n",
  sep = ""
)
RS

total_objects=$(( $(wc -l < "${PRE}") - 1 ))

printf \
'dataset\tobject_name\tsource_file\tsize_bytes\tlocal_relpath\tsha256\n' \
> "${FINAL}.tmp"

object_index=0

while IFS=$'\t' read -r \
    dataset object_name source_file size_bytes local_relpath
do
    object_index=$((object_index + 1))

    printf '[%d/%d] SHA256 %s\n' \
        "${object_index}" \
        "${total_objects}" \
        "${object_name}"

    file_hash=$(
        sha256sum -- "${source_file}" |
        awk '{print $1}'
    )

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${dataset}" \
        "${object_name}" \
        "${source_file}" \
        "${size_bytes}" \
        "${local_relpath}" \
        "${file_hash}" \
        >> "${FINAL}.tmp"

done < <(tail -n +2 "${PRE}")

mv "${FINAL}.tmp" "${FINAL}"

: > "${SUPPORT_LIST}"

for item in \
    AGENTS.md \
    PROJECT_STATUS.md \
    RESULTS_MANIFEST.tsv \
    scripts_v1 \
    phase4d_95_object_compatibility_audit_v1 \
    phase4e_clinical_metadata_harmonization_v1 \
    phase4e_manual_clinical_overrides.tsv
do
    if [[ -e "${ROOT}/${item}" ]]; then
        printf '%s\n' "${item}" >> "${SUPPORT_LIST}"
    fi
done

if [[ -f \
"${ROOT}/inputs/phase4e_manual_clinical_overrides.tsv" ]]
then
    printf '%s\n' \
        "inputs/phase4e_manual_clinical_overrides.tsv" \
        >> "${SUPPORT_LIST}"
fi

tar \
    -C "${ROOT}" \
    -czf "${SUPPORT_ARCHIVE}" \
    -T "${SUPPORT_LIST}"

rds_objects=$(
    awk 'END {print NR - 1}' "${FINAL}"
)

rds_bytes=$(
    awk -F'\t' '
    NR > 1 {
        total += $4
    }
    END {
        printf "%.0f", total
    }
    ' "${FINAL}"
)

max_rds_bytes=$(
    awk -F'\t' '
    NR > 1 && $4 > maximum {
        maximum = $4
    }
    END {
        printf "%.0f", maximum
    }
    ' "${FINAL}"
)

support_bytes=$(
    stat -c '%s' "${SUPPORT_ARCHIVE}"
)

rds_human=$(
    numfmt --to=iec-i --suffix=B "${rds_bytes}"
)

max_rds_human=$(
    numfmt --to=iec-i --suffix=B "${max_rds_bytes}"
)

support_human=$(
    numfmt --to=iec-i --suffix=B "${support_bytes}"
)

cat > "${OUT}/TRANSFER_SUMMARY.txt" <<SUMMARY
transfer_root=${ROOT}
objects=${rds_objects}
cells=623431
rds_total_bytes=${rds_bytes}
rds_total_human=${rds_human}
largest_rds_bytes=${max_rds_bytes}
largest_rds_human=${max_rds_human}
support_archive_bytes=${support_bytes}
support_archive_human=${support_human}
local_target=D:\\project26210_external_joint_v1
SUMMARY

(
    cd "${OUT}"

    sha256sum \
        RDS_FILES.tsv \
        project_support_v1.tar.gz \
        > CONTROL_SHA256SUMS.txt
)

rm -f "${PRE}"

printf '\n===== TRANSFER PREPARATION COMPLETE =====\n'
cat "${OUT}/TRANSFER_SUMMARY.txt"

printf '\n===== OUTPUT FILES =====\n'
ls -lh \
    "${OUT}/RDS_FILES.tsv" \
    "${OUT}/project_support_v1.tar.gz" \
    "${OUT}/CONTROL_SHA256SUMS.txt" \
    "${OUT}/TRANSFER_SUMMARY.txt"
