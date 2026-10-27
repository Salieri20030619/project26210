#!/usr/bin/env bash
set -Eeuo pipefail
umask 027

EXPECTED_HOST="osv01"
SOURCE_ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_singlecell_downloads_v1"
WORK_ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
INPUT_ROOT="${WORK_ROOT}/external_inputs_v1"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${WORK_ROOT}/phase2_prepare_${RUN_ID}"

mkdir -p "${INPUT_ROOT}" "${OUT_ROOT}"

LOG="${OUT_ROOT}/phase2_prepare.log"
DATASET_SUMMARY="${OUT_ROOT}/dataset_summary.tsv"
FILE_INVENTORY="${OUT_ROOT}/file_inventory.tsv"
COMPRESSION_QC="${OUT_ROOT}/compression_qc.tsv"
TENX_QC="${OUT_ROOT}/tenx_triplet_qc.tsv"
SPATIAL_QC="${OUT_ROOT}/spatial_bundle_qc.tsv"
HDF5_QC="${OUT_ROOT}/hdf5_inventory.tsv"
NESTED_ARCHIVES="${OUT_ROOT}/nested_archives.tsv"
METADATA_CANDIDATES="${OUT_ROOT}/metadata_candidates.tsv"

exec > >(tee -a "${LOG}") 2>&1

printf 'dataset\tarchive\tarchive_bytes\tarchive_sha256\textract_status\textract_path\tn_files\tn_gz\tn_10x\tn_10x_pass\tn_spatial_full\tn_hdf5\tn_nested_archives\n' > "${DATASET_SUMMARY}"
printf 'dataset\trelative_path\tbytes\tmtime\tfile_type\n' > "${FILE_INVENTORY}"
printf 'dataset\trelative_path\tbytes\tstatus\n' > "${COMPRESSION_QC}"
printf 'dataset\tsample_key\tmatrix\tbarcodes\tfeatures\tmatrix_genes\tmatrix_barcodes\tnnz\tfeature_rows\tbarcode_rows\tstatus\n' > "${TENX_QC}"
printf 'dataset\tsample_key\tpositions\tscalefactors\ttissue_image\tstatus\n' > "${SPATIAL_QC}"
printf 'dataset\trelative_path\tbytes\tfile_type\th5ls_status\th5ls_preview\n' > "${HDF5_QC}"
printf 'dataset\trelative_path\tbytes\ttype\n' > "${NESTED_ARCHIVES}"
printf 'dataset\trelative_path\tbytes\n' > "${METADATA_CANDIDATES}"

declare -A ARCHIVES=(
  [GSE178481]="${SOURCE_ROOT}/GSE178481/GSE178481_RAW.tar"
  [GSE207493]="${SOURCE_ROOT}/GSE207493/GSE207493_scRNA-seq.tar.gz"
  [GSE294006]="${SOURCE_ROOT}/GSE294006/GSE294006_RAW.tar"
  [GSE294109]="${SOURCE_ROOT}/GSE294109/GSE294109_RAW.tar"
  [GSE264586]="${SOURCE_ROOT}/GSE264586/GSE264586_RAW.tar"
  [GSE202374]="${SOURCE_ROOT}/GSE202374/GSE202374_RAW.tar"
)

DATASETS=(
  GSE178481
  GSE207493
  GSE294006
  GSE294109
  GSE264586
  GSE202374
)

relpath() {
  local path="$1"
  local root="$2"
  printf '%s' "${path#"${root}"/}"
}

stream_text() {
  local path="$1"

  if [[ "${path}" == *.gz ]]; then
    gzip -cd -- "${path}"
  else
    cat -- "${path}"
  fi
}

count_lines() {
  local path="$1"
  stream_text "${path}" | wc -l | awk '{print $1}'
}

matrix_dims() {
  local path="$1"
  local result

  set +o pipefail
  result="$(
    stream_text "${path}" 2>/dev/null |
      awk '!/^%/ && NF >= 3 {print $1, $2, $3; exit}'
  )"
  set -o pipefail

  printf '%s' "${result}"
}

find_first_file() {
  local dir="$1"
  shift
  local candidate

  for candidate in "$@"; do
    if [[ -f "${dir}/${candidate}" ]]; then
      printf '%s' "${dir}/${candidate}"
      return 0
    fi
  done

  return 1
}

echo "PHASE2_EXTERNAL_INPUT_PREPARATION"
echo "START=$(date --iso-8601=seconds)"
echo "HOST=$(hostname)"
echo "USER=$(whoami)"
echo "SOURCE_ROOT=${SOURCE_ROOT}"
echo "INPUT_ROOT=${INPUT_ROOT}"
echo "OUT_ROOT=${OUT_ROOT}"
df -h /mnt/BioDisk

if [[ "$(hostname -s)" != "${EXPECTED_HOST}" ]]; then
  echo "ERROR: expected host ${EXPECTED_HOST}; current host is $(hostname -s)."
  echo "Run this phase on sgs10. No archive has been extracted by this invocation."
  exit 20
fi

for command_name in \
  tar gzip find stat sha256sum file awk sed wc sort tee
do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "ERROR: missing required command: ${command_name}"
    exit 21
  fi
done

echo
echo "[EXTRACT_ARCHIVES]"

for dataset in "${DATASETS[@]}"; do
  archive="${ARCHIVES[${dataset}]}"
  final_dir="${INPUT_ROOT}/${dataset}"
  extract_status=""

  echo "Preparing ${dataset}: ${archive}"

  if [[ ! -s "${archive}" ]]; then
    echo "ERROR: archive missing or empty: ${archive}"

    printf '%s\t%s\tNA\tNA\tMISSING_ARCHIVE\t%s\t0\t0\t0\t0\t0\t0\t0\n' \
      "${dataset}" \
      "${archive}" \
      "${final_dir}" \
      >> "${DATASET_SUMMARY}"

    continue
  fi

  archive_bytes="$(stat -c '%s' "${archive}")"
  archive_sha="$(sha256sum "${archive}" | awk '{print $1}')"

  if [[ -f "${final_dir}/.extract_complete" ]]; then
    extract_status="REUSED_COMPLETE"
    echo "${dataset}: reusing completed extraction at ${final_dir}"

  elif [[ -e "${final_dir}" ]]; then
    extract_status="BLOCKED_EXISTING_UNVERIFIED"

    echo "ERROR: ${final_dir} exists without .extract_complete."
    echo "The directory will not be overwritten."

    printf '%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' \
      "${dataset}" \
      "${archive}" \
      "${archive_bytes}" \
      "${archive_sha}" \
      "${extract_status}" \
      "${final_dir}" \
      >> "${DATASET_SUMMARY}"

    continue
  else
    staging_dir="${INPUT_ROOT}/.${dataset}.incoming_${RUN_ID}_$$"
    mkdir -p "${staging_dir}"

    if ! tar -tf "${archive}" | awk '
      BEGIN {
        bad=0
      }

      /^\// {
        bad=1
      }

      /(^|\/)\.\.($|\/)/ {
        bad=1
      }

      END {
        exit(bad ? 1 : 0)
      }
    '; then
      echo "ERROR: unsafe archive member detected in ${archive}."
      echo "Staging directory retained: ${staging_dir}"

      printf '%s\t%s\t%s\t%s\tUNSAFE_MEMBER\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' \
        "${dataset}" \
        "${archive}" \
        "${archive_bytes}" \
        "${archive_sha}" \
        "${staging_dir}" \
        >> "${DATASET_SUMMARY}"

      continue
    fi

    if ! tar \
      --extract \
      --file "${archive}" \
      --directory "${staging_dir}" \
      --no-same-owner \
      --no-same-permissions
    then
      echo "ERROR: extraction failed for ${dataset}."
      echo "Incomplete staging directory retained: ${staging_dir}"

      printf '%s\t%s\t%s\t%s\tEXTRACT_FAILED\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' \
        "${dataset}" \
        "${archive}" \
        "${archive_bytes}" \
        "${archive_sha}" \
        "${staging_dir}" \
        >> "${DATASET_SUMMARY}"

      continue
    fi

    printf \
      'dataset=%s\narchive=%s\narchive_bytes=%s\narchive_sha256=%s\nextracted_at=%s\nhost=%s\n' \
      "${dataset}" \
      "${archive}" \
      "${archive_bytes}" \
      "${archive_sha}" \
      "$(date --iso-8601=seconds)" \
      "$(hostname -s)" \
      > "${staging_dir}/.extract_complete"

    mv -- "${staging_dir}" "${final_dir}"

    extract_status="EXTRACTED"
    echo "${dataset}: extracted to ${final_dir}"
  fi

  echo "Inventorying ${dataset} ..."

  while IFS= read -r -d '' input_file; do
    relative="$(relpath "${input_file}" "${final_dir}")"
    bytes="$(stat -c '%s' "${input_file}")"
    mtime="$(stat -c '%y' "${input_file}" | sed 's/\t/ /g')"
    type_text="$(
      file -b -- "${input_file}" |
        tr '\t\n' '  '
    )"

    printf '%s\t%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${relative}" \
      "${bytes}" \
      "${mtime}" \
      "${type_text}" \
      >> "${FILE_INVENTORY}"

  done < <(
    find "${final_dir}" \
      -type f \
      ! -name '.extract_complete' \
      -print0 |
      sort -z
  )

  echo "Validating extracted gzip members for ${dataset} ..."

  while IFS= read -r -d '' gzip_file; do
    relative="$(relpath "${gzip_file}" "${final_dir}")"
    bytes="$(stat -c '%s' "${gzip_file}")"

    if gzip -t -- "${gzip_file}"; then
      gzip_status="PASS"
    else
      gzip_status="FAIL"
    fi

    printf '%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${relative}" \
      "${bytes}" \
      "${gzip_status}" \
      >> "${COMPRESSION_QC}"

  done < <(
    find "${final_dir}" \
      -type f \
      -name '*.gz' \
      -print0 |
      sort -z
  )

  while IFS= read -r -d '' nested_file; do
    relative="$(relpath "${nested_file}" "${final_dir}")"
    bytes="$(stat -c '%s' "${nested_file}")"
    nested_type="$(
      file -b -- "${nested_file}" |
        tr '\t\n' '  '
    )"

    printf '%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${relative}" \
      "${bytes}" \
      "${nested_type}" \
      >> "${NESTED_ARCHIVES}"

  done < <(
    find "${final_dir}" \
      -type f \
      \( \
        -name '*.tar' \
        -o -name '*.tar.gz' \
        -o -name '*.tgz' \
        -o -name '*.zip' \
        -o -name '*.bz2' \
        -o -name '*.xz' \
      \) \
      -print0 |
      sort -z
  )

  while IFS= read -r -d '' metadata_file; do
    relative="$(relpath "${metadata_file}" "${final_dir}")"
    bytes="$(stat -c '%s' "${metadata_file}")"

    printf '%s\t%s\t%s\n' \
      "${dataset}" \
      "${relative}" \
      "${bytes}" \
      >> "${METADATA_CANDIDATES}"

  done < <(
    find "${final_dir}" \
      -type f \
      \( \
        -iname '*meta*.csv' \
        -o -iname '*meta*.tsv' \
        -o -iname '*meta*.txt' \
        -o -iname '*annot*.csv' \
        -o -iname '*annot*.tsv' \
        -o -iname '*annot*.txt' \
        -o -iname '*clinical*.csv' \
        -o -iname '*clinical*.tsv' \
        -o -iname '*clinical*.txt' \
        -o -iname '*sample*.csv' \
        -o -iname '*sample*.tsv' \
        -o -iname '*sample*.txt' \
        -o -iname '*.xlsx' \
        -o -iname '*.xls' \
      \) \
      -print0 |
      sort -z
  )

  echo "Checking 10x triplets for ${dataset} ..."

  while IFS= read -r -d '' matrix_file; do
    matrix_dir="$(dirname "${matrix_file}")"
    matrix_base="$(basename "${matrix_file}")"

    if [[ "${matrix_base}" == *matrix.mtx.gz ]]; then
      prefix="${matrix_base%matrix.mtx.gz}"

      barcodes_file="$(
        find_first_file "${matrix_dir}" \
          "${prefix}barcodes.tsv.gz" \
          "${prefix}barcodes.tsv" \
          || true
      )"

      features_file="$(
        find_first_file "${matrix_dir}" \
          "${prefix}features.tsv.gz" \
          "${prefix}features.tsv" \
          "${prefix}genes.tsv.gz" \
          "${prefix}genes.tsv" \
          || true
      )"
    else
      prefix="${matrix_base%matrix.mtx}"

      barcodes_file="$(
        find_first_file "${matrix_dir}" \
          "${prefix}barcodes.tsv" \
          "${prefix}barcodes.tsv.gz" \
          || true
      )"

      features_file="$(
        find_first_file "${matrix_dir}" \
          "${prefix}features.tsv" \
          "${prefix}features.tsv.gz" \
          "${prefix}genes.tsv" \
          "${prefix}genes.tsv.gz" \
          || true
      )"
    fi

    sample_key="${prefix%_}"

    if [[ -z "${sample_key}" ]]; then
      sample_key="$(basename "${matrix_dir}")"
    fi

    dims="$(matrix_dims "${matrix_file}")"

    matrix_genes="NA"
    matrix_barcodes="NA"
    nnz="NA"

    read -r matrix_genes matrix_barcodes nnz <<< "${dims}" || true

    matrix_genes="${matrix_genes:-NA}"
    matrix_barcodes="${matrix_barcodes:-NA}"
    nnz="${nnz:-NA}"

    feature_rows="NA"
    barcode_rows="NA"
    triplet_status="PASS"

    if [[ -z "${barcodes_file}" || -z "${features_file}" ]]; then
      triplet_status="MISSING_COMPANION"

    elif [[ "${matrix_genes}" == "NA" || "${matrix_barcodes}" == "NA" ]]; then
      triplet_status="HEADER_UNREADABLE"
    else
      feature_rows="$(count_lines "${features_file}")"
      barcode_rows="$(count_lines "${barcodes_file}")"

      if [[ "${matrix_genes}" != "${feature_rows}" ||
            "${matrix_barcodes}" != "${barcode_rows}" ]]
      then
        triplet_status="DIMENSION_MISMATCH"
      fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${sample_key}" \
      "$(relpath "${matrix_file}" "${final_dir}")" \
      "${barcodes_file:+$(relpath "${barcodes_file}" "${final_dir}")}" \
      "${features_file:+$(relpath "${features_file}" "${final_dir}")}" \
      "${matrix_genes}" \
      "${matrix_barcodes}" \
      "${nnz}" \
      "${feature_rows}" \
      "${barcode_rows}" \
      "${triplet_status}" \
      >> "${TENX_QC}"

    positions_file="$(
      find_first_file "${matrix_dir}" \
        "${prefix}tissue_positions_list.csv.gz" \
        "${prefix}tissue_positions_list.csv" \
        "${prefix}tissue_positions.csv.gz" \
        "${prefix}tissue_positions.csv" \
        || true
    )"

    scalefactors_file="$(
      find_first_file "${matrix_dir}" \
        "${prefix}scalefactors_json.json.gz" \
        "${prefix}scalefactors_json.json" \
        || true
    )"

    set +o pipefail
    tissue_image="$(
      find "${matrix_dir}" \
        -maxdepth 1 \
        -type f \
        \( \
          -iname "${prefix}*tissue*image*.png" \
          -o -iname "${prefix}*tissue*image*.png.gz" \
          -o -iname "${prefix}*hires*.png" \
          -o -iname "${prefix}*hires*.png.gz" \
        \) \
        -print \
        -quit \
        2>/dev/null
    )"
    set -o pipefail

    if [[ -n "${positions_file}" &&
          -n "${scalefactors_file}" &&
          -n "${tissue_image}" ]]
    then
      spatial_status="SPATIAL_FULL"

    elif [[ -n "${positions_file}" ||
            -n "${scalefactors_file}" ||
            -n "${tissue_image}" ]]
    then
      spatial_status="SPATIAL_INCOMPLETE"
    else
      spatial_status="NO_SPATIAL_FILES"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${sample_key}" \
      "${positions_file:+$(relpath "${positions_file}" "${final_dir}")}" \
      "${scalefactors_file:+$(relpath "${scalefactors_file}" "${final_dir}")}" \
      "${tissue_image:+$(relpath "${tissue_image}" "${final_dir}")}" \
      "${spatial_status}" \
      >> "${SPATIAL_QC}"

  done < <(
    find "${final_dir}" \
      -type f \
      \( \
        -name '*matrix.mtx.gz' \
        -o -name '*matrix.mtx' \
      \) \
      -print0 |
      sort -z
  )

  echo "Inventorying HDF5-family files for ${dataset} ..."

  while IFS= read -r -d '' hdf5_file; do
    relative="$(relpath "${hdf5_file}" "${final_dir}")"
    bytes="$(stat -c '%s' "${hdf5_file}")"
    type_text="$(
      file -b -- "${hdf5_file}" |
        tr '\t\n' '  '
    )"

    h5ls_status="NOT_AVAILABLE"
    h5ls_preview=""

    if command -v h5ls >/dev/null 2>&1; then
      if h5_output="$(h5ls "${hdf5_file}" 2>&1 | head -n 8)"; then
        h5ls_status="PASS"
      else
        h5ls_status="FAIL"
      fi

      h5ls_preview="$(
        printf '%s' "${h5_output}" |
          tr '\t\n' '  ' |
          sed 's/[[:space:]][[:space:]]*/ /g'
      )"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${dataset}" \
      "${relative}" \
      "${bytes}" \
      "${type_text}" \
      "${h5ls_status}" \
      "${h5ls_preview}" \
      >> "${HDF5_QC}"

  done < <(
    find "${final_dir}" \
      -type f \
      \( \
        -iname '*.h5' \
        -o -iname '*.h5ad' \
        -o -iname '*.hdf5' \
      \) \
      -print0 |
      sort -z
  )

  n_files="$(
    find "${final_dir}" \
      -type f \
      ! -name '.extract_complete' |
      wc -l |
      awk '{print $1}'
  )"

  n_gz="$(
    find "${final_dir}" \
      -type f \
      -name '*.gz' |
      wc -l |
      awk '{print $1}'
  )"

  n_10x="$(
    awk -F '\t' -v d="${dataset}" '
      NR>1 && $1==d {
        n++
      }

      END {
        print n+0
      }
    ' "${TENX_QC}"
  )"

  n_10x_pass="$(
    awk -F '\t' -v d="${dataset}" '
      NR>1 && $1==d && $11=="PASS" {
        n++
      }

      END {
        print n+0
      }
    ' "${TENX_QC}"
  )"

  n_spatial_full="$(
    awk -F '\t' -v d="${dataset}" '
      NR>1 && $1==d && $6=="SPATIAL_FULL" {
        n++
      }

      END {
        print n+0
      }
    ' "${SPATIAL_QC}"
  )"

  n_hdf5="$(
    awk -F '\t' -v d="${dataset}" '
      NR>1 && $1==d {
        n++
      }

      END {
        print n+0
      }
    ' "${HDF5_QC}"
  )"

  n_nested="$(
    awk -F '\t' -v d="${dataset}" '
      NR>1 && $1==d {
        n++
      }

      END {
        print n+0
      }
    ' "${NESTED_ARCHIVES}"
  )"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${dataset}" \
    "${archive}" \
    "${archive_bytes}" \
    "${archive_sha}" \
    "${extract_status}" \
    "${final_dir}" \
    "${n_files}" \
    "${n_gz}" \
    "${n_10x}" \
    "${n_10x_pass}" \
    "${n_spatial_full}" \
    "${n_hdf5}" \
    "${n_nested}" \
    >> "${DATASET_SUMMARY}"
done

gzip_failures="$(
  awk -F '\t' '
    NR>1 && $4!="PASS" {
      n++
    }

    END {
      print n+0
    }
  ' "${COMPRESSION_QC}"
)"

tenx_failures="$(
  awk -F '\t' '
    NR>1 && $11!="PASS" {
      n++
    }

    END {
      print n+0
    }
  ' "${TENX_QC}"
)"

blocked_datasets="$(
  awk -F '\t' '
    NR>1 && $5 ~ /MISSING|BLOCKED|FAILED|UNSAFE/ {
      n++
    }

    END {
      print n+0
    }
  ' "${DATASET_SUMMARY}"
)"

echo
echo "[SUMMARY]"

column -s $'\t' -t "${DATASET_SUMMARY}" 2>/dev/null ||
  cat "${DATASET_SUMMARY}"

echo
echo "10X_QC_STATUS_COUNTS"

awk -F '\t' '
  NR>1 {
    count[$11]++
  }

  END {
    for (x in count) {
      print x, count[x]
    }
  }
' "${TENX_QC}" |
  sort

echo
echo "SPATIAL_QC_STATUS_COUNTS"

awk -F '\t' '
  NR>1 {
    count[$6]++
  }

  END {
    for (x in count) {
      print x, count[x]
    }
  }
' "${SPATIAL_QC}" |
  sort

echo
echo "GZIP_FAILURES=${gzip_failures}"
echo "TENX_FAILURES=${tenx_failures}"
echo "BLOCKED_DATASETS=${blocked_datasets}"
echo "END=$(date --iso-8601=seconds)"

echo
echo "PHASE2_COMPLETE"
echo "OUT_ROOT=${OUT_ROOT}"
echo "INPUT_ROOT=${INPUT_ROOT}"
echo "DATASET_SUMMARY=${DATASET_SUMMARY}"
echo "TENX_QC=${TENX_QC}"
echo "SPATIAL_QC=${SPATIAL_QC}"
echo "HDF5_QC=${HDF5_QC}"
echo "NESTED_ARCHIVES=${NESTED_ARCHIVES}"
echo "METADATA_CANDIDATES=${METADATA_CANDIDATES}"
echo "LOG=${LOG}"

if (( gzip_failures > 0 ||
      tenx_failures > 0 ||
      blocked_datasets > 0 ))
then
  echo "QC_ATTENTION_REQUIRED"
  exit 2
fi

echo "QC_PASS_FOR_DETECTED_10X_INPUTS"
