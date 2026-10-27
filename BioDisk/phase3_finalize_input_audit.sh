#!/usr/bin/env bash
set -Eeuo pipefail

umask 027

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
INPUT_ROOT="${ROOT}/external_inputs_v1"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${ROOT}/phase3_input_audit_${STAMP}"

GSE207_DIR="${INPUT_ROOT}/GSE207493"
GSE207_ZIP="${GSE207_DIR}/scRNA-seq.zip"
GSE207_FINAL="${GSE207_DIR}/scRNA-seq_extracted_v1"
GSE178_DIR="${INPUT_ROOT}/GSE178481"
GSE294006_DIR="${INPUT_ROOT}/GSE294006"

mkdir -p "${OUT_ROOT}"

exec > >(tee -a "${OUT_ROOT}/phase3_audit.log") 2>&1

echo "PHASE3_FINALIZE_INPUT_AUDIT"
echo "START=$(date --iso-8601=seconds)"
echo "HOST=$(hostname -s)"
echo "USER=$(id -un)"
echo "CONDA_ENV=${CONDA_DEFAULT_ENV:-NONE}"
echo "ROOT=${ROOT}"
echo "INPUT_ROOT=${INPUT_ROOT}"
echo "OUT_ROOT=${OUT_ROOT}"

if [[ "${CONDA_DEFAULT_ENV:-}" != "batch_bench" ]]; then
  echo "ERROR: activate conda environment batch_bench before running."
  exit 1
fi

for cmd in python unzip zipinfo gzip curl sha256sum find sort awk stat; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: required command not found: ${cmd}"
    exit 1
  fi
done

for required in "${GSE207_ZIP}" "${GSE178_DIR}" "${GSE294006_DIR}"; do
  if [[ ! -e "${required}" ]]; then
    echo "ERROR: required input missing: ${required}"
    exit 1
  fi
done

echo
echo "[1/4] GSE207493 ZIP SAFETY AND EXTRACTION"

ZIP_MEMBERS="${OUT_ROOT}/gse207493_zip_members.txt"
ZIP_SAFETY="${OUT_ROOT}/gse207493_zip_safety.tsv"

zipinfo -1 "${GSE207_ZIP}" > "${ZIP_MEMBERS}"

printf 'check\tstatus\tdetail\n' > "${ZIP_SAFETY}"

if awk '
  BEGIN { bad=0 }
  /^\// { bad=1 }
  /^[A-Za-z]:[\\/]/ { bad=1 }
  /(^|[\\/])\.\.([\\/]|$)/ { bad=1 }
  END { exit(bad ? 1 : 0) }
' "${ZIP_MEMBERS}"; then
  printf 'member_paths\tPASS\tno_absolute_or_parent_paths\n' >> "${ZIP_SAFETY}"
else
  printf 'member_paths\tFAIL\tabsolute_or_parent_path_detected\n' >> "${ZIP_SAFETY}"
  echo "ERROR: unsafe path detected in GSE207493 ZIP."
  exit 1
fi

duplicate_members="$(sort "${ZIP_MEMBERS}" | uniq -d | wc -l | awk '{print $1}')"

if [[ "${duplicate_members}" == "0" ]]; then
  printf 'duplicate_members\tPASS\t0\n' >> "${ZIP_SAFETY}"
else
  printf 'duplicate_members\tFAIL\t%s\n' "${duplicate_members}" >> "${ZIP_SAFETY}"
  echo "ERROR: duplicate member names detected in GSE207493 ZIP."
  exit 1
fi

symlink_members="$({ zipinfo -l "${GSE207_ZIP}" || true; } | awk '$1 ~ /^l/ {n++} END {print n+0}')"

if [[ "${symlink_members}" == "0" ]]; then
  printf 'symlink_members\tPASS\t0\n' >> "${ZIP_SAFETY}"
else
  printf 'symlink_members\tFAIL\t%s\n' "${symlink_members}" >> "${ZIP_SAFETY}"
  echo "ERROR: symbolic-link members detected in GSE207493 ZIP."
  exit 1
fi

if unzip -t "${GSE207_ZIP}" >/dev/null; then
  printf 'zip_integrity\tPASS\tunzip_test_passed\n' >> "${ZIP_SAFETY}"
else
  printf 'zip_integrity\tFAIL\tunzip_test_failed\n' >> "${ZIP_SAFETY}"
  echo "ERROR: GSE207493 ZIP integrity test failed."
  exit 1
fi

if [[ -f "${GSE207_FINAL}/.extract_complete" ]]; then
  echo "Reusing validated extraction: ${GSE207_FINAL}"
  GSE207_WORK="${GSE207_FINAL}"
  GSE207_MODE="REUSED"
elif [[ -e "${GSE207_FINAL}" ]]; then
  echo "ERROR: extraction target exists without completion marker: ${GSE207_FINAL}"
  exit 1
else
  GSE207_STAGE="${GSE207_DIR}/.scRNA-seq_extracted_v1.staging_${STAMP}_$$"
  mkdir "${GSE207_STAGE}"
  echo "Extracting to staging directory: ${GSE207_STAGE}"
  unzip -q "${GSE207_ZIP}" -d "${GSE207_STAGE}"
  GSE207_WORK="${GSE207_STAGE}"
  GSE207_MODE="EXTRACTED"
fi

while IFS= read -r -d '' gz_file; do
  gzip -t -- "${gz_file}"
done < <(find "${GSE207_WORK}" -type f -name '*.gz' -print0 | sort -z)

export GSE207_WORK OUT_ROOT

python - <<'PY'
import gzip
import os
from pathlib import Path

root = Path(os.environ["GSE207_WORK"])
out = Path(os.environ["OUT_ROOT"]) / "gse207493_10x_qc.tsv"

def count_lines(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    count = 0
    last = b""
    with opener(path, "rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            count += block.count(b"\n")
            last = block[-1:]
    if last and last != b"\n":
        count += 1
    return count

def matrix_header(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.strip().split()
            if len(fields) != 3:
                raise ValueError(f"invalid MatrixMarket dimensions: {path}")
            return tuple(map(int, fields))
    raise ValueError(f"missing MatrixMarket dimensions: {path}")

rows = []
for matrix in sorted(root.rglob("matrix.mtx.gz")):
    sample_dir = matrix.parent
    barcodes = sample_dir / "barcodes.tsv.gz"
    features = sample_dir / "features.tsv.gz"
    status = "PASS"
    reason = "OK"

    if not barcodes.is_file() or not features.is_file():
        genes = cells = nnz = feature_rows = barcode_rows = -1
        status = "FAIL"
        reason = "MISSING_COMPANION"
    else:
        genes, cells, nnz = matrix_header(matrix)
        feature_rows = count_lines(features)
        barcode_rows = count_lines(barcodes)
        if genes != feature_rows or cells != barcode_rows:
            status = "FAIL"
            reason = "DIMENSION_MISMATCH"

    rows.append((sample_dir.name, matrix.relative_to(root), genes, cells, nnz,
                 feature_rows, barcode_rows, status, reason))

with out.open("w", encoding="utf-8", newline="") as handle:
    handle.write("sample\tmatrix\tgenes\tbarcodes\tnnz\tfeature_rows\tbarcode_rows\tstatus\treason\n")
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")

if len(rows) != 19:
    raise SystemExit(f"expected 19 matrices, detected {len(rows)}")
if any(row[7] != "PASS" for row in rows):
    raise SystemExit("one or more GSE207493 10x triplets failed")

print(f"GSE207493_10X_MATRICES={len(rows)}")
print("GSE207493_10X_QC=PASS")
PY

if [[ "${GSE207_MODE}" == "EXTRACTED" ]]; then
  printf 'archive=%s\narchive_sha256=%s\nextracted_at=%s\nhost=%s\n' \
    "${GSE207_ZIP}" \
    "$(sha256sum "${GSE207_ZIP}" | awk '{print $1}')" \
    "$(date --iso-8601=seconds)" \
    "$(hostname -s)" \
    > "${GSE207_WORK}/.extract_complete"

  mv -- "${GSE207_WORK}" "${GSE207_FINAL}"
  echo "GSE207493 extracted to: ${GSE207_FINAL}"
fi

echo
echo "[2/4] GSE178481 DENSE CSV STREAMING AUDIT"

export GSE178_DIR

python - <<'PY'
import csv
import gzip
import os
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)
root = Path(os.environ["GSE178_DIR"])
out = Path(os.environ["OUT_ROOT"]) / "gse178481_dense_csv_qc.tsv"

rows_out = []
for path in sorted(root.glob("*.count.csv.gz")):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) < 2:
            rows_out.append((path.name, -1, -1, -1, -1, -1, "FAIL", "INVALID_HEADER"))
            continue

        barcodes = header[1:]
        n_cells = len(barcodes)
        duplicate_barcodes = n_cells - len(set(barcodes))
        genes_seen = set()
        duplicate_genes = 0
        width_mismatch = 0
        n_genes = 0

        for row in reader:
            n_genes += 1
            if len(row) != len(header):
                width_mismatch += 1
                continue
            gene = row[0]
            if gene in genes_seen:
                duplicate_genes += 1
            else:
                genes_seen.add(gene)

        status = "PASS"
        reason = "OK"
        if duplicate_barcodes or duplicate_genes or width_mismatch or n_genes == 0:
            status = "FAIL"
            reason = "DUPLICATE_OR_WIDTH_ERROR"

        rows_out.append((path.name, n_genes, n_cells, duplicate_genes,
                         duplicate_barcodes, width_mismatch, status, reason))

with out.open("w", encoding="utf-8", newline="") as handle:
    handle.write("file\tgenes\tcells\tduplicate_genes\tduplicate_barcodes\twidth_mismatch_rows\tstatus\treason\n")
    for row in rows_out:
        handle.write("\t".join(map(str, row)) + "\n")

if len(rows_out) != 26:
    raise SystemExit(f"expected 26 CSV matrices, detected {len(rows_out)}")
if any(row[-2] != "PASS" for row in rows_out):
    raise SystemExit("one or more GSE178481 CSV matrices failed")

print(f"GSE178481_CSV_MATRICES={len(rows_out)}")
print("GSE178481_CSV_QC=PASS")
PY

echo
echo "[3/4] OFFICIAL GEO SOFT METADATA"

GEO_DIR="${OUT_ROOT}/geo_soft"
mkdir -p "${GEO_DIR}"

while read -r gse bucket; do
  url="https://ftp.ncbi.nlm.nih.gov/geo/series/${bucket}/${gse}/soft/${gse}_family.soft.gz"
  target="${GEO_DIR}/${gse}_family.soft.gz"

  echo "Downloading ${gse} official SOFT metadata ..."
  curl -fL \
    --retry 3 \
    --retry-delay 5 \
    --connect-timeout 30 \
    --output "${target}.part" \
    "${url}"
  mv -- "${target}.part" "${target}"
  gzip -t "${target}"
done <<'EOF'
GSE178481 GSE178nnn
GSE207493 GSE207nnn
GSE294006 GSE294nnn
GSE294109 GSE294nnn
GSE264586 GSE264nnn
GSE202374 GSE202nnn
EOF

export GEO_DIR

python - <<'PY'
import gzip
import os
from pathlib import Path

root = Path(os.environ["GEO_DIR"])
out_long = Path(os.environ["OUT_ROOT"]) / "geo_sample_metadata_long.tsv"
out_summary = Path(os.environ["OUT_ROOT"]) / "geo_sample_metadata_summary.tsv"

long_rows = []
summary_rows = []

for path in sorted(root.glob("GSE*_family.soft.gz")):
    gse = path.name.split("_", 1)[0]
    current = None
    fields = {}
    occurrences = {}

    def flush():
        if current is None:
            return
        title = " | ".join(fields.get("title", []))
        source = " | ".join(fields.get("source_name_ch1", []))
        organism = " | ".join(fields.get("organism_ch1", []))
        characteristics = " | ".join(fields.get("characteristics_ch1", []))
        summary_rows.append((gse, current, title, source, organism, characteristics))

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                flush()
                current = line.split("=", 1)[1].strip()
                fields = {}
                occurrences = {}
                continue
            if current is None or not line.startswith("!Sample_") or " = " not in line:
                continue
            key, value = line[8:].split(" = ", 1)
            key = key.strip()
            value = value.strip()
            occurrences[key] = occurrences.get(key, 0) + 1
            fields.setdefault(key, []).append(value)
            long_rows.append((gse, current, key, occurrences[key], value))
        flush()

with out_long.open("w", encoding="utf-8", newline="") as handle:
    handle.write("dataset\tgsm\tfield\toccurrence\tvalue\n")
    for row in long_rows:
        handle.write("\t".join(map(str, row)) + "\n")

with out_summary.open("w", encoding="utf-8", newline="") as handle:
    handle.write("dataset\tgsm\ttitle\tsource\torganism\tcharacteristics_all\n")
    for row in summary_rows:
        handle.write("\t".join(str(x).replace("\t", " ") for x in row) + "\n")

counts = {}
for row in summary_rows:
    counts[row[0]] = counts.get(row[0], 0) + 1

for gse in sorted(counts):
    print(f"GEO_SAMPLES_{gse}={counts[gse]}")
PY

echo
echo "[4/4] GSE294006 SPATIAL BUNDLE DEEP AUDIT"

export GSE294006_DIR

python - <<'PY'
import csv
import gzip
import json
import os
import struct
from pathlib import Path

root = Path(os.environ["GSE294006_DIR"])
out = Path(os.environ["OUT_ROOT"]) / "gse294006_spatial_deep_qc.tsv"

def read_lines(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        return [line.rstrip("\r\n").split("\t")[0] for line in handle]

def matrix_header(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line.startswith("%"):
                return tuple(map(int, line.strip().split()))
    raise ValueError(f"missing dimensions: {path}")

def read_json(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)

def png_size(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"invalid PNG: {path}")
    return struct.unpack(">II", header[16:24])

rows = []
for matrix in sorted(root.glob("*matrix.mtx.gz")):
    prefix = matrix.name[:-len("matrix.mtx.gz")]
    sample = prefix.rstrip("._-")
    barcodes = root / f"{prefix}barcodes.tsv.gz"
    features = root / f"{prefix}features.tsv.gz"
    positions = root / f"{prefix}tissue_positions_list.csv.gz"
    scalefactors = root / f"{prefix}scalefactors_json.json.gz"
    images = sorted(root.glob(f"{prefix}*tissue*image*.png.gz"))

    genes, matrix_barcodes, nnz = matrix_header(matrix)
    feature_ids = read_lines(features)
    barcode_ids = read_lines(barcodes)

    with gzip.open(positions, "rt", encoding="utf-8", newline="") as handle:
        position_rows = list(csv.reader(handle))
    if position_rows and position_rows[0] and position_rows[0][0].lower() == "barcode":
        position_rows = position_rows[1:]

    position_map = {row[0]: row for row in position_rows if len(row) >= 6}
    missing_positions = sum(b not in position_map for b in barcode_ids)
    matrix_out_of_tissue = sum(
        b in position_map and position_map[b][1] != "1" for b in barcode_ids
    )
    n_in_tissue = sum(row[1] == "1" for row in position_map.values())

    scale = read_json(scalefactors)
    image = images[0] if images else None
    width, height = png_size(image) if image else (-1, -1)
    image_kind = "hires" if image and "hires" in image.name else "lowres" if image else "missing"

    status = "PASS"
    reason = "OK"
    if (
        genes != len(feature_ids)
        or matrix_barcodes != len(barcode_ids)
        or len(barcode_ids) != len(set(barcode_ids))
        or missing_positions > 0
        or matrix_out_of_tissue > 0
        or not scale
        or image is None
    ):
        status = "FAIL"
        reason = "DIMENSION_POSITION_JSON_OR_IMAGE_ERROR"

    rows.append((sample, genes, matrix_barcodes, nnz, len(position_map),
                 n_in_tissue, missing_positions, matrix_out_of_tissue,
                 image_kind, width, height, ",".join(sorted(scale)), status, reason))

with out.open("w", encoding="utf-8", newline="") as handle:
    handle.write(
        "sample\tgenes\tmatrix_barcodes\tnnz\tposition_rows\tin_tissue_rows\t"
        "missing_positions\tmatrix_out_of_tissue\timage_kind\timage_width\t"
        "image_height\tscalefactor_keys\tstatus\treason\n"
    )
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")

if len(rows) != 3:
    raise SystemExit(f"expected 3 spatial matrices, detected {len(rows)}")
if any(row[-2] != "PASS" for row in rows):
    raise SystemExit("one or more GSE294006 spatial bundles failed")

print("GSE294006_SPATIAL_BUNDLES=3")
print("GSE294006_SPATIAL_DEEP_QC=PASS")
PY

echo
echo "[SUMMARY]"
echo "GSE207493_MODE=${GSE207_MODE}"
echo "GSE207493_FINAL=${GSE207_FINAL}"
echo "OUT_ROOT=${OUT_ROOT}"
echo "END=$(date --iso-8601=seconds)"
echo "PHASE3_COMPLETE"
