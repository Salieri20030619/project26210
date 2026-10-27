#!/usr/bin/env bash
set -uo pipefail

# Project 26210: phase-1 read-only input audit.
# Run on sgs10. This script does not extract archives or modify project files.

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
PERSONAL_ROOT="/mnt/BioDisk/wangmingyuan"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${PERSONAL_ROOT}/project26210_external_joint_v1/audit_${RUN_ID}"
LIST_DIR="${OUT_ROOT}/archive_lists"

mkdir -p "${LIST_DIR}"

REPORT="${OUT_ROOT}/audit_report.txt"
ARCHIVE_TSV="${OUT_ROOT}/archive_inventory.tsv"
STRUCTURE_TSV="${OUT_ROOT}/archive_structure.tsv"
PATHOLOGY_TSV="${OUT_ROOT}/pathology_json_audit.tsv"
FOUND_TSV="${OUT_ROOT}/located_inputs.tsv"

exec > >(tee -a "${REPORT}") 2>&1

tsv_safe() {
    printf '%s' "$1" | tr '\t\r\n' '   '
}

human_size() {
    numfmt --to=iec-i --suffix=B "$1" 2>/dev/null || printf '%sB' "$1"
}

archive_kind() {
    case "$1" in
        *.tar.gz|*.tgz) printf 'tar.gz' ;;
        *.tar) printf 'tar' ;;
        *.zip) printf 'zip' ;;
        *.gz) printf 'gz' ;;
        *) printf 'other' ;;
    esac
}

list_archive() {
    local path="$1"
    local list_file="$2"
    local kind
    kind="$(archive_kind "${path}")"
    case "${kind}" in
        tar.gz) tar -tzf "${path}" > "${list_file}" ;;
        tar) tar -tf "${path}" > "${list_file}" ;;
        zip) unzip -Z1 "${path}" > "${list_file}" ;;
        gz) gzip -t "${path}" && printf '%s\n' "$(basename "${path}" .gz)" > "${list_file}" ;;
        *) return 2 ;;
    esac
}

archive_names=(
    "GSE178481_RAW.tar"
    "GSE207493_scRNA-seq.tar.gz"
    "GSE294006_RAW.tar"
    "GSE294109_RAW.tar"
    "GSE264586_RAW.tar"
    "GSE202374_RAW.tar"
)

declare -A dataset_of=(
    [GSE178481_RAW.tar]="GSE178481"
    [GSE207493_scRNA-seq.tar.gz]="GSE207493"
    [GSE294006_RAW.tar]="GSE294006"
    [GSE294109_RAW.tar]="GSE294109"
    [GSE264586_RAW.tar]="GSE264586"
    [GSE202374_RAW.tar]="GSE202374"
)

printf 'PROJECT26210 PHASE-1 INPUT AUDIT\n'
printf 'run_id=%s\n' "${RUN_ID}"
printf 'host=%s\n' "$(hostname)"
printf 'user=%s\n' "$(whoami)"
printf 'date=%s\n' "$(date --iso-8601=seconds)"
printf 'project_root=%s\n' "${PROJECT_ROOT}"
printf 'output_root=%s\n' "${OUT_ROOT}"
printf '\n[DISK]\n'
df -h "${PERSONAL_ROOT}" "${PROJECT_ROOT}" 2>&1 || true

if [[ "$(hostname -s)" != "sgs10" ]]; then
    printf '\nWARNING: expected host sgs10, current host is %s.\n' "$(hostname -s)"
fi

printf 'dataset\tinput_name\tpath\tbytes\tmtime\n' > "${FOUND_TSV}"
printf 'dataset\tpath\tbytes\thuman_size\tmtime\tarchive_kind\tintegrity\tmember_count\tzero_member_count\tpartial_marker_count\tlist_file\n' > "${ARCHIVE_TSV}"
printf 'dataset\tpath\textension\tmember_count\n' > "${STRUCTURE_TSV}"

printf '\n[LOCATE_AND_VALIDATE_ARCHIVES]\n'
for input_name in "${archive_names[@]}"; do
    dataset="${dataset_of[${input_name}]}"
    mapfile -t matches < <(
        find "${PERSONAL_ROOT}" "${PROJECT_ROOT}" -xdev \
            -type f -name "${input_name}" \
            -not -path "${PERSONAL_ROOT}/project26210_external_joint_v1/*" \
            -print 2>/dev/null | sort -u
    )

    if (( ${#matches[@]} == 0 )); then
        printf '%s\t%s\tMISSING\t0\tNA\n' "${dataset}" "${input_name}" >> "${FOUND_TSV}"
        printf '%s: MISSING (%s)\n' "${dataset}" "${input_name}"
        continue
    fi

    if (( ${#matches[@]} > 1 )); then
        printf '%s: WARNING multiple copies found: %d\n' "${dataset}" "${#matches[@]}"
    fi

    copy_index=0
    for path in "${matches[@]}"; do
        copy_index=$((copy_index + 1))
        bytes="$(stat -c '%s' "${path}" 2>/dev/null || printf '0')"
        mtime="$(stat -c '%y' "${path}" 2>/dev/null || printf 'NA')"
        kind="$(archive_kind "${path}")"
        safe_base="${dataset}_copy${copy_index}_$(basename "${path}" | tr -c 'A-Za-z0-9._-' '_')"
        list_file="${LIST_DIR}/${safe_base}.members.txt"

        printf '%s\t%s\t%s\t%s\t%s\n' \
            "${dataset}" "${input_name}" "$(tsv_safe "${path}")" "${bytes}" "$(tsv_safe "${mtime}")" \
            >> "${FOUND_TSV}"

        printf 'Validating %s (%s, %s) ...\n' "${path}" "$(human_size "${bytes}")" "${kind}"
        if list_archive "${path}" "${list_file}"; then
            integrity="PASS"
            member_count="$(wc -l < "${list_file}" | awk '{print $1}')"
            zero_member_count="0"
            if [[ "${kind}" == "tar" || "${kind}" == "tar.gz" ]]; then
                if [[ "${kind}" == "tar.gz" ]]; then
                    zero_member_count="$(tar -tvzf "${path}" 2>/dev/null | awk '$3==0 && $1 !~ /^d/ {n++} END{print n+0}')"
                else
                    zero_member_count="$(tar -tvf "${path}" 2>/dev/null | awk '$3==0 && $1 !~ /^d/ {n++} END{print n+0}')"
                fi
            fi
            partial_marker_count="$(awk 'BEGIN{IGNORECASE=1} /(^|\/)([^\/]*\.(part|partial|aria2|tmp|download))$/ {n++} END{print n+0}' "${list_file}")"

            awk -v ds="${dataset}" -v p="${path}" '
                function ext(path, base,n,a) {
                    n=split(path,a,"/"); base=a[n]
                    if (base !~ /\./) return "[no_ext]"
                    sub(/^.*\./,"",base); return tolower(base)
                }
                { c[ext($0)]++ }
                END { for (e in c) print ds "\t" p "\t" e "\t" c[e] }
            ' "${list_file}" | sort -t $'\t' -k3,3 >> "${STRUCTURE_TSV}"

            printf '%s: PASS; members=%s; zero-size-members=%s; partial-markers=%s\n' \
                "${dataset}" "${member_count}" "${zero_member_count}" "${partial_marker_count}"
        else
            integrity="FAIL"
            member_count="0"
            zero_member_count="NA"
            partial_marker_count="NA"
            : > "${list_file}"
            printf '%s: FAIL archive listing/integrity test\n' "${dataset}"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${dataset}" "$(tsv_safe "${path}")" "${bytes}" "$(human_size "${bytes}")" \
            "$(tsv_safe "${mtime}")" "${kind}" "${integrity}" "${member_count}" \
            "${zero_member_count}" "${partial_marker_count}" "$(tsv_safe "${list_file}")" \
            >> "${ARCHIVE_TSV}"

        if [[ "${integrity}" == "PASS" ]]; then
            printf '%s candidate metadata/matrix files:\n' "${dataset}"
            awk 'BEGIN{IGNORECASE=1}
                /(^|\/)([^\/]*(meta|annot|clinical|sample|patient|matrix|barcodes|features|genes|spatial|scalefactor|tissue_positions)[^\/]*)$/ {print}
            ' "${list_file}" | head -n 40
        fi
    done
done

printf '\n[PATHOLOGY_JSON_AUDIT]\n'
printf 'sample\tpath\tbytes\tmtime\tjson_valid\ttop_type\ttop_keys\tpolygon_like_arrays\tcoordinate_pairs\tfinite_pairs\tmin_x\tmax_x\tmin_y\tmax_y\tnote\n' > "${PATHOLOGY_TSV}"

PYTHON_BIN=""
if command -v conda >/dev/null 2>&1; then
    PYTHON_BIN="$(conda run -n batch_bench which python 2>/dev/null || true)"
fi
if [[ -z "${PYTHON_BIN}" ]] && command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3)"
fi

for json_name in GLMF1_pathology_polygons.json GLMF2_pathology_polygons.json; do
    sample="${json_name%%_*}"
    mapfile -t json_matches < <(find "${PROJECT_ROOT}" -xdev -type f -name "${json_name}" -print 2>/dev/null | sort -u)
    if (( ${#json_matches[@]} == 0 )); then
        printf '%s\tMISSING\t0\tNA\tNO\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\tNA\tfile_not_found\n' "${sample}" >> "${PATHOLOGY_TSV}"
        printf '%s: MISSING\n' "${json_name}"
        continue
    fi
    for json_path in "${json_matches[@]}"; do
        if [[ -z "${PYTHON_BIN}" ]]; then
            printf '%s\t%s\t%s\t%s\tUNKNOWN\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tpython_not_found\n' \
                "${sample}" "$(tsv_safe "${json_path}")" "$(stat -c '%s' "${json_path}")" \
                "$(tsv_safe "$(stat -c '%y' "${json_path}")")" >> "${PATHOLOGY_TSV}"
            continue
        fi
        "${PYTHON_BIN}" - "${sample}" "${json_path}" >> "${PATHOLOGY_TSV}" <<'PY'
import json, math, os, sys

sample, path = sys.argv[1], sys.argv[2]
st = os.stat(path)

def clean(x):
    return str(x).replace("\t", " ").replace("\n", " ").replace("\r", " ")

try:
    with open(path, "r", encoding="utf-8-sig") as fh:
        obj = json.load(fh)
except Exception as exc:
    row = [sample, path, st.st_size, st.st_mtime, "NO", "NA", "NA", 0, 0, 0,
           "NA", "NA", "NA", "NA", f"parse_error:{type(exc).__name__}:{exc}"]
    print("\t".join(clean(x) for x in row))
    raise SystemExit(0)

pairs = []
polygon_like = 0

def walk(x):
    global polygon_like
    if isinstance(x, dict):
        for v in x.values():
            walk(v)
    elif isinstance(x, list):
        local_pairs = []
        for item in x:
            if (isinstance(item, (list, tuple)) and len(item) >= 2 and
                    isinstance(item[0], (int, float)) and isinstance(item[1], (int, float))):
                local_pairs.append((item[0], item[1]))
        if len(local_pairs) >= 3:
            polygon_like += 1
        pairs.extend(local_pairs)
        for item in x:
            walk(item)

walk(obj)
finite = [(float(x), float(y)) for x, y in pairs
          if math.isfinite(float(x)) and math.isfinite(float(y))]
keys = ",".join(map(str, obj.keys())) if isinstance(obj, dict) else "NA"
if finite:
    xs, ys = zip(*finite)
    bounds = [min(xs), max(xs), min(ys), max(ys)]
else:
    bounds = ["NA"] * 4
note = "AI_assisted_regions;requires_image_overlay_QC_and_pathologist_limitation_statement"
row = [sample, path, st.st_size, st.st_mtime, "YES", type(obj).__name__, keys,
       polygon_like, len(pairs), len(finite), *bounds, note]
print("\t".join(clean(x) for x in row))
PY
        printf '%s: JSON parsed; see %s\n' "${json_name}" "${PATHOLOGY_TSV}"
    done
done

printf '\n[ACTIVE_RELEVANT_PROCESSES]\n'
ps -u "$(whoami)" -o pid=,etime=,%cpu=,%mem=,cmd= | \
    awk 'BEGIN{IGNORECASE=1} /GSE178481|GSE207493|GSE294006|GSE294109|GSE264586|GSE202374|SRZ190804|cell2loc|pathology|signed.?distance/ {print}' || true

printf '\n[OUTPUTS]\n'
printf '%s\t%s\t%s\n' "${REPORT}" "$(stat -c '%s' "${REPORT}" 2>/dev/null || printf '0')" "LIVE_LOG_HASH_AFTER_COMPLETION"
for output in "${FOUND_TSV}" "${ARCHIVE_TSV}" "${STRUCTURE_TSV}" "${PATHOLOGY_TSV}"; do
    printf '%s\t%s\t%s\n' "${output}" "$(stat -c '%s' "${output}" 2>/dev/null || printf '0')" \
        "$(sha256sum "${output}" 2>/dev/null | awk '{print $1}')"
done

printf '\nAUDIT_COMPLETE\n'
printf 'OUT_ROOT=%s\n' "${OUT_ROOT}"
printf 'REPORT=%s\n' "${REPORT}"
printf 'ARCHIVE_TSV=%s\n' "${ARCHIVE_TSV}"
printf 'STRUCTURE_TSV=%s\n' "${STRUCTURE_TSV}"
printf 'PATHOLOGY_TSV=%s\n' "${PATHOLOGY_TSV}"
printf '\nNo archives were extracted and no project files were modified.\n'
