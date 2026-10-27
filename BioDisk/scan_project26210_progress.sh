#!/usr/bin/env bash
# Read-only evidence scan for project 26210 on sgs10.
# The script never writes inside PROJECT_ROOT. Reports are written to OUT_DIR.

set -uo pipefail

PROJECT_ROOT="${1:-/mnt/BioDisk/ludingyu/project26210}"
OUT_DIR="${2:-/mnt/BioDisk/wangmingyuan/project26210_audit}"
ARCHIVE_MODE="${ARCHIVE_MODE:-quick}"       # none | quick | full
QUICK_ARCHIVE_MAX_BYTES="${QUICK_ARCHIVE_MAX_BYTES:-1073741824}"  # 1 GiB
ARCHIVE_TIMEOUT_SEC="${ARCHIVE_TIMEOUT_SEC:-300}"
MAX_DEPTH="${MAX_DEPTH:-12}"
TEXT_MAX_BYTES="${TEXT_MAX_BYTES:-20971520}" # 20 MiB
NOW="$(date +%Y%m%d_%H%M%S)"

if [[ ! -d "$PROJECT_ROOT" ]]; then
  echo "ERROR: project root does not exist: $PROJECT_ROOT" >&2
  exit 2
fi

PROJECT_ROOT="$(cd "$PROJECT_ROOT" && pwd -P)"
case "$OUT_DIR/" in
  "$PROJECT_ROOT"/*)
    echo "ERROR: OUT_DIR must be outside PROJECT_ROOT to keep the project read-only." >&2
    exit 2
    ;;
esac

mkdir -p "$OUT_DIR" || {
  echo "ERROR: cannot create output directory: $OUT_DIR" >&2
  exit 2
}
OUT_DIR="$(cd "$OUT_DIR" && pwd -P)"

REPORT="$OUT_DIR/project26210_progress_scan_${NOW}.md"
FILE_INDEX="$OUT_DIR/project26210_file_index_${NOW}.tsv"
DATASET_INDEX="$OUT_DIR/project26210_dataset_index_${NOW}.tsv"
COMMAND_LOG="$OUT_DIR/project26210_scan_commands_${NOW}.log"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/project26210_scan.XXXXXX")"
trap 'rm -rf -- "$TMP_DIR"' EXIT

exec 3>>"$COMMAND_LOG"

say() { printf '%s\n' "$*" | tee -a "$REPORT" >/dev/null; }
section() {
  say ""
  say "## $*"
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}
subsection() { say ""; say "### $*"; }
run_report() {
  printf '$ %q' "$1" >&3
  shift
  printf ' %q' "$@" >&3
  printf '\n' >&3
  "$@" 2>&1 | sed 's/^/    /' | tee -a "$REPORT" >/dev/null
  local rc=${PIPESTATUS[0]}
  printf 'exit_code=%s\n' "$rc" >&3
  return "$rc"
}
have() { command -v "$1" >/dev/null 2>&1; }
human_size() {
  local n="${1:-0}"
  if have numfmt; then numfmt --to=iec-i --suffix=B "$n" 2>/dev/null || printf '%s bytes' "$n"
  else printf '%s bytes' "$n"
  fi
}
stat_line() {
  local p="$1"
  if [[ -e "$p" || -L "$p" ]]; then
    stat -c '%n | type=%F | bytes=%s | mtime=%y | owner=%U:%G | mode=%A' -- "$p" 2>/dev/null || true
  else
    printf '%s | MISSING\n' "$p"
  fi
}
du_limited() {
  local p="$1"
  if have timeout; then timeout 120 du -sh --apparent-size -- "$p" 2>/dev/null | awk '{print $1}'
  else du -sh --apparent-size -- "$p" 2>/dev/null | awk '{print $1}'
  fi
}
count_files_limited() {
  local p="$1"
  if [[ -f "$p" ]]; then printf '1\n'; return; fi
  if have timeout; then timeout 120 find "$p" -type f -printf '.' 2>/dev/null | wc -c
  else find "$p" -type f -printf '.' 2>/dev/null | wc -c
  fi
}
first_n() { sed -n "1,${1}p"; }

printf '# Project 26210 read-only progress scan\n' >"$REPORT"
say ""
say "- Scan time: $(date --iso-8601=seconds)"
say "- Host: $(hostname 2>/dev/null || printf unknown)"
say "- User: $(whoami 2>/dev/null || printf unknown)"
say "- Project root: $PROJECT_ROOT"
say "- Output directory: $OUT_DIR"
say "- Archive mode: $ARCHIVE_MODE"
say "- Scope: metadata, lightweight text/log evidence, process state and optional archive integrity checks"
say "- Safety: no file under the project root is created, edited, renamed or deleted"

section "1. Runtime and storage"
say '```text'
{
  echo "date=$(date --iso-8601=seconds)"
  echo "hostname=$(hostname 2>/dev/null || true)"
  echo "user=$(whoami 2>/dev/null || true)"
  echo "pwd=$(pwd -P)"
  echo "project_root=$PROJECT_ROOT"
  echo "output_dir=$OUT_DIR"
  echo "shell=${SHELL:-unknown}"
  echo "conda_default_env=${CONDA_DEFAULT_ENV:-unset}"
  echo "http_proxy=${http_proxy:+set}"
  echo "https_proxy=${https_proxy:+set}"
  echo "HTTP_PROXY=${HTTP_PROXY:+set}"
  echo "HTTPS_PROXY=${HTTPS_PROXY:+set}"
  echo "NO_PROXY=${NO_PROXY:+set}"
  echo "rg=$(command -v rg 2>/dev/null || echo unavailable)"
  echo "h5ls=$(command -v h5ls 2>/dev/null || echo unavailable)"
  echo "timeout=$(command -v timeout 2>/dev/null || echo unavailable)"
  echo
  echo "[disk space]"
  df -h "$PROJECT_ROOT" "$OUT_DIR" 2>&1
  echo
  echo "[inode space]"
  df -ih "$PROJECT_ROOT" "$OUT_DIR" 2>&1
} | tee -a "$REPORT" >&3
say '```'

section "2. Project file index"
printf 'path\tbytes\tmtime\ttype\n' >"$FILE_INDEX"
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" \
  \( -type f -o -type l \) \
  -printf '%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS%Tz\t%y\n' 2>>"$COMMAND_LOG" \
  | sort >>"$FILE_INDEX"
INDEX_ROWS="$(awk 'END{print NR-1}' "$FILE_INDEX")"
say "Indexed files/symlinks (max depth $MAX_DEPTH): $INDEX_ROWS"
say "Full index: $FILE_INDEX"
say ""
say 'Top 30 largest indexed files:'
say '```text'
{ tail -n +2 "$FILE_INDEX" | sort -t$'\t' -k2,2nr | head -n 30 | awk -F'\t' '{printf "%12s  %s  %s\n", $2, $3, $1}'; } | tee -a "$REPORT" >&3
say '```'

subsection "Project directory overview"
say '```text'
find "$PROJECT_ROOT" -maxdepth 3 -type d -printf '%p\n' 2>/dev/null | sort | first_n 300 | tee -a "$REPORT" >&3
say '```'

section "3. Governing files and handoffs"
GOVERNING=(
  "$PROJECT_ROOT/AGENTS.md"
  "$PROJECT_ROOT/PROJECT_STATUS.md"
  "$PROJECT_ROOT/RESULTS_MANIFEST.tsv"
)
say '```text'
for p in "${GOVERNING[@]}"; do stat_line "$p"; done | tee -a "$REPORT" >&3
say '```'

for p in "${GOVERNING[@]}"; do
  if [[ -f "$p" ]]; then
    subsection "$(basename "$p") selected evidence"
    say '```text'
    {
      echo "[first 80 lines]"
      sed -n '1,80p' "$p" 2>/dev/null
      echo
      echo "[matched status/evidence lines; maximum 250]"
      grep -Ein 'EGA|EGAD00001008781|TJ-RCC|PRJCA014547|24599295|14620362|SCP3155|GSE178481|GSE207493|GSE294006|GSE294109|GSE264586|GSE202374|SCP1288|SRZ190804|GLMF1|GLMF2|cell2location|BayesSpace|SPARK|Hotspot|MISTy|COMMOT|tumor.?NAT|34,?609|9,?999|912|38/38|19条|TYMP|SPP1|APOE|C1QC|APP|CD74|TREM2' "$p" 2>/dev/null | head -n 250 || true
    } | tee -a "$REPORT" >&3
    say '```'
  fi
done

subsection "All AGENTS.md files within the project"
say '```text'
find "$PROJECT_ROOT" -maxdepth 8 -type f -name AGENTS.md \
  -printf '%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null \
  | sort | head -n 100 | tee "$TMP_DIR/agents_files.tsv" | tee -a "$REPORT" >&3
say '```'
if [[ -s "$TMP_DIR/agents_files.tsv" ]]; then
  cut -f3- "$TMP_DIR/agents_files.tsv" | head -n 20 | while IFS= read -r p; do
    subsection "Rules excerpt: $p"
    say '```text'
    sed -n '1,120p' "$p" 2>/dev/null | tee -a "$REPORT" >&3
    say '```'
  done
fi

subsection "Latest handoff files"
say '```text'
find "$PROJECT_ROOT" -maxdepth 5 -type f \
  \( -name '26210_CHAT_HANDOFF_*.md' -o -name 'CODEX_HANDOFF.md' \) \
  -printf '%T@\t%TY-%Tm-%Td %TH:%TM:%TS\t%s\t%p\n' 2>/dev/null \
  | sort -nr | head -n 3 | cut -f2- | tee -a "$REPORT" >&3
say '```'

section "4. Manifest structural audit"
MANIFEST="$PROJECT_ROOT/RESULTS_MANIFEST.tsv"
if [[ -f "$MANIFEST" ]]; then
  say '```text'
  awk -F'\t' '
    NR==1 {expected=NF; print "header_columns=" NF; print "header=" $0; next}
    NF!=expected {bad++; print "BAD_COLUMN_COUNT line=" NR " fields=" NF}
    END {print "data_rows=" NR-1; print "bad_column_rows=" bad+0}
  ' "$MANIFEST" | tee -a "$REPORT" >&3

  awk -F'\t' '
    NR==1 {
      for (i=1;i<=NF;i++) {
        h=tolower($i)
        if (h=="result_id" || h=="id") idcol=i
        if (h=="qc_status" || h=="status") statuscol=i
        if (h=="path" || h=="result_path" || h=="output_path" || h=="absolute_path") pathcol=i
      }
      next
    }
    idcol && $idcol!="" {seen[$idcol]++}
    statuscol && $statuscol!="" {status[$statuscol]++}
    pathcol && $pathcol!="" {paths[$pathcol]++}
    END {
      print "id_column=" (idcol?idcol:"not_detected")
      print "status_column=" (statuscol?statuscol:"not_detected")
      print "path_column=" (pathcol?pathcol:"not_detected")
      print "[duplicate result IDs]"
      n=0; for (k in seen) if (seen[k]>1) {print k, seen[k]; n++}
      if (!n) print "none_detected"
      print "[status distribution]"
      for (k in status) print k, status[k]
      if (pathcol) {
        print "[manifest paths]"
        for (k in paths) print k
      }
    }
  ' "$MANIFEST" >"$TMP_DIR/manifest_audit.txt"
  cat "$TMP_DIR/manifest_audit.txt" | tee -a "$REPORT" >&3
  say '```'

  subsection "Manifest path existence"
  say '```text'
  awk -F'\t' '
    NR==1 {for(i=1;i<=NF;i++){h=tolower($i); if(h=="path"||h=="result_path"||h=="output_path"||h=="absolute_path") p=i}; next}
    p && $p!="" {print $p}
  ' "$MANIFEST" | sort -u >"$TMP_DIR/manifest_paths.txt"
  if [[ -s "$TMP_DIR/manifest_paths.txt" ]]; then
    while IFS= read -r p; do
      if [[ "$p" = /* ]]; then
        if [[ -e "$p" || -L "$p" ]]; then printf 'EXISTS\t%s\n' "$p"; else printf 'MISSING\t%s\n' "$p"; fi
      else
        printf 'NON_ABSOLUTE\t%s\n' "$p"
      fi
    done <"$TMP_DIR/manifest_paths.txt" | tee -a "$REPORT" >&3
  else
    echo 'No recognized path column or no paths present.' | tee -a "$REPORT" >&3
  fi
  say '```'
else
  say "RESULTS_MANIFEST.tsv not found."
fi

section "5. Expected reports and formal objects"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f \
  \( -iname 'part4_compact_report_v13r1.pdf' \
     -o -iname 'project26210_external_singlecell_dataset_client_summary_20260805.docx' \
     -o -iname '*.h5ad' -o -iname '*.rds' -o -iname '*.pdf' -o -iname '*.docx' \) \
  -printf '%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null \
  | sort -t$'\t' -k2,2r | head -n 250 | tee -a "$REPORT" >&3
say '```'

section "6. Part 4 module evidence by path/name"
declare -a MODULES=(
  'GLMF1|GLMF1'
  'GLMF2|GLMF2'
  'Reference/regression|reference|regression|signature'
  'cell2location|cell2location|cell2loc|cell_abundance'
  'TJ-RCC/12-patient|TJ.?RCC|PRJCA014547|24599295|patient|patients'
  'BayesSpace|BayesSpace'
  'SPARK-X/Hotspot|SPARK.?X|Hotspot'
  'MISTy|MISTy'
  'COMMOT|COMMOT|commot'
  'TYMP-SPP1/MTA|TYMP|SPP1|MTA'
  'Vascular niche|vascular|endothelial|angiogen|vessel'
  'T-cell balance|T.?cell|cytotoxic|exhaust|balance'
  'Directed communication|directed|directional|19.?axis|38.?38|912|9999'
  'POSTN CAF-APOE TAM|POSTN|APOE|CAF|TAM'
  'Tumor-NAT annotation/distance|tumor.?NAT|signed.?distance|polygon|annotation.*json|pathology'
  'Niche/ecotype|niche|ecotype|microenvironment'
)

for item in "${MODULES[@]}"; do
  label="${item%%|*}"
  regex="${item#*|}"
  subsection "$label"
  say '```text'
  { tail -n +2 "$FILE_INDEX" | cut -f1 | grep -Ei "$regex" | head -n 120 || true; } | tee "$TMP_DIR/module_hits.txt" | tee -a "$REPORT" >&3
  if [[ ! -s "$TMP_DIR/module_hits.txt" ]]; then
    echo 'NO_PATH_HIT_WITHIN_SCAN_DEPTH' | tee -a "$REPORT" >&3
  fi
  say '```'
done

subsection "Numeric/result claims found in lightweight text files"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -size -"${TEXT_MAX_BYTES}"c \
  \( -iname '*.md' -o -iname '*.txt' -o -iname '*.tsv' -o -iname '*.csv' \
     -o -iname '*.log' -o -iname '*.out' -o -iname '*.err' -o -iname '*.json' \
     -o -iname '*.yaml' -o -iname '*.yml' -o -iname '*.R' -o -iname '*.py' -o -iname '*.sh' \) \
  -print0 2>/dev/null \
  | xargs -0 -r grep -IEn '34,?609|9,?999|38/38|912([^0-9]|$)|19条|19 (axes|axis)|12/12|14类|14 (cell|reference)|tumor.?NAT|signed.?distance' 2>/dev/null \
  | head -n 500 | tee -a "$REPORT" >&3 || true
say '```'

section "7. External single-cell/spatial dataset inventory"
SEARCH_ROOTS=("$PROJECT_ROOT")
if [[ -d /mnt/BioDisk/wangmingyuan && /mnt/BioDisk/wangmingyuan != "$PROJECT_ROOT" ]]; then
  SEARCH_ROOTS+=("/mnt/BioDisk/wangmingyuan")
fi
DATASETS=(GSE178481 GSE207493 GSE294006 GSE294109 GSE264586 GSE202374 SCP1288 SRZ190804 PRJCA014547 24599295 14620362 SCP3155 EGAD00001008781)
printf 'dataset\tmatch_type\tpath\tbytes_if_file\tmtime\n' >"$DATASET_INDEX"

for ds in "${DATASETS[@]}"; do
  : >"$TMP_DIR/${ds}.matches"
  for root in "${SEARCH_ROOTS[@]}"; do
    find "$root" -maxdepth 7 \( -type f -o -type d -o -type l \) \
      -iname "*${ds}*" -printf '%y\t%p\t%s\t%TY-%Tm-%TdT%TH:%TM:%TS%Tz\n' 2>/dev/null \
      >>"$TMP_DIR/${ds}.matches"
  done
  sort -u "$TMP_DIR/${ds}.matches" -o "$TMP_DIR/${ds}.matches"
  while IFS=$'\t' read -r typ p bytes mtime; do
    [[ -n "${p:-}" ]] || continue
    printf '%s\t%s\t%s\t%s\t%s\n' "$ds" "$typ" "$p" "$bytes" "$mtime" >>"$DATASET_INDEX"
  done <"$TMP_DIR/${ds}.matches"

  subsection "$ds"
  say '```text'
  if [[ ! -s "$TMP_DIR/${ds}.matches" ]]; then
    echo 'NO_MATCH_IN_SEARCH_ROOTS' | tee -a "$REPORT" >&3
  else
    echo '[matched paths; maximum 80]' | tee -a "$REPORT" >&3
    head -n 80 "$TMP_DIR/${ds}.matches" | tee -a "$REPORT" >&3
    cut -f2 "$TMP_DIR/${ds}.matches" | while IFS= read -r p; do
      printf '%s\n' "$p"
    done | sort -u | head -n 12 >"$TMP_DIR/${ds}.candidate_dirs"
    echo '[candidate directory summaries; maximum 12]' | tee -a "$REPORT" >&3
    while IFS= read -r d; do
      [[ -e "$d" ]] || continue
      size="$(du_limited "$d")"; [[ -n "$size" ]] || size='timeout_or_error'
      count="$(count_files_limited "$d")"; [[ -n "$count" ]] || count='timeout_or_error'
      partial="$(find "$d" -maxdepth 5 -type f \( -name '*.part' -o -name '*.aria2' -o -name '*.tmp' -o -name '*.download' -o -name '*.crdownload' \) -printf '.' 2>/dev/null | wc -c)"
      zero="$(find "$d" -maxdepth 5 -type f -size 0 -printf '.' 2>/dev/null | wc -c)"
      printf 'dir=%s | apparent_size=%s | files=%s | partial_markers=%s | zero_byte=%s\n' "$d" "$size" "$count" "$partial" "$zero" | tee -a "$REPORT" >&3
      if [[ -d "$d" ]]; then
        find "$d" -maxdepth 5 -type f \
          \( -iname '*matrix*' -o -iname '*barcode*' -o -iname '*feature*' \
             -o -iname '*tissue_position*' -o -iname '*scalefactor*' -o -iname '*spatial*' \
             -o -iname '*.tif' -o -iname '*.tiff' -o -iname '*.png' -o -iname '*.jpg' \
             -o -iname '*.h5ad' -o -iname '*metadata*' -o -iname '*sample*' -o -iname '*patient*' \) \
          -printf '  component=%s bytes | %p\n' 2>/dev/null | head -n 80 | tee -a "$REPORT" >&3
      fi
    done <"$TMP_DIR/${ds}.candidate_dirs"
  fi
  say '```'
done

section "8. Download completeness and obvious filesystem anomalies"
subsection "Partial-download markers"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f \
  \( -name '*.part' -o -name '*.aria2' -o -name '*.tmp' -o -name '*.download' -o -name '*.crdownload' \) \
  -printf '%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null | sort -k2,3r | head -n 500 | tee -a "$REPORT" >&3
say '```'

subsection "Zero-byte files"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -size 0 \
  -printf '%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null | sort -r | head -n 500 | tee -a "$REPORT" >&3
say '```'

subsection "Broken symlinks"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -xtype l -printf '%p -> %l\n' 2>/dev/null | head -n 500 | tee -a "$REPORT" >&3
say '```'

subsection "Archive candidates and optional integrity tests"
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f \
  \( -iname '*.zip' -o -iname '*.gz' -o -iname '*.tgz' -o -iname '*.tar' -o -iname '*.tar.gz' \) \
  -printf '%s\t%p\n' 2>/dev/null | sort -nr >"$TMP_DIR/archives.tsv"
say "Archive candidates: $(wc -l <"$TMP_DIR/archives.tsv")"
say '```text'
head -n 300 "$TMP_DIR/archives.tsv" | tee -a "$REPORT" >&3
say '```'

if [[ "$ARCHIVE_MODE" == 'none' ]]; then
  say 'Integrity tests skipped because ARCHIVE_MODE=none.'
else
  say '```text'
  while IFS=$'\t' read -r bytes p; do
    [[ -n "${p:-}" ]] || continue
    if [[ "$ARCHIVE_MODE" == 'quick' && "$bytes" -gt "$QUICK_ARCHIVE_MAX_BYTES" ]]; then
      printf 'SKIPPED_TOO_LARGE\t%s\t%s\n' "$bytes" "$p"
      continue
    fi
    rc=99
    if [[ "$p" =~ \.zip$ ]] && have unzip; then
      if have timeout; then timeout "$ARCHIVE_TIMEOUT_SEC" unzip -tqq "$p" >/dev/null 2>&1; rc=$?; else unzip -tqq "$p" >/dev/null 2>&1; rc=$?; fi
    elif [[ "$p" =~ \.(gz|tgz)$ ]] && have gzip; then
      if have timeout; then timeout "$ARCHIVE_TIMEOUT_SEC" gzip -t "$p" >/dev/null 2>&1; rc=$?; else gzip -t "$p" >/dev/null 2>&1; rc=$?; fi
    elif [[ "$p" =~ \.tar$ ]] && have tar; then
      if have timeout; then timeout "$ARCHIVE_TIMEOUT_SEC" tar -tf "$p" >/dev/null 2>&1; rc=$?; else tar -tf "$p" >/dev/null 2>&1; rc=$?; fi
    else
      printf 'NO_TEST_TOOL_OR_UNSUPPORTED\t%s\t%s\n' "$bytes" "$p"
      continue
    fi
    if [[ "$rc" -eq 0 ]]; then state=PASSED; elif [[ "$rc" -eq 124 ]]; then state=TIMEOUT; else state=FAILED; fi
    printf '%s\trc=%s\t%s\t%s\n' "$state" "$rc" "$bytes" "$p"
  done <"$TMP_DIR/archives.tsv" | tee -a "$REPORT" >&3
  say '```'
fi

section "9. HDF5/H5AD lightweight metadata"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f \( -iname '*.h5' -o -iname '*.h5ad' -o -iname '*.hdf5' \) \
  -printf '%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null | sort -nr | head -n 120 | tee -a "$REPORT" >&3
say '```'
if have h5ls; then
  say 'Top-level h5ls output (maximum 40 files; 30-second timeout each):'
  say '```text'
  find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f \( -iname '*.h5' -o -iname '*.h5ad' -o -iname '*.hdf5' \) -print0 2>/dev/null \
    | head -z -n 40 \
    | while IFS= read -r -d '' p; do
        echo "[$p]"
        if have timeout; then timeout 30 h5ls "$p" 2>&1 | head -n 80; else h5ls "$p" 2>&1 | head -n 80; fi
      done | tee -a "$REPORT" >&3
  say '```'
else
  say 'h5ls is unavailable; no environment was modified.'
fi

section "10. Signature and joint-analysis preparation evidence"
declare -a SIGNATURES=(
  'TYMP/F1 malignant|TYMP|F1.*malignant|malignant.*F1'
  'SPP1 myeloid|SPP1.*myeloid|SPP1.*TAM|SPP1.*macro'
  'APOE/C1QC TAM|APOE|C1QC|TREM2'
  'Hypoxia/vascular|hypoxia|VEGFA|vascular|angiogen'
  'EMT/invasion|EMT|invasion|epithelial.?mesenchymal'
  'T-cell programs|cytotoxic|exhaust|T.?cell.*balance'
  'Metabolic signatures|metabolic|metabolism|pathway.*signature'
  'Pseudobulk/paired model|pseudobulk|paired|patient.?level|mixed.?effect'
)
for item in "${SIGNATURES[@]}"; do
  label="${item%%|*}"; regex="${item#*|}"
  subsection "$label"
  say '```text'
  {
    tail -n +2 "$FILE_INDEX" | cut -f1 | grep -Ei "$regex" | head -n 100 || true
    find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -size -"${TEXT_MAX_BYTES}"c \
      \( -iname '*.md' -o -iname '*.txt' -o -iname '*.tsv' -o -iname '*.csv' -o -iname '*.json' -o -iname '*.R' -o -iname '*.py' -o -iname '*.sh' \) \
      -print0 2>/dev/null | xargs -0 -r grep -IEn "$regex" 2>/dev/null | head -n 120 || true
  } | awk '!seen[$0]++' | tee -a "$REPORT" >&3
  say '```'
done

section "11. Current processes"
say 'Only the current Linux user is queried; no process is stopped or signalled.'
say '```text'
ps -u "$(id -u)" -o pid,ppid,lstart,etime,%cpu,%mem,stat,args --sort=-etime 2>/dev/null \
  | awk 'NR==1 || tolower($0) ~ /(aria2|wget|curl|aspera|python|rscript|cell2loc|cell2location|tar |gzip |unzip|rsync|scp |prefetch|fasterq|download)/' \
  | head -n 300 | tee -a "$REPORT" >&3
say '```'

section "12. Logs and recent failures"
subsection "Recently modified logs"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -mtime -30 \
  \( -iname '*.log' -o -iname '*.out' -o -iname '*.err' -o -iname 'nohup.out' \) \
  -printf '%T@\t%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null \
  | sort -nr | head -n 250 | cut -f2- | tee -a "$REPORT" >&3
say '```'

subsection "Error-like lines in recent lightweight logs"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -mtime -30 -size -"${TEXT_MAX_BYTES}"c \
  \( -iname '*.log' -o -iname '*.out' -o -iname '*.err' -o -iname 'nohup.out' \) -print0 2>/dev/null \
  | xargs -0 -r grep -IEin 'error|failed|failure|exception|traceback|killed|out of memory|oom|segmentation|corrupt|unexpected end|no space left|permission denied|connection reset|timed out' 2>/dev/null \
  | tail -n 600 | tee -a "$REPORT" >&3 || true
say '```'

section "13. Files modified in the last 7 days"
say '```text'
find "$PROJECT_ROOT" -maxdepth "$MAX_DEPTH" -type f -mtime -7 \
  -printf '%T@\t%s\t%TY-%Tm-%Td %TH:%TM:%TS\t%p\n' 2>/dev/null \
  | sort -nr | head -n 1000 | cut -f2- | tee -a "$REPORT" >&3
say '```'

section "14. Git state (read-only)"
say '```text'
mapfile -t GIT_DIRS < <(find "$PROJECT_ROOT" -maxdepth 5 -type d -name .git -printf '%h\n' 2>/dev/null | head -n 20)
if [[ "${#GIT_DIRS[@]}" -eq 0 ]]; then
  echo 'No Git worktree detected within depth 5.' | tee -a "$REPORT" >&3
else
  for repo in "${GIT_DIRS[@]}"; do
    echo "[$repo]" | tee -a "$REPORT" >&3
    if git -C "$repo" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
      git -C "$repo" status --short 2>&1 | head -n 300 | tee -a "$REPORT" >&3
      git -C "$repo" log -n 5 --date=iso --pretty='format:%h %ad %an %s' 2>&1 | tee -a "$REPORT" >&3
    else
      echo 'NOT_A_VALID_GIT_WORKTREE' | tee -a "$REPORT" >&3
    fi
    echo | tee -a "$REPORT" >&3
  done
fi
say '```'

section "15. Evidence summary and interpretation limits"
say "- This script reports evidence presence, file integrity indicators and inconsistencies; it does not declare a scientific module complete merely because a path exists."
say "- Patient counts, slice-to-patient mapping, 34,609 spots, 912 tests, 19 axes and 38/38 QC require confirmation from the listed result tables/logs or object metadata."
say "- External single-cell and spatial cohorts must not be described as matched multi-omics unless patient identity is directly documented."
say "- SCP/login-restricted and EGA policy-excluded datasets may have no local file evidence; absence in the scan is not evidence of data quality."
if [[ "$ARCHIVE_MODE" == 'quick' ]]; then
  say "- ARCHIVE_MODE=quick. Archives larger than $(human_size "$QUICK_ARCHIVE_MAX_BYTES") are listed but not fully read."
elif [[ "$ARCHIVE_MODE" == 'full' ]]; then
  say "- ARCHIVE_MODE=full. All recognized archives are tested, subject to a ${ARCHIVE_TIMEOUT_SEC}-second timeout per file."
else
  say "- ARCHIVE_MODE=none. Archive files are listed, but integrity tests are skipped."
fi
say "- The project root remained read-only; all generated audit artifacts are outside it."

section "16. Generated audit artifacts"
say '```text'
for p in "$REPORT" "$FILE_INDEX" "$DATASET_INDEX" "$COMMAND_LOG"; do stat_line "$p"; done | tee -a "$REPORT" >&3
say '```'

if have sha256sum; then
  sha256sum "$REPORT" "$FILE_INDEX" "$DATASET_INDEX" "$COMMAND_LOG" >"$OUT_DIR/project26210_scan_sha256_${NOW}.txt"
fi

printf '\nSCAN_COMPLETE\nREPORT=%s\nFILE_INDEX=%s\nDATASET_INDEX=%s\nCOMMAND_LOG=%s\n' \
  "$REPORT" "$FILE_INDEX" "$DATASET_INDEX" "$COMMAND_LOG"
if [[ -f "$OUT_DIR/project26210_scan_sha256_${NOW}.txt" ]]; then
  printf 'SHA256=%s\n' "$OUT_DIR/project26210_scan_sha256_${NOW}.txt"
fi
