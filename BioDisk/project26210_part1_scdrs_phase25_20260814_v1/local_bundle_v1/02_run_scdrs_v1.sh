#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 SCDRS_EXEC H5AD COV GS OUTPUT_DIR" >&2
  exit 2
fi

scdrs_exec="$1"
h5ad="$2"
cov="$3"
gs="$4"
output="$5"

[[ -x "$scdrs_exec" ]] || { echo "scdrs executable missing: $scdrs_exec" >&2; exit 3; }
for source in "$h5ad" "$cov" "$gs"; do
  [[ -f "$source" ]] || { echo "input missing: $source" >&2; exit 4; }
done
[[ ! -e "$output" ]] || { echo "refusing existing output: $output" >&2; exit 5; }
part="${output}.part"
[[ ! -e "$part" ]] || { echo "refusing stale partial output: $part" >&2; exit 6; }
mkdir -p "$part/scores"

set +e
"$scdrs_exec" compute-score \
  --h5ad-file "$h5ad" \
  --h5ad-species human \
  --cov-file "$cov" \
  --gs-file "$gs" \
  --gs-species human \
  --out-folder "$part/scores" \
  --flag-filter-data True \
  --flag-raw-count True \
  --n-ctrl 1000 \
  > "$part/scdrs_compute.log" 2>&1
exit_code=$?
set -e
if [[ $exit_code -ne 0 ]]; then
  printf '{"status":"FAIL_SCDRS_COMPUTE","exit_code":%d}\n' "$exit_code" > "$part/status.json"
  exit "$exit_code"
fi

score_count=$(find "$part/scores" -maxdepth 1 -type f -name '*.score.gz' | wc -l)
full_count=$(find "$part/scores" -maxdepth 1 -type f -name '*.full_score.gz' | wc -l)
[[ "$score_count" -eq 1 ]] || { echo "expected one score.gz, found $score_count" >&2; exit 7; }
[[ "$full_count" -eq 1 ]] || { echo "expected one full_score.gz, found $full_count" >&2; exit 8; }
sha256sum "$part"/scores/* "$part/scdrs_compute.log" > "$part/SHA256SUMS.txt"
printf '{"status":"PASS_SCDRS_COMPUTE","n_ctrl":1000,"raw_count":true}\n' > "$part/status.json"
mv "$part" "$output"

