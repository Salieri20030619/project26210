#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: bash 07_download_pmc_xlsx_v1.sh <url> <output_xlsx> <proxy_url> <pow_solver.py>" >&2
  exit 2
fi

url="$1"
output="$2"
proxy_url="$3"
solver="$4"
work_prefix="${output}.pmc_pow"

if [[ -e "$output" ]]; then
  echo "Refusing to overwrite existing output: $output" >&2
  exit 3
fi
mkdir -p "$(dirname "$output")"

for attempt in 1 2 3; do
  challenge_page="${work_prefix}.attempt${attempt}.html"
  candidate="${output}.attempt${attempt}.part"
  curl --http1.1 --compressed --proxy "$proxy_url" --fail --location \
    --retry 5 --retry-all-errors --user-agent "Mozilla/5.0" \
    --output "$challenge_page" "$url"
  cookie="$(python3 "$solver" "$challenge_page")"
  curl --http1.1 --compressed --proxy "$proxy_url" --fail --location \
    --retry 5 --retry-all-errors --user-agent "Mozilla/5.0" \
    --cookie "$cookie" --output "$candidate" "$url"
  if unzip -t "$candidate" >/dev/null 2>&1; then
    mv "$candidate" "$output"
    sha256sum "$output" > "${output}.sha256"
    echo "status=PASS"
    echo "attempt=$attempt"
    echo "bytes=$(stat -c '%s' "$output")"
    exit 0
  fi
  mv "$candidate" "${candidate}.not_xlsx"
done

echo "PMC PoW succeeded computationally but no valid XLSX was returned" >&2
exit 4
