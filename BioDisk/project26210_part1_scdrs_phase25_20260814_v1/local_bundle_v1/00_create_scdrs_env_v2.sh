#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 BASE_PYTHON SCDRS_SOURCE ENV_DIR RECORD_DIR" >&2
  exit 2
fi

base_python="$1"
source_dir="$2"
env_dir="$3"
record_dir="$4"

[[ -x "$base_python" ]] || { echo "base Python is not executable: $base_python" >&2; exit 3; }
[[ -f "$source_dir/setup.py" && -f "$source_dir/scdrs/version.py" ]] || {
  echo "scDRS source is incomplete: $source_dir" >&2; exit 4;
}
[[ ! -e "$env_dir" ]] || { echo "refusing existing environment: $env_dir" >&2; exit 5; }
[[ ! -e "$record_dir" ]] || { echo "refusing existing record directory: $record_dir" >&2; exit 6; }

record_part="${record_dir}.part"
[[ ! -e "$record_part" ]] || { echo "refusing stale partial record: $record_part" >&2; exit 7; }
mkdir -p "$record_part"

cleanup() {
  status=$?
  if [[ $status -ne 0 ]]; then
    printf '{"status":"FAIL","exit_code":%d}\n' "$status" > "$record_part/status.json"
  fi
}
trap cleanup EXIT

git -c "safe.directory=$source_dir" -C "$source_dir" status --short > "$record_part/source_git_status.txt"
[[ ! -s "$record_part/source_git_status.txt" ]] || {
  echo "scDRS source worktree is not clean" >&2; exit 8;
}
git -c "safe.directory=$source_dir" -C "$source_dir" rev-parse HEAD > "$record_part/source_git_commit.txt"
git -c "safe.directory=$source_dir" -C "$source_dir" describe --tags --always \
  > "$record_part/source_git_describe.txt"
find "$source_dir" -path "$source_dir/.git" -prune -o -type f -print0 \
  | sort -z | xargs -0 sha256sum > "$record_part/source_files_sha256.txt"

"$base_python" -m venv --system-site-packages "$env_dir"
"$env_dir/bin/python" -m pip install --no-deps --no-build-isolation "$source_dir" \
  > "$record_part/pip_install.log" 2>&1

"$env_dir/bin/python" - <<'PY' > "$record_part/self_test.txt" 2>&1
from importlib.metadata import version
import scdrs
print("scdrs_import=PASS")
print("scdrs_version=" + version("scdrs"))
PY
"$env_dir/bin/scdrs" --help > "$record_part/scdrs_help.txt" 2>&1
"$env_dir/bin/python" -m pip freeze --all > "$record_part/pip_freeze.txt"
"$env_dir/bin/python" -m pip show scdrs > "$record_part/pip_show_scdrs.txt"

sha256sum "$record_part"/*.txt "$record_part"/*.log > "$record_part/SHA256SUMS.txt"
commit=$(tr -d '\r\n' < "$record_part/source_git_commit.txt")
installed=$("$env_dir/bin/python" -c 'from importlib.metadata import version; print(version("scdrs"))')
printf '{"status":"PASS_ENVIRONMENT_COMPLETE","scdrs_version":"%s","source_commit":"%s"}\n' \
  "$installed" "$commit" > "$record_part/status.json"
mv "$record_part" "$record_dir"
trap - EXIT

