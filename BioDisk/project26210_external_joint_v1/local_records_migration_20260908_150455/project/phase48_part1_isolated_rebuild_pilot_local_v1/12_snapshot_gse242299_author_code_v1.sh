#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 12_snapshot_gse242299_author_code_v1.sh REBUILD_ROOT PROXY_URL}
proxy_url=${2:?Usage: 12_snapshot_gse242299_author_code_v1.sh REBUILD_ROOT PROXY_URL}
url=https://github.com/zvirblyte/2023_ccRCC.git
final_dir="$rebuild_root/01_input_inventory/gse242299_author_code_snapshot_v1"
part_dir="$final_dir.part"
log_file="$rebuild_root/logs/12_snapshot_gse242299_author_code_v1.log"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
[[ ! -e "$final_dir" ]] || { echo "refusing existing final directory: $final_dir" >&2; exit 1; }
[[ ! -e "$part_dir" ]] || { echo "refusing existing part directory: $part_dir" >&2; exit 1; }
git --version
git -c http.proxy="$proxy_url" clone --depth 1 "$url" "$part_dir"
commit=$(git -C "$part_dir" rev-parse HEAD)
tree=$(git -C "$part_dir" rev-parse HEAD^{tree})
status=$(git -C "$part_dir" status --porcelain)
[[ -z "$status" ]] || { echo "snapshot worktree is not clean" >&2; exit 1; }
{
  echo -e "repository\tcommit\ttree\tclone_depth\tstatus"
  echo -e "$url\t$commit\t$tree\t1\tPASS"
} > "$part_dir/SNAPSHOT_MANIFEST_v1.tsv"
mv "$part_dir" "$final_dir"
echo "commit=$commit"
echo "tree=$tree"
echo "finished_at=$(date --iso-8601=seconds)"
