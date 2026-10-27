#!/usr/bin/env bash
set -euo pipefail

root=/mnt/BioDisk/wangmingyuan/project26210_part4_composition_niche_phase18_20260813_v1
python=/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python
phase17=/mnt/BioDisk/wangmingyuan/project26210_part4_tymp_reset_phase17_20260812_v1
output="$root/independent_validation_v2"

if [[ -e "$output" ]]; then
    echo "Refusing to overwrite validation directory: $output" >&2
    exit 2
fi

"$python" "$root/provenance/02_validate_composition_niches_v1.py" \
    --analysis-root "$root/results_v1" \
    --phase17-audit "$phase17/mapping_queue_audit_v1/mapping_queue_audit_v1.tsv" \
    --output-root "$output" \
    > "$root/logs/validation_v2.log" 2>&1

cat "$root/logs/validation_v2.log"
