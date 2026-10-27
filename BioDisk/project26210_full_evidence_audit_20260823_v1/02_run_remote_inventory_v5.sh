#!/usr/bin/env bash
set -euo pipefail
AUDIT_ROOT=/mnt/BioDisk/wangmingyuan/project26210_full_evidence_audit_20260823_v1
PROJECT_ROOT=/mnt/BioDisk/ludingyu/project26210
OUTDIR="$AUDIT_ROOT/inventory_remote_v4"
mkdir -p "$OUTDIR" "$AUDIT_ROOT/logs_v1"
"$PROJECT_ROOT/check_memory_gate.sh" | tee "$AUDIT_ROOT/logs_v1/memory_gate_inventory_v4.log"
args=(--root "LUDINGYU_PROJECT=$PROJECT_ROOT")
for d in /mnt/BioDisk/wangmingyuan/project26210*; do
  [[ -d "$d" ]] || continue
  [[ "$d" == "$AUDIT_ROOT" ]] && continue
  [[ "$(basename "$d")" == *env* ]] && continue
  args+=(--root "$(basename "$d")=$d")
done
python3 "$AUDIT_ROOT/01_inventory_and_evidence_scan_v4.py" "${args[@]}" --outdir "$OUTDIR" 2>&1 | tee "$AUDIT_ROOT/logs_v1/inventory_remote_v4.log"
