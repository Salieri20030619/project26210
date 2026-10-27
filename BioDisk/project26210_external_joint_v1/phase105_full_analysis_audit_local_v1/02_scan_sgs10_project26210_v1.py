from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime
from pathlib import Path


OUT = Path(r"C:\Users\31937\project26210_phase105_full_analysis_audit_v1")
HOST = "sgs10-codex"


def ssh(command: str) -> str:
    completed = subprocess.run(
        ["ssh", "-o", "ClearAllForwardings=yes", "-o", "BatchMode=yes", HOST, command],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return completed.stdout


directory_command = r'''bash -lc '
for d in /mnt/BioDisk/wangmingyuan/*26210* /mnt/BioDisk/ludingyu/project26210; do
  [ -d "$d" ] || continue
  shallow=$(find "$d" -maxdepth 4 -type f 2>/dev/null | wc -l)
  markers=$(find "$d" -maxdepth 4 -type f \( -iname "*PASS*" -o -iname "*COMPLETE*" -o -iname "*STATUS*" -o -iname "*VALIDATION*" \) 2>/dev/null | wc -l)
  latest=$(find "$d" -maxdepth 4 -type f -printf "%T@\n" 2>/dev/null | sort -nr | head -n 1)
  printf "%s\t%s\t%s\t%s\n" "$d" "$shallow" "$markers" "$latest"
done
' '''

lines = [line for line in ssh(directory_command).splitlines() if line.strip()]
rows: list[dict[str, object]] = []
for line in lines:
    path, shallow, markers, latest = line.split("\t")
    rows.append(
        {
            "remote_path": path,
            "shallow_file_count_maxdepth4": int(shallow),
            "completion_marker_count_maxdepth4": int(markers),
            "latest_file_epoch_maxdepth4": latest,
        }
    )

with (OUT / "sgs10_project26210_directory_inventory_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda row: row["remote_path"]))

required_command = r'''bash -lc '
files="
/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/PROJECT_STATUS.md
/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/DECISIONS.md
/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/RESULTS_MANIFEST.tsv
/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/CODEX_HANDOFF.md
/mnt/BioDisk/wangmingyuan/project26210_monocyte_spp1_phase21_20260814_v1/results_v2/PHASE21_MONOCYTE_SPP1_TRAJECTORY_STATUS_v1.json
/mnt/BioDisk/wangmingyuan/project26210_monocyte_spp1_phase21_20260814_v1/independent_validation_v1/PHASE21_INDEPENDENT_VALIDATION_SUMMARY_v1.json
/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/analysis_phase99_validation_v1/validation_status_phase99_v1.json
"
for f in $files; do
  if [ -f "$f" ]; then
    size=$(stat -c %s "$f")
    hash=$(sha256sum "$f" | cut -d" " -f1)
    printf "PASS\t%s\t%s\t%s\n" "$size" "$hash" "$f"
  else
    printf "MISSING\t\t\t%s\n" "$f"
  fi
done
' '''
required_lines = [line for line in ssh(required_command).splitlines() if line.strip()]
required_rows = []
for line in required_lines:
    status, size, digest, path = line.split("\t", 3)
    required_rows.append({"status": status, "size_bytes": size, "sha256": digest, "remote_path": path})
with (OUT / "sgs10_required_evidence_audit_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(required_rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(required_rows)

processes = ssh("pgrep -af 'project26210' || true")
(OUT / "sgs10_project26210_process_snapshot_v1.txt").write_text(processes, encoding="utf-8")

summary = {
    "status": "PASS_SGS10_SCAN_COMPLETE",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "remote_project_directories": len(rows),
    "required_evidence_files": len(required_rows),
    "required_evidence_missing": sum(row["status"] != "PASS" for row in required_rows),
    "process_snapshot_nonempty": bool(processes.strip()),
}
(OUT / "sgs10_scan_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
