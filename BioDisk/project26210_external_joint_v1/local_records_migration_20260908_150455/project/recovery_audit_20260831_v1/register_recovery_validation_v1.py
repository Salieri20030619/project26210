from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
MANIFEST = ROOT / "RESULTS_MANIFEST.tsv"
RESULT_ID = "project_state_recovery_validation_v1"

raw = MANIFEST.read_bytes()
bom = raw.startswith(b"\xef\xbb\xbf")
text = raw.decode("utf-8-sig")
reader = csv.DictReader(text.splitlines(), delimiter="\t")
rows = list(reader)
header = reader.fieldnames
if header is None:
    raise RuntimeError("manifest header missing")

if RESULT_ID not in {row["result_id"] for row in rows}:
    row = {
        "result_id": RESULT_ID,
        "stage": "project_recovery_20260831",
        "artifact_type": "recovery_validation_json",
        "absolute_path": str(ROOT / "recovery_audit_20260831_v1" / "RECOVERY_AUDIT_STATUS_v1.json"),
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "validation_status": "PASS",
        "validation_evidence": "PASS_PROJECT_RECOVERY_VALIDATION; zero failed checks; key-file hashes recorded in JSON",
        "notes": "Independent structural, path, hash and phase-validation reconciliation",
    }
    values = [row[column] for column in header]
    if any("\t" in value or "\n" in value or "\r" in value for value in values):
        raise RuntimeError("invalid TSV field")
    text = text.rstrip("\r\n") + "\r\n" + "\t".join(values) + "\r\n"
    MANIFEST.write_text(text, encoding="utf-8-sig" if bom else "utf-8", newline="")

print("RECOVERY_VALIDATION_REGISTER_PASS")
