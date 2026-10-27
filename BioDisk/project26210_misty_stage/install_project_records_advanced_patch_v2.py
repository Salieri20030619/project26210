#!/usr/bin/env python3
"""Install final project records after the one-paragraph report-date correction."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
QC = ROOT / "results/part4_spatial/advanced_patch_report_qc_v1"
FILES = {
    "PROJECT_STATUS.md": {
        "stage": "PROJECT_STATUS_v2.md",
        "before": "5ce86df76bc12780a1cc54455e127bdc4c5bac0a440da14b9f0c2e7789e6c697",
        "after": "92f28ab8032f5194f86da1dd4d019e17c9d608bfa03cd29f4127eefed358b3fe",
    },
    "CODEX_HANDOFF.md": {
        "stage": "CODEX_HANDOFF_v2.md",
        "before": "13b7ef794333e9754ecf0fe1f9b0a743da6eeb8dd2125a48a8ada7e2a84a9eff",
        "after": "fdab7725753db9565953b6d32cfbcde2fb530db81b07e2a92764899afab0e006",
    },
    "RESULTS_MANIFEST.tsv": {
        "stage": "RESULTS_MANIFEST_v2.tsv",
        "before": "9efba51d852ec063abf09b2ed9e0a5ea1a5863cbf3ac7f183dbd93a3e00d0589",
        "after": "7574a67bea0ce372ea64eb1fe8e9079462daf42393fd8a7407904cb114bae937",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


for name, info in FILES.items():
    official = ROOT / name
    staged = STAGE / info["stage"]
    if sha256(official) != info["before"] or sha256(staged) != info["after"]:
        raise RuntimeError(f"Record hash gate failed: {name}")
with (STAGE / "RESULTS_MANIFEST_v2.tsv").open(newline="") as handle:
    rows = list(csv.reader(handle, delimiter="\t"))
ids = [row[0] for row in rows[1:]]
if len(rows) != 286 or any(len(row) != 18 for row in rows) or len(ids) != len(set(ids)):
    raise RuntimeError("Final staged manifest validation failed")

records = []
for name, info in FILES.items():
    official = ROOT / name
    staged = STAGE / info["stage"]
    backup = QC / f"{name}.before_datefix_record_update"
    if backup.exists() and sha256(backup) != info["before"]:
        raise RuntimeError(f"Unexpected backup hash: {backup}")
    if not backup.exists():
        shutil.copy2(official, backup)
    temp = ROOT / f".{name}.datefix_record_v2.tmp"
    if temp.exists():
        raise RuntimeError(f"Unexpected temporary record: {temp}")
    shutil.copy2(staged, temp)
    if sha256(temp) != info["after"]:
        temp.unlink(missing_ok=True)
        raise RuntimeError(f"Temporary record hash failed: {name}")
    os.replace(temp, official)
    if sha256(official) != info["after"]:
        raise RuntimeError(f"Official record hash failed: {name}")
    records.append({
        "path": str(official), "size": official.stat().st_size,
        "before_sha256": info["before"], "after_sha256": info["after"],
        "backup": str(backup), "backup_sha256": sha256(backup),
    })

audit = {
    "status": "PASS",
    "installed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "records": records,
    "manifest_rows_excluding_header": 285,
    "manifest_columns": 18,
    "manifest_duplicate_result_ids": 0,
    "official_report_sha256_recorded": "f3b61f3334c8bdb5dbbd6cb5b674f5ea13acbb3bb7fa0ebc863c4d3b2cb86bf4",
}
(QC / "project_record_update_audit_v2.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
shutil.copy2(Path(__file__), QC / Path(__file__).name)
print(json.dumps(audit, ensure_ascii=False, indent=2))
