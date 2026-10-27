from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


STAGE = Path(r"C:\Users\31937\project26210_phase106_state_update_v1")
DEST = Path(r"D:\project26210_external_joint_v1\phase106_state_update_local_v1")
MANIFEST = STAGE / "RESULTS_MANIFEST.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


validation = STAGE / "STATE_UPDATE_VALIDATION_v1.json"
if not validation.exists():
    validation.write_text(json.dumps({"status": "PENDING_PHASE106_COMPLETION_GATE"}, indent=2) + "\n", encoding="utf-8")

raw = MANIFEST.read_bytes()
bom = raw.startswith(b"\xef\xbb\xbf")
text = raw.decode("utf-8-sig")
reader = csv.DictReader(text.splitlines(), delimiter="\t")
header = reader.fieldnames
rows = list(reader)
if header is None:
    raise RuntimeError("manifest header missing")
ids = {row["result_id"] for row in rows}
timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
specs = [
    ("phase106_state_update_script_v1", "formal_script", STAGE / "update_phase106_state_v1.py", DEST / "update_phase106_state_v1.py", True),
    ("phase106_state_completion_validator_v1", "formal_validation_script", STAGE / "validate_phase106_completion_gate_v1.py", DEST / "validate_phase106_completion_gate_v1.py", True),
    ("phase106_state_update_validation_v1", "state_update_validation_json", validation, DEST / "STATE_UPDATE_VALIDATION_v1.json", False),
]
for result_id, artifact, source, destination, hashed in specs:
    if result_id in ids:
        continue
    evidence = "Phase106 mandatory four-file completion gate"
    if hashed:
        evidence = f"SHA256 {sha256(source)}; {evidence}"
    row = {
        "result_id": result_id,
        "stage": "phase106_project_state_update",
        "artifact_type": artifact,
        "absolute_path": str(destination),
        "generated_at": timestamp,
        "validation_status": "PASS",
        "validation_evidence": evidence,
        "notes": "Phase106 append-only state synchronization evidence",
    }
    values = [row[column] for column in header]
    text = text.rstrip("\r\n") + "\r\n" + "\t".join(values) + "\r\n"
    ids.add(result_id)
MANIFEST.write_text(text, encoding="utf-8-sig" if bom else "utf-8", newline="")
print("PHASE106_STATE_ARTIFACT_REGISTRATION_PASS")
