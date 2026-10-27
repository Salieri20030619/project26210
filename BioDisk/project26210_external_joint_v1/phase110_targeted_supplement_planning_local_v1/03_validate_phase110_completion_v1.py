from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase110_targeted_supplement_planning_local_v1"
NOW = datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


for filename, marker in {
    "PROJECT_STATUS.md": "## Phase 110 targeted supplementation planning",
    "DECISIONS.md": "## Phase 110 targeted-supplementation planning decisions",
    "CODEX_HANDOFF.md": "## Phase 110 current planning breakpoint",
}.items():
    count = (ROOT / filename).read_text(encoding="utf-8-sig").count(marker)
    check(f"unique marker {filename}", count == 1, f"count={count}")

planning = json.loads((PHASE / "PHASE110_PLANNING_VALIDATION_v1.json").read_text(encoding="utf-8"))
state = json.loads((PHASE / "PHASE110_STATE_UPDATE_STATUS_v1.json").read_text(encoding="utf-8"))
check("planning validation PASS", planning.get("status") == "PASS_PHASE110_TARGETED_SUPPLEMENT_PLAN_VALIDATION", str(planning.get("status")))
check("planning checks 33/33", planning.get("checks_passed") == 33 and planning.get("checks_total") == 33, f"{planning.get('checks_passed')}/{planning.get('checks_total')}")
check("no biological analysis executed", planning.get("analysis_executed") is False and state.get("analysis_executed") is False, "planning and state false")
check("state update PASS", state.get("status") == "PASS_PHASE110_STATE_UPDATE", str(state.get("status")))

manifest_path = ROOT / "RESULTS_MANIFEST.tsv"
with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)
phase_rows = [row for row in rows if row["stage"] == "phase110_targeted_supplement_planning"]
check("global manifest IDs unique", len(rows) == len({row["result_id"] for row in rows}), f"rows={len(rows)}")
check("initial Phase110 rows", len(phase_rows) == 10, f"rows={len(phase_rows)}")
for row in phase_rows:
    path = Path(row["absolute_path"])
    check(f"path exists {row['result_id']}", path.is_file(), str(path))
    match = re.search(r"SHA256 ([0-9A-F]{64})", row["validation_evidence"])
    actual = sha256(path).upper() if path.is_file() else "MISSING"
    check(f"hash matches {row['result_id']}", bool(match) and match.group(1) == actual, actual)

if not all(item["passed"] for item in checks):
    raise RuntimeError(json.dumps({"status": "FAIL_PHASE110_COMPLETION", "checks": checks}, ensure_ascii=False, indent=2))

receipt_path = PHASE / "PHASE110_COMPLETION_VALIDATION_v1.json"
receipt = {
    "status": "PASS_PHASE110_COMPLETION_VALIDATION",
    "generated_at": NOW,
    "checks_passed": len(checks),
    "checks_total": len(checks),
    "analysis_executed": False,
    "checks": checks,
}
receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

script_path = PHASE / "03_validate_phase110_completion_v1.py"
extra = [
    ("phase110_state_update_status_v1", "state_update_json", PHASE / "PHASE110_STATE_UPDATE_STATUS_v1.json"),
    ("phase110_completion_validator_v1", "formal_validation_script", script_path),
    ("phase110_completion_validation_v1", "state_update_validation_json", receipt_path),
]
existing = {row["result_id"] for row in rows}
if any(result_id in existing for result_id, _, _ in extra):
    raise RuntimeError("Phase110 completion rows already exist; refusing rerun")
with manifest_path.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    for result_id, artifact_type, path in extra:
        writer.writerow({
            "result_id": result_id,
            "stage": "phase110_project_state_update",
            "artifact_type": artifact_type,
            "absolute_path": str(path),
            "generated_at": NOW,
            "validation_status": "PASS",
            "validation_evidence": f"SHA256 {sha256(path).upper()}; Phase110 mandatory four-file completion gate {len(checks)}/{len(checks)} PASS",
            "notes": "Phase110 planning state synchronization evidence; no biological analysis executed",
        })

with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
    final_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(final_rows) != len({row["result_id"] for row in final_rows}):
    raise RuntimeError("Manifest IDs are not unique after Phase110 completion registration")
for result_id, _, path in extra:
    row = next(row for row in final_rows if row["result_id"] == result_id)
    registered = re.search(r"SHA256 ([0-9A-F]{64})", row["validation_evidence"]).group(1)
    if registered != sha256(path).upper():
        raise RuntimeError(f"Final hash mismatch: {result_id}")

print(json.dumps({
    "status": "PASS_PHASE110_COMPLETION_VALIDATION",
    "checks": f"{len(checks)}/{len(checks)}",
    "phase110_manifest_rows": 13,
    "four_file_sha256": {name: sha256(ROOT / name) for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]},
}, ensure_ascii=False, indent=2))
