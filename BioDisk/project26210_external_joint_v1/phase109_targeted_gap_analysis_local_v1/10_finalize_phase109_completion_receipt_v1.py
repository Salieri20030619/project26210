from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = ROOT / "phase109_targeted_gap_analysis_local_v1"
NOW = datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


manifest = ROOT / "RESULTS_MANIFEST.tsv"
with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)
phase_rows = [r for r in rows if r["stage"] in {"phase109_targeted_gap_analysis", "phase109_project_state_update"}]

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


check("global result IDs unique", len(rows) == len({r["result_id"] for r in rows}), f"rows={len(rows)}")
check("Phase109 registered rows before receipt", len(phase_rows) == 24, f"rows={len(phase_rows)}")
for filename, marker in {
    "PROJECT_STATUS.md": "## Phase 109 targeted supplementation",
    "DECISIONS.md": "## Phase 109 targeted-gap decisions",
    "CODEX_HANDOFF.md": "## Phase 109 current breakpoint",
}.items():
    count = (ROOT / filename).read_text(encoding="utf-8-sig").count(marker)
    check(f"unique state marker {filename}", count == 1, f"count={count}")
for row in phase_rows:
    path = Path(row["absolute_path"])
    check(f"registered path exists {row['result_id']}", path.is_file(), str(path))
    match = re.search(r"SHA256 ([0-9A-F]{64})", row["validation_evidence"])
    actual = sha256(path).upper() if path.is_file() else "MISSING"
    check(f"registered hash matches {row['result_id']}", bool(match) and match.group(1) == actual, actual)

if not all(item["passed"] for item in checks):
    raise RuntimeError(json.dumps({"status": "FAIL_PHASE109_FINAL_COMPLETION_RECEIPT", "checks": checks}, ensure_ascii=False, indent=2))

receipt_path = PHASE / "PHASE109_FINAL_COMPLETION_RECEIPT_v1.json"
receipt = {
    "status": "PASS_PHASE109_FINAL_COMPLETION_RECEIPT",
    "generated_at": NOW,
    "checks_passed": len(checks),
    "checks_total": len(checks),
    "phase109_registered_rows_validated": len(phase_rows),
    "four_file_sha256_before_receipt_registration": {
        name: sha256(ROOT / name)
        for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
    },
    "checks": checks,
}
receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

script_path = PHASE / "10_finalize_phase109_completion_receipt_v1.py"
extra = [
    ("phase109_final_completion_receipt_script_v1", "formal_validation_script", script_path),
    ("phase109_final_completion_receipt_v1", "state_update_validation_json", receipt_path),
]
existing = {r["result_id"] for r in rows}
if any(result_id in existing for result_id, _, _ in extra):
    raise RuntimeError("Final receipt rows already exist; refusing rerun")
with manifest.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    for result_id, artifact_type, path in extra:
        writer.writerow(
            {
                "result_id": result_id,
                "stage": "phase109_project_state_update",
                "artifact_type": artifact_type,
                "absolute_path": str(path),
                "generated_at": NOW,
                "validation_status": "PASS",
                "validation_evidence": f"SHA256 {sha256(path).upper()}; persisted Phase109 final completion receipt {len(checks)}/{len(checks)} PASS",
                "notes": "Phase109 persisted four-file completion and hash-validation evidence",
            }
        )

with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
    final_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(final_rows) != len({r["result_id"] for r in final_rows}):
    raise RuntimeError("Manifest IDs are not unique after receipt registration")
for result_id, _, path in extra:
    row = next(r for r in final_rows if r["result_id"] == result_id)
    registered = re.search(r"SHA256 ([0-9A-F]{64})", row["validation_evidence"]).group(1)
    if registered != sha256(path).upper():
        raise RuntimeError(f"Receipt registration hash mismatch: {result_id}")

print(json.dumps({
    "status": "PASS_PHASE109_FINAL_COMPLETION_RECEIPT",
    "persisted_checks": f"{len(checks)}/{len(checks)}",
    "phase109_manifest_rows_after_receipt": len([r for r in final_rows if r["stage"] in {"phase109_targeted_gap_analysis", "phase109_project_state_update"}]),
    "four_file_sha256_final": {name: sha256(ROOT / name) for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]},
}, ensure_ascii=False, indent=2))
