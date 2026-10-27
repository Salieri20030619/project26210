from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(r"D:\project26210_external_joint_v1")
AUDIT = ROOT / "recovery_audit_20260831_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: object = "") -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


required = ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
for name in required:
    path = ROOT / name
    check(f"required_{name}", path.is_file() and path.stat().st_size > 0, path)

status = (ROOT / "PROJECT_STATUS.md").read_text(encoding="utf-8-sig")
decisions = (ROOT / "DECISIONS.md").read_text(encoding="utf-8-sig")
handoff = (ROOT / "CODEX_HANDOFF.md").read_text(encoding="utf-8-sig")

for phase in range(100, 105):
    check(f"status_phase{phase}_present", f"Phase {phase}" in status)
check("status_recovery_section_once", status.count("## Recovery audit and authoritative breakpoint (2026-08-31)") == 1)
check("decisions_recovery_section_once", decisions.count("## 2026-08-31 - project recovery and authoritative-state decision") == 1)
check("handoff_phase104_breakpoint", "Phase 104" in handoff and "DIRECT_MTA_NOT_EVALUABLE" in handoff)
check("handoff_786_boundary", "786-0" in handoff and "匹配对照" in handoff)

with (ROOT / "RESULTS_MANIFEST.tsv").open("r", encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

expected_columns = [
    "result_id",
    "stage",
    "artifact_type",
    "absolute_path",
    "generated_at",
    "validation_status",
    "validation_evidence",
    "notes",
]
check("manifest_columns", bool(rows) and list(rows[0].keys()) == expected_columns, list(rows[0].keys()) if rows else [])
ids = [row["result_id"] for row in rows]
check("manifest_unique_ids", len(ids) == len(set(ids)), f"rows={len(ids)} unique={len(set(ids))}")

minimum_counts = {
    "phase100_client_mta_immune_brief": 6,
    "phase101_public_cellline_data": 7,
    "phase102_nci60_two_level_preprocessing": 7,
    "phase103_nci60_two_level_analysis": 11,
    "phase104_ccl180_ccle_adjacent_metabolites": 13,
    "project_recovery_20260831": 2,
}
for stage, minimum in minimum_counts.items():
    stage_rows = [row for row in rows if row["stage"] == stage]
    check(f"manifest_{stage}_count", len(stage_rows) >= minimum, len(stage_rows))
    check(f"manifest_{stage}_all_pass", bool(stage_rows) and all(row["validation_status"] == "PASS" for row in stage_rows), len(stage_rows))

recent_rows = [row for row in rows if row["stage"] in minimum_counts]
for row in recent_rows:
    path = Path(row["absolute_path"])
    check(f"path_{row['result_id']}", path.is_file(), path)
    match = re.search(r"SHA256[ =]([0-9A-Fa-f]{64})", row["validation_evidence"])
    if match and path.is_file():
        check(f"hash_{row['result_id']}", sha256(path).lower() == match.group(1).lower(), path)

validation_specs = [
    (ROOT / "phase102_nci60_two_level_preprocessing_local_v1" / "qc_v1" / "independent_validation_status_v1.json", 15),
    (ROOT / "phase103_nci60_two_level_analysis_local_v1" / "qc_v1" / "independent_validation_status_v1.json", 20),
    (ROOT / "phase104_ccl180_ccle_adjacent_metabolites_local_v1" / "qc_v1" / "independent_validation_status_v1.json", 18),
]
for path, expected in validation_specs:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    label = path.parent.parent.name
    check(f"validation_status_{label}", payload.get("status") == "PASS_INDEPENDENT_VALIDATION")
    check(
        f"validation_count_{label}",
        payload.get("checks_passed") == expected and payload.get("checks_total") == expected,
        payload,
    )

failed = [item for item in checks if not item["passed"]]
output = {
    "status": "PASS_PROJECT_RECOVERY_VALIDATION" if not failed else "FAIL_PROJECT_RECOVERY_VALIDATION",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "failed_checks": failed,
    "authoritative_root": str(ROOT),
    "latest_completed_phase": 104,
    "key_file_sha256": {name: sha256(ROOT / name) for name in required if (ROOT / name).is_file()},
}
AUDIT.mkdir(parents=True, exist_ok=True)
(AUDIT / "RECOVERY_AUDIT_STATUS_v1.json").write_text(
    json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(output, ensure_ascii=False, indent=2))
raise SystemExit(0 if not failed else 1)
