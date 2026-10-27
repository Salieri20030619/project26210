from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


STAGE = Path(r"C:\Users\31937\project26210_phase105_state_update_v1")
BACKUP = STAGE / "preupdate_backup"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig")


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: object = "") -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]:
    check(f"required_{name}", (STAGE / name).is_file() and (STAGE / name).stat().st_size > 0)

markers = {
    "PROJECT_STATUS.md": "## Phase 105 full-project analysis audit and client-scope re-lock (2026-08-31)",
    "DECISIONS.md": "## 2026-08-31 - Phase 105 internal-single-cell redesign decision",
    "CODEX_HANDOFF.md": "## Phase 105 full-project audit and current client task (2026-08-31)",
}
for name, marker in markers.items():
    old = read(BACKUP / name)
    new = read(STAGE / name)
    check(f"append_only_{name}", new.startswith(old), f"old={len(old)} new={len(new)}")
    check(f"marker_once_{name}", new.count(marker) == 1)

old_manifest = tsv(BACKUP / "RESULTS_MANIFEST.tsv")
new_manifest = tsv(STAGE / "RESULTS_MANIFEST.tsv")
check("manifest_append_count_23", len(new_manifest) == len(old_manifest) + 23, f"old={len(old_manifest)} new={len(new_manifest)}")
check("manifest_old_rows_preserved", new_manifest[: len(old_manifest)] == old_manifest)
ids = [row["result_id"] for row in new_manifest]
check("manifest_unique_ids", len(ids) == len(set(ids)), f"rows={len(ids)} unique={len(set(ids))}")

audit_rows = [row for row in new_manifest if row["stage"] == "phase105_full_analysis_audit"]
state_rows = [row for row in new_manifest if row["stage"] == "phase105_project_state_update"]
check("phase105_audit_rows_20", len(audit_rows) == 20, len(audit_rows))
check("phase105_state_rows_3", len(state_rows) == 3, len(state_rows))
check("phase105_rows_all_pass", all(row["validation_status"] == "PASS" for row in audit_rows + state_rows))

for row in audit_rows + state_rows:
    path = Path(row["absolute_path"])
    if row["result_id"] == "phase105_state_update_validation_v1":
        path = STAGE / "STATE_UPDATE_VALIDATION_v1.json"
    elif row["stage"] == "phase105_project_state_update":
        path = STAGE / Path(row["absolute_path"]).name
    check(f"path_{row['result_id']}", path.is_file(), path)
    match = re.search(r"SHA256 ([0-9A-Fa-f]{64})", row["validation_evidence"])
    if match and path.is_file():
        check(f"hash_{row['result_id']}", sha256(path).lower() == match.group(1).lower(), path)

audit_validation = json.loads(Path(r"D:\project26210_external_joint_v1\phase105_full_analysis_audit_local_v1\PHASE105_AUDIT_VALIDATION_v1.json").read_text(encoding="utf-8"))
check("audit_validation_56_of_56", audit_validation.get("status") == "PASS_PHASE105_AUDIT_VALIDATION" and audit_validation.get("checks_passed") == 56 and audit_validation.get("checks_total") == 56)

status = read(STAGE / "PROJECT_STATUS.md")
decisions = read(STAGE / "DECISIONS.md")
handoff = read(STAGE / "CODEX_HANDOFF.md")
check("status_internal_five", all(dataset in status for dataset in ["GSE156632", "GSE159115", "GSE178481", "GSE207493", "GSE242299"]))
check("status_gap_explicit", "No formal internal-five-study" in status)
check("decision_client_data_boundary", "Do not request, import, reanalyse, redraw" in decisions)
check("handoff_resume_gap", "Unresolved gap" in handoff and "fine-state-specific communication" in handoff)
check("handoff_no_active_phase105_job", "No current Phase 105 background job remains" in handoff)

failed = [item for item in checks if not item["passed"]]
payload = {
    "status": "PASS_PHASE105_COMPLETION_GATE" if not failed else "FAIL_PHASE105_COMPLETION_GATE",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "failed_checks": failed,
    "key_file_sha256": {
        name: sha256(STAGE / name)
        for name in ["PROJECT_STATUS.md", "DECISIONS.md", "RESULTS_MANIFEST.tsv", "CODEX_HANDOFF.md"]
    },
}
(STAGE / "STATE_UPDATE_VALIDATION_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(payload, ensure_ascii=False, indent=2))
raise SystemExit(0 if not failed else 1)
