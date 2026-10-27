from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(r"C:\Users\31937\project26210_phase107_evidence_escalation_audit_v1")


def tsv(name: str) -> list[dict[str, str]]:
    with (ROOT / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: object = "") -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


required = [
    "PHASE107_EVIDENCE_ESCALATION_LOGIC_v1.md",
    "FOLLOWUP_DATASET_RATIONALE_v1.tsv",
    "CHAIN_EDGE_EVIDENCE_MATRIX_v1.tsv",
    "CLIENT_RESPONSE_BLUEPRINT_v3.md",
    "PHASE107_EVIDENCE_ESCALATION_AUDIT_SUMMARY_v1.md",
]
for name in required:
    path = ROOT / name
    check(f"required_{name}", path.is_file() and path.stat().st_size > 0)

rationale = tsv("FOLLOWUP_DATASET_RATIONALE_v1.tsv")
check("rationale_rows_7", len(rationale) == 7, len(rationale))
check("rationale_ids_unique", len({row["step_id"] for row in rationale}) == 7)
check("immune_added_for_gap", any(row["step_id"] == "S01" and "Rebuilt Part1" in row["starting_gap"] for row in rationale))
check("spatial_upstream_downstream", any(row["step_id"] == "S03" and row["gap_status"] == "DOWNSTREAM_ONLY" for row in rationale))
check("metabolomics_added_for_direct_mta", any(row["step_id"] == "S04" and "measure MTA directly" in row["starting_gap"] for row in rationale))
check("client_pending", any(row["step_id"] == "S07" and row["gap_status"] == "PENDING_CLIENT_OWNED" for row in rationale))

edges = tsv("CHAIN_EDGE_EVIDENCE_MATRIX_v1.tsv")
check("edge_rows_9", len(edges) == 9, len(edges))
check("edge_ids_unique", len({row["edge_id"] for row in edges}) == 9)
check("part1_local_support", any(row["edge_id"] == "E01" and row["current_edge_decision"] == "SUPPORTED_ASSOCIATION_ONLY" for row in edges))
check("direct_mta_not_supported", any(row["edge_id"] == "E03" and row["current_edge_decision"] == "CURRENT_DATA_DO_NOT_SUPPORT" for row in edges))
check("mta_spp1_pending", any(row["edge_id"] == "E04" and row["current_edge_decision"] == "NOT_CLOSED_PENDING_CLIENT" for row in edges))
check("spatial_downstream_specific", any(row["edge_id"] == "E06" and row["current_edge_decision"] == "DATASET_SPECIFIC_SPATIAL_SUPPORT" for row in edges))
check("complete_chain_not_supported", any(row["edge_id"] == "E08" and row["current_edge_decision"] == "NOT_SUPPORTED" for row in edges))
check("clinical_not_supported", any(row["edge_id"] == "E09" and row["current_edge_decision"] == "CURRENT_DATA_DO_NOT_SUPPORT_CLINICAL_TRANSLATION" for row in edges))

logic = (ROOT / "PHASE107_EVIDENCE_ESCALATION_LOGIC_v1.md").read_text(encoding="utf-8")
blueprint = (ROOT / "CLIENT_RESPONSE_BLUEPRINT_v3.md").read_text(encoding="utf-8")
summary = (ROOT / "PHASE107_EVIDENCE_ESCALATION_AUDIT_SUMMARY_v1.md").read_text(encoding="utf-8")
check("logic_connected_not_pooled", "connected to gaps exposed by Part 1" in logic and "cannot be pooled" in logic)
check("logic_part1_only_downgrade", "full chain must be downgraded" in logic)
check("blueprint_explains_history", "重做后的第一部分只支持了F1–TYMP等局部关联" in blueprint)
check("blueprint_no_combined_statistics", "不是可以合并统计的同一队列" in blueprint)
check("blueprint_two_options", "Option A" in blueprint and "Option B" in blueprint)
check("summary_complete_chain_not_established", "The complete chain is therefore not established" in summary)
check("summary_part1_only_consequence", "conclusion must be reduced" in summary)

failed = [row for row in checks if not row["passed"]]
payload = {
    "status": "PASS_PHASE107_EVIDENCE_ESCALATION_VALIDATION" if not failed else "FAIL_PHASE107_EVIDENCE_ESCALATION_VALIDATION",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "failed_checks": failed,
}
(ROOT / "PHASE107_EVIDENCE_ESCALATION_VALIDATION_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

manifest_rows = []
for path in sorted(item for item in ROOT.iterdir() if item.is_file() and item.name != "FILE_MANIFEST_SHA256_v1.tsv"):
    manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
with (ROOT / "FILE_MANIFEST_SHA256_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

print(json.dumps(payload, ensure_ascii=False, indent=2))
raise SystemExit(0 if not failed else 1)
