from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(r"C:\Users\31937\project26210_phase105_full_analysis_audit_v1")


def read_tsv(name: str) -> list[dict[str, str]]:
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
    "PHASE105_CLIENT_REQUIREMENT_LOCK_v1.md",
    "CLIENT_RESPONSE_BLUEPRINT_v1.md",
    "PHASE105_FULL_ANALYSIS_AUDIT_SUMMARY_v1.md",
    "ANALYSIS_MODULE_LEDGER_v1.tsv",
    "CLIENT_INTERNAL_MAINLINE_EVIDENCE_v1.tsv",
    "scan_summary_v1.json",
    "phase_presence_matrix_v1.tsv",
    "top_level_phase_directory_inventory_v1.tsv",
    "manifest_stage_summary_v1.tsv",
    "manifest_path_audit_v1.tsv",
    "sgs10_scan_summary_v1.json",
    "sgs10_project26210_directory_inventory_v1.tsv",
    "sgs10_required_evidence_audit_v1.tsv",
    "sgs10_project26210_process_snapshot_v1.txt",
]
for name in required:
    path = ROOT / name
    check(f"required_{name}", path.is_file() and path.stat().st_size > 0, path)

scan = json.loads((ROOT / "scan_summary_v1.json").read_text(encoding="utf-8"))
check("scan_status", scan.get("status") == "PASS_SCAN_COMPLETE", scan.get("status"))
check("manifest_rows_574", scan.get("manifest_rows") == 574, scan.get("manifest_rows"))
check("manifest_ids_unique", scan.get("manifest_unique_ids") == 574, scan.get("manifest_unique_ids"))
check("local_phase_dirs_91", scan.get("top_level_phase_directories") == 91, scan.get("top_level_phase_directories"))
check("manifest_no_missing_local_paths", scan.get("missing_local_manifest_paths") == 0)
check("unregistered_directory_count_16", scan.get("directory_without_manifest_count") == 16)

presence = read_tsv("phase_presence_matrix_v1.tsv")
check("presence_rows_104", len(presence) == 104, len(presence))
check("phase21_registration_gap_captured", any(row["phase"] == "21" and "DIRECTORY_WITHOUT_MANIFEST_STAGE" in row["audit_flags"] for row in presence))
check("phase57_registration_gap_captured", any(row["phase"] == "57" and "DIRECTORY_WITHOUT_MANIFEST_STAGE" in row["audit_flags"] for row in presence))
check("phase58_registration_gap_captured", any(row["phase"] == "58" and "DIRECTORY_WITHOUT_MANIFEST_STAGE" in row["audit_flags"] for row in presence))
check("phase92_registration_gap_captured", any(row["phase"] == "92" and "DIRECTORY_WITHOUT_MANIFEST_STAGE" in row["audit_flags"] for row in presence))

remote = json.loads((ROOT / "sgs10_scan_summary_v1.json").read_text(encoding="utf-8"))
check("remote_scan_status", remote.get("status") == "PASS_SGS10_SCAN_COMPLETE")
check("remote_directories_81", remote.get("remote_project_directories") == 81, remote.get("remote_project_directories"))
check("remote_required_missing_zero", remote.get("required_evidence_missing") == 0)
required_remote = read_tsv("sgs10_required_evidence_audit_v1.tsv")
check("remote_required_7_of_7", len(required_remote) == 7 and all(row["status"] == "PASS" for row in required_remote), len(required_remote))

ledger = read_tsv("ANALYSIS_MODULE_LEDGER_v1.tsv")
ledger_ids = [row["module_id"] for row in ledger]
check("ledger_rows_31", len(ledger) == 31, len(ledger))
check("ledger_ids_unique", len(ledger_ids) == len(set(ledger_ids)))
check("ledger_phase21_main", any(row["phase_scope"] == "Phase 21" and row["client_mainline_use"] == "MAIN_INTERNAL" for row in ledger))
check("ledger_phase33_main", any(row["phase_scope"] == "Phase 33" and row["client_mainline_use"] == "MAIN_INTERNAL" for row in ledger))
check("ledger_phase34_main", any(row["phase_scope"] == "Phase 34" and row["client_mainline_use"] == "MAIN_INTERNAL" for row in ledger))
check("ledger_phase35_gap", any(row["phase_scope"] == "Phase 35" and "fine-state rerun" in row["reason_or_boundary"] for row in ledger))
check("ledger_external_reports_superseded", any(row["phase_scope"] == "Phase 89-94" and row["client_mainline_use"] == "SUPERSEDED_FOR_CURRENT_REQUEST" for row in ledger))

evidence = read_tsv("CLIENT_INTERNAL_MAINLINE_EVIDENCE_v1.tsv")
evidence_ids = [row["evidence_id"] for row in evidence]
check("internal_evidence_rows_10", len(evidence) == 10, len(evidence))
check("internal_evidence_ids_unique", len(evidence_ids) == len(set(evidence_ids)))
check("internal_evidence_no_external_cohorts", all("GSE264586" not in row["internal_scope"] and "GSE294109" not in row["internal_scope"] for row in evidence))
check("macrophage_subdivision_counts", any(row["evidence_id"] == "I01" and "SPP1-state 15322" in row["result"] for row in evidence))
check("internal_fine_state_gap_explicit", any(row["evidence_id"] == "I10" and row["evidence_class"] == "UNRESOLVED_REQUIRED_GAP" for row in evidence))

lock = (ROOT / "PHASE105_CLIENT_REQUIREMENT_LOCK_v1.md").read_text(encoding="utf-8")
for dataset in ["GSE156632", "GSE159115", "GSE178481", "GSE207493", "GSE242299"]:
    check(f"lock_internal_{dataset}", dataset in lock)
check("lock_external_not_main", "External validation cohorts are not required" in lock)
check("lock_spp1_after_subdivision", "SPP1 must be introduced only after" in lock)
check("lock_client_data_excluded", "must not request, import, reanalyse, redraw" in lock)

blueprint = (ROOT / "CLIENT_RESPONSE_BLUEPRINT_v1.md").read_text(encoding="utf-8")
check("blueprint_client_reply", "收到。后续内容会按您的意见全部重构" in blueprint)
check("blueprint_order", "免疫细胞图谱 → 髓系细胞提取 → 巨噬细胞细分 → SPP1" in blueprint)
check("blueprint_gap", "internal-five-study, macrophage-fine-state-specific communication analysis" in blueprint)

summary = (ROOT / "PHASE105_FULL_ANALYSIS_AUDIT_SUMMARY_v1.md").read_text(encoding="utf-8")
check("summary_not_absence", "not absence of analysis" in summary)
check("summary_mainline_counts", "Other Macrophage: 40,429" in summary and "SPP1-state Macrophage: 15,322" in summary)
check("summary_phase21_scope", "50,761 myeloid cells from 40 patients" in summary)
check("summary_client_data_boundary", "client's unpublished validation experiments" in summary)

failed = [row for row in checks if not row["passed"]]
payload = {
    "status": "PASS_PHASE105_AUDIT_VALIDATION" if not failed else "FAIL_PHASE105_AUDIT_VALIDATION",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "failed_checks": failed,
}
(ROOT / "PHASE105_AUDIT_VALIDATION_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

manifest_rows = []
for path in sorted(item for item in ROOT.iterdir() if item.is_file() and item.name != "FILE_MANIFEST_SHA256_v1.tsv"):
    manifest_rows.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
with (ROOT / "FILE_MANIFEST_SHA256_v1.tsv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["file", "size_bytes", "sha256"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(manifest_rows)

print(json.dumps(payload, ensure_ascii=False, indent=2))
raise SystemExit(0 if not failed else 1)
