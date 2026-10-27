#!/usr/bin/env python
"""Validate the Phase 24 completion audit without modifying prior outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
AUDIT = ROOT / "ANALYSIS_COMPLETION_AUDIT_v1.tsv"
SOURCE = ROOT / "reference" / "分析待做（新）_reference.docx"
EXPECTED_SOURCE_SHA256 = "a32766f2a477d6112b6b57a51a647c70dc74961830f320e698361dbe74da7753"
EXPECTED_COUNTS = {
    "COMPLETE_VALIDATED": 11,
    "ANALYSIS_COMPLETE_PRESENTATION_PENDING": 4,
    "PARTIAL_CURRENT": 4,
    "LEGACY_ONLY_REDO_REQUIRED": 8,
    "NOT_COMPLETED_CURRENT": 8,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


with AUDIT.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, detail: object) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


check("source_exists", SOURCE.is_file(), str(SOURCE))
check("source_sha256", sha256(SOURCE) == EXPECTED_SOURCE_SHA256, sha256(SOURCE))
check("row_count", len(rows) == 35, len(rows))
ids = [row["audit_id"] for row in rows]
check("unique_ids", len(ids) == len(set(ids)), len(set(ids)))
counts = Counter(row["status"] for row in rows)
check("status_counts", dict(counts) == EXPECTED_COUNTS, dict(counts))
check("all_requirements_present", all(row["checklist_requirement"].strip() for row in rows), len(rows))
check("all_pages_present", all(row["checklist_page"].strip() for row in rows), len(rows))
check("all_evidence_present", all(row["authoritative_evidence"].strip() for row in rows), len(rows))

local_evidence_missing: list[str] = []
for row in rows:
    for item in row["authoritative_evidence"].split(";"):
        path_text = item.strip()
        if not path_text or path_text == "PROJECT_STATUS.md":
            candidate = PROJECT / path_text if path_text else None
        else:
            candidate = PROJECT / path_text
        if candidate is not None and not candidate.exists():
            local_evidence_missing.append(path_text)
check("local_evidence_paths_exist", not local_evidence_missing, sorted(set(local_evidence_missing)))

reportable = counts["COMPLETE_VALIDATED"] + counts["ANALYSIS_COMPLETE_PRESENTATION_PENDING"]
check("reportable_count", reportable == 15, reportable)

status = "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL"
result = {
    "schema_version": "phase24_analysis_completion_audit_validation_v1",
    "status": status,
    "checks_passed": sum(item["status"] == "PASS" for item in checks),
    "checks_total": len(checks),
    "audit_rows": len(rows),
    "status_counts": dict(counts),
    "source_sha256": sha256(SOURCE),
    "checks": checks,
}
(ROOT / "ANALYSIS_COMPLETION_AUDIT_VALIDATION_v1.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(result, ensure_ascii=False, indent=2))
raise SystemExit(0 if status == "PASS" else 1)
