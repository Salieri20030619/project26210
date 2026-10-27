#!/usr/bin/env python3
"""Read-only post-install verification for the two Part 4 v3 reports."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3")
DOCS = PROJECT / "docs/client_reports/part4_spatial"
QC = PROJECT / "results/part4_spatial/tymp_vascular_reports_v3_qc"
FULL = DOCS / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3.docx"
CURATED = DOCS / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3.docx"
FULL_V2 = DOCS / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
CURATED_V2 = DOCS / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
FULL_PDF = QC / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3.pdf"
CURATED_PDF = QC / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3.pdf"
OUT_TSV = STAGE / "post_install_verification_v3.tsv"
OUT_JSON = STAGE / "post_install_verification_summary_v3.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


checks: list[dict[str, str]] = []


def add(name: str, passed: bool, observed: str, expected: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})


expected_hashes = {
    FULL: "d9ed4172e3a9b675946cdf0fa0b26c4d44154c33e40c72048001c2dee2b19940",
    CURATED: "919b0af905bd86f5742649699b428c24d612b6fa8ae149f9f75f2bb6791cde5d",
    FULL_V2: "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b",
    CURATED_V2: "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39",
    FULL_PDF: "52b484ead1ccc10e635bbfe5e374e6510a37624deefd8fb879aa4a713db2e22c",
    CURATED_PDF: "f17221f2ff0ff9c6a6b5f16f05bc0508bfb53602414cf5ccdfe3823be3532e37",
}
for path, expected in expected_hashes.items():
    observed = sha256(path) if path.is_file() else "missing"
    add(f"hash_{path.name}", observed == expected, observed, expected)

required = {
    "README.md",
    "build_part4_vascular_reports_v3.py",
    "verify_part4_vascular_reports_v3.py",
    "verify_rendered_reports_v3.py",
    "install_part4_vascular_reports_v3.py",
    "report_generation_metadata_v3.json",
    "report_install_metadata_v3.json",
    "independent_verification_v3.tsv",
    "independent_verification_summary_v3.json",
    "render_verification_v3.tsv",
    "render_verification_summary_v3.json",
    "report_render_manual_review_v3.tsv",
}
present = {path.name for path in QC.iterdir()} if QC.is_dir() else set()
missing = sorted(required - present)
add("qc_required_files", not missing, ",".join(missing) or "none missing", "none missing")

metadata = json.loads((QC / "report_install_metadata_v3.json").read_text(encoding="utf-8"))
add("install_metadata_candidate", metadata.get("status") == "candidate", str(metadata.get("status")), "candidate")
add("install_metadata_qc", metadata.get("qc") == {"docx": "50/50 PASS", "pdf": "29/29 PASS", "manual": "11/11 page groups PASS"}, json.dumps(metadata.get("qc"), ensure_ascii=False), "50/50;29/29;11/11 PASS")

status = (PROJECT / "PROJECT_STATUS.md").read_text(encoding="utf-8")
status_ok = "Part 4完整报告与TYMP精选报告v3已生成" in status and FULL.name in status and CURATED.name in status and "candidate" in status
add("project_status_record", status_ok, str(status_ok), "True")

handoff = (PROJECT / "CODEX_HANDOFF.md").read_text(encoding="utf-8")
handoff_ok = "Part 4完整与精选报告v3正式候选" in handoff and str(FULL) in handoff and str(CURATED) in handoff
add("handoff_record", handoff_ok, str(handoff_ok), "True")

manifest_path = PROJECT / "RESULTS_MANIFEST.tsv"
with manifest_path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    rows = list(reader)
field_count = len(reader.fieldnames or [])
add("manifest_18_columns", field_count == 18 and all(len(row) == 18 for row in rows), f"header={field_count};rows={len(rows)}", "header=18;all rows=18")
ids = [row["result_id"] for row in rows]
add("manifest_ids_unique", len(ids) == len(set(ids)), f"rows={len(ids)};unique={len(set(ids))}", "equal")
expected_ids = {
    "p4_tymp_vascular_full_report_v3_20260804",
    "p4_tymp_vascular_full_report_pdf_v3_20260804",
    "p4_tymp_vascular_curated_report_v3_20260804",
    "p4_tymp_vascular_curated_report_pdf_v3_20260804",
    "p4_tymp_vascular_report_independent_qc_v3_20260804",
    "p4_tymp_vascular_report_render_qc_v3_20260804",
    "p4_tymp_vascular_report_manual_qc_v3_20260804",
}
found_ids = expected_ids & set(ids)
add("manifest_v3_ids", found_ids == expected_ids, ",".join(sorted(found_ids)), ",".join(sorted(expected_ids)))
selected = {row["result_id"]: row for row in rows if row["result_id"] in expected_ids}
paths_ok = all(Path(row["path"]).exists() for row in selected.values()) and len(selected) == len(expected_ids)
add("manifest_v3_paths_exist", paths_ok, str(paths_ok), "True")

with (QC / "report_render_manual_review_v3.tsv").open(encoding="utf-8", newline="") as handle:
    manual = list(csv.DictReader(handle, delimiter="\t"))
add("manual_review_all_pass", len(manual) == 11 and all(row["status"] == "PASS" for row in manual), f"rows={len(manual)};statuses={','.join(row['status'] for row in manual)}", "11 rows;all PASS")

with OUT_TSV.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected"], delimiter="\t")
    writer.writeheader()
    writer.writerows(checks)

failed = [row for row in checks if row["status"] != "PASS"]
summary = {
    "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
    "status": "passed" if not failed else "failed",
    "checks_total": len(checks),
    "checks_passed": len(checks) - len(failed),
    "checks_failed": len(failed),
    "full_docx": str(FULL),
    "curated_docx": str(CURATED),
    "qc_directory": str(QC),
    "failed_checks": failed,
}
OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
raise SystemExit(1 if failed else 0)
