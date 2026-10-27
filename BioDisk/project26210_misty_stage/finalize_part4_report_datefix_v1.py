#!/usr/bin/env python3
"""Finalize the one-paragraph report-date correction with a recoverable backup."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
REPORT_DIR = ROOT / "docs/client_reports/part4_spatial"
REPORT = REPORT_DIR / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
META = REPORT_DIR / "project26210_part4_advanced_patch_v1_report_metadata.json"
QC = ROOT / "results/part4_spatial/advanced_patch_report_qc_v1"
BACKUP = QC / "source_before_report_date_fix_v1.docx"
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
CANDIDATE = STAGE / "project26210_part4_advanced_patch_v1_datefix_candidate.docx"
DATE_AUDIT = STAGE / "project26210_part4_advanced_patch_v1_datefix_audit.json"
DATE_VISUAL = STAGE / "project26210_part4_advanced_patch_v1_datefix_visual_review.tsv"
DATE_PDF = STAGE / "render/project26210_part4_advanced_patch_v1_datefix_candidate.pdf"
DATE_PAGE = STAGE / "render/datefix_page-01.png"
SCRIPT_DIR = Path("/mnt/BioDisk/wangmingyuan/project26210_misty_stage")
FIX_SCRIPT = SCRIPT_DIR / "fix_part4_report_date_v1.py"
FINALIZE_SCRIPT = SCRIPT_DIR / "finalize_part4_report_datefix_v1.py"
EXPECTED_BEFORE = "646b1abf7d7580e9770e60021bc26634bace6055f2b969e0cbb0b6c13919ec0d"
EXPECTED_AFTER = "f3b61f3334c8bdb5dbbd6cb5b674f5ea13acbb3bb7fa0ebc863c4d3b2cb86bf4"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


for path in [REPORT, META, CANDIDATE, DATE_AUDIT, DATE_VISUAL, DATE_PDF, DATE_PAGE, FIX_SCRIPT, FINALIZE_SCRIPT]:
    if not path.is_file():
        raise FileNotFoundError(path)
if sha256(REPORT) != EXPECTED_BEFORE or sha256(CANDIDATE) != EXPECTED_AFTER:
    raise RuntimeError("Report or date-fix candidate hash mismatch")
audit = json.loads(DATE_AUDIT.read_text())
if audit.get("status") != "PASS" or audit.get("changed_members") != ["word/document.xml"] or audit.get("changed_paragraphs") != 1:
    raise RuntimeError("Date-fix audit did not pass")
if "\tPASS\t" not in DATE_VISUAL.read_text():
    raise RuntimeError("Date-fix visual review missing")
with ZipFile(CANDIDATE) as archive:
    if archive.testzip() is not None:
        raise RuntimeError("Date-fix candidate failed ZIP test")

QC.mkdir(parents=True, exist_ok=True)
if BACKUP.exists() and sha256(BACKUP) != EXPECTED_BEFORE:
    raise RuntimeError("Existing date-fix backup hash mismatch")
if not BACKUP.exists():
    shutil.copy2(REPORT, BACKUP)
if sha256(BACKUP) != EXPECTED_BEFORE:
    raise RuntimeError("Date-fix backup verification failed")

temp = REPORT_DIR / ".project26210_part4_datefix_v1.tmp.docx"
if temp.exists():
    raise RuntimeError(f"Unexpected temporary file: {temp}")
shutil.copy2(CANDIDATE, temp)
if sha256(temp) != EXPECTED_AFTER:
    temp.unlink(missing_ok=True)
    raise RuntimeError("Temporary date-fix copy verification failed")
os.replace(temp, REPORT)
if sha256(REPORT) != EXPECTED_AFTER:
    raise RuntimeError("Official date-fixed report hash mismatch")

copy_map = {
    FIX_SCRIPT: QC / FIX_SCRIPT.name,
    FINALIZE_SCRIPT: QC / FINALIZE_SCRIPT.name,
    DATE_AUDIT: QC / DATE_AUDIT.name,
    DATE_VISUAL: QC / DATE_VISUAL.name,
    DATE_PDF: QC / "project26210_part4_advanced_patch_v1_datefix_rendered.pdf",
    DATE_PAGE: QC / "project26210_part4_advanced_patch_v1_datefix_page_1.png",
}
for source, target in copy_map.items():
    shutil.copy2(source, target)
    if sha256(source) != sha256(target):
        raise RuntimeError(f"Date-fix QC copy failed: {source}")

metadata = json.loads(META.read_text())
metadata.update({
    "status": "finalized_verified_and_date_corrected",
    "report_date": "2026-08-02",
    "report_date_label": "2026年8月2日（高级分析补丁更新）",
    "date_corrected_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "pre_date_fix_sha256": EXPECTED_BEFORE,
    "official_size": REPORT.stat().st_size,
    "official_sha256": sha256(REPORT),
    "date_fix_changed_members": ["word/document.xml"],
    "date_fix_changed_paragraphs": 1,
    "date_fix_visual_review": "2/2 PASS; all 53 pages retained; page 1 inspected",
})
META.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
shutil.copy2(META, QC / META.name)

final_audit = {
    "status": "PASS",
    "finalized_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "before": {"path": str(BACKUP), "size": BACKUP.stat().st_size, "sha256": sha256(BACKUP)},
    "after": {"path": str(REPORT), "size": REPORT.stat().st_size, "sha256": sha256(REPORT)},
    "date_text": "2026年8月2日（高级分析补丁更新）",
    "changed_members": ["word/document.xml"],
    "changed_paragraphs": 1,
    "zip_test": "PASS",
    "rendered_pdf_pages": 53,
    "visual_review": "2/2 PASS",
    "metadata": {"path": str(META), "size": META.stat().st_size, "sha256": sha256(META)},
    "qc_files": [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in copy_map.values()
    ],
}
(QC / "report_datefix_finalization_audit_v1.json").write_text(json.dumps(final_audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(final_audit, ensure_ascii=False, indent=2))
