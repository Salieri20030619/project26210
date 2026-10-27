#!/usr/bin/env python3
"""Back up and atomically replace the single Part 4 client-report candidate."""

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
BASELINE = REPORT_DIR / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版.docx"
QC_DIR = ROOT / "results/part4_spatial/advanced_patch_report_qc_v1"
BACKUP = QC_DIR / "source_before_advanced_patch_v1.docx"
FINAL_META = REPORT_DIR / "project26210_part4_advanced_patch_v1_report_metadata.json"

STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
CANDIDATE = STAGE / "project26210_part4_advanced_patch_v1_candidate.docx"
STAGED_META = STAGE / "project26210_part4_advanced_patch_v1_report_metadata.json"
STAGED_VERIFY_JSON = STAGE / "project26210_part4_advanced_patch_v1_verification.json"
STAGED_VERIFY_TSV = STAGE / "project26210_part4_advanced_patch_v1_verification.tsv"
STAGED_VISUAL = STAGE / "project26210_part4_advanced_patch_v1_manual_visual_review.tsv"
STAGED_PDF = STAGE / "render/project26210_part4_advanced_patch_v1_candidate.pdf"
STAGED_CONTACT = STAGE / "render/advanced_patch_pages_contact_v2.png"
SCRIPT_DIR = Path("/mnt/BioDisk/wangmingyuan/project26210_misty_stage")
UPDATER = SCRIPT_DIR / "update_part4_report_advanced_patch_v1.py"
VERIFIER = SCRIPT_DIR / "verify_part4_report_advanced_patch_v1.py"
FINALIZER = SCRIPT_DIR / "finalize_part4_report_advanced_patch_v1.py"

EXPECTED_BEFORE = "06332668ac51b1513eed799ac802f747be937695275cd90846f9677b4ab1315a"
EXPECTED_AFTER = "646b1abf7d7580e9770e60021bc26634bace6055f2b969e0cbb0b6c13919ec0d"
EXPECTED_BASELINE = "ecdd12f18f290229593a531e05e37d27677ff192043ef9abc2e336e677faf1c0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


required = [
    REPORT, BASELINE, CANDIDATE, STAGED_META, STAGED_VERIFY_JSON, STAGED_VERIFY_TSV,
    STAGED_VISUAL, STAGED_PDF, STAGED_CONTACT, UPDATER, VERIFIER, FINALIZER,
]
for path in required:
    if not path.is_file():
        raise FileNotFoundError(path)
if sha256(REPORT) != EXPECTED_BEFORE:
    raise RuntimeError("Official report changed since staging; refusing to overwrite")
if sha256(CANDIDATE) != EXPECTED_AFTER:
    raise RuntimeError("Candidate hash differs from the verified candidate")
if sha256(BASELINE) != EXPECTED_BASELINE:
    raise RuntimeError("Read-only baseline report changed unexpectedly")
verify = json.loads(STAGED_VERIFY_JSON.read_text())
if verify.get("status") != "PASS" or verify.get("checks_failed") != 0 or verify.get("checks_total") != 56:
    raise RuntimeError("Independent staged verification did not pass 56/56")
visual_lines = STAGED_VISUAL.read_text().strip().splitlines()
if len(visual_lines) != 14 or any("\tPASS\t" not in line for line in visual_lines[1:]):
    raise RuntimeError("Manual visual review is incomplete")
with ZipFile(CANDIDATE) as archive:
    if archive.testzip() is not None:
        raise RuntimeError("Candidate DOCX ZIP test failed")

QC_DIR.mkdir(parents=True, exist_ok=True)
if BACKUP.exists() and sha256(BACKUP) != EXPECTED_BEFORE:
    raise RuntimeError("Existing pre-update backup has an unexpected hash")
if not BACKUP.exists():
    shutil.copy2(REPORT, BACKUP)
if sha256(BACKUP) != EXPECTED_BEFORE:
    raise RuntimeError("Pre-update backup verification failed")

temp_report = REPORT_DIR / ".project26210_part4_advanced_patch_v1.tmp.docx"
if temp_report.exists():
    raise RuntimeError(f"Unexpected temporary report exists: {temp_report}")
shutil.copy2(CANDIDATE, temp_report)
if sha256(temp_report) != EXPECTED_AFTER:
    temp_report.unlink(missing_ok=True)
    raise RuntimeError("Temporary official-report copy failed hash verification")
os.replace(temp_report, REPORT)
if sha256(REPORT) != EXPECTED_AFTER:
    raise RuntimeError("Official report hash mismatch after atomic replacement")
if sha256(BASELINE) != EXPECTED_BASELINE:
    raise RuntimeError("Read-only baseline report changed during finalization")

copy_map = {
    UPDATER: QC_DIR / UPDATER.name,
    VERIFIER: QC_DIR / VERIFIER.name,
    FINALIZER: QC_DIR / FINALIZER.name,
    STAGED_META: QC_DIR / STAGED_META.name,
    STAGED_VERIFY_JSON: QC_DIR / STAGED_VERIFY_JSON.name,
    STAGED_VERIFY_TSV: QC_DIR / STAGED_VERIFY_TSV.name,
    STAGED_VISUAL: QC_DIR / STAGED_VISUAL.name,
    STAGED_PDF: QC_DIR / "project26210_part4_advanced_patch_v1_rendered.pdf",
    STAGED_CONTACT: QC_DIR / "project26210_part4_advanced_patch_v1_pages_39_49_contact.png",
}
for source, target in copy_map.items():
    shutil.copy2(source, target)
    if sha256(source) != sha256(target):
        raise RuntimeError(f"QC artifact copy failed: {source}")

now = datetime.now().astimezone().isoformat(timespec="seconds")
staged_meta = json.loads(STAGED_META.read_text())
final_meta = {
    **staged_meta,
    "status": "finalized_and_verified",
    "finalized_at": now,
    "official_report": str(REPORT),
    "official_size": REPORT.stat().st_size,
    "official_sha256": sha256(REPORT),
    "pre_update_backup": str(BACKUP),
    "pre_update_size": BACKUP.stat().st_size,
    "pre_update_sha256": sha256(BACKUP),
    "baseline_report": str(BASELINE),
    "baseline_sha256_unchanged": sha256(BASELINE),
    "independent_verification": "56/56 PASS",
    "manual_visual_review": "13/13 scopes PASS; rendered pages 39-49 inspected",
    "rendered_pdf_pages": 53,
}
FINAL_META.write_text(json.dumps(final_meta, ensure_ascii=False, indent=2) + "\n")
shutil.copy2(FINAL_META, QC_DIR / FINAL_META.name)

audit_files = [BACKUP, REPORT, FINAL_META, *copy_map.values()]
audit = {
    "status": "PASS",
    "finalized_at": now,
    "before": {"path": str(BACKUP), "size": BACKUP.stat().st_size, "sha256": sha256(BACKUP)},
    "after": {"path": str(REPORT), "size": REPORT.stat().st_size, "sha256": sha256(REPORT)},
    "baseline_unchanged": {"path": str(BASELINE), "size": BASELINE.stat().st_size, "sha256": sha256(BASELINE)},
    "structural_and_content_verification": "56/56 PASS",
    "manual_visual_review": "13/13 PASS",
    "rendered_pdf_pages": 53,
    "audit_files": [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in audit_files
    ],
}
(QC_DIR / "report_finalization_audit_v1.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(audit, ensure_ascii=False, indent=2))
