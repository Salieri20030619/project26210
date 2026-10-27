#!/usr/bin/env python3
from __future__ import annotations

import csv
import fcntl
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile

ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
REPORT = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage")
CANDIDATE = STAGE / "project26210_part4_full_advanced_tymp_v2_candidate.docx"
PDF = STAGE / "rendered_v2/project26210_part4_full_advanced_tymp_v2_candidate.pdf"
META = STAGE / "project26210_part4_full_advanced_tymp_v2_metadata.json"
VERIFY = STAGE / "independent_verification_v2.tsv"
VERIFY_SUMMARY = STAGE / "independent_verification_summary_v2.json"
MANUAL = STAGE / "report_render_manual_review_v2.tsv"
RENDER_SUMMARY = STAGE / "report_render_qc_summary_v2.json"
QC = ROOT / "results/part4_spatial/full_advanced_tymp_report_qc_v1"
EXPECTED_SOURCE = "f3b61f3334c8bdb5dbbd6cb5b674f5ea13acbb3bb7fa0ebc863c4d3b2cb86bf4"
EXPECTED_CANDIDATE = "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b"
EXPECTED_PDF = "95e19e9d592cd1b3e0ef147e728fec061da16161d4b5c5379be4f71039888159"
MARKER = "2026-08-03 Part 4全12例高级分析与TYMP状态报告整合终态"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def atomic_write(path: Path, data: bytes) -> None:
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
        temp = Path(handle.name)
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)

def copy_artifact(source: Path, destination: Path) -> None:
    shutil.copy2(source, destination)
    if sha256(source) != sha256(destination):
        raise RuntimeError(f"Copy hash mismatch: {source}")

for path in [REPORT, CANDIDATE, PDF, META, VERIFY, VERIFY_SUMMARY, MANUAL, RENDER_SUMMARY]:
    if not path.is_file():
        raise FileNotFoundError(path)
if QC.exists():
    raise FileExistsError(f"QC destination already exists: {QC}")
if sha256(REPORT) != EXPECTED_SOURCE:
    raise RuntimeError("Official report changed since staging; refusing replacement")
if sha256(CANDIDATE) != EXPECTED_CANDIDATE or sha256(PDF) != EXPECTED_PDF:
    raise RuntimeError("Staged candidate or PDF hash mismatch")
with ZipFile(CANDIDATE) as archive:
    if archive.testzip() is not None:
        raise RuntimeError("Candidate ZIP CRC failed")
verification = json.loads(VERIFY_SUMMARY.read_text(encoding="utf-8"))
render = json.loads(RENDER_SUMMARY.read_text(encoding="utf-8"))
if verification.get("status") != "passed" or verification.get("passed") != verification.get("total") or verification.get("total") != 60:
    raise RuntimeError("Independent verification is not 60/60 PASS")
if render.get("status") != "passed" or render.get("manual_review_records_passed") != 7:
    raise RuntimeError("Rendered report review is not 7/7 PASS")
with MANUAL.open(encoding="utf-8", newline="") as handle:
    manual_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(manual_rows) != 7 or {row["status"] for row in manual_rows} != {"PASS"}:
    raise RuntimeError("Manual review table is not 7/7 PASS")

status_path = ROOT / "PROJECT_STATUS.md"
handoff_path = ROOT / "CODEX_HANDOFF.md"
manifest_path = ROOT / "RESULTS_MANIFEST.tsv"
lock_path = ROOT / ".part4_full_advanced_tymp_report_install_v1.lock"
with lock_path.open("w") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    status_before = status_path.read_bytes()
    handoff_before = handoff_path.read_bytes()
    manifest_before = manifest_path.read_bytes()
    if MARKER.encode() in status_before or MARKER.encode() in handoff_before:
        raise RuntimeError("Project records already contain install marker")

    QC.mkdir(parents=True)
    copy_artifact(REPORT, QC / "source_before_full_advanced_tymp_v1.docx")
    (QC / "project_records_before_install_v1").mkdir()
    atomic_write(QC / "project_records_before_install_v1/PROJECT_STATUS.before.md", status_before)
    atomic_write(QC / "project_records_before_install_v1/CODEX_HANDOFF.before.md", handoff_before)
    atomic_write(QC / "project_records_before_install_v1/RESULTS_MANIFEST.before.tsv", manifest_before)
    artifacts = [
        (CANDIDATE, "final_candidate_v2.docx"),
        (PDF, "rendered_report_v2.pdf"),
        (META, "report_generation_metadata_v2.json"),
        (VERIFY, "independent_verification_v2.tsv"),
        (VERIFY_SUMMARY, "independent_verification_summary_v2.json"),
        (MANUAL, "report_render_manual_review_v2.tsv"),
        (RENDER_SUMMARY, "report_render_qc_summary_v2.json"),
        (STAGE.parent / "project26210_part4_report_update_stage/ooxml_report_utils_v1.py", "ooxml_report_utils_v1.py"),
        (STAGE.parent / "project26210_part4_report_update_stage/update_part4_report_full_advanced_tymp_v1.py", "update_part4_report_full_advanced_tymp_v1.py"),
        (STAGE.parent / "project26210_part4_report_update_stage/verify_part4_report_full_advanced_tymp_v1.py", "verify_part4_report_full_advanced_tymp_v1.py"),
        (STAGE.parent / "project26210_part4_report_update_stage/make_report_page_contact_sheets_v1.py", "make_report_page_contact_sheets_v1.py"),
        (STAGE.parent / "project26210_part4_report_update_stage/install_part4_full_advanced_tymp_report_v1.py", "install_part4_full_advanced_tymp_report_v1.py"),
    ]
    for source, name in artifacts:
        copy_artifact(source, QC / name)
    shutil.copytree(STAGE / "contact_sheets_v2", QC / "contact_sheets_v2")

    with tempfile.NamedTemporaryFile(dir=REPORT.parent, prefix=".part4_full_advanced_tymp_v1.", suffix=".docx", delete=False) as handle:
        report_temp = Path(handle.name)
    shutil.copy2(CANDIDATE, report_temp)
    if sha256(report_temp) != EXPECTED_CANDIDATE:
        raise RuntimeError("Temporary report hash mismatch")
    os.replace(report_temp, REPORT)
    if sha256(REPORT) != EXPECTED_CANDIDATE:
        raise RuntimeError("Installed report hash mismatch")

    now = datetime.now().astimezone()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S %z")
    status_addition = f"""

## 15. {MARKER}

- 唯一Part 4增强版主报告已从SHA-256 `{EXPECTED_SOURCE}` 更新为 `{EXPECTED_CANDIDATE}`，大小由31,463,480 bytes变为33,702,763 bytes；更新前文件已备份到 `results/part4_spatial/full_advanced_tymp_report_qc_v1/source_before_full_advanced_tymp_v1.docx`。原26 MB基线报告未修改。
- 新增第4.13节（TJ-RCC 12例BayesSpace、SPARK-X/Hotspot、MISTy、靶向COMMOT全套高级分析）、第4.14节（TYMP–MTA-context–SPP1巨噬细胞状态与探索性中介/下游表型）、图4.38–4.49、6张汇总表，并同步更新封面、摘要、第5章结论/限制/状态表和最终总述。
- 高级队列已完成12/12患者、34,609 spots、733/733独立检查和56/56原始图检；状态分析完成211/211独立检查和16/16原始图检。报告DOCX另行独立复核60/60通过；LibreOffice渲染71页，封面、目录和新增页52–71共7/7版面审查记录通过。
- 结论边界保持：客户湿实验作为TYMP–5′-MTA–SPP1机制的因果锚点；空间数据支持人组织相关性与患者异质性。参考轻依赖无SPP1巨噬细胞状态两尺度10/12正向但全局FDR=0.154；SPP1相对C1QC对比主要为负，MTA转录代理中介未获支持，不能升级为独立代谢中介或促恶性因果证明。
- 图4.45初版预览出现左侧裁切，已通过缩窄图片修复；修正版第61页图题和全部纵轴标签完整。全部生成脚本、输入哈希、DOCX核验、PDF、人工版式记录和安装前项目记录均归档于 `results/part4_spatial/full_advanced_tymp_report_qc_v1`。
- 报告仍为待交付候选；服务器没有外发或客户签收证据，`report_status`保持 `candidate`。

## 16. 最近更新时间

{timestamp}（服务器实际时间）
"""
    handoff_addition = f"""

## {MARKER}

- 正式主报告：`{REPORT}`；最终大小33,702,763 bytes，SHA-256=`{EXPECTED_CANDIDATE}`。更新前版本SHA-256=`{EXPECTED_SOURCE}`已备份，26 MB基线报告保持不变。
- 正文新增4.13、4.14节及图4.38–4.49，覆盖12例四类高级分析、TYMP–MTA-context–SPP1状态路线、探索性中介和下游表型；封面、摘要、第5章结论、限制、状态表和最终总述同步修订。
- 数值/结构核验60/60通过；PDF 71页，7/7版式检查范围通过。v1图4.45左溢出已在v2修正，最终DOCX/PDF及全部审计位于`{QC}`。
- 对客户方向的最终证据层级：湿实验是因果锚点；空间结果支持人组织相关性与患者异质性，但不独立证明5′-MTA中介、SPP1特异性、直接通讯或恶性进展因果。报告状态仍为candidate，未外发。

## 本次交接更新时间

{timestamp}（Asia/Shanghai）
"""
    atomic_write(status_path, status_before + status_addition.encode("utf-8"))
    atomic_write(handoff_path, handoff_before + handoff_addition.encode("utf-8"))

    rows = [
        ["p4_full_advanced_tymp_report_v1_20260803", "Part4", "GLMF1;GLMF2;12 TJ-RCC tumor cases", "full advanced and TYMP state integrated report", "client_report", str(REPORT), str(QC / "update_part4_report_full_advanced_tymp_v1.py"), "advanced_full_v2 cohort synthesis;tymp_spp1_state_v1", "batch_bench;LibreOffice7.3", "sections4.13-4.14;figures4.38-4.49;wet-lab causal anchor", "v1", "2026-08-03", str(REPORT.stat().st_size), "passed", "DOCX60/60;render7/7;PDF71pages", "candidate", "", f"SHA-256={EXPECTED_CANDIDATE};not delivered"],
        ["p4_full_advanced_tymp_report_pdf_v1_20260803", "Part4", "GLMF1;GLMF2;12 TJ-RCC tumor cases", "integrated report rendered preview", "rendered_report", str(QC / "rendered_report_v2.pdf"), "LibreOffice 7.3 writer_pdf_Export", str(REPORT), "LibreOffice7.3", "71 pages;visual review cover TOC and pages52-71", "v2", "2026-08-03", str((QC / "rendered_report_v2.pdf").stat().st_size), "passed", "report_render_manual_review_v2.tsv 7/7", "not_used", "", f"SHA-256={EXPECTED_PDF}"],
        ["p4_full_advanced_tymp_report_verification_v2_20260803", "Part4", "integrated Part4 report", "DOCX independent verification", "qc_table", str(QC / "independent_verification_v2.tsv"), str(QC / "verify_part4_report_full_advanced_tymp_v1.py"), "final_candidate_v2.docx;formal cohort tables", "batch_bench", "ZIP;relationships;hashes;tables;figures;source recomputation", "v2", "2026-08-03", str((QC / "independent_verification_v2.tsv").stat().st_size), "passed", "60/60", "not_used", "", f"SHA-256={sha256(QC / 'independent_verification_v2.tsv')}"],
        ["p4_full_advanced_tymp_report_manual_v2_20260803", "Part4", "cover;TOC;pages52-71", "integrated report visual review", "manual_qc", str(QC / "report_render_manual_review_v2.tsv"), "Codex visual inspection", str(QC / "rendered_report_v2.pdf"), "manual visual review", "titles;tables;captions;axes;clipping;page flow", "v2", "2026-08-03", str((QC / "report_render_manual_review_v2.tsv").stat().st_size), "passed", "7/7;figure4.45 v1 overflow corrected in v2", "not_used", "", f"SHA-256={sha256(QC / 'report_render_manual_review_v2.tsv')}"],
        ["p4_full_advanced_tymp_report_metadata_v2_20260803", "Part4", "integrated Part4 report", "report source and embedded-image provenance", "qc_metadata", str(QC / "report_generation_metadata_v2.json"), str(QC / "update_part4_report_full_advanced_tymp_v1.py"), "12 embedded cohort images;9 formal source tables", "batch_bench", "source and image SHA-256;report hash;inference boundary", "v2", "2026-08-03", str((QC / "report_generation_metadata_v2.json").stat().st_size), "passed", "candidate hash and 12 image hashes", "not_used", "", f"SHA-256={sha256(QC / 'report_generation_metadata_v2.json')}"],
    ]
    manifest_text = manifest_before.decode("utf-8")
    if not manifest_text.endswith("\n"):
        manifest_text += "\n"
    manifest_text += "".join("\t".join(row) + "\n" for row in rows)
    atomic_write(manifest_path, manifest_text.encode("utf-8"))

    audit = {
        "status": "installed_and_verified",
        "installed_at": timestamp,
        "official_report": str(REPORT),
        "official_size": REPORT.stat().st_size,
        "official_sha256": sha256(REPORT),
        "source_before_sha256": EXPECTED_SOURCE,
        "backup": str(QC / "source_before_full_advanced_tymp_v1.docx"),
        "backup_sha256": sha256(QC / "source_before_full_advanced_tymp_v1.docx"),
        "independent_verification": "60/60",
        "render_manual_review": "7/7",
        "rendered_pdf_pages": 71,
        "project_record_hashes": {"PROJECT_STATUS.md": sha256(status_path), "CODEX_HANDOFF.md": sha256(handoff_path), "RESULTS_MANIFEST.tsv": sha256(manifest_path)},
    }
    (QC / "final_install_audit_v1.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fcntl.flock(lock, fcntl.LOCK_UN)

print(json.dumps(audit, ensure_ascii=False, indent=2))
