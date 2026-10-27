#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_curated_report_stage_v2")
REPORT_SOURCE = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2_candidate.docx"
PDF_SOURCE = STAGE / "rendered_v2/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2_candidate.pdf"
REPORT_TARGET = ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
RESULT_TARGET = ROOT / "results/part4_spatial/tymp_curated_validation_report_v2"
FIGURE1_PNG = RESULT_TARGET / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png"
FIGURE1_PDF = RESULT_TARGET / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.pdf"
FIGURE2_PNG = RESULT_TARGET / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png"
FIGURE2_PDF = RESULT_TARGET / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.pdf"
REPORT_PDF = RESULT_TARGET / "project26210_TYMP_5MTA_SPP1_mechanism_spatial_atlas_curated_report_v2.pdf"

PROJECT_STATUS = ROOT / "PROJECT_STATUS.md"
HANDOFF = ROOT / "CODEX_HANDOFF.md"
MANIFEST = ROOT / "RESULTS_MANIFEST.tsv"

EXPECTED = {
    REPORT_SOURCE: "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39",
    PDF_SOURCE: "0edf5100e08495e02ac152d12f6c05cd181d4e391117e7dad56999c809ad806c",
    STAGE / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png": "7c70c5dfb4692e0efe5be1201be7f96849e78ec64bae78f82151e08ce890bf76",
    STAGE / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.pdf": "59421b9dd38deba6cabfb161403b0899c192aa0110eaf7c1fcda93f66e899f1d",
    STAGE / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png": "8d3180730ff099d94db0b5d328e99bc55506a9eb80ffbda852e731e31ec4794e",
    STAGE / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.pdf": "c3380dc543603053ad85ed88d1d09b46de3aef68efd96c3b558e3827bc1d2630",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if sha256(source) != sha256(target):
        raise RuntimeError(f"Copy hash mismatch: {source} -> {target}")


def append_text(path: Path, text: str) -> None:
    with path.open("a", encoding="utf-8", newline="") as handle:
        handle.write(text)


for path, expected_hash in EXPECTED.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    actual_hash = sha256(path)
    if actual_hash != expected_hash:
        raise RuntimeError(f"Unexpected hash for {path}: {actual_hash}")
for path in [PROJECT_STATUS, HANDOFF, MANIFEST]:
    if not path.is_file():
        raise FileNotFoundError(path)
if REPORT_TARGET.exists() or RESULT_TARGET.exists():
    raise FileExistsError(f"Refusing to overwrite existing curated delivery: {REPORT_TARGET} or {RESULT_TARGET}")

verification = json.loads((STAGE / "curated_report_independent_verification_summary_v2.json").read_text(encoding="utf-8"))
if verification.get("status") != "passed" or verification.get("passed") != 58 or verification.get("total") != 58:
    raise RuntimeError(f"Independent verification not complete: {verification}")
render_qc = json.loads((STAGE / "curated_report_render_qc_v2.json").read_text(encoding="utf-8"))
if render_qc.get("status") != "passed" or render_qc.get("pdf", {}).get("pages") != 7:
    raise RuntimeError(f"Render QC not complete: {render_qc}")

backup_dir = STAGE / "project_records_before_curated_install_v2"
if backup_dir.exists():
    raise FileExistsError(backup_dir)
backup_dir.mkdir(parents=True)
for source in [PROJECT_STATUS, HANDOFF, MANIFEST]:
    copy(source, backup_dir / source.name)

RESULT_TARGET.mkdir(parents=True)
copy(REPORT_SOURCE, REPORT_TARGET)
copy(PDF_SOURCE, REPORT_PDF)
copy(STAGE / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png", FIGURE1_PNG)
copy(STAGE / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.pdf", FIGURE1_PDF)
copy(STAGE / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png", FIGURE2_PNG)
copy(STAGE / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.pdf", FIGURE2_PDF)

artifact_sources = {
    "panel_source_manifest_v2.tsv": STAGE / "panel_source_manifest_v2.tsv",
    "selected_cohort_evidence_v2.tsv": STAGE / "selected_cohort_evidence_v2.tsv",
    "selected_patient_effect_matrix_v2.tsv": STAGE / "selected_patient_effect_matrix_v2.tsv",
    "representative_patient_selection_v2.tsv": STAGE / "representative_patient_selection_v2.tsv",
    "figure_build_metadata_v2.json": STAGE / "figure_build_metadata_v2.json",
    "curated_report_generation_metadata_v2.json": STAGE / "curated_report_generation_metadata_v2.json",
    "curated_report_independent_verification_v2.tsv": STAGE / "curated_report_independent_verification_v2.tsv",
    "curated_report_independent_verification_summary_v2.json": STAGE / "curated_report_independent_verification_summary_v2.json",
    "curated_report_manual_page_review_v2.tsv": STAGE / "curated_report_manual_page_review_v2.tsv",
    "curated_report_render_qc_v2.json": STAGE / "curated_report_render_qc_v2.json",
    "contact_sheet_pages_01_04_v2.jpg": STAGE / "contact_sheets_v2/pages_01_04.jpg",
    "contact_sheet_pages_05_07_v2.jpg": STAGE / "contact_sheets_v2/pages_05_07.jpg",
    "build_tymp_curated_figures_v2.py": STAGE / "build_tymp_curated_figures_v2.py",
    "build_tymp_curated_report_v2.py": STAGE / "build_tymp_curated_report_v2.py",
    "verify_tymp_curated_report_v2.py": STAGE / "verify_tymp_curated_report_v2.py",
    "make_curated_report_contact_sheets_v2.py": STAGE / "make_curated_report_contact_sheets_v2.py",
    "pagebreak.lua": STAGE / "pagebreak.lua",
    "project26210_tymp_curated_report_v2.md": STAGE / "project26210_tymp_curated_report_v2.md",
}
for name, source in artifact_sources.items():
    if not source.is_file():
        raise FileNotFoundError(source)
    copy(source, RESULT_TARGET / name)
for source in backup_dir.iterdir():
    copy(source, RESULT_TARGET / "project_records_before_install" / source.name)

timestamp = datetime.now().astimezone()
timestamp_text = timestamp.strftime("%Y-%m-%d %H:%M:%S %z")
date_text = timestamp.strftime("%Y-%m-%d")
status_block = f"""

## 19. {date_text} Part 4 TYMP–5′-MTA–SPP1机制空间图谱增强精选版

- 新增独立精选版报告：`{REPORT_TARGET}`，仅保留与客户湿实验机制方向直接相容的空间结果，不覆盖现有完整Part 4主报告。
- Figure 1为8-panel组织空间图，汇总GLMF1/GLMF2的H&E底图、TYMP-high、Polyamine–MTA轮廓、SPP1+ TAM丰度和联合高状态生态位；Figure 2为6-panel患者图谱，包含代表患者4张空间状态图、12例完整效应矩阵和12例共激活空间图谱。
- 关键证据层级保持：GLMF1/GLMF2提供样本内空间分块支持；TJ-RCC 12例提供患者方向重复趋势。湿实验为因果锚点，空间数据不单独证明5′-MTA实测蓄积、统计中介、直接通讯或恶性演进因果。
- DOCX独立核验58/58通过；LibreOffice渲染为A4共7页，逐页人工检查7/7通过。表格边框、两张内嵌主图哈希、14个panel来源、12例效应矩阵及正式数值均已核验。
- 报告状态为`candidate`，未发现服务器外发或客户签收证据。

## 20. 最近更新时间

{timestamp_text}（服务器实际时间）
"""
append_text(PROJECT_STATUS, status_block)

handoff_block = f"""

## {date_text} Part 4 TYMP–5′-MTA–SPP1机制空间图谱增强精选版

- 精选版DOCX：`{REPORT_TARGET}`；SHA-256=`{sha256(REPORT_TARGET)}`。
- 正式主图：8-panel组织空间生态位图`{FIGURE1_PNG}`与6-panel患者空间图谱`{FIGURE2_PNG}`；完整PNG/PDF、14个panel来源、12例效应矩阵、代表患者选择表、生成脚本、独立核验、PDF预览及人工逐页检查均归档于`{RESULT_TARGET}`。
- 独立核验58/58通过；A4 PDF共7页，人工逐页检查7/7通过。v1和完整Part 4主报告均未覆盖。
- 对客户方向的表述边界：空间分析支持人体ccRCC组织中的机制相容性和患者异质性补充，不替代5′-MTA定量、分子中介和干预因果证据。当前状态`candidate`，未外发。

## 本次交接更新时间

{timestamp_text}（Asia/Shanghai）
"""
append_text(HANDOFF, handoff_block)

manifest_rows = [
    ["p4_tymp_curated_report_v2_20260803", "Part4", "GLMF1;GLMF2;12 TJ-RCC tumor cases", "TYMP-5MTA-SPP1 spatial-atlas curated validation", "client_report", str(REPORT_TARGET), str(RESULT_TARGET / "build_tymp_curated_report_v2.py"), "formal TYMP immune ecology tables;12-patient effect matrix;8-panel and 6-panel figures", "batch_bench;Pandoc;LibreOffice7.3", "positive mechanism-aligned spatial evidence only;causal boundary retained", "v2", date_text, str(REPORT_TARGET.stat().st_size), "passed", str(RESULT_TARGET / "curated_report_independent_verification_v2.tsv") + ";" + str(RESULT_TARGET / "curated_report_manual_page_review_v2.tsv"), "candidate", "", f"SHA-256={sha256(REPORT_TARGET)};58/58;7 pages;not delivered"],
    ["p4_tymp_curated_figure1_v2_20260803", "Part4", "GLMF1;GLMF2", "TYMP-MTA-context-SPP1 TAM spatial niche atlas", "main_figure", str(FIGURE1_PNG), str(RESULT_TARGET / "build_tymp_curated_figures_v2.py"), "spot scores;H5AD H&E;formal spatial-block tests", "batch_bench", "8 panels A-H;spatial contours and joint niche overlays", "v2", date_text, str(FIGURE1_PNG.stat().st_size), "passed", str(RESULT_TARGET / "panel_source_manifest_v2.tsv"), "not_used", "", f"SHA-256={sha256(FIGURE1_PNG)};8 panels;PNG and PDF"],
    ["p4_tymp_curated_figure2_v2_20260803", "Part4", "12 TJ-RCC tumor cases", "patient-resolved TYMP-SPP1 spatial atlas", "main_figure", str(FIGURE2_PNG), str(RESULT_TARGET / "build_tymp_curated_figures_v2.py"), "patient spot-state scores;patient effect matrix;formal patient tables", "batch_bench", "6 panels A-F;all 12 patients included;patient is biological replicate", "v2", date_text, str(FIGURE2_PNG.stat().st_size), "passed", str(RESULT_TARGET / "panel_source_manifest_v2.tsv"), "not_used", "", f"SHA-256={sha256(FIGURE2_PNG)};6 panels;PNG and PDF"],
    ["p4_tymp_curated_report_verification_v2_atlas_20260803", "Part4", "curated mechanism spatial-atlas report", "DOCX and source independent verification", "qc_table", str(RESULT_TARGET / "curated_report_independent_verification_v2.tsv"), str(RESULT_TARGET / "verify_tymp_curated_report_v2.py"), "curated DOCX;formal source tables;figure metadata;12-patient matrix", "batch_bench", "ZIP;media hashes;tables;A4 layout;source recomputation", "v2", date_text, str((RESULT_TARGET / "curated_report_independent_verification_v2.tsv").stat().st_size), "passed", "58/58", "not_used", "", f"SHA-256={sha256(RESULT_TARGET / 'curated_report_independent_verification_v2.tsv')}"],
    ["p4_tymp_curated_report_manual_v2_atlas_20260803", "Part4", "all 7 report pages", "spatial-atlas curated report visual review", "manual_qc", str(RESULT_TARGET / "curated_report_manual_page_review_v2.tsv"), "Codex visual inspection", str(REPORT_PDF), "manual visual review", "tables;captions;spatial figures;labels;clipping;page flow", "v2", date_text, str((RESULT_TARGET / "curated_report_manual_page_review_v2.tsv").stat().st_size), "passed", "7/7", "not_used", "", f"SHA-256={sha256(RESULT_TARGET / 'curated_report_manual_page_review_v2.tsv')}"],
]
with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerows(manifest_rows)

installed = [REPORT_TARGET, REPORT_PDF, FIGURE1_PNG, FIGURE1_PDF, FIGURE2_PNG, FIGURE2_PDF]
installation_metadata = {
    "status": "installed",
    "installed_at": timestamp_text,
    "report_status": "candidate",
    "artifacts": [
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)}
        for path in installed
    ],
    "project_records_after_install": {
        str(path): {"size": path.stat().st_size, "sha256": sha256(path)}
        for path in [PROJECT_STATUS, HANDOFF, MANIFEST]
    },
    "project_records_before_install": str(RESULT_TARGET / "project_records_before_install"),
    "independent_verification": "58/58 passed",
    "manual_render_review": "7/7 pages passed",
}
(RESULT_TARGET / "installation_metadata_v2.json").write_text(json.dumps(installation_metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(installation_metadata, ensure_ascii=False, indent=2))
