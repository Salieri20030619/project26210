#!/usr/bin/env python3
"""Install the two verified Part 4 v3 reports without overwriting v2."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3")
DOCS = PROJECT / "docs/client_reports/part4_spatial"
RESULTS = PROJECT / "results/part4_spatial"
QC_DEST = RESULTS / "tymp_vascular_reports_v3_qc"

FULL_V2 = DOCS / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
CURATED_V2 = DOCS / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
FULL_SRC = STAGE / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3_candidate.docx"
CURATED_SRC = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3_candidate.docx"
FULL_DEST = DOCS / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3.docx"
CURATED_DEST = DOCS / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3.docx"
FULL_PDF_SRC = STAGE / "rendered_full_v3/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3_candidate.pdf"
CURATED_PDF_SRC = STAGE / "rendered_curated_v3/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3_candidate.pdf"
FULL_PDF_DEST = QC_DEST / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3.pdf"
CURATED_PDF_DEST = QC_DEST / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3.pdf"

RECORDS = {
    PROJECT / "CODEX_HANDOFF.md": "094f3c7ae41a45cb1963fd89b7871ed132dc3e31d025b0be51cc8695cf6fcca9",
    PROJECT / "PROJECT_STATUS.md": "7eb359aa37212a3d52aba0ac7396c079d1a33da1b58e338c7ac01afd74cad4fc",
    PROJECT / "RESULTS_MANIFEST.tsv": "e7b257b6fedd33e4e9e7f3f928bdca40c3e7abcf676be4b1bec53b9a5a8ddd76",
}
EXPECTED = {
    FULL_V2: "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b",
    CURATED_V2: "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39",
    FULL_SRC: "d9ed4172e3a9b675946cdf0fa0b26c4d44154c33e40c72048001c2dee2b19940",
    CURATED_SRC: "919b0af905bd86f5742649699b428c24d612b6fa8ae149f9f75f2bb6791cde5d",
    FULL_PDF_SRC: "52b484ead1ccc10e635bbfe5e374e6510a37624deefd8fb879aa4a713db2e22c",
    CURATED_PDF_SRC: "f17221f2ff0ff9c6a6b5f16f05bc0508bfb53602414cf5ccdfe3823be3532e37",
    STAGE / "independent_verification_v3.tsv": "69325b8242f486de57aeaca241f99adf5b8907106fb52ba2a656b17a05cc9a43",
    STAGE / "render_verification_v3.tsv": "a90e99dd7cdfd9ac2890fc8059104384b08dc7f31dd601b4b1dae9494cd84c58",
    STAGE / "report_render_manual_review_v3.tsv": "952f8c2ba3647d0529d62fb5a8779b4daec8a44ab03b0dad4240aae6cbab4459",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_hashes(items: dict[Path, str]) -> None:
    for path, expected in items.items():
        observed = sha256(path)
        if observed != expected:
            raise RuntimeError(f"hash changed: {path}\nobserved={observed}\nexpected={expected}")


def atomic_write(path: Path, content: str) -> None:
    mode = path.stat().st_mode
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.chmod(temporary, mode)
    os.replace(temporary, path)


def manifest_line(values: list[str]) -> str:
    if len(values) != 18 or any("\t" in value or "\n" in value for value in values):
        raise ValueError("invalid manifest row")
    return "\t".join(values)


subprocess.run([str(PROJECT / "check_memory_gate.sh")], check=True)
assert_hashes(EXPECTED)
assert_hashes(RECORDS)
for destination in (FULL_DEST, CURATED_DEST, QC_DEST):
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite: {destination}")

independent = json.loads((STAGE / "independent_verification_summary_v3.json").read_text(encoding="utf-8"))
rendered = json.loads((STAGE / "render_verification_summary_v3.json").read_text(encoding="utf-8"))
if independent["status"] != "passed" or independent["checks_passed"] != 50 or independent["checks_failed"] != 0:
    raise RuntimeError("independent DOCX verification is not 50/50 PASS")
if rendered["status"] != "passed" or rendered["checks_passed"] != 29 or rendered["checks_failed"] != 0:
    raise RuntimeError("render verification is not 29/29 PASS")

with (STAGE / "report_render_manual_review_v3.tsv").open(encoding="utf-8", newline="") as handle:
    manual_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(manual_rows) != 11 or any(row["status"] != "PASS" for row in manual_rows):
    raise RuntimeError("manual review is not 11/11 PASS")

manifest_path = PROJECT / "RESULTS_MANIFEST.tsv"
manifest = manifest_path.read_text(encoding="utf-8")
reader = csv.DictReader(manifest.splitlines(), delimiter="\t")
if reader.fieldnames is None or len(reader.fieldnames) != 18:
    raise RuntimeError("manifest schema is not 18 columns")
existing_ids = {row["result_id"] for row in reader}
new_ids = {
    "p4_tymp_vascular_full_report_v3_20260804",
    "p4_tymp_vascular_full_report_pdf_v3_20260804",
    "p4_tymp_vascular_curated_report_v3_20260804",
    "p4_tymp_vascular_curated_report_pdf_v3_20260804",
    "p4_tymp_vascular_report_independent_qc_v3_20260804",
    "p4_tymp_vascular_report_render_qc_v3_20260804",
    "p4_tymp_vascular_report_manual_qc_v3_20260804",
}
if existing_ids & new_ids:
    raise RuntimeError(f"manifest IDs already exist: {sorted(existing_ids & new_ids)}")

status_path = PROJECT / "PROJECT_STATUS.md"
status = status_path.read_text(encoding="utf-8")
if "Part 4完整报告与TYMP精选报告v3已生成" in status:
    raise RuntimeError("status block already exists")
status_match = re.search(r"^## 0\. 当前权威状态（[^\n]+）\n", status, flags=re.MULTILINE)
if not status_match:
    raise RuntimeError("authoritative status heading not found")

handoff_path = PROJECT / "CODEX_HANDOFF.md"
handoff = handoff_path.read_text(encoding="utf-8")
if "Part 4完整与精选报告v3正式候选" in handoff:
    raise RuntimeError("handoff block already exists")

now = datetime.now().astimezone()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S %z")
iso_timestamp = now.isoformat(timespec="seconds")
status_block = f"""

- **Part 4完整报告与TYMP精选报告v3已生成、渲染并通过验收，原v2未覆盖。** 完整版为`docs/client_reports/part4_spatial/{FULL_DEST.name}`（{FULL_SRC.stat().st_size:,} bytes；SHA-256=`{sha256(FULL_SRC)}`；76页），精选版为`docs/client_reports/part4_spatial/{CURATED_DEST.name}`（{CURATED_SRC.stat().st_size:,} bytes；SHA-256=`{sha256(CURATED_SRC)}`；11页）。两份报告均为`candidate`，尚无客户外发或签收证据。
- **新增内容已纳入两种报告深度。** 完整版新增4.15节、图4.50–4.53、9条5×5患者级空间关系、两尺度关键重复、7类DPT评分及第五章同步结论；精选版保留原两张组合图并新增3.2节、图3–4，只突出两条全局FDR支持的血管关系和直接限制机制解释的髓系状态结果。空间相关不作因果/耐药结论，DPT不作RNA velocity、谱系转化或纵向进展解释。
- **报告QC全部通过。** DOCX结构、旧媒体保真、新媒体源哈希、表格边框、禁用项目符号和源数值复算共50/50项通过；LibreOffice 7.3实际渲染完整版76页、精选版11页，PDF级29/29项通过，11/11页组人工目检通过。正式审计目录为`results/part4_spatial/tymp_vascular_reports_v3_qc`；最初一次受沙箱dconf只读限制的渲染失败属于环境性非科学失败，随后使用独立临时LibreOffice用户配置成功完成渲染。
"""
status = status[: status_match.end()] + status_block + status[status_match.end() :]
status = re.sub(r"^## 0\. 当前权威状态（[^\n]+）$", f"## 0. 当前权威状态（{now.strftime('%Y-%m-%d %H:%M')} +0800）", status, count=1, flags=re.MULTILINE)

handoff_block = f"""

## 2026-08-04 Part 4完整与精选报告v3正式候选

- 用户明确要求同时修改完整报告和精选报告。本轮以只读v2为基础新建v3，未覆盖原文件。完整v3正式路径为`{FULL_DEST}`；精选v3正式路径为`{CURATED_DEST}`。
- 完整版新增4.15节与图4.50–4.53，完整列出9条5×5患者级空间关系和7类DPT评分；精选版新增3.2节与图3–4，仅突出两条两尺度均12/12且全局BH FDR=0.002197的血管关系，以及C1QC-state、SPP1-state和二者差值的横断面DPT结果。
- 完整版DOCX为{FULL_SRC.stat().st_size:,} bytes，SHA-256=`{sha256(FULL_SRC)}`，LibreOffice渲染76页；精选版DOCX为{CURATED_SRC.stat().st_size:,} bytes，SHA-256=`{sha256(CURATED_SRC)}`，渲染11页。DOCX独立检查50/50、PDF级检查29/29、人工页组检查11/11通过。
- QC目录：`{QC_DEST}`，含两份渲染PDF、生成/独立核验/渲染核验/安装脚本、元数据、逐页PNG、联系表和人工复核记录。v2源文件哈希仍分别为`{sha256(FULL_V2)}`与`{sha256(CURATED_V2)}`。
- 两份v3当前状态均为`candidate`、未外发。空间结果不能证明直接通讯、免疫治疗耐药或治疗反应；DPT只表示依赖既定单核细胞根和邻居图的横断面状态排序，不能证明单核细胞向SPP1+ TAM的谱系转化。

## 本次交接更新时间

{iso_timestamp}（Asia/Shanghai）
"""
handoff = handoff.rstrip() + handoff_block + "\n"

qc_relative = "results/part4_spatial/tymp_vascular_reports_v3_qc"
rows = [
    ["p4_tymp_vascular_full_report_v3_20260804", "Part4", "GLMF1/GLMF2 and TJ-RCC 12 patients; 53794 myeloid cells", "complete Part 4 report with vascular niche and cross-sectional DPT", "client_report", str(FULL_DEST), str(QC_DEST / "build_part4_vascular_reports_v3.py"), "full v2 report;formal vascular and myeloid synthesis outputs", "batch_bench;Python3.11.15;python-docx", "versioned v3;section 4.15;figures 4.50-4.53;causal and lineage limits", "v3", "2026-08-04", str(FULL_SRC.stat().st_size), "passed", "50/50 DOCX;29/29 PDF;11/11 manual page groups", "candidate", "", f"SHA-256={sha256(FULL_SRC)};76 rendered pages;v2 preserved"],
    ["p4_tymp_vascular_full_report_pdf_v3_20260804", "Part4", "76 pages", "rendered complete Part 4 v3 report", "report_pdf", str(FULL_PDF_DEST), "LibreOffice 7.3.7.2 headless", str(FULL_DEST), "LibreOffice 7.3.7.2", "PDF rendering and page-level text/layout verification", "v3", "2026-08-04", str(FULL_PDF_SRC.stat().st_size), "passed", "29/29 PDF checks;manual visual PASS", "not_used", "", f"SHA-256={sha256(FULL_PDF_SRC)};Letter;preview only"],
    ["p4_tymp_vascular_curated_report_v3_20260804", "Part4", "TJ-RCC 12 patients; 63 eligible myeloid samples", "curated TYMP-5MTA-SPP1 spatial report with vascular and DPT summaries", "client_report", str(CURATED_DEST), str(QC_DEST / "build_part4_vascular_reports_v3.py"), "curated v2 report;formal vascular and myeloid synthesis outputs", "batch_bench;Python3.11.15;python-docx", "versioned v3;section 3.2;figures 3-4;selected supported results", "v3", "2026-08-04", str(CURATED_SRC.stat().st_size), "passed", "50/50 DOCX;29/29 PDF;11/11 manual page groups", "candidate", "", f"SHA-256={sha256(CURATED_SRC)};11 rendered pages;v2 preserved"],
    ["p4_tymp_vascular_curated_report_pdf_v3_20260804", "Part4", "11 pages", "rendered curated TYMP-5MTA-SPP1 v3 report", "report_pdf", str(CURATED_PDF_DEST), "LibreOffice 7.3.7.2 headless", str(CURATED_DEST), "LibreOffice 7.3.7.2", "PDF rendering and page-level text/layout verification", "v3", "2026-08-04", str(CURATED_PDF_SRC.stat().st_size), "passed", "29/29 PDF checks;manual visual PASS", "not_used", "", f"SHA-256={sha256(CURATED_PDF_SRC)};A4;preview only"],
    ["p4_tymp_vascular_report_independent_qc_v3_20260804", "Part4", "2 DOCX reports", "independent DOCX structure media and numeric verification", "qc_table", f"{QC_DEST}/independent_verification_v3.tsv", f"{QC_DEST}/verify_part4_vascular_reports_v3.py", "v2 reports;v3 candidates;formal source tables and figures", "batch_bench;Python3.11.15", "OOXML;CRC;relationships;media hashes;borders;markers;numeric replay", "v3", "2026-08-04", str((STAGE / "independent_verification_v3.tsv").stat().st_size), "passed", "50/50 checks", "not_used", "", f"SHA-256={sha256(STAGE / 'independent_verification_v3.tsv')};old media byte-identical"],
    ["p4_tymp_vascular_report_render_qc_v3_20260804", "Part4", "76-page full and 11-page curated PDFs", "rendered PDF verification", "qc_table", f"{QC_DEST}/render_verification_v3.tsv", f"{QC_DEST}/verify_rendered_reports_v3.py", "rendered full and curated v3 PDFs", "batch_bench;Python3.11.15;poppler", "page count;page size;text markers;manual review linkage", "v3", "2026-08-04", str((STAGE / "render_verification_v3.tsv").stat().st_size), "passed", "29/29 checks", "not_used", "", f"SHA-256={sha256(STAGE / 'render_verification_v3.tsv')};no clipping or missing pages"],
    ["p4_tymp_vascular_report_manual_qc_v3_20260804", "Part4", "11 reviewed page groups", "manual rendered report visual review", "manual_qc", f"{QC_DEST}/report_render_manual_review_v3.tsv", "Codex visual inspection", "rendered pages and contact sheets", "manual visual review", "original-resolution pages;layout;clipping;tables;figures;semantic limits", "v3", "2026-08-04", str((STAGE / "report_render_manual_review_v3.tsv").stat().st_size), "passed", "11/11 page groups PASS", "not_used", "", f"SHA-256={sha256(STAGE / 'report_render_manual_review_v3.tsv')};full and curated reports reviewed"],
]
manifest = manifest.rstrip("\n") + "\n" + "\n".join(manifest_line(row) for row in rows) + "\n"

qc_tmp = RESULTS / f".tymp_vascular_reports_v3_qc.installing.{os.getpid()}"
installed_paths: list[Path] = []
try:
    qc_tmp.mkdir(parents=False)
    for filename in (
        "build_part4_vascular_reports_v3.py",
        "verify_part4_vascular_reports_v3.py",
        "verify_rendered_reports_v3.py",
        "install_part4_vascular_reports_v3.py",
        "report_generation_metadata_v3.json",
        "independent_verification_v3.tsv",
        "independent_verification_summary_v3.json",
        "render_verification_v3.tsv",
        "render_verification_summary_v3.json",
        "report_render_manual_review_v3.tsv",
    ):
        shutil.copy2(STAGE / filename, qc_tmp / filename)
    shutil.copytree(STAGE / "contact_sheets_v3", qc_tmp / "contact_sheets_v3")
    shutil.copytree(STAGE / "full_review_pages_v3", qc_tmp / "full_review_pages_v3")
    shutil.copytree(STAGE / "curated_review_pages_v3", qc_tmp / "curated_review_pages_v3")
    shutil.copy2(FULL_PDF_SRC, qc_tmp / FULL_PDF_DEST.name)
    shutil.copy2(CURATED_PDF_SRC, qc_tmp / CURATED_PDF_DEST.name)
    readme = f"""# Part 4 vascular/DPT report v3 QC\n\nInstalled: {timestamp}\n\n- Full v3 DOCX: `{FULL_DEST}`\n- Curated v3 DOCX: `{CURATED_DEST}`\n- DOCX verification: 50/50 PASS\n- PDF verification: 29/29 PASS\n- Manual review: 11/11 page groups PASS\n- Full PDF: 76 pages, SHA-256 `{sha256(FULL_PDF_SRC)}`\n- Curated PDF: 11 pages, SHA-256 `{sha256(CURATED_PDF_SRC)}`\n- v2 sources were verified unchanged and were not overwritten.\n- Report status: candidate; no external delivery evidence.\n"""
    (qc_tmp / "README.md").write_text(readme, encoding="utf-8")
    metadata = {
        "installed_at": iso_timestamp,
        "status": "candidate",
        "full_v2": {"path": str(FULL_V2), "sha256": sha256(FULL_V2)},
        "curated_v2": {"path": str(CURATED_V2), "sha256": sha256(CURATED_V2)},
        "full_v3": {"path": str(FULL_DEST), "sha256": sha256(FULL_SRC), "bytes": FULL_SRC.stat().st_size, "pages": 76},
        "curated_v3": {"path": str(CURATED_DEST), "sha256": sha256(CURATED_SRC), "bytes": CURATED_SRC.stat().st_size, "pages": 11},
        "qc": {"docx": "50/50 PASS", "pdf": "29/29 PASS", "manual": "11/11 page groups PASS"},
        "environment": {"generation": "batch_bench Python 3.11.15", "rendering": "LibreOffice 7.3.7.2 headless", "validation": "batch_bench Python 3.11.15 and poppler"},
    }
    (qc_tmp / "report_install_metadata_v3.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    for source, destination in ((FULL_SRC, FULL_DEST), (CURATED_SRC, CURATED_DEST)):
        temporary = destination.with_name(f".{destination.name}.installing.{os.getpid()}")
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
        installed_paths.append(destination)
    os.replace(qc_tmp, QC_DEST)
    installed_paths.append(QC_DEST)

    atomic_write(status_path, status)
    atomic_write(handoff_path, handoff)
    atomic_write(manifest_path, manifest)
except Exception:
    if qc_tmp.exists():
        shutil.rmtree(qc_tmp)
    for path in reversed(installed_paths):
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    raise

assert sha256(FULL_DEST) == EXPECTED[FULL_SRC]
assert sha256(CURATED_DEST) == EXPECTED[CURATED_SRC]
assert sha256(FULL_PDF_DEST) == EXPECTED[FULL_PDF_SRC]
assert sha256(CURATED_PDF_DEST) == EXPECTED[CURATED_PDF_SRC]
print(f"PART4_VASCULAR_REPORTS_V3_INSTALL_PASS\nfull={FULL_DEST}\ncurated={CURATED_DEST}\nqc={QC_DEST}")
