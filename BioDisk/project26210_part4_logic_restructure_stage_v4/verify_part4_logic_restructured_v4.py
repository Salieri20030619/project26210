#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta
from pathlib import Path
from zipfile import ZipFile


STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v4")
DOCX = STAGE / "project26210_第四部分空间转录组分析报告_主线逻辑重构版_v4_candidate.docx"
PDF = STAGE / "rendered_v4/part4_logic_v4.pdf"
FLOWCHART = STAGE / "part4_mainline_framework_v4.png"
SUMMARY = STAGE / "independent_verification_summary_v4.json"
SOURCE = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3/"
    "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_"
    "POSTN_APOE论文对标增强版_v3_candidate.docx"
)
EXPECTED_SOURCE_SHA = "d9ed4172e3a9b675946cdf0fa0b26c4d44154c33e40c72048001c2dee2b19940"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
R = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check(condition: bool, label: str, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(label)
    checks.append(label)


def main() -> None:
    checks: list[str] = []
    for path in (SOURCE, DOCX, PDF, FLOWCHART):
        check(path.is_file() and path.stat().st_size > 0, f"file exists: {path.name}", checks)
    check(sha256(SOURCE) == EXPECTED_SOURCE_SHA, "source fingerprint unchanged", checks)

    with ZipFile(SOURCE) as source_zip:
        source_media = {name for name in source_zip.namelist() if name.startswith("word/media/")}
        source_root = ET.fromstring(source_zip.read("word/document.xml"))
    with ZipFile(DOCX) as archive:
        check(archive.testzip() is None, "DOCX ZIP CRC", checks)
        names = set(archive.namelist())
        media = {name for name in names if name.startswith("word/media/")}
        root = ET.fromstring(archive.read("word/document.xml"))
        rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))

    check(source_media.issubset(media), "all source media retained", checks)
    check(len(source_media) == 57 and len(media) == 58, "57 source media + 1 flowchart", checks)
    check(len(list(source_root.iter(W + "drawing"))) == 57, "source has 57 drawing objects", checks)
    check(len(list(root.iter(W + "drawing"))) == 58, "output has 58 drawing objects", checks)
    check(len(list(root.iter(W + "tbl"))) == 29, "output has 29 tables", checks)
    check(len(list(root.iter(W + "p"))) == 1279, "output paragraph count", checks)

    image_targets = []
    for rel in rels.findall(R + "Relationship"):
        if rel.get("Type", "").endswith("/image"):
            image_targets.append("word/" + rel.get("Target", ""))
    check(len(image_targets) == 58, "58 image relationships", checks)
    check(all(target in media for target in image_targets), "all image relationships resolve", checks)
    check(len(set(image_targets)) == 58, "image targets are unique", checks)

    text = "\n".join("".join(t.text or "" for t in p.iter(W + "t")) for p in root.iter(W + "p"))
    required_headings = [
        "执行摘要",
        "1. 研究问题与证据策略",
        "2. 数据基础与空间图谱",
        "3. GLMF1/GLMF2发现层：候选关系与方法深化",
        "4. TJ-RCC 12名患者验证层：从空间重复到严格收窄",
        "5. 机制收敛层：TYMP–MTA–SPP1、血管生态位与平行对标",
        "6. 证据分级、主要结论与下一步",
    ]
    for heading in required_headings:
        check(text.count(heading) >= 1, f"heading present: {heading}", checks)
    for forbidden in ("[FLOWCHART_IMAGE_PLACEHOLDER]", "图FLOW", "4.12 【高级分析补丁】", "4.13 【TJ-RCC", "4.14 【TYMP", "4.15 【血管"):
        check(forbidden not in text, f"forbidden token absent: {forbidden}", checks)

    captions = []
    for line in text.splitlines():
        match = re.match(r"图\s*4\.(\d+)([AB]?)", line.strip())
        if match:
            captions.append(match.group(1) + match.group(2))
    expected_captions = ["1", "2", "3", "4A", "4B", *[str(i) for i in range(5, 55)]]
    check(captions == expected_captions, "55 figure captions are complete and ordered", checks)
    check("2026年8月5日（主线逻辑重构版）" in text, "cover date updated", checks)
    check("空间定位 → GLMF1/GLMF2发现 → TJ-RCC患者验证 → 机制收敛 → 证据分级" in text,
          "cover framework updated", checks)

    with tempfile.TemporaryDirectory(prefix="part4_v4_verify_") as temp_dir:
        text_path = Path(temp_dir) / "report.txt"
        subprocess.run(["pdftotext", "-layout", str(PDF), str(text_path)], check=True)
        pdf_pages = text_path.read_text(errors="replace").split("\f")[:-1]
        pdfinfo = subprocess.run(["pdfinfo", str(PDF)], check=True, text=True, capture_output=True).stdout
    check("Pages:           73" in pdfinfo, "PDF has 73 pages", checks)
    check(len(pdf_pages) == 73, "PDF text has 73 pages", checks)
    check(all(len("".join(page.split())) >= 30 for page in pdf_pages), "PDF has no blank/near-blank pages", checks)
    normalized_first_page = "".join(pdf_pages[0].split())
    check("空间定位→GLMF1/GLMF2发现→TJ-RCC患者验证→机制收敛→证据分级" in normalized_first_page,
          "rendered cover framework visible", checks)
    check("2026年8月5日（主线逻辑重构版）" in normalized_first_page,
          "rendered cover date visible", checks)
    check("图4.1第四部分空间分析主线" in "".join(pdf_pages[1].split()),
          "flowchart caption visible in rendered PDF", checks)
    for page_no, heading in [(4, "1.研究问题与证据策略"), (7, "2.数据基础与空间图谱"),
                             (12, "3.GLMF1/GLMF2发现层"), (28, "4.TJ-RCC12名患者验证层"),
                             (50, "5.机制收敛层"), (72, "6.证据分级、主要结论与下一步")]:
        check(heading in "".join(pdf_pages[page_no - 1].split()),
              f"rendered chapter starts on page {page_no}", checks)

    summary = {
        "status": "passed",
        "verified_at": datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds"),
        "checks_passed": len(checks),
        "checks_failed": 0,
        "docx": str(DOCX),
        "docx_size": DOCX.stat().st_size,
        "docx_sha256": sha256(DOCX),
        "pdf": str(PDF),
        "pdf_size": PDF.stat().st_size,
        "pdf_sha256": sha256(PDF),
        "flowchart_sha256": sha256(FLOWCHART),
        "pages": 73,
        "tables": 29,
        "drawings": 58,
        "captioned_figures": 55,
        "checks": checks,
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
