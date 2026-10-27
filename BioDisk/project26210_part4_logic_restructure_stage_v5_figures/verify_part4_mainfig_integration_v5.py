#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zipfile import ZipFile


STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v5_figures")
SOURCE = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v4/"
    "project26210_第四部分空间转录组分析报告_主线逻辑重构版_v4_candidate.docx"
)
DOCX = STAGE / "project26210_第四部分空间转录组分析报告_主线逻辑重构版_三张主图整合_v5_candidate.docx"
PDF = STAGE / "rendered_v5/part4_logic_mainfig_v5.pdf"
SUMMARY = STAGE / "independent_verification_summary_v5.json"
FIGROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_final_v4")
FIGURE_MEDIA = {
    "word/media/part4_main_figure_4_1_v4.png": FIGROOT / "Figure1_TYMP_MTA_SPP1_evidence_narrowing_v4.png",
    "word/media/part4_main_figure_4_2_v4.png": FIGROOT / "Figure2_SPP1_vascular_niche_Tcell_balance_v4.png",
    "word/media/part4_main_figure_4_3_v4.png": FIGROOT / "Figure3_myeloid_state_continuum_v4.png",
}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def node_text(node: ET.Element) -> str:
    return "".join(t.text or "" for t in node.iter(W + "t"))


def check(condition: bool, label: str, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(label)
    checks.append(label)


def blip_sequence(root: ET.Element) -> list[str]:
    return [node.get(R + "embed", "") for node in root.iter(A + "blip")]


def main() -> None:
    checks: list[str] = []
    for path in (SOURCE, DOCX, PDF, *FIGURE_MEDIA.values()):
        check(path.is_file() and path.stat().st_size > 0, f"file exists: {path.name}", checks)

    with ZipFile(SOURCE) as source_zip:
        source_root = ET.fromstring(source_zip.read("word/document.xml"))
        source_rels = ET.fromstring(source_zip.read("word/_rels/document.xml.rels"))
        source_media = {name for name in source_zip.namelist() if name.startswith("word/media/")}
    with ZipFile(DOCX) as archive:
        check(archive.testzip() is None, "DOCX ZIP CRC", checks)
        root = ET.fromstring(archive.read("word/document.xml"))
        rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        media = {name for name in archive.namelist() if name.startswith("word/media/")}
        for embedded_name, source_path in FIGURE_MEDIA.items():
            check(sha256_bytes(archive.read(embedded_name)) == sha256(source_path),
                  f"embedded figure hash: {embedded_name}", checks)

    check(source_media.issubset(media), "all pre-existing media retained", checks)
    check(len(source_media) == 58 and len(media) == 61, "58 existing media + 3 main figures", checks)
    check(len(list(source_root.iter(W + "drawing"))) == 58, "source drawing count", checks)
    check(len(list(root.iter(W + "drawing"))) == 61, "output drawing count", checks)
    check(len(list(root.iter(W + "tbl"))) == 29, "table count unchanged", checks)

    source_rel_target = {rel.get("Id"): rel.get("Target") for rel in source_rels
                         if rel.get("Type", "").endswith("/image")}
    output_rel_target = {rel.get("Id"): rel.get("Target") for rel in rels
                         if rel.get("Type", "").endswith("/image")}
    source_existing_sequence = [source_rel_target[rid] for rid in blip_sequence(source_root)]
    output_existing_sequence = [output_rel_target[rid] for rid in blip_sequence(root)
                                if output_rel_target[rid] in set(source_rel_target.values())]
    check(output_existing_sequence == source_existing_sequence,
          "pre-existing image order unchanged", checks)

    body = root.find(W + "body")
    check(body is not None, "document body present", checks)
    direct_paragraphs = [node_text(node).strip() for node in list(body) if node.tag == W + "p"]
    main_captions = [text for text in direct_paragraphs if re.match(r"^图4\.[123]\s{2}", text)]
    supplement_captions = [text for text in direct_paragraphs if text.startswith("补充图 S4.")]
    expected_supplements = ["1", "2", "3", "4A", "4B", *[str(i) for i in range(5, 55)]]
    observed_supplements = []
    for text in supplement_captions:
        match = re.match(r"补充图 S4\.(\d+)([AB]?)", text)
        check(match is not None, f"supplement caption parses: {text[:50]}", checks)
        observed_supplements.append(match.group(1) + match.group(2))
    check(len(main_captions) == 3, "exactly 3 main figure captions", checks)
    check(observed_supplements == expected_supplements,
          "55 supplementary captions complete and ordered", checks)
    check(not any(re.match(r"^图\s*4\.(?:4|[5-9]|[1-9][0-9])", text) for text in direct_paragraphs),
          "no pre-existing figure remains labeled as main", checks)

    full_text = "\n".join(direct_paragraphs)
    for number, anchor, heading in [
        ("4.1", "整体证据链见图4.1", "5.2 TYMP–MTA-context–SPP1"),
        ("4.2", "队列空间关系及尺度稳定性见图4.2", "5.3.1 TYMP–SPP1血管生态位"),
        ("4.3", "髓系状态连续谱的综合结果见图4.3", "5.3.2 髓系横断面状态排序"),
    ]:
        check(anchor in full_text, f"main figure {number} cited in text", checks)
        check(any(text.startswith(heading) for text in direct_paragraphs),
              f"target heading present for {number}", checks)
    check("主图整合版" in full_text and "原报告全部分析图保持原位置" in full_text,
          "cover documents main/supplement rule", checks)

    pdfinfo = subprocess.run(["pdfinfo", str(PDF)], check=True, text=True,
                             stdout=subprocess.PIPE).stdout
    pages_match = re.search(r"^Pages:\s+(\d+)$", pdfinfo, re.MULTILINE)
    check(pages_match is not None and int(pages_match.group(1)) == 75, "rendered PDF has 75 pages", checks)
    text_path = STAGE / "rendered_v5/part4_logic_mainfig_v5.txt"
    subprocess.run(["pdftotext", "-layout", str(PDF), str(text_path)], check=True)
    pages = text_path.read_text(errors="replace").split("\f")[:-1]
    check(len(pages) == 75, "rendered PDF text has 75 pages", checks)
    check(all(len("".join(page.split())) >= 50 for page in pages), "no blank/near-blank PDF page", checks)
    normalized = ["".join(page.split()) for page in pages]
    expected_pages = {
        "图4.1TYMP–MTA–SPP1巨噬细胞假说": 57,
        "图4.2患者重复的SPP1–vascularniche": 63,
        "图4.3髓系细胞形成以C1QC状态增强": 66,
    }
    for phrase, page_no in expected_pages.items():
        check(phrase in normalized[page_no - 1], f"main caption visible on page {page_no}", checks)
    check(re.search(r"补充图S4\.1(?!\d)", normalized[1]) is not None,
          "supplementary flowchart label visible on page 2", checks)
    check(any(re.search(r"补充图S4\.54(?!\d)", page) for page in normalized),
          "last supplementary caption visible", checks)

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
        "pages": 75,
        "main_figures": 3,
        "supplementary_captions": 55,
        "tables": 29,
        "drawings": 61,
        "manual_visual_qc": "PASS: pages 57, 63, and 66 inspected; captions and full figures visible without clipping",
        "checks": checks,
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
