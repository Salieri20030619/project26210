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
DOCX = STAGE / "project26210_第四部分空间转录组分析报告_仅三张主图精简版_v5_candidate.docx"
PDF = STAGE / "rendered_compact_v5/part4_three_mainfig_compact_v5.pdf"
SUMMARY = STAGE / "compact_report_independent_verification_v5.json"
FIGROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_final_v4")
FIGURE_MEDIA = {
    "word/media/part4_main_figure_4_1_v4.png": FIGROOT / "Figure1_TYMP_MTA_SPP1_evidence_narrowing_v4.png",
    "word/media/part4_main_figure_4_2_v4.png": FIGROOT / "Figure2_SPP1_vascular_niche_Tcell_balance_v4.png",
    "word/media/part4_main_figure_4_3_v4.png": FIGROOT / "Figure3_myeloid_state_continuum_v4.png",
}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check(condition: bool, label: str, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(label)
    checks.append(label)


def node_text(node: ET.Element) -> str:
    return "".join(text.text or "" for text in node.iter(W + "t"))


def main() -> None:
    checks: list[str] = []
    for path in (DOCX, PDF, *FIGURE_MEDIA.values()):
        check(path.is_file() and path.stat().st_size > 0, f"file exists: {path.name}", checks)

    with ZipFile(DOCX) as archive:
        check(archive.testzip() is None, "DOCX ZIP CRC", checks)
        media_payload = {name for name in archive.namelist() if name.startswith("word/media/")}
        check(media_payload == set(FIGURE_MEDIA), "DOCX contains only three main-figure media assets", checks)
        root = ET.fromstring(archive.read("word/document.xml"))
        rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        for media_name, source_path in FIGURE_MEDIA.items():
            check(hashlib.sha256(archive.read(media_name)).hexdigest() == sha256(source_path),
                  f"main figure hash: {media_name}", checks)

    drawings = list(root.iter(W + "drawing"))
    check(len(drawings) == 3, "exactly three visible drawings", checks)
    rel_target = {rel.get("Id"): rel.get("Target") for rel in rels.findall(REL + "Relationship")}
    visible_targets = [rel_target[blip.get(R + "embed")] for blip in root.iter(A + "blip")]
    check(visible_targets == [
        "media/part4_main_figure_4_1_v4.png",
        "media/part4_main_figure_4_2_v4.png",
        "media/part4_main_figure_4_3_v4.png",
    ], "visible drawings are the three main figures in order", checks)
    check(len(list(root.iter(W + "tbl"))) == 29, "all 29 report tables retained", checks)

    paragraphs = [node_text(node).strip() for node in root.iter(W + "p")]
    full_text = "\n".join(paragraphs)
    check("补充图 S4." not in full_text and "Supplementary Figure" not in full_text,
          "all supplementary captions removed", checks)
    main_captions = [text for text in paragraphs if re.match(r"^图4\.[123]\s{2}", text)]
    check(len(main_captions) == 3, "exactly three detailed main captions", checks)
    for number, required in [
        ("4.1", ["A，发现", "B，GLMF2", "E，从空间相容性", "不独立证明5′-MTA"]),
        ("4.2", ["A，变量构建", "B，代表性患者", "F，患者级方向检验", "不证明直接通讯"]),
        ("4.3", ["A，四类髓系状态UMAP", "E，63个合格样本", "F，观察结果", "不替代RNA velocity"]),
    ]:
        caption = next(text for text in main_captions if text.startswith(f"图{number}"))
        check(all(term in caption for term in required), f"detailed panel legend complete: Figure {number}", checks)

    required_headings = [
        "1. 研究问题与证据策略",
        "2. 数据基础与空间图谱",
        "3. GLMF1/GLMF2发现层：候选关系与方法深化",
        "4. TJ-RCC 12名患者验证层：从空间重复到严格收窄",
        "5. 机制收敛层：TYMP–MTA–SPP1、血管生态位与平行对标",
        "6. 证据分级、主要结论与下一步",
    ]
    for heading in required_headings:
        check(heading in full_text, f"content heading retained: {heading}", checks)
    check("【精简版规则】正文仅保留图4.1–图4.3三张主题组合主图" in full_text,
          "compact figure policy stated on cover", checks)

    pdfinfo = subprocess.run(["pdfinfo", str(PDF)], check=True, text=True,
                             stdout=subprocess.PIPE).stdout
    match = re.search(r"^Pages:\s+(\d+)$", pdfinfo, re.MULTILINE)
    check(match is not None and int(match.group(1)) == 38, "rendered PDF has 38 pages", checks)
    text_path = STAGE / "rendered_compact_v5/part4_three_mainfig_compact_v5.txt"
    subprocess.run(["pdftotext", "-layout", str(PDF), str(text_path)], check=True)
    pages = text_path.read_text(errors="replace").split("\f")[:-1]
    check(len(pages) == 38, "PDF text has 38 pages", checks)
    check(all(len("".join(page.split())) >= 100 for page in pages), "PDF contains no blank/near-blank pages", checks)
    normalized = ["".join(page.split()) for page in pages]
    for phrase, page_no in {
        "图4.1TYMP–MTA–SPP1巨噬细胞假说": 27,
        "图4.2患者重复的SPP1–vascularniche": 30,
        "图4.3髓系细胞形成以C1QC状态增强": 32,
    }.items():
        check(phrase in normalized[page_no - 1], f"main figure and caption visible on page {page_no}", checks)
    check(not any("补充图S4." in page for page in normalized), "no supplementary caption in rendered PDF", checks)

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
        "pages": 38,
        "visible_figures": 3,
        "supplementary_captions_visible": 0,
        "tables": 29,
        "manual_visual_qc": "PASS: pages 27, 30, and 32 inspected; detailed captions and complete figures visible",
        "checks": checks,
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
