#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


OUTDIR = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_compact_revision_v6")
SOURCE = OUTDIR / "project26210_第四部分空间转录组分析报告_仅三张主图精简修订版_v8.docx"
OUTPUT = OUTDIR / "project26210_第四部分空间转录组分析报告_仅三张主图精简修订版_v9_compat.docx"
METADATA = OUTDIR / "compact_report_generation_metadata_v9_compat.json"
EXPECTED_SOURCE_SHA256 = "f1aa2a49728598dfd948c2f5e9bee199a736dcf63b83b50c6fba2c21611fcf03"

MEDIA = {
    "word/media/part4_main_figure_4_1_v4.png": OUTDIR / "figures/Figure1_TYMP_MTA_SPP1_real_data_v9.png",
    "word/media/part4_main_figure_4_2_v4.png": OUTDIR / "figures/Figure2_SPP1_vascular_niche_real_data_v9.png",
    "word/media/part4_main_figure_4_3_v4.png": OUTDIR / "figures/Figure3_myeloid_continuum_real_data_v9.png",
}
FLOWCHART = OUTDIR / "article_logic_flowchart_v2.png"
FLOWCHART_MEDIA = "word/media/article_logic_flowchart_v2.png"
LEAD = (
    "本报告围绕一个主问题展开：TYMP相关恶性细胞状态是否对应可重复的髓系与血管空间生态，并能否在独立患者中收敛为可信的人组织证据。"
    "分析不再按工具或补丁时间排列，而按‘发现—外部验证—机制收敛—证据分级’组织。"
)
FLOWCHART_CAPTION = (
    "报告逻辑主线图。报告按‘数据基础与统计单位—分析与验证策略—三条生物学主线—证据分级与边界—结论表达’展开；"
    "图4.1、图4.2和图4.3分别聚焦巨噬细胞状态、血管生态位/T细胞平衡和髓系横断面连续谱。"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def replace_once(value: str, old: str, new: str, label: str) -> str:
    count = value.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return value.replace(old, new)


def drawing_paragraph(document: str, rid: str) -> tuple[int, int, str]:
    marker = f'embed="{rid}"'
    hit = document.index(marker)
    start = document.rfind("<ns0:p>", 0, hit)
    end = document.index("</ns0:p>", hit) + len("</ns0:p>")
    if start < 0:
        raise RuntimeError(f"Cannot locate drawing paragraph for {rid}")
    return start, end, document[start:end]


def main() -> None:
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("Source fingerprint changed")
    for path in (*MEDIA.values(), FLOWCHART):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in (OUTPUT, METADATA):
        if path.exists():
            raise FileExistsError(path)

    with ZipFile(SOURCE) as archive:
        names = archive.namelist()
        document = archive.read("word/document.xml").decode("utf-8")
        rels = archive.read("word/_rels/document.xml.rels").decode("utf-8")
        core = archive.read("docProps/core.xml").decode("utf-8")

    document = replace_once(
        document,
        "逻辑重构版 v8｜仅三张主图精简修订版",
        "逻辑重构版 v9｜仅三张主图精简修订版",
        "version",
    )
    document = replace_once(
        document,
        "【精简修订版规则】正文内容与表格完整保留，仅嵌入图4.1–图4.3三张主题组合主图；每张主图含7–8个真实数据面板，不含流程、网络或解释框类示意面板。原始单项图未在本版本嵌入，完整图集见主图整合完整版。",
        "【精简修订版规则】正文内容与表格完整保留；前置1张报告逻辑主线图，正文仅嵌入图4.1–图4.3三张主题组合主图。每张主图含7–8个真实数据面板，不含流程、网络或解释框类示意面板。原始单项图未在本版本嵌入，完整图集见主图整合完整版。",
        "compact rule",
    )
    document = replace_once(document, "图内四色键联合表示", "图下方四项横排图例联合表示", "Figure 4.1 caption")
    document = replace_once(document, "图内四色键标示", "图下方四项横排图例标示", "Figure 4.2 caption")

    # Preserve the source namespace/prefix layout: clone the known-good first drawing as raw OOXML.
    _, _, first_drawing = drawing_paragraph(document, "rId82")
    flow_drawing = first_drawing
    flow_drawing = flow_drawing.replace('cx="6250000" cy="3810976"', 'cx="6250000" cy="2520984"')
    flow_drawing = flow_drawing.replace('id="60" name="part4_main_figure_4_1_v4.png"', 'id="90" name="article_logic_flowchart_v2.png"')
    flow_drawing = flow_drawing.replace('id="0" name="part4_main_figure_4_1_v4.png"', 'id="0" name="article_logic_flowchart_v2.png"')
    flow_drawing = flow_drawing.replace('embed="rId82"', 'embed="rId85"')
    flow_caption = (
        '<ns0:p><ns0:pPr><ns0:keepLines /><ns0:jc ns0:val="center" />'
        '<ns0:spacing ns0:before="60" ns0:after="120" /></ns0:pPr><ns0:r><ns0:rPr>'
        f'<ns0:b /></ns0:rPr><ns0:t>{escape(FLOWCHART_CAPTION)}</ns0:t></ns0:r></ns0:p>'
    )
    lead_hit = document.index(LEAD)
    lead_end = document.index("</ns0:p>", lead_hit) + len("</ns0:p>")
    document = document[:lead_end] + flow_drawing + flow_caption + document[lead_end:]

    # Correct figure heights to native aspect ratios without touching other OOXML.
    for rid, cy in (("rId82", "3664053"), ("rId83", "3739095"), ("rId84", "3837939")):
        start, end, paragraph = drawing_paragraph(document, rid)
        paragraph = paragraph.replace('cx="6250000" cy="3810976"', f'cx="6250000" cy="{cy}"')
        document = document[:start] + paragraph + document[end:]

    relationship = (
        '<Relationship Id="rId85" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" '
        'Target="media/article_logic_flowchart_v2.png" />'
    )
    rels = replace_once(rels, "</Relationships>", relationship + "</Relationships>", "relationship insertion")
    core = core.replace("仅三张主图精简修订版v8", "仅三张主图精简修订版v9")

    replacements = {
        "word/document.xml": document.encode("utf-8"),
        "word/_rels/document.xml.rels": rels.encode("utf-8"),
        "docProps/core.xml": core.encode("utf-8"),
        **{name: path.read_bytes() for name, path in MEDIA.items()},
    }
    with ZipFile(SOURCE) as source_zip, ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
        for name in names:
            if name in replacements:
                output_zip.writestr(name, replacements[name])
            else:
                output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
        output_zip.writestr(FLOWCHART_MEDIA, FLOWCHART.read_bytes())

    with ZipFile(OUTPUT) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("ZIP CRC failure")
        if len([n for n in archive.namelist() if n.startswith("word/media/")]) != 4:
            raise RuntimeError("Unexpected media count")
    now = datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds")
    metadata = {
        "status": "generated_pending_pdf_render_qc",
        "generated_at": now,
        "source": str(SOURCE),
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
        "output_size": OUTPUT.stat().st_size,
        "logic_flowchart_inserted": 1,
        "main_figures_replaced": 3,
        "tables_preserved": 29,
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
