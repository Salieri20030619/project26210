from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


NAVY = "23364D"
BLUE = "2E74B5"
DARK_BLUE = "1F4D78"
TEAL = "0D8C82"
MAGENTA = "B72C63"
GOLD = "E99A2E"
INK = "20252B"
MUTED = "66717C"
LIGHT = "F2F4F7"
CALLOUT = "F4F6F9"
SOFT_TEAL = "E8F5F3"
SOFT_MAGENTA = "F7E8EF"
WHITE = "FFFFFF"
CONTENT_DXA = 9360


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def set_run_font(run, size=None, bold=None, italic=None, color=INK, east_asia="Microsoft YaHei"):
    run.font.name = "Calibri"
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), "Calibri")
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), "Calibri")
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), east_asia)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic
    if color:
        run.font.color.rgb = RGBColor.from_string(color)


def set_style_font(style, size, color, bold=False, east_asia="Microsoft YaHei"):
    style.font.name = "Calibri"
    style._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), "Calibri")
    style._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), "Calibri")
    style._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), east_asia)
    style.font.size = Pt(size)
    style.font.color.rgb = RGBColor.from_string(color)
    style.font.bold = bold


def set_cell_shading(cell, fill):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=120, bottom=80, end=120):
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for tag, value in [("top", top), ("start", start), ("bottom", bottom), ("end", end)]:
        node = tc_mar.find(qn(f"w:{tag}"))
        if node is None:
            node = OxmlElement(f"w:{tag}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_table_borders(table, color="D7DDE3", size="5", inside=True):
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.find(qn("w:tblBorders"))
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    names = ["top", "left", "bottom", "right"] + (["insideH", "insideV"] if inside else [])
    for name in names:
        node = borders.find(qn(f"w:{name}"))
        if node is None:
            node = OxmlElement(f"w:{name}")
            borders.append(node)
        node.set(qn("w:val"), "single")
        node.set(qn("w:sz"), size)
        node.set(qn("w:space"), "0")
        node.set(qn("w:color"), color)


def set_table_geometry(table, widths_dxa, indent=120):
    if sum(widths_dxa) != CONTENT_DXA:
        raise ValueError(f"Table widths must sum to {CONTENT_DXA}: {widths_dxa}")
    table.autofit = False
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    tbl = table._tbl
    tbl_pr = tbl.tblPr
    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), str(CONTENT_DXA))
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), str(indent))
    tbl_ind.set(qn("w:type"), "dxa")
    layout = tbl_pr.find(qn("w:tblLayout"))
    if layout is None:
        layout = OxmlElement("w:tblLayout")
        tbl_pr.append(layout)
    layout.set(qn("w:type"), "fixed")
    grid = tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths_dxa:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)
    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            width = widths_dxa[min(idx, len(widths_dxa) - 1)]
            tc_pr = cell._tc.get_or_add_tcPr()
            tc_w = tc_pr.find(qn("w:tcW"))
            if tc_w is None:
                tc_w = OxmlElement("w:tcW")
                tc_pr.append(tc_w)
            tc_w.set(qn("w:w"), str(width))
            tc_w.set(qn("w:type"), "dxa")
            set_cell_margins(cell)


def set_repeat_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def add_table(doc, headers, rows, widths, font_size=9.0, header_fill=LIGHT):
    table = doc.add_table(rows=1, cols=len(headers))
    set_table_geometry(table, widths)
    set_table_borders(table)
    set_repeat_header(table.rows[0])
    for idx, header in enumerate(headers):
        cell = table.rows[0].cells[idx]
        set_cell_shading(cell, header_fill)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        p = cell.paragraphs[0]
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        run = p.add_run(str(header))
        set_run_font(run, size=9, bold=True, color=NAVY)
    for row_data in rows:
        row = table.add_row()
        for idx, value in enumerate(row_data):
            cell = row.cells[idx]
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            p.paragraph_format.line_spacing = 1.05
            run = p.add_run(str(value))
            set_run_font(run, size=font_size, color=INK)
    set_table_geometry(table, widths)
    return table


def add_callout(doc, title, text, fill=CALLOUT, accent=BLUE):
    table = doc.add_table(rows=1, cols=1)
    set_table_geometry(table, [CONTENT_DXA])
    set_table_borders(table, color=fill, size="2", inside=False)
    cell = table.cell(0, 0)
    set_cell_shading(cell, fill)
    set_cell_margins(cell, top=130, bottom=130, start=180, end=180)
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(3)
    r = p.add_run(title)
    set_run_font(r, size=10.5, bold=True, color=accent)
    p2 = cell.add_paragraph()
    p2.paragraph_format.space_after = Pt(0)
    p2.paragraph_format.line_spacing = 1.08
    r2 = p2.add_run(text)
    set_run_font(r2, size=10, color=INK)
    return table


def add_paragraph(doc, text="", style=None, bold_prefix=None, color=INK, size=None, italic=False, align=None, after=None):
    p = doc.add_paragraph(style=style)
    if align is not None:
        p.alignment = align
    if after is not None:
        p.paragraph_format.space_after = Pt(after)
    if bold_prefix and text.startswith(bold_prefix):
        r1 = p.add_run(bold_prefix)
        set_run_font(r1, size=size, bold=True, color=color)
        r2 = p.add_run(text[len(bold_prefix):])
        set_run_font(r2, size=size, color=color, italic=italic)
    else:
        r = p.add_run(text)
        set_run_font(r, size=size, color=color, italic=italic)
    return p


def add_bullet(doc, text, level=0):
    p = doc.add_paragraph(style="List Bullet" if level == 0 else "List Bullet 2")
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.line_spacing = 1.15
    r = p.add_run(text)
    set_run_font(r, size=10.5, color=INK)
    return p


def add_number(doc, text):
    p = doc.add_paragraph(style="List Number")
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.line_spacing = 1.15
    r = p.add_run(text)
    set_run_font(r, size=10.5, color=INK)
    return p


def set_alt_text(inline_shape, title, description):
    doc_pr = inline_shape._inline.docPr
    doc_pr.set("title", title)
    doc_pr.set("descr", description)


def add_figure(doc, path, width_in, caption, alt_title, alt_description):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(5)
    p.paragraph_format.space_after = Pt(4)
    shape = p.add_run().add_picture(str(path), width=Inches(width_in))
    set_alt_text(shape, alt_title, alt_description)
    cap = doc.add_paragraph(style="Caption")
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.paragraph_format.space_after = Pt(6)
    r = cap.add_run(caption)
    set_run_font(r, size=8.5, color=MUTED, italic=True)
    return shape


def add_page_number(paragraph):
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    fld_sep = OxmlElement("w:fldChar")
    fld_sep.set(qn("w:fldCharType"), "separate")
    text = OxmlElement("w:t")
    text.text = "1"
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    run._r.extend([fld_begin, instr, fld_sep, text, fld_end])
    set_run_font(run, size=8.5, color=MUTED)


def configure_styles(doc):
    styles = doc.styles
    normal = styles["Normal"]
    set_style_font(normal, 10.5, INK)
    normal.paragraph_format.space_before = Pt(0)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.10
    normal.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

    h1 = styles["Heading 1"]
    set_style_font(h1, 16, BLUE, True)
    h1.paragraph_format.space_before = Pt(16)
    h1.paragraph_format.space_after = Pt(8)
    h1.paragraph_format.keep_with_next = True
    h2 = styles["Heading 2"]
    set_style_font(h2, 13, BLUE, True)
    h2.paragraph_format.space_before = Pt(12)
    h2.paragraph_format.space_after = Pt(6)
    h2.paragraph_format.keep_with_next = True
    h3 = styles["Heading 3"]
    set_style_font(h3, 12, DARK_BLUE, True)
    h3.paragraph_format.space_before = Pt(8)
    h3.paragraph_format.space_after = Pt(4)
    h3.paragraph_format.keep_with_next = True
    cap = styles["Caption"]
    set_style_font(cap, 8.5, MUTED, False)
    cap.paragraph_format.space_before = Pt(0)
    cap.paragraph_format.space_after = Pt(6)
    cap.paragraph_format.keep_with_next = False
    for list_name in ["List Bullet", "List Bullet 2", "List Number"]:
        style = styles[list_name]
        set_style_font(style, 10.5, INK)
        style.paragraph_format.space_after = Pt(4)
        style.paragraph_format.line_spacing = 1.15
    if "Source Note" not in [style.name for style in styles]:
        style = styles.add_style("Source Note", WD_STYLE_TYPE.PARAGRAPH)
    source_note = styles["Source Note"]
    set_style_font(source_note, 8, MUTED)
    source_note.paragraph_format.space_before = Pt(4)
    source_note.paragraph_format.space_after = Pt(4)


def configure_section(section):
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(1.0)
    section.header_distance = Inches(0.492)
    section.footer_distance = Inches(0.492)
    section.different_first_page_header_footer = True
    header = section.header
    p = header.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT
    r = p.add_run("肾透明细胞癌外部免疫多组学整合分析")
    set_run_font(r, size=8.5, bold=True, color=MUTED)
    p.add_run(" " * 6)
    r2 = p.add_run("客户报告")
    set_run_font(r2, size=8.5, color=MUTED)
    footer = section.footer
    fp = footer.paragraphs[0]
    add_page_number(fp)
    first_header = section.first_page_header
    first_header.paragraphs[0].text = ""
    first_footer = section.first_page_footer
    first_footer.paragraphs[0].text = ""


def add_heading(doc, text, level=1):
    p = doc.add_paragraph(style=f"Heading {level}")
    r = p.add_run(text)
    set_run_font(r, size={1: 16, 2: 13, 3: 12}[level], bold=True, color={1: BLUE, 2: BLUE, 3: DARK_BLUE}[level])
    return p


def page_break(doc):
    p = doc.add_paragraph()
    p.add_run().add_break(WD_BREAK.PAGE)


def assert_source_values(project_root: Path):
    feature = pd.read_csv(project_root / "phase65_part1_immune_mechanistic_convergence_local_v1" / "inputs_v1" / "gse207493_celltype_feature_summary_v1.tsv", sep="\t")
    retained = pd.read_csv(project_root / "phase67_part1_direct_program_replication_local_v1" / "results_v5" / "phase67_cross_dataset_retained_v1.tsv", sep="\t")
    concordance = pd.read_csv(project_root / "phase87_three_cohort_myeloid_trajectory_local_v1" / "results_v1" / "phase87_dataset_method_concordance_summary_v1.tsv", sep="\t")
    communication = pd.read_csv(project_root / "phase88_three_cohort_fine_state_cellchat_local_v1" / "results_v1" / "phase88_dataset_endpoint_summary_v1.tsv", sep="\t")
    spatial = pd.read_csv(project_root / "phase69_85_server_status_sync_local_v1" / "phase78" / "analysis_v1" / "phase78_preparation_meta_v1.tsv", sep="\t")
    myeloid_r2 = retained[(retained.endpoint == "myeloid_dc_tymp") & (retained.replication_grade == "R2_TWO_EXTERNAL_FDR")]
    treg = retained[(retained.endpoint == "treg_tymp") & (retained.datasets_fdr_replication >= 1)]
    assert len(myeloid_r2) == 18
    assert len(treg) == 13
    assert set(concordance.dataset) == {"GSE207493", "GSE264586", "GSE294109"}
    assert len(communication[communication.tier.str.contains("PRIMARY", na=False)]) == 18
    assert len(spatial) == 10
    return feature, retained, concordance, communication, spatial


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit("Usage: 02_build_report_v1.py <project_root> <assets_dir> <main_figure> <output_docx>")
    project_root = Path(sys.argv[1])
    assets = Path(sys.argv[2])
    main_figure = Path(sys.argv[3])
    output = Path(sys.argv[4])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    feature, retained, concordance, communication, spatial = assert_source_values(project_root)
    now = datetime.now().astimezone()

    doc = Document()
    configure_styles(doc)
    configure_section(doc.sections[0])
    core = doc.core_properties
    core.title = "肾透明细胞癌外部免疫多组学整合分析报告"
    core.subject = "TYMP相关髓系程序、巨噬细胞状态、空间生态与细胞通讯"
    core.author = "项目分析团队"
    core.keywords = "ccRCC; TYMP; myeloid; macrophage; spatial transcriptomics; CellChat"
    core.comments = "从零重构的客户报告；正文不使用内部阶段编号。"

    # Cover: editorial_cover pattern with a restrained scientific accent.
    add_paragraph(doc, "外部免疫多组学整合分析", size=11, color=GOLD, align=WD_ALIGN_PARAGRAPH.CENTER, after=16)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(10)
    r = p.add_run("肾透明细胞癌")
    set_run_font(r, size=28, bold=True, color=NAVY)
    p2 = doc.add_paragraph()
    p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p2.paragraph_format.space_after = Pt(8)
    r2 = p2.add_run("外部免疫多组学整合分析报告")
    set_run_font(r2, size=22, bold=True, color=NAVY)
    add_paragraph(doc, "TYMP相关髓系程序、巨噬细胞状态、空间生态与临床转化边界", size=13, color=TEAL, align=WD_ALIGN_PARAGRAPH.CENTER, after=20)
    add_figure(
        doc, assets / "panel_A_analysis_architecture.png", 6.1,
        "分析架构：发现队列、外部单细胞验证、空间验证与临床转化检验。",
        "免疫多组学分析架构", "展示单细胞发现、外部验证、空间切片和临床队列之间的分析关系。",
    )
    add_paragraph(doc, f"报告日期：{now.strftime('%Y年%m月%d日')}  |  客户交付版", size=10, color=MUTED, align=WD_ALIGN_PARAGRAPH.CENTER, after=6)
    add_paragraph(doc, "本报告为全新逻辑重构，不沿用既有队列顺序或内部工作编号。", size=9.5, color=MUTED, italic=True, align=WD_ALIGN_PARAGRAPH.CENTER)
    page_break(doc)

    # Reading guide and executive summary.
    add_heading(doc, "阅读导航", 1)
    toc_items = [
        "执行摘要", "研究问题与证据路线", "队列与分析架构", "免疫图谱与TYMP细胞来源",
        "TYMP相关免疫程序的跨队列复现", "髓系状态排序", "空间免疫生态",
        "细胞通讯特异性", "bulk可辨识性与临床转化", "综合结论与验证建议", "统计方法与质量控制",
    ]
    for item in toc_items:
        add_number(doc, item)
    add_heading(doc, "执行摘要", 1)
    add_callout(
        doc,
        "核心结论",
        "TYMP与广义髓系免疫代谢/调控程序及巨噬细胞富集生态存在可重复关联；但现有证据不支持将其解释为选择性SPP1巨噬细胞极化、5-MTA介导的完整因果链，或可用于response/PFS/OS预测的临床标志物。",
        fill=SOFT_TEAL,
        accent=TEAL,
    )
    add_bullet(doc, "三个独立单细胞队列共覆盖45位患者、389,477个细胞；每个队列独立处理和展示，不进行跨队列细胞池化推断。")
    add_bullet(doc, "18个TYMP去锚定的Myeloid/DC程序在两个外部单细胞队列中同时复现，形成稳定的缺氧、呼吸链、ROS、糖代谢和炎症转录调控核心。")
    add_bullet(doc, "Monocyte到巨噬细胞状态的横断面排序在三个队列和两种轨迹方法中重复出现，但不等同于真实谱系转换或治疗诱导。")
    add_bullet(doc, "24张空间切片未支持肿瘤TYMP到邻近SPP1状态的双制备关联；三队列细胞通讯也未支持SPP1相对C1QC/Other的特异性。")
    add_bullet(doc, "临床复合评分、校正后的PFS模型和bulk精细状态投影均未完成临床升级；后续需要匹配代谢组、状态读出和干预/救援设计。")
    page_break(doc)

    add_heading(doc, "1 研究问题与证据路线", 1)
    add_paragraph(doc, "本报告围绕一个中心问题展开：肿瘤或微环境中的TYMP升高，是否对应可重复的髓系免疫程序，并进一步形成选择性SPP1巨噬细胞状态、空间免疫生态和免疫治疗结局？")
    add_paragraph(doc, "为避免将相关性直接升级为机制或临床结论，证据按由近及远的顺序组织：先确定细胞来源和谱系内程序，再做外部单细胞复现；随后检验髓系状态排序、空间邻域和细胞通讯；最后进入bulk可辨识性和临床结局压力测试。")
    route_rows = [
        ("细胞来源定位", "TYMP与F1分别主要出现在哪些细胞谱系？", "GSE207493单细胞图谱"),
        ("程序复现", "TYMP相关代谢与TF调控程序能否跨队列复现？", "GSE264586、GSE294109"),
        ("状态排序", "Monocyte与三类巨噬细胞状态是否呈稳定顺序？", "Slingshot与PCA-kNN测地拟时序"),
        ("空间生态", "肿瘤TYMP是否邻近SPP1状态或血管信号？", "12张FFPE与12张冰冻切片"),
        ("细胞通讯", "SPP1状态是否具有相对C1QC/Other的通讯特异性？", "三队列患者级CellChat"),
        ("临床转化", "程序或状态是否预测response/PFS/OS？", "CheckMate与JAVELIN"),
    ]
    add_table(doc, ["分析名称", "核心问题", "证据层"], route_rows, [1900, 4200, 3260], font_size=8.8, header_fill=SOFT_TEAL)
    add_callout(doc, "解释原则", "单细胞复现用于支持机制背景；空间、通讯和临床模型用于检验特异性与转化。任一远端证据失败，均不能由更近端的表达相关性替代。", fill=CALLOUT, accent=BLUE)

    add_heading(doc, "2 队列与分析架构", 1)
    cohort_rows = [
        ("发现型单细胞", "GSE207493", "19位患者；160,710个细胞", "细胞来源、代谢程序、TF调控和髓系状态"),
        ("外部单细胞验证", "GSE264586", "12位患者；88,152个细胞", "Myeloid/DC程序复现、状态排序和通讯"),
        ("外部单细胞验证", "GSE294109", "15位患者；140,615个细胞", "Myeloid/DC与Treg程序复现、状态排序和通讯"),
        ("空间转录组", "GSE175540", "24张切片；FFPE与冰冻各12张", "细胞状态映射与跨区室邻域检验"),
        ("临床转录组", "CheckMate", "311份RNA样本", "校正后的PFS与治疗交互检验"),
        ("临床转录组", "JAVELIN Renal 101", "726份RNA样本", "PFS与治疗交互检验"),
    ]
    add_table(doc, ["数据层", "队列", "规模", "分析角色"], cohort_rows, [1600, 1900, 2400, 3460], font_size=8.5)
    add_paragraph(doc, "三个单细胞队列分别保留自身的归一化、降维、注释和统计单位。跨队列输出是患者层面的方向与FDR复现，不是合并后的细胞级检验。")
    add_figure(doc, main_figure, 5.55, "图1  外部免疫多组学整合主图。各面板依次展示分析架构、三队列图谱、细胞来源、程序复现、髓系状态排序、空间生态和通讯特异性。", "免疫多组学整合主图", "六个面板展示从单细胞发现到空间、通讯和临床压力测试的完整逻辑。")
    page_break(doc)

    add_heading(doc, "3 免疫图谱与TYMP细胞来源", 1)
    add_heading(doc, "3.1 三个队列的独立细胞图谱", 2)
    add_paragraph(doc, "三个队列均可分辨肿瘤/上皮、髓系、T/NK、内皮、B/浆细胞和基质等主要区室。图谱仅用于队列内结构展示；图形配色一致，但坐标不做跨队列对齐。")
    add_figure(doc, assets / "panel_B_three_cohort_atlas.png", 6.15, "图2  三个队列的独立UMAP图谱。", "三队列单细胞图谱", "分别展示GSE207493、GSE264586和GSE294109的主要细胞区室。")
    add_heading(doc, "3.2 TYMP与F1的细胞来源", 2)
    add_paragraph(doc, "在发现队列中，F1投影主要集中于肿瘤细胞；TYMP总体上在Myeloid/DC中最强，并在Treg中呈现较弱但可评估的富集。Myeloid细胞的TYMP平均表达为1.764，阳性比例为87.96%；Dendritic细胞分别为1.597和88.08%。")
    add_paragraph(doc, "F1与Treg之间的bulk关联更符合组成/生态背景，而不是Treg内源F1程序。报告因此将F1肿瘤程序、Myeloid/DC-TYMP程序和Treg-TYMP程序分开分析。")
    add_figure(doc, assets / "panel_C_cellular_source_programs.png", 6.15, "图3  TYMP与F1的细胞来源、谱系特征矩阵及程序覆盖。", "TYMP和F1细胞来源", "展示TYMP表达、F1投影、谱系特征气泡矩阵及三类程序覆盖。")
    add_callout(doc, "本节结论", "TYMP不是单一肿瘤细胞标志。其最稳定的转录背景位于Myeloid/DC；Treg结果需要保持为单外部队列证据。", fill=SOFT_MAGENTA, accent=MAGENTA)
    page_break(doc)

    add_heading(doc, "4 TYMP相关免疫程序的跨队列复现", 1)
    add_heading(doc, "4.1 发现队列中的谱系内程序", 2)
    add_paragraph(doc, "在控制细胞周期后，F1-high肿瘤细胞有22/31条目标代谢通路通过FDR；TYMP-high Myeloid/DC为17/31条，TYMP-high Treg为10/31条。Myeloid/DC程序集中于缺氧、呼吸链、ROS解毒、糖酵解/糖异生、磷酸戊糖途径和嘌呤代谢。")
    add_paragraph(doc, "de novo SCENIC与CollecTRI-VIPER在发现队列中形成176个同方向双方法FDR结果，其中肿瘤F1为95、Myeloid/DC-TYMP为67、Treg-TYMP为14。该一致性代表单队列方法收敛，不替代外部复现。")
    add_heading(doc, "4.2 外部单细胞直接复现", 2)
    add_paragraph(doc, "使用冻结且排除TYMP本身的程序进行患者内残差相关分析。Myeloid/DC共有18个签名在GSE264586与GSE294109中同时达到方向一致和FDR：9个代谢程序与9个TF调控子。")
    metabolic_names = "缺氧反应、呼吸链、ROS解毒、糖酵解、糖异生、磷酸戊糖途径、嘌呤分解等"
    regulon_names = "IRF1、SPI1、RELB、STAT6、IRF7、STAT2、ZBTB7A、CEBPB和RELA"
    add_bullet(doc, f"代谢核心：{metabolic_names}。")
    add_bullet(doc, f"调控核心：{regulon_names}。")
    add_bullet(doc, "Treg共有13个签名在GSE294109达到FDR；GSE264586缺少可用于该终点的Treg注释，因此不能称为双外部队列复现。")
    add_figure(doc, assets / "panel_D_cross_cohort_replication.png", 6.30, "图4  TYMP相关免疫程序的跨队列复现。左侧为患者级效应桥，右侧为双外部队列Myeloid/DC代谢程序-TYMP核心-TF调控子双轨图。", "跨队列免疫程序复现", "展示18个双队列Myeloid/DC签名和13个单外部Treg签名。")
    add_callout(doc, "本节结论", "外部数据支持的是广义Myeloid/DC免疫代谢与转录调控核心，而不是TYMP到SPP1状态的特异因果边。", fill=SOFT_TEAL, accent=TEAL)
    page_break(doc)

    add_heading(doc, "5 髓系状态排序", 1)
    add_paragraph(doc, "在三个队列中分别提取Monocyte、SPP1-state macrophage、C1QC-state macrophage和Other macrophage，共68,711个细胞、45位队列内患者和62个患者-样本单元。")
    add_paragraph(doc, "Slingshot与独立PCA-kNN测地拟时序使用同一Monocyte根规则。两种方法的患者内一致性中位数分别为0.841、0.559和0.430；三个队列均通过方向一致性门槛。")
    add_bullet(doc, "SPP1-state在三个队列、两种方法中均晚于Monocyte。")
    add_bullet(doc, "Monocyte程序在三个队列的两种方法中均随排序下降。")
    add_bullet(doc, "SPP1去锚定程序在GSE207493和GSE294109双方法支持；GSE264586仅测地方法支持，Slingshot未通过。")
    add_figure(doc, assets / "panel_E_myeloid_ordering.png", 5.85, "图5  三队列髓系状态排序及程序动态。", "髓系状态排序", "展示三个队列的状态分布、轨迹方向和沿拟时序的程序变化。")
    add_callout(doc, "解释边界", "该结果支持横断面状态顺序，不证明Monocyte真实转化为SPP1巨噬细胞，也不证明治疗诱导、时间先后、5-MTA介导或因果极化。", fill=CALLOUT, accent=BLUE)
    page_break(doc)

    add_heading(doc, "6 空间免疫生态", 1)
    add_paragraph(doc, "使用同一13状态参考将24张GSE175540切片分别映射，12张为FFPE、12张为冰冻制备；24/24完成且无失败。肿瘤TYMP、SPP1/C1QC巨噬细胞状态和内皮信号均可在切片内形成局部空间结构。")
    add_paragraph(doc, "核心检验不是观察空间图是否“看起来相邻”，而是用肿瘤区室校正后的连续TYMP检验邻域状态丰度，并要求FFPE和冰冻制备同时支持。")
    spatial_rows = [
        ("肿瘤TYMP - 邻近SPP1状态", "FFPE", "rho=0.0089；q=0.3928", "不支持"),
        ("肿瘤TYMP - 邻近SPP1状态", "冰冻", "rho=0.0099；q=0.3802", "不支持"),
        ("肿瘤TYMP - 邻近C1QC状态", "FFPE", "rho=0.0268；q=0.0252", "单制备支持"),
        ("肿瘤TYMP - 邻近C1QC状态", "冰冻", "rho=0.0304；q=0.1471", "不支持"),
    ]
    add_table(doc, ["空间关系", "制备", "汇总效应", "判定"], spatial_rows, [3000, 1200, 2600, 2560], font_size=8.7, header_fill=SOFT_MAGENTA)
    add_figure(doc, assets / "panel_F_spatial_maps.png", 6.15, "图6  代表性FFPE与冰冻切片中的肿瘤TYMP、SPP1/C1QC巨噬细胞状态和内皮丰度。", "空间免疫生态", "展示两种制备的代表性切片和四类连续空间信号。")
    add_callout(doc, "本节结论", "空间数据允许报告广义巨噬细胞/血管生态背景，但不支持双制备的肿瘤TYMP到邻近SPP1状态桥。", fill=SOFT_MAGENTA, accent=MAGENTA)
    page_break(doc)

    add_heading(doc, "7 细胞通讯特异性", 1)
    add_paragraph(doc, "CellChat在三个单细胞队列中分别运行，保留36个患者-样本单元、28位队列内患者、46,964个确定性抽样细胞和161,865条配体-受体记录。所有推断使用冻结的细胞级状态标签，不使用bulk投影状态。")
    add_paragraph(doc, "主检验关注SPP1-state相对于C1QC-state与Other macrophage平均水平的通讯特异性，包括肿瘤到SPP1、SPP1到肿瘤、SPP1到T细胞、SPP1到内皮以及总流出/流入六个终点。")
    comm_rows = [
        ("GSE207493", "16", "六个特异性终点均未通过FDR", "绝对肿瘤-SPP1双向强度与TYMP相关，但状态特异性消失"),
        ("GSE264586", "4", "描述性", "患者数小且TYMP变异不足"),
        ("GSE294109", "6-8", "未达到正式样本量与FDR门槛", "方向多为正，但统计证据不足"),
    ]
    add_table(doc, ["队列", "有效患者数", "特异性判定", "补充解释"], comm_rows, [1700, 1500, 2900, 3260], font_size=8.5)
    add_figure(doc, assets / "panel_F_cellchat_chord_specificity.png", 6.15, "图7  GSE207493患者中位绝对通讯弦图及三队列SPP1状态特异性矩阵。弦带宽度按患者中位通讯概率的平方根缩放，颜色代表发送细胞。", "CellChat弦图与特异性矩阵", "展示六个细胞群之间的绝对通讯以及三队列SPP1特异性终点。")
    add_callout(doc, "本节结论", "绝对通讯强度可以随广义细胞组成或整体活性升高；目前没有证据证明SPP1状态相对C1QC/Other具有可重复的特异通讯优势。", fill=CALLOUT, accent=BLUE)
    page_break(doc)

    add_heading(doc, "8 bulk可辨识性与临床转化", 1)
    add_heading(doc, "8.1 精细巨噬细胞状态能否从bulk中恢复", 2)
    add_paragraph(doc, "使用与CheckMate/JAVELIN共同平台一致的478基因面板，在38位患者的细胞级状态表达上构建1,824个确定性伪bulk混合物，并比较冻结参考、去锚定参考、留一队列参考和三状态参考。")
    add_paragraph(doc, "所有预注册门槛均失败。即使在纯巨噬细胞混合物中，Other macrophage仍被系统性估计为零，SPP1状态出现补偿性高估。组织情境中的总巨噬细胞负荷仍可部分恢复（rho约0.745-0.894），但SPP1/C1QC/Other不可作为可靠bulk细分状态。")
    add_heading(doc, "8.2 临床复合评分与校正后的PFS模型", 2)
    add_paragraph(doc, "18个双外部队列Myeloid/DC签名被等家族加权形成不含TYMP的连续复合评分，并在E-MTAB-3218、CheckMate和JAVELIN中独立测试。七个主要临床模型没有一个通过FDR。")
    clinical_rows = [
        ("Myeloid/DC复合评分", "7个预设response/PFS模型", "0/7通过FDR", "机制复现未转化为临床预测"),
        ("完整链代理", "CheckMate PFS", "HR=1.017；q=0.844", "不支持"),
        ("完整链代理", "JAVELIN PFS", "HR=1.114；q=0.546", "方向一致但不支持"),
        ("TYMP-SPP1联合特征", "CheckMate PFS", "HR=0.928；q=0.546", "方向冲突"),
        ("TYMP-SPP1联合特征", "JAVELIN PFS", "HR=1.082；q=0.546", "方向冲突"),
        ("投影SPP1状态比例", "CheckMate/JAVELIN", "HR=0.970 / 1.130", "跨队列方向冲突"),
    ]
    add_table(doc, ["特征", "终点/队列", "结果", "解释"], clinical_rows, [2300, 2200, 2350, 2510], font_size=8.5)
    add_heading(doc, "8.3 因果闭环与湿实验数据边界", 2)
    add_paragraph(doc, "虚拟扰动、网络富集、MTAP自然扰动和跨情境MTA响应签名均未形成两外部队列一致的TYMP-5-MTA-SPP1闭环。现有客户文件也未提供可重算的样本级5-MTA丰度矩阵、重复结构和明确注释，因此不能把客户报告的倍数变化视为本项目独立复算结论。")
    add_callout(doc, "临床与因果结论", "直接F1/TYMP临床预测、选择性SPP1状态bulk定量和5-MTA介导的完整因果链均不成立为当前报告结论。", fill=SOFT_MAGENTA, accent=MAGENTA)
    page_break(doc)

    add_heading(doc, "9 综合结论与验证建议", 1)
    conclusion_rows = [
        ("外部可重复", "TYMP相关Myeloid/DC免疫代谢与TF调控核心；广义巨噬细胞富集生态；Monocyte到巨噬细胞的横断面状态排序"),
        ("背景性支持", "部分空间血管/巨噬细胞共现；绝对通讯强度随生态组成变化"),
        ("当前不支持", "选择性SPP1极化、SPP1特异通讯、可靠bulk三状态分解、response/PFS/OS预测、5-MTA因果介导"),
    ]
    add_table(doc, ["证据层级", "当前可报告内容"], conclusion_rows, [1900, 7460], font_size=9.2, header_fill=SOFT_TEAL)
    add_paragraph(doc, "建议后续验证按照能直接改变证据等级的顺序推进：")
    recommendations = [
        "在同一批样本中同步获取TYMP、5-MTA靶向定量、巨噬细胞状态和空间邻域信息，建立可审计的样本清单与重复结构。",
        "采用肿瘤细胞TYMP干预，并设置5-MTA补充/清除或代谢救援条件，区分TYMP表达效应与5-MTA介导效应。",
        "在单核/巨噬细胞体系中同步检测SPP1、C1QC、抗原呈递、补体-溶酶体和免疫抑制读出，而不是只检测SPP1单基因。",
        "空间验证使用配对FFPE/冰冻或相同制备的独立样本，预先冻结跨区室距离和邻域统计规则。",
        "临床验证只使用可在独立队列复现且不依赖不可辨识bulk细状态的粗粒度Myeloid/DC程序。",
    ]
    for item in recommendations:
        add_number(doc, item)
    add_callout(doc, "最终报告口径", "TYMP-5-MTA-SPP1应保留为具有广义髓系背景支持的待验证生物学假说，而不是已经建立的选择性SPP1机制、完整因果链或临床预测标志物。", fill=SOFT_TEAL, accent=TEAL)

    add_heading(doc, "10 统计方法与质量控制", 1)
    method_rows = [
        ("单细胞队列", "每个队列独立处理；患者为推断单位", "不进行跨队列细胞池化检验"),
        ("程序复现", "患者内TYMP残差相关；队列内FDR；外部队列方向/FDR门槛", "冻结且排除锚定基因的签名"),
        ("状态排序", "Slingshot与PCA-kNN测地拟时序", "双方法一致性；不作谱系因果解释"),
        ("空间分析", "切片内Spearman与置换；制备内汇总；跨制备门槛", "连续信号而非主观选图"),
        ("细胞通讯", "患者/样本独立CellChat；患者中位聚合；六终点FDR", "区分绝对强度与状态特异性"),
        ("临床模型", "预设协变量Cox/Logistic/交互模型；队列内终点校正", "校正事件编码；不复用旧方向"),
        ("bulk可辨识性", "确定性伪bulk、留一队列和去锚定敏感性", "细状态不可辨识时禁止临床使用"),
    ]
    add_table(doc, ["分析名称", "统计框架", "关键边界"], method_rows, [1900, 4100, 3360], font_size=8.4)
    add_heading(doc, "质量控制摘要", 2)
    add_bullet(doc, "主图使用389,477个细胞来源对象中的115,605个确定性绘图点；分析结果使用完整权威表，不以抽样图替代统计。")
    add_bullet(doc, "髓系状态排序验证33/33通过；精细状态通讯验证97/97通过；主图数据、SVG、PDF与渲染一致性验证35/35通过。")
    add_bullet(doc, "空间代表切片按接近制备层总体效应选择，未挑选极端显著切片。")
    add_bullet(doc, "所有报告结论均按分析名称表达，未使用内部工作编号作为证据。")

    add_heading(doc, "附录 关键数字索引", 1)
    appendix_rows = [
        ("单细胞规模", "45位患者；389,477个细胞", "三个队列独立"),
        ("Myeloid/DC直接复现", "18个双外部队列签名", "9代谢 + 9调控子"),
        ("Treg直接复现", "13个单外部队列签名", "另一队列不可评估"),
        ("髓系状态排序", "68,711个细胞；45位患者；62个样本单元", "三队列、两方法"),
        ("空间映射", "24/24切片完成", "12 FFPE + 12冰冻"),
        ("CellChat", "36个样本单元；28位患者；161,865条LR记录", "六个状态特异终点"),
        ("伪bulk可辨识性", "1,824个混合物；478基因", "三状态不可辨识"),
        ("临床主模型", "7个复合评分模型；0个通过FDR", "无临床等级升级"),
    ]
    add_table(doc, ["项目", "规模/结果", "说明"], appendix_rows, [2300, 3900, 3160], font_size=8.7, header_fill=LIGHT)
    add_paragraph(doc, "报告结束", size=9, color=MUTED, align=WD_ALIGN_PARAGRAPH.CENTER, after=0)

    doc.save(output)
    visible_text = "\n".join(p.text for p in doc.paragraphs)
    visible_text += "\n" + "\n".join(cell.text for table in doc.tables for row in table.rows for cell in row.cells)
    forbidden = [r"Phase\s*\d+", r"phase\s*\d+", "五队列", "5队列"]
    hits = {pattern: re.findall(pattern, visible_text) for pattern in forbidden}
    if any(hits.values()):
        raise SystemExit(f"Forbidden visible narrative found: {hits}")
    status = {
        "status": "GENERATED_PENDING_RENDER_QA",
        "output": str(output),
        "bytes": output.stat().st_size,
        "sha256": sha256(output),
        "paragraphs": len(doc.paragraphs),
        "tables": len(doc.tables),
        "inline_shapes": len(doc.inline_shapes),
        "forbidden_visible_hits": hits,
        "design_preset": "standard_business_brief",
        "header_pattern": "editorial_cover",
    }
    (output.parent / "BUILD_STATUS_v1.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
