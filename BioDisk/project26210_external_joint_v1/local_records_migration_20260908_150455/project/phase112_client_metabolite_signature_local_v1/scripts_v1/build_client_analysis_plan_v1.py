from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_ROW_HEIGHT_RULE, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(r"D:\project26210_external_joint_v1\phase112_client_metabolite_signature_local_v1")
OUT = ROOT / "output_v1"
ASSETS = ROOT / "assets_v1"
OUT.mkdir(parents=True, exist_ok=True)
ASSETS.mkdir(parents=True, exist_ok=True)

DOCX_PATH = OUT / "项目26210_无权重七特征公共数据迁移分析方案_v1.docx"
FIG_PATH = ASSETS / "seven_feature_analysis_workflow_v1.png"

BLUE = "2E74B5"
DARK_BLUE = "17365D"
NAVY = "0B2545"
MID_BLUE = "3F6F91"
PALE_BLUE = "E8EEF5"
LIGHT_GRAY = "F2F4F7"
CALLOUT = "F4F6F9"
MUTED = "626B75"
GOLD = "A97800"
WHITE = "FFFFFF"
BLACK = "1F2328"


def rgb(hex_color: str) -> RGBColor:
    return RGBColor.from_string(hex_color)


def set_run_font(run, size=None, color=None, bold=None, italic=None, latin="Calibri", east_asia="Microsoft YaHei"):
    run.font.name = latin
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), latin)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), latin)
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), east_asia)
    if size is not None:
        run.font.size = Pt(size)
    if color is not None:
        run.font.color.rgb = rgb(color)
    if bold is not None:
        run.bold = bold
    if italic is not None:
        run.italic = italic


def set_cell_shading(cell, fill: str):
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
    for name, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{name}"))
        if node is None:
            node = OxmlElement(f"w:{name}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_cell_width(cell, width_dxa: int):
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(width_dxa))
    tc_w.set(qn("w:type"), "dxa")


def set_table_geometry(table, widths_dxa, indent_dxa=120):
    table.alignment = WD_TABLE_ALIGNMENT.LEFT
    table.autofit = False
    tbl_pr = table._tbl.tblPr
    tbl_w = tbl_pr.find(qn("w:tblW"))
    if tbl_w is None:
        tbl_w = OxmlElement("w:tblW")
        tbl_pr.append(tbl_w)
    tbl_w.set(qn("w:w"), str(sum(widths_dxa)))
    tbl_w.set(qn("w:type"), "dxa")
    tbl_ind = tbl_pr.find(qn("w:tblInd"))
    if tbl_ind is None:
        tbl_ind = OxmlElement("w:tblInd")
        tbl_pr.append(tbl_ind)
    tbl_ind.set(qn("w:w"), str(indent_dxa))
    tbl_ind.set(qn("w:type"), "dxa")
    grid = table._tbl.tblGrid
    for child in list(grid):
        grid.remove(child)
    for width in widths_dxa:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(width))
        grid.append(col)
    for row in table.rows:
        for idx, cell in enumerate(row.cells):
            set_cell_width(cell, widths_dxa[idx])
            set_cell_margins(cell)


def set_table_borders(table, color="D8DDE3", size=4):
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.find(qn("w:tblBorders"))
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        tag = borders.find(qn(f"w:{edge}"))
        if tag is None:
            tag = OxmlElement(f"w:{edge}")
            borders.append(tag)
        tag.set(qn("w:val"), "single")
        tag.set(qn("w:sz"), str(size))
        tag.set(qn("w:color"), color)


def mark_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    header = OxmlElement("w:tblHeader")
    header.set(qn("w:val"), "true")
    tr_pr.append(header)


def add_page_number(paragraph):
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("第 ")
    set_run_font(run, size=9, color=MUTED)
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    r = OxmlElement("w:r")
    r_pr = OxmlElement("w:rPr")
    r_fonts = OxmlElement("w:rFonts")
    r_fonts.set(qn("w:ascii"), "Calibri")
    r_fonts.set(qn("w:hAnsi"), "Calibri")
    r_pr.append(r_fonts)
    sz = OxmlElement("w:sz")
    sz.set(qn("w:val"), "18")
    r_pr.append(sz)
    r.append(r_pr)
    t = OxmlElement("w:t")
    t.text = "1"
    r.append(t)
    fld.append(r)
    paragraph._p.append(fld)
    run2 = paragraph.add_run(" 页")
    set_run_font(run2, size=9, color=MUTED)


def add_bullet_num_definition(doc, num_fmt="bullet", glyph="•"):
    numbering = doc.part.numbering_part.element
    existing_abs = [int(x.get(qn("w:abstractNumId"))) for x in numbering.findall(qn("w:abstractNum"))]
    existing_num = [int(x.get(qn("w:numId"))) for x in numbering.findall(qn("w:num"))]
    abstract_id = max(existing_abs + [0]) + 1
    num_id = max(existing_num + [0]) + 1
    abstract = OxmlElement("w:abstractNum")
    abstract.set(qn("w:abstractNumId"), str(abstract_id))
    multi = OxmlElement("w:multiLevelType")
    multi.set(qn("w:val"), "singleLevel")
    abstract.append(multi)
    lvl = OxmlElement("w:lvl")
    lvl.set(qn("w:ilvl"), "0")
    start = OxmlElement("w:start")
    start.set(qn("w:val"), "1")
    lvl.append(start)
    fmt = OxmlElement("w:numFmt")
    fmt.set(qn("w:val"), num_fmt)
    lvl.append(fmt)
    text = OxmlElement("w:lvlText")
    text.set(qn("w:val"), glyph if num_fmt == "bullet" else "%1.")
    lvl.append(text)
    jc = OxmlElement("w:lvlJc")
    jc.set(qn("w:val"), "left")
    lvl.append(jc)
    ppr = OxmlElement("w:pPr")
    tabs = OxmlElement("w:tabs")
    tab = OxmlElement("w:tab")
    tab.set(qn("w:val"), "num")
    tab.set(qn("w:pos"), "720")
    tabs.append(tab)
    ppr.append(tabs)
    ind = OxmlElement("w:ind")
    ind.set(qn("w:left"), "720")
    ind.set(qn("w:hanging"), "360")
    ppr.append(ind)
    lvl.append(ppr)
    rpr = OxmlElement("w:rPr")
    rfonts = OxmlElement("w:rFonts")
    rfonts.set(qn("w:ascii"), "Calibri")
    rfonts.set(qn("w:hAnsi"), "Calibri")
    rfonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    rpr.append(rfonts)
    lvl.append(rpr)
    abstract.append(lvl)
    numbering.append(abstract)
    num = OxmlElement("w:num")
    num.set(qn("w:numId"), str(num_id))
    abs_id = OxmlElement("w:abstractNumId")
    abs_id.set(qn("w:val"), str(abstract_id))
    num.append(abs_id)
    numbering.append(num)
    return num_id


def add_list_item(doc, text, num_id, bold_prefix=None, after=8):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(after)
    p.paragraph_format.line_spacing = 1.167
    p.paragraph_format.left_indent = Inches(0.18)
    ppr = p._p.get_or_add_pPr()
    pbdr = OxmlElement("w:pBdr")
    left = OxmlElement("w:left")
    left.set(qn("w:val"), "single")
    left.set(qn("w:sz"), "12")
    left.set(qn("w:space"), "7")
    left.set(qn("w:color"), "B8CEE1")
    pbdr.append(left)
    ppr.append(pbdr)
    if bold_prefix and text.startswith(bold_prefix):
        r1 = p.add_run(bold_prefix)
        set_run_font(r1, size=11, color=BLACK, bold=True)
        r2 = p.add_run(text[len(bold_prefix):])
        set_run_font(r2, size=11, color=BLACK)
    else:
        r = p.add_run(text)
        set_run_font(r, size=11, color=BLACK)
    return p


def add_callout(doc, label, text, fill=CALLOUT, accent=BLUE):
    table = doc.add_table(rows=1, cols=1)
    set_table_geometry(table, [9360], indent_dxa=180)
    set_table_borders(table, color="D8DDE3", size=4)
    cell = table.cell(0, 0)
    set_cell_shading(cell, fill)
    set_cell_margins(cell, top=150, bottom=150, start=180, end=180)
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.line_spacing = 1.15
    r1 = p.add_run(label + "  ")
    set_run_font(r1, size=10.5, color=accent, bold=True)
    r2 = p.add_run(text)
    set_run_font(r2, size=10.5, color=BLACK)
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def add_heading(doc, text, level=1):
    p = doc.add_paragraph(text, style=f"Heading {level}")
    p.paragraph_format.keep_with_next = True
    return p


def add_para_callout(doc, label, text, fill=CALLOUT, accent=BLUE):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.08)
    p.paragraph_format.right_indent = Inches(0.08)
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(8)
    p.paragraph_format.line_spacing = 1.12
    ppr = p._p.get_or_add_pPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:fill'), fill)
    ppr.append(shd)
    pbdr = OxmlElement('w:pBdr')
    for edge in ('top', 'left', 'bottom', 'right'):
        node = OxmlElement(f'w:{edge}')
        node.set(qn('w:val'), 'single')
        node.set(qn('w:sz'), '4')
        node.set(qn('w:space'), '6')
        node.set(qn('w:color'), 'D8DDE3')
        pbdr.append(node)
    ppr.append(pbdr)
    r1 = p.add_run(label + '  ')
    set_run_font(r1, size=10.2, color=accent, bold=True)
    r2 = p.add_run(text)
    set_run_font(r2, size=10.2, color=BLACK)
    return p


def add_body(doc, text, bold_prefix=None, align=WD_ALIGN_PARAGRAPH.LEFT, after=6):
    p = doc.add_paragraph()
    p.alignment = align
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(after)
    p.paragraph_format.line_spacing = 1.10
    if bold_prefix and text.startswith(bold_prefix):
        r1 = p.add_run(bold_prefix)
        set_run_font(r1, size=11, color=BLACK, bold=True)
        r2 = p.add_run(text[len(bold_prefix):])
        set_run_font(r2, size=11, color=BLACK)
    else:
        r = p.add_run(text)
        set_run_font(r, size=11, color=BLACK)
    return p


def add_page_break(doc):
    p = doc.paragraphs[-1]
    p.add_run().add_break(WD_BREAK.PAGE)


def add_table_text(cell, text, bold=False, color=BLACK, size=9.5, align=WD_ALIGN_PARAGRAPH.LEFT):
    p = cell.paragraphs[0]
    p.alignment = align
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.line_spacing = 1.08
    r = p.add_run(text)
    set_run_font(r, size=size, color=color, bold=bold)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def make_workflow_figure(path: Path):
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    chosen = "Microsoft YaHei" if "Microsoft YaHei" in fonts else ("SimHei" if "SimHei" in fonts else "DejaVu Sans")
    plt.rcParams["font.family"] = chosen
    plt.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(figsize=(12.2, 4.2))
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 4.2)
    ax.axis("off")
    boxes = [
        (0.25, 2.35, 2.15, 1.2, "01 数据与元数据核查", "确认样本独立性\n冻结24张切片输入"),
        (2.7, 2.35, 2.15, 1.2, "02 空间生态位构建", "发现SPP1巨噬细胞\n与内皮细胞生态位"),
        (5.15, 2.35, 2.15, 1.2, "03 通讯空间验证", "验证限定的双向\n配体-受体关系"),
        (7.6, 2.35, 2.15, 1.2, "04 TLS与免疫梯度", "解析生态位与TLS\n及免疫功能的关系"),
        (10.05, 2.35, 1.9, 1.2, "05 机制整合", "生态位转录程序\n更新报告主线"),
    ]
    for idx, (x, y, w, h, title, sub) in enumerate(boxes):
        face = "#E8EEF5" if idx < 4 else "#EAF2EA"
        edge = "#3F6F91" if idx < 4 else "#447A4B"
        rect = plt.Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, linewidth=1.6)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center", fontsize=10.2, fontweight="bold", color="#17365D")
        ax.text(x + w / 2, y + h * 0.28, sub, ha="center", va="center", fontsize=8.8, color="#3C4650", linespacing=1.35)
        if idx < len(boxes) - 1:
            next_x = boxes[idx + 1][0]
            ax.annotate("", xy=(next_x - 0.08, y + h / 2), xytext=(x + w + 0.08, y + h / 2), arrowprops=dict(arrowstyle="->", color="#7E8994", lw=1.5))
    ax.text(0.35, 1.1, "执行原则", fontsize=10.5, fontweight="bold", color="#A97800")
    ax.text(1.55, 1.1, "复用既有Cell2location结果", fontsize=9.4, color="#3C4650")
    ax.text(4.1, 1.1, "FFPE发现 / frozen验证", fontsize=9.4, color="#3C4650")
    ax.text(6.6, 1.1, "以切片而非spot作为统计样本", fontsize=9.4, color="#3C4650")
    ax.text(9.65, 1.1, "不以算法替代因果实验", fontsize=9.4, color="#3C4650")
    ax.plot([0.35, 11.85], [0.82, 0.82], color="#D5DBE2", lw=1)
    ax.text(0.35, 0.43, "计划输出：生态位图谱 · 空间通讯证据 · TLS/免疫梯度 · 候选通路与报告更新", fontsize=10.1, color="#17365D", fontweight="bold")
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_evidence_route_figure(path: Path):
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    chosen = 'Microsoft YaHei' if 'Microsoft YaHei' in fonts else ('SimHei' if 'SimHei' in fonts else 'DejaVu Sans')
    plt.rcParams['font.family'] = chosen
    plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(figsize=(12.4, 6.0))
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 6.0)
    ax.axis('off')

    def box(x, y, w, h, title, sub, face, edge):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, linewidth=1.7))
        ax.text(x + w / 2, y + h * 0.66, title, ha='center', va='center', fontsize=10.2, fontweight='bold', color='#17365D')
        ax.text(x + w / 2, y + h * 0.28, sub, ha='center', va='center', fontsize=8.6, color='#3C4650', linespacing=1.25)

    box(0.25, 4.15, 2.25, 1.15, 'F1 / TYMP肿瘤状态', '17/20患者方向为正\n关联有异质性', '#E8EEF5', '#2E74B5')
    box(3.35, 4.15, 2.25, 1.15, '直接测量MTA', 'UnitedMet肿瘤n=255：无支持\nNCI60全谱系n=58：无支持', '#FCEBEC', '#B94A48')
    box(6.45, 4.15, 2.25, 1.15, 'SPP1巨噬细胞状态', '分群和状态顺序已建立\nMTA诱导SPP1尚未直接证明', '#FFF5D9', '#A97800')
    ax.annotate('', xy=(3.22, 4.73), xytext=(2.58, 4.73), arrowprops=dict(arrowstyle='->', color='#B94A48', lw=2.0, linestyle='--'))
    ax.text(2.90, 4.94, '未闭合', ha='center', fontsize=8.4, color='#9B2C2C', fontweight='bold')
    ax.annotate('', xy=(6.32, 4.73), xytext=(5.68, 4.73), arrowprops=dict(arrowstyle='->', color='#A97800', lw=2.0, linestyle='--'))
    ax.text(6.00, 4.94, '待扰动验证', ha='center', fontsize=8.4, color='#8A6200', fontweight='bold')
    box(9.55, 4.55, 2.45, 1.05, '路线A：空间结构', 'SPP1—内皮通讯与血管关联\n已有数据集特异性支持', '#EAF4EA', '#447A4B')
    box(9.55, 2.75, 2.45, 1.05, '路线B：预后 / 治疗', '响应、PFS、OS与治疗交互\n当前均未形成稳健验证', '#FCEBEC', '#B94A48')
    ax.annotate('', xy=(9.42, 5.02), xytext=(8.78, 4.85), arrowprops=dict(arrowstyle='->', color='#447A4B', lw=2.0))
    ax.annotate('', xy=(9.42, 3.30), xytext=(8.78, 4.48), arrowprops=dict(arrowstyle='->', color='#B94A48', lw=2.0, linestyle='--'))
    ax.text(0.35, 2.80, '下一步只补真实断点', fontsize=11.2, fontweight='bold', color='#17365D')
    ax.text(0.35, 2.30, '空间线：固定生态位 → 限定LR空间验证 → TLS/免疫梯度 → 生态位伪bulk', fontsize=9.7, color='#3C4650')
    ax.text(0.35, 1.82, '临床线：先冻结空间表型 → 再投射独立临床队列 → 响应/PFS/OS与治疗交互', fontsize=9.7, color='#3C4650')
    ax.text(0.35, 1.34, '因果线：TYMP扰动 + 直接MTA + SPP1 RNA/蛋白 + A2B/STAT3阻断/救援', fontsize=9.7, color='#3C4650')
    ax.plot([0.35, 12.0], [0.92, 0.92], color='#D5DBE2', lw=1)
    ax.text(0.35, 0.48, '颜色：蓝=已有支持  绿=下游数据集特异性支持  黄=待直接验证  红=当前数据不支持/未闭合', fontsize=9.2, color='#626B75')
    fig.savefig(path, dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)


make_evidence_route_figure(FIG_PATH)

doc = Document()
section = doc.sections[0]
section.page_width = Inches(8.5)
section.page_height = Inches(11)
section.top_margin = Inches(1.0)
section.bottom_margin = Inches(1.0)
section.left_margin = Inches(1.0)
section.right_margin = Inches(1.0)
section.header_distance = Inches(0.492)
section.footer_distance = Inches(0.492)

styles = doc.styles
normal = styles["Normal"]
normal.font.name = "Calibri"
normal._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
normal._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
normal.font.size = Pt(11)
normal.font.color.rgb = rgb(BLACK)
normal.paragraph_format.space_before = Pt(0)
normal.paragraph_format.space_after = Pt(6)
normal.paragraph_format.line_spacing = 1.10

for name, size, color, before, after in [
    ("Heading 1", 16, BLUE, 16, 8),
    ("Heading 2", 13, BLUE, 12, 6),
    ("Heading 3", 12, DARK_BLUE, 8, 4),
]:
    style = styles[name]
    style.font.name = "Calibri"
    style._element.rPr.rFonts.set(qn("w:ascii"), "Calibri")
    style._element.rPr.rFonts.set(qn("w:hAnsi"), "Calibri")
    style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    style.font.size = Pt(size)
    style.font.bold = True
    style.font.color.rgb = rgb(color)
    style.paragraph_format.space_before = Pt(before)
    style.paragraph_format.space_after = Pt(after)
    style.paragraph_format.keep_with_next = True

header = section.header
hp = header.paragraphs[0]
hp.alignment = WD_ALIGN_PARAGRAPH.LEFT
hp.paragraph_format.space_after = Pt(0)
hr1 = hp.add_run("项目26210  |  下一阶段分析计划")
set_run_font(hr1, size=9, color=MUTED, bold=True)
footer = section.footer
fp = footer.paragraphs[0]
add_page_number(fp)

bullet_id = add_bullet_num_definition(doc, "bullet", "•")
number_id = add_bullet_num_definition(doc, "decimal")

# Cover
doc.add_paragraph().paragraph_format.space_after = Pt(38)
kicker = doc.add_paragraph()
kicker.paragraph_format.space_after = Pt(4)
kr = kicker.add_run("NEXT-STAGE RESEARCH PLAN")
set_run_font(kr, size=9.5, color=GOLD, bold=True)

title = doc.add_paragraph()
title.paragraph_format.space_before = Pt(0)
title.paragraph_format.space_after = Pt(8)
tr = title.add_run("项目26210下一阶段\n生物信息分析计划")
set_run_font(tr, size=29, color=NAVY, bold=True)

subtitle = doc.add_paragraph()
subtitle.paragraph_format.space_after = Pt(22)
sr = subtitle.add_run("围绕TYMP—代谢—SPP1巨噬细胞—血管/免疫空间生态的定向补充")
set_run_font(sr, size=13.5, color=MID_BLUE)

meta = doc.add_table(rows=2, cols=2)
set_table_geometry(meta, [1800, 7560], indent_dxa=120)
set_table_borders(meta, color="D8DDE3", size=4)
for row in meta.rows:
    row.height_rule = WD_ROW_HEIGHT_RULE.AT_LEAST
for i, (label, value) in enumerate([
    ("文件用途", "客户沟通与下一阶段工作确认"),
    ("版本日期", "V1.0  |  2026年9月1日"),
]):
    set_cell_shading(meta.cell(i, 0), LIGHT_GRAY)
    add_table_text(meta.cell(i, 0), label, bold=True, color=DARK_BLUE, size=10)
    add_table_text(meta.cell(i, 1), value, size=10)

doc.add_paragraph().paragraph_format.space_after = Pt(12)
add_callout(
    doc,
    "计划核心",
    "下一阶段不再重复已有去卷积和常规相关性分析，而是优先深度利用已完成预处理的24张空间转录组切片，构建可跨制备方式验证的SPP1巨噬细胞—内皮细胞生态位，并进一步解析其与TLS及免疫功能的空间关系。",
    fill="EEF4F8",
    accent=DARK_BLUE,
)

lead = doc.add_paragraph()
lead.paragraph_format.space_before = Pt(16)
lead.paragraph_format.space_after = Pt(10)
lead.paragraph_format.line_spacing = 1.2
lr = lead.add_run("本计划以现有证据链中的真实缺口为导向，优先选择能够提高证据等级的算法和数据，而不是单纯增加分析数量。")
set_run_font(lr, size=12, color=BLACK, bold=True)

add_body(doc, "现有结果已经为肿瘤细胞TYMP/F1状态、SPP1巨噬细胞分群以及SPP1—血管空间关联建立了基础；下一阶段将重点补充空间组织方式、细胞通讯定位和免疫生态后果。")

# Page 2
add_page_break(doc)
add_heading(doc, "一、现有基础与下一阶段重点", 1)
add_body(doc, "本项目已完成多个单细胞、空间转录组、代谢组和临床队列的分层分析。当前最适合继续深化的部分，是已经获得多层证据但仍缺少空间机制解释的SPP1巨噬细胞—内皮/血管轴。")

comparison = doc.add_table(rows=4, cols=2)
mark_table_header(comparison.rows[0])
set_table_geometry(comparison, [4680, 4680], indent_dxa=120)
set_table_borders(comparison, color="D8DDE3", size=4)
headers = ["已建立的分析基础", "下一阶段需要补充的内容"]
for j, text in enumerate(headers):
    set_cell_shading(comparison.cell(0, j), PALE_BLUE)
    add_table_text(comparison.cell(0, j), text, bold=True, color=DARK_BLUE, size=10.2, align=WD_ALIGN_PARAGRAPH.CENTER)
rows = [
    ("已完成巨噬细胞亚群划分，并形成SPP1、C1QC及其他巨噬细胞状态。", "将细胞状态从“丰度相关”提升为空间生态位和组织结构。"),
    ("单细胞通讯分析显示SPP1巨噬细胞与内皮细胞存在稳定的双向通讯相容性。", "定位具体配体—受体对，并检验其是否在空间邻域中重复出现。"),
    ("24张空间切片支持SPP1状态与血管重塑的空间关联。", "进一步解析该生态位与TLS、细胞毒性及T细胞功能状态的关系。"),
]
for i, pair in enumerate(rows, 1):
    for j, text in enumerate(pair):
        add_table_text(comparison.cell(i, j), text, size=9.5)
        if i % 2 == 0:
            set_cell_shading(comparison.cell(i, j), "FAFBFC")

add_heading(doc, "总体目标", 2)
for text in [
    "系统利用24张空间转录组切片，建立可重复的SPP1巨噬细胞—内皮细胞空间生态位图谱。",
    "把单细胞通讯结果与空间邻域中的具体配体—受体表达相连接。",
    "判断该生态位更接近免疫排斥、TLS邻近免疫组织化，或不同样本中的情境依赖状态。",
    "形成可直接纳入最终报告的空间主图、统计表和结论性文字。",
]:
    add_list_item(doc, text, bullet_id)

add_callout(doc, "数据优势", "24张切片均已有可复用的原始计数、空间坐标、组织图像和细状态后验结果，其中12张为FFPE、12张为frozen；无需重新运行耗时的Cell2location。", fill=CALLOUT, accent=GOLD)

# Page 3
add_page_break(doc)
add_heading(doc, "二、重点工作：24张空间转录组切片深度挖掘", 1)

add_heading(doc, "1. 构建细状态空间生态位图谱", 2)
add_body(doc, "基于每张切片的空间邻接图和细状态丰度，识别由肿瘤细胞、SPP1巨噬细胞、C1QC巨噬细胞、内皮细胞、单核细胞及淋巴细胞共同构成的局部微环境。生态位发现阶段不使用TYMP、TLS或临床结局信息，以降低结果导向偏倚。")
for text in [
    "在FFPE切片中建立生态位中心，并将其保持不变地转移至frozen切片验证。",
    "比较保守丰度估计与后验均值，检验生态位对参数选择的稳定性。",
    "重点判断SPP1巨噬细胞—内皮细胞生态位是否在两类制备方式中均可重复识别。",
]:
    add_list_item(doc, text, bullet_id, after=5)

add_heading(doc, "2. 验证双向配体—受体的空间定位", 2)
add_body(doc, "从既有单细胞患者级通讯结果中预先筛选有限数量的SPP1巨噬细胞→内皮细胞及内皮细胞→SPP1巨噬细胞候选关系，再检验发送细胞、配体表达、受体表达和接收细胞是否在相邻spot中形成方向一致的空间组合。")
for text in [
    "主要检验家族限制为每个方向最多5对，避免事后扩大终点。",
    "分别在FFPE和frozen切片中汇总，不把两类切片简单合并。",
    "结果定位为空间通讯相容性，不表述为直接结合或因果诱导。",
]:
    add_list_item(doc, text, bullet_id, after=5)

add_heading(doc, "3. 分析TLS边界与免疫功能梯度", 2)
add_body(doc, "利用官方提供的23张TLS区域标注，计算生态位到TLS内部、边界和远端肿瘤区域的空间距离，评估SPP1—内皮生态位与细胞毒性、T细胞功能障碍及HLA-II/APC支持之间的连续梯度。")
add_callout(doc, "该分析回答的问题", "SPP1—血管生态位究竟对应免疫排斥，还是位于免疫组织化/TLS相关区域，或呈现样本依赖的双重状态。两种方向均为可接受结果。", fill="FFF9EB", accent=GOLD)

# Page 4
add_page_break(doc)
add_heading(doc, "三、机制深化与外部验证策略", 1)

add_heading(doc, "4. 生态位特异转录程序与通路", 2)
add_body(doc, "在空间生态位标签冻结后，对每张切片分别构建生态位区域与匹配非生态位区域的伪bulk计数，比较其差异基因和通路。匹配过程中控制肿瘤细胞丰度、总巨噬细胞丰度、测序深度和组织面积。")
for text in [
    "统计单位为独立切片，而不是spot数量。",
    "FFPE与frozen分别建模，仅将方向一致的结果提升为主要候选。",
    "将重复基因与单细胞来源表达及空间配体—受体结果交叉核对。",
]:
    add_list_item(doc, text, bullet_id, after=5)

add_heading(doc, "5. 条件性TYMP—SPP1空间关系", 2)
add_body(doc, "既有结果不支持普遍性的肿瘤TYMP→SPP1巨噬细胞关系。因此，本项只作为次级分析：在不使用TYMP信息形成的生态位标签中，检验低氧/血管生成和POSTN/ECM环境是否构成可重复的条件因素。即使出现阳性结果，也仅表述为情境依赖，不改变总体结论。")

add_heading(doc, "6. 外部数据集的条件性补充", 2)
add_body(doc, "新增数据集不作为“寻找阳性结果”的手段，而只用于验证已经在24张切片中冻结的空间表型。候选数据包括：")
for text in [
    "GSE210038/GSE210041：用于固定空间生态位的定位和迁移验证；空间样本量较小，主要承担定位而非总体推断。",
    "PRJNA705464：用于多区域肿瘤、邻近肾组织、外周血和淋巴结中的免疫背景验证。",
    "SCP1288：仅在形成固定临床表型后，用于晚期RCC免疫治疗语境的敏感性分析。",
]:
    add_list_item(doc, text, bullet_id, after=5)

add_heading(doc, "7. 与客户验证实验的接口", 2)
add_body(doc, "直接TYMP→MTA以及MTA→SPP1仍需扰动和代谢物/蛋白测量完成因果闭环。本阶段将为客户实验保留明确的结果接口和图位，但不提前介入、预测或重绘客户数据。")

# Page 5
add_page_break(doc)
add_heading(doc, "四、实施顺序与阶段性产出", 1)
pic = doc.add_picture(str(FIG_PATH), width=Inches(6.5))
pic._inline.docPr.set("title", "下一阶段分析流程")
pic._inline.docPr.set("descr", "从数据核查、空间生态位构建、通讯空间验证、TLS与免疫梯度到机制整合的五阶段分析流程")
pic_par = doc.paragraphs[-1]
pic_par.alignment = WD_ALIGN_PARAGRAPH.CENTER
pic_par.paragraph_format.space_after = Pt(5)
cap = doc.add_paragraph()
cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
cap.paragraph_format.space_after = Pt(10)
cr = cap.add_run("图1  下一阶段分析流程与主要输出")
set_run_font(cr, size=9, color=MUTED, italic=True)

schedule = doc.add_table(rows=6, cols=4)
mark_table_header(schedule.rows[0])
set_table_geometry(schedule, [1050, 2550, 3510, 2250], indent_dxa=120)
set_table_borders(schedule, color="D8DDE3", size=4)
for j, text in enumerate(["阶段", "工作重点", "核心任务", "阶段产出"]):
    set_cell_shading(schedule.cell(0, j), PALE_BLUE)
    add_table_text(schedule.cell(0, j), text, bold=True, color=DARK_BLUE, size=9.3, align=WD_ALIGN_PARAGRAPH.CENTER)
schedule_rows = [
    ("阶段1", "输入与元数据核查", "确认24张切片完整性、样本独立性和TLS/临床元数据来源", "冻结输入与分析方案"),
    ("阶段2", "空间生态位构建", "FFPE发现、frozen验证，形成全切片生态位分布", "生态位图谱与稳定性结果"),
    ("阶段3", "空间通讯与TLS梯度", "验证限定配体—受体并分析TLS/免疫距离关系", "通讯和免疫生态主结果"),
    ("阶段4", "生态位转录机制", "配对伪bulk、差异基因、通路和细胞来源核对", "候选机制与补充表"),
    ("阶段5", "整合与交付", "合并空间、单细胞及客户验证接口，更新报告主线", "正式图表、统计表和报告内容"),
]
for i, row in enumerate(schedule_rows, 1):
    for j, text in enumerate(row):
        add_table_text(schedule.cell(i, j), text, size=8.8, align=WD_ALIGN_PARAGRAPH.CENTER if j == 0 else WD_ALIGN_PARAGRAPH.LEFT)
        if i % 2 == 0:
            set_cell_shading(schedule.cell(i, j), "FAFBFC")

add_heading(doc, "预期交付内容", 2)
for text in [
    "24张切片的空间生态位全景图、样本一致性图和代表性组织图。",
    "SPP1巨噬细胞—内皮细胞双向通讯的空间配体—受体结果。",
    "TLS边界、细胞毒性及T细胞功能状态的距离梯度图。",
    "生态位特异基因、通路和细胞来源整合表。",
    "更新后的项目主线图、客户报告结果章节及完整统计补充材料。",
]:
    add_list_item(doc, text, bullet_id, after=5)

# Page 6
add_page_break(doc)
add_heading(doc, "五、统计与质量控制", 1)
for text in [
    "独立生物学单位：跨样本统计以切片/患者为单位，spot仅用于描述空间结构。",
    "制备方式控制：FFPE与frozen分别分析，只有方向在两类制备中重复时才称为双制备支持。",
    "候选冻结：生态位发现不使用TYMP、TLS或临床结局；配体—受体候选在查看空间结果前限定。",
    "多重检验：在预先定义的同类假设家族内进行FDR校正，不混合基因、通路、通讯和临床终点。",
    "稳健性与结论边界：比较不同丰度估计、邻域半径和leave-one-section-out结果；空间共现和通讯相容性不写成因果。",
]:
    add_list_item(doc, text, bullet_id, bold_prefix=text.split("：", 1)[0] + "：", after=6)

add_heading(doc, "风险与应对", 2)
risk = doc.add_table(rows=4, cols=3)
mark_table_header(risk.rows[0])
set_table_geometry(risk, [2100, 3780, 3480], indent_dxa=120)
set_table_borders(risk, color="D8DDE3", size=4)
for j, text in enumerate(["潜在情况", "处理方式", "报告口径"]):
    set_cell_shading(risk.cell(0, j), LIGHT_GRAY)
    add_table_text(risk.cell(0, j), text, bold=True, color=DARK_BLUE, size=9.2, align=WD_ALIGN_PARAGRAPH.CENTER)
risk_rows = [
    ("生态位仅在一种制备中稳定", "保留制备差异并开展参数敏感性核查", "定位为制备依赖结果"),
    ("TLS与免疫梯度方向不一致", "报告样本异质性，不扩大亚组筛选", "当前数据不支持一般性方向"),
    ("外部数据缺少完整元数据", "限制为定位或描述，不进入推断性统计", "不作为独立临床验证"),
]
for i, row in enumerate(risk_rows, 1):
    for j, text in enumerate(row):
        add_table_text(risk.cell(i, j), text, size=8.9)
        set_cell_margins(risk.cell(i, j), top=45, bottom=45, start=120, end=120)

add_heading(doc, "六、预期对最终报告的提升", 1)
add_body(doc, "完成后，报告将不再把单细胞、空间和临床结果拆成彼此独立的部分，而是围绕一条证据等级清晰的连续主线组织：")
add_callout(doc, "整合主线", "肿瘤TYMP/F1状态 → 巨噬细胞状态划分 → SPP1巨噬细胞—内皮空间生态位 → 血管/TLS/免疫功能 → 客户扰动验证接口", fill="EEF4F8", accent=DARK_BLUE)
add_body(doc, "其中已获得支持、情境依赖、当前数据不支持以及需要客户实验验证的关系将分别标识，既保证结果主线完整，也避免把不同数据集或不同证据等级直接混用。")

doc.core_properties.title = "项目26210下一阶段生物信息分析计划"
doc.core_properties.subject = "现有数据深度利用与关键证据链定向补充"
doc.core_properties.author = "项目分析团队"
doc.core_properties.keywords = "26210, 空间转录组, SPP1, TYMP, 内皮细胞, TLS"

doc.save(DOCX_PATH)
print(DOCX_PATH)


def build_v2():
    vdoc = Document()
    sec = vdoc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(0.82)
    sec.bottom_margin = Inches(0.82)
    sec.left_margin = Inches(0.82)
    sec.right_margin = Inches(0.82)
    sec.header_distance = Inches(0.38)
    sec.footer_distance = Inches(0.38)

    base = vdoc.styles['Normal']
    base.font.name = 'Calibri'
    base._element.rPr.rFonts.set(qn('w:ascii'), 'Calibri')
    base._element.rPr.rFonts.set(qn('w:hAnsi'), 'Calibri')
    base._element.rPr.rFonts.set(qn('w:eastAsia'), 'Microsoft YaHei')
    base.font.size = Pt(10.5)
    base.font.color.rgb = rgb(BLACK)
    base.paragraph_format.space_before = Pt(0)
    base.paragraph_format.space_after = Pt(5)
    base.paragraph_format.line_spacing = 1.08
    for name, size, color, before, after in [
        ('Heading 1', 16, BLUE, 14, 7),
        ('Heading 2', 13, BLUE, 10, 5),
        ('Heading 3', 11.5, DARK_BLUE, 7, 3),
    ]:
        style = vdoc.styles[name]
        style.font.name = 'Calibri'
        style._element.rPr.rFonts.set(qn('w:ascii'), 'Calibri')
        style._element.rPr.rFonts.set(qn('w:hAnsi'), 'Calibri')
        style._element.rPr.rFonts.set(qn('w:eastAsia'), 'Microsoft YaHei')
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = rgb(color)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True

    hp = sec.header.paragraphs[0]
    hp.alignment = WD_ALIGN_PARAGRAPH.LEFT
    hr = hp.add_run('项目26210  |  双路线证据与补充分析计划')
    set_run_font(hr, size=8.6, color=MUTED, bold=True)
    add_page_number(sec.footer.paragraphs[0])
    bullets = add_bullet_num_definition(vdoc, 'bullet', '•')

    def evidence_table(rows):
        table = vdoc.add_table(rows=len(rows) + 1, cols=4)
        mark_table_header(table.rows[0])
        set_table_geometry(table, [1800, 3200, 1660, 2700], indent_dxa=120)
        set_table_borders(table, color='D8DDE3', size=4)
        for j, text in enumerate(['证据节点', '现有证据', '当前判断', '需要补全']):
            set_cell_shading(table.cell(0, j), PALE_BLUE)
            add_table_text(table.cell(0, j), text, bold=True, color=DARK_BLUE, size=8.8, align=WD_ALIGN_PARAGRAPH.CENTER)
        for i, row in enumerate(rows, 1):
            for j, text in enumerate(row):
                add_table_text(table.cell(i, j), text, size=8.2, align=WD_ALIGN_PARAGRAPH.CENTER if j == 2 else WD_ALIGN_PARAGRAPH.LEFT)
                set_cell_margins(table.cell(i, j), top=55, bottom=55, start=105, end=105)
                if i % 2 == 0:
                    set_cell_shading(table.cell(i, j), 'FAFBFC')
        return table

    def plan_table(rows):
        table = vdoc.add_table(rows=len(rows) + 1, cols=4)
        mark_table_header(table.rows[0])
        set_table_geometry(table, [1050, 2480, 3880, 1950], indent_dxa=120)
        set_table_borders(table, color='D8DDE3', size=4)
        for j, text in enumerate(['优先级', '拟补缺口', '可能分析', '成功 / 停止标准']):
            set_cell_shading(table.cell(0, j), LIGHT_GRAY)
            add_table_text(table.cell(0, j), text, bold=True, color=DARK_BLUE, size=8.7, align=WD_ALIGN_PARAGRAPH.CENTER)
        for i, row in enumerate(rows, 1):
            for j, text in enumerate(row):
                add_table_text(table.cell(i, j), text, size=8.15, align=WD_ALIGN_PARAGRAPH.CENTER if j == 0 else WD_ALIGN_PARAGRAPH.LEFT)
                set_cell_margins(table.cell(i, j), top=55, bottom=55, start=105, end=105)
                if i % 2 == 0:
                    set_cell_shading(table.cell(i, j), 'FAFBFC')
        return table

    vdoc.add_paragraph().paragraph_format.space_after = Pt(24)
    p = vdoc.add_paragraph()
    p.paragraph_format.space_after = Pt(4)
    r = p.add_run('EVIDENCE-LED NEXT-STAGE PLAN')
    set_run_font(r, size=9.2, color=GOLD, bold=True)
    p = vdoc.add_paragraph()
    p.paragraph_format.space_after = Pt(8)
    r = p.add_run('项目26210双路线证据与补充分析计划')
    set_run_font(r, size=25, color=NAVY, bold=True)
    p = vdoc.add_paragraph()
    p.paragraph_format.space_after = Pt(15)
    r = p.add_run('F1/TYMP → MTA → SPP1 → 空间结构  |  预后与治疗')
    set_run_font(r, size=13.2, color=MID_BLUE, bold=True)

    add_body(vdoc, '文件定位：现有证据盘点、关键断点与可执行补充分析', bold_prefix='文件定位：', after=2)
    add_body(vdoc, '版本日期：V2.0  |  2026年9月1日', bold_prefix='版本日期：', after=8)
    add_para_callout(
        vdoc,
        '核心判断',
        'F1/TYMP肿瘤状态与下游SPP1—内皮/血管结构分别已有可复核证据，但TYMP→直接MTA和MTA→SPP1两个中央环节尚未闭合。空间线可继续提高机制深度；预后/治疗线应等待固定空间表型后再进行独立临床验证，不能重复空泛相关性分析。',
        fill='EEF4F8',
        accent=DARK_BLUE,
    )
    pic = vdoc.add_picture(str(FIG_PATH), width=Inches(6.75))
    pic._inline.docPr.set('title', 'F1/TYMP-MTA-SPP1双路线证据状态图')
    pic._inline.docPr.set('descr', '展示F1/TYMP到MTA到SPP1后分为空间结构和预后治疗两条路线，并按证据状态着色')
    vdoc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap = vdoc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cap.paragraph_format.space_after = Pt(0)
    cr = cap.add_run('图1  两条路线的现有证据状态与下一步补缺逻辑')
    set_run_font(cr, size=8.8, color=MUTED, italic=True)

    add_page_break(vdoc)
    add_heading(vdoc, '一、两条路线共享的上游证据：先明确已经知道什么', 1)
    add_body(vdoc, '空间结构和预后/治疗两条路线共享同一上游链条。现有项目结果并非“没有证据”，而是证据集中在链条两端，中央代谢与诱导环节仍缺少直接闭环。')
    evidence_table([
        ('F1与TYMP肿瘤状态', '原始单细胞恶性细胞分析中，20名合格患者有17名方向为正；GSE207493为10/10正向，但整体强度受该研究权重影响。', '支持关联，有异质性', '无需再做通用相关性；保留患者一致性和研究权重说明。'),
        ('TYMP特异功能', '严格TYMP功能轴0/4端点通过；更广泛的TYMP/嘧啶代谢敏感性3/4为次级支持，但缺少轨迹收敛。', '部分/未闭合', '只有细胞来源明确或扰动级数据才能升级，空间签名不能证明酶功能。'),
        ('TYMP到直接测量MTA', 'UnitedMet 255个肿瘤样本：meta rho=0.0395，P=0.538；NCI60全谱系n=58同样无支持；肾癌n=8仅为正向探索线索；CCL180/CCLE不可估计直接MTA。', '当前数据不支持', '需要TYMP扰动、匹配对照、独立生物学重复和直接MTA定量。'),
        ('MTA到SPP1', 'UnitedMet中MTA与bulk SPP1仅名义正向，BH q=0.0976，校正TYMP/MTAP后消失；公开MTA扰动签名未通过跨情境资格门；直接MTA诱导SPP1未被测量。', '待直接验证', 'MTA刺激后测SPP1 RNA/蛋白，并加入A2B或STAT3阻断/救援。'),
        ('F1/TYMP到SPP1巨噬细胞', '原始单细胞20名患者的4个丰度/极化端点为0/4支持；14名患者的恶性细胞↔SPP1细状态通讯端点为0/2支持。', '不支持一般性直接边', '只允许在结局盲的固定生态位中检验条件异质性，不能推翻总体阴性。'),
    ])
    add_para_callout(vdoc, '共同断点', '两条路线都不能绕过TYMP→MTA和MTA→SPP1。新增空间算法可以深化SPP1下游结构，但不能替代直接代谢物与扰动实验。', fill='FFF9EB', accent=GOLD)

    add_page_break(vdoc)
    add_heading(vdoc, '二、路线A：F1/TYMP → MTA → SPP1 → 空间结构', 1)
    add_body(vdoc, '这条路线的中央因果桥尚未闭合，但SPP1巨噬细胞到内皮/血管结构的下游证据已经形成两种独立分析层面的三角验证，因此最适合先做“机制深化”，而不是重新证明是否存在SPP1状态。')
    add_heading(vdoc, '2.1 现有证据', 2)
    evidence_table([
        ('SPP1细状态定义', '71,045个巨噬细胞池已分为Other 40,429、SPP1-state 15,322和C1QC-state 15,294；50,761个髓系细胞、40名患者的研究分开轨迹中，SPP1-state均晚于Monocyte。', '操作性状态已支持', 'SPP1-versus-C1QC并未在5个研究中一致上升，不能写成C1QC替代或谱系转换。'),
        ('SPP1→内皮通讯', '14名患者中14/14方向为正，中位细状态特异性=0.21065，BH q=0.00305。', '支持计算相容性', '需要在空间邻域中验证预先限定的配体—受体；不能写成结合或诱导。'),
        ('内皮→SPP1通讯', '14名患者中13/14方向为正，中位细状态特异性=0.14642，BH q=0.00305。', '支持计算相容性', '同上，需空间方向性与制备方式重复。'),
        ('SPP1→血管重塑', 'GSE175540中，FFPE n=12：中位rho=0.2121，P=0.000488，q=0.00220；frozen n=12：rho=0.2083，P=0.000488，q=0.00440。', '数据集特异性空间支持', '需要生态位结构、TLS位置和转录程序；不能外推到所有空间队列。'),
        ('肿瘤TYMP→邻近SPP1', '同一24张切片中，FFPE rho=0.00894，q=0.3928；frozen rho=0.00985，q=0.3802。', '两种制备均不支持', '后续只能做预先限定的条件异质性，不得把下游阳性当作上游救援。'),
        ('血管结构→免疫排斥', '既有GLMF/TJ-RCC和后续空间汇总结果方向异质或无支持。', '尚未闭合', '需利用官方TLS几何和免疫功能连续梯度直接检验。'),
    ])

    add_heading(vdoc, '2.2 需要补全的空间断点', 2)
    for text in [
        '空间组织方式：目前知道SPP1状态与血管信号相关，但尚不知道它们是否形成可重复的多细胞生态位。',
        '通讯定位：单细胞CellChat支持双向相容性，但尚未证明候选配体、受体、发送细胞和接收细胞在相邻空间邻域共同出现。',
        'TLS与免疫后果：尚不清楚SPP1—血管生态位位于TLS边界、免疫排斥区还是组织化免疫区。',
        '生态位转录机制：尚缺按独立切片汇总的生态位特异基因与通路，现有spot数不能替代生物学样本量。',
    ]:
        add_list_item(vdoc, text, bullets, after=5)

    add_page_break(vdoc)
    add_heading(vdoc, '2.3 可执行的空间补充分析', 2)
    plan_table([
        ('P0', '输入与统计单位', '核对24张已处理空间H5AD、24套细状态q05/均值后验、组织图像及23张官方TLS图；冻结12 FFPE与12 frozen，不重跑Cell2location。', '24/24对象与坐标/元数据一致；任一关键输入不一致则先停止下游。'),
        ('P1-A', '多细胞生态位', '在FFPE中用空间邻接与细状态丰度无监督发现生态位；不使用TYMP、TLS或结局。冻结中心后原样转移到frozen，反向转移仅作敏感性。', '两种制备均识别同方向SPP1—内皮生态位；若只在一种制备稳定，则定位为制备依赖。'),
        ('P1-B', '限定LR空间验证', '从独立的14患者CellChat表预先选择每方向最多5对LR，检验发送细胞/配体与接收细胞/受体的邻域共现；FFPE、frozen分开形成FDR家族。', '限定家族在两类制备中同方向并通过切片级FDR；否则保留通讯相容性但不升级。'),
        ('P1-C', 'TLS与免疫梯度', '用官方TLS多边形计算内部、边界及远端距离；检验细胞毒性、T细胞功能障碍、HLA-II/APC支持和SPP1—内皮生态位的切片级连续梯度。', '免疫排斥或TLS组织化任一方向均可接受，但须在两类制备中重复；方向不一致则明确异质性。'),
        ('P1-D', '生态位转录程序', '生态位标签冻结后，对每张切片构建生态位与匹配非生态位伪bulk；控制肿瘤/巨噬细胞丰度、测序深度和组织面积，做DE与通路并核对细胞来源。', '以独立切片为单位，跨制备方向一致；若生态位未冻结则不得开展。'),
        ('P2', '条件TYMP异质性', '仅在结局盲生态位冻结后，检验低氧、血管生成、POSTN/ECM是否修饰TYMP—SPP1空间关系；预先限定交互与同方向门槛。', '只允许“情境依赖”结论；无论结果如何都不能逆转总体TYMP→SPP1阴性。'),
    ])
    add_para_callout(vdoc, '空间线优先输出', '生态位全景图、跨制备转移矩阵、限定LR空间图、TLS/免疫距离梯度、切片级伪bulk通路表，以及一张标注证据等级的空间机制主图。', fill='EAF4EA', accent='447A4B')

    add_page_break(vdoc)
    add_heading(vdoc, '三、路线B：F1/TYMP → MTA → SPP1 → 预后与治疗', 1)
    add_body(vdoc, '这条路线目前不能称为临床转化链。既有队列已经充分检验了通用F1/TYMP、SPP1及组合代理，但没有形成跨队列稳健的响应、PFS、OS或治疗交互证据。下一步只有在空间线先冻结一个可重复表型后，才值得开展新的临床投射。')
    add_heading(vdoc, '3.1 现有证据', 2)
    evidence_table([
        ('单细胞治疗语境', 'GSE264586/GSE294109的纵向、响应分层及桥接分析未形成稳健F1/TYMP临床标志；最终CNV门控后的响应/纵向比较均未通过预设FDR。', '当前数据不支持', '如重用，只能投射已冻结的新表型，不能重新定义SPP1或寻找最佳阈值。'),
        ('三队列临床复合分数', 'E-MTAB-3218、Braun CheckMate和JAVELIN中，预注册的7个R2 Myeloid/DC—TYMP主要模型为0/7通过FDR。', '无临床等级升级', '无需重复通用复合分数；新分析必须来自独立冻结的空间表型。'),
        ('校正后的ICI生存', 'CheckMate n=311与JAVELIN n=726采用正确事件编码后，FULL_CHAIN_PROXY与TYMP—SPP1联合PFS均无稳健支持；POSTN-CAF仅在JAVELIN出现，未在CheckMate复制。', 'PFS链未验证', '任何新Cox模型必须锁定事件编码、协变量、终点和阈值。'),
        ('SPP1细状态投射', '同一13状态参考投射至两队列后，SPP1分数、TYMP×SPP1、OS及治疗交互无支持或方向冲突；Other Macrophage在两队列100%为零。', '细状态不可可靠识别', '不得继续把bulk去卷积权重当SPP1细状态真值；优先使用固定基因程序而非三状态比例。'),
        ('广义巨噬生态', 'TYMP与原始SPP1权重及总巨噬细胞权重在两队列同向，但与巨噬细胞内SPP1比例和SPP1/C1QC比无关。', '仅支持广义共扩增', '可作为协变量/背景，不得写成选择性SPP1极化或疗效预测。'),
        ('TCGA下游受体/生态轴', '部分受体联合分组与OS/PFS有关，但连续TYMP×受体交互多数无支持，不能证明TYMP—MTA—SPP1临床链。', '下游情境线索', '只有与冻结空间表型同源且跨ICI队列重复，才可提升为候选转化证据。'),
    ])
    add_para_callout(vdoc, '临床线当前结论', '当前数据不支持F1/TYMP—MTA—SPP1作为ICI响应、PFS或OS预测标志。阴性结果是已完成的证据，不应通过新增阈值、重复分组或重新组合代理来“寻找阳性”。', fill='FCEBEC', accent='B94A48')

    add_heading(vdoc, '3.2 需要补全的临床断点', 2)
    for text in [
        '缺少可独立复现的空间表型：现有临床模型直接使用F1/TYMP或bulk SPP1代理，无法代表SPP1—内皮/TLS结构。',
        '缺少从空间表型到bulk/单细胞临床队列的冻结映射：评分基因、方向、阈值和协变量尚未在结局盲状态下固定。',
        '缺少跨队列重复：单队列预后或治疗交互不足以升级，必须在独立队列中同方向并按同一模型验证。',
        '缺少因果上游：即使空间表型与疗效相关，也不能反推TYMP通过MTA诱导SPP1。',
    ]:
        add_list_item(vdoc, text, bullets, after=5)

    add_page_break(vdoc)
    add_heading(vdoc, '3.3 可执行的预后/治疗补充分析', 2)
    plan_table([
        ('前置门', '固定空间表型', '只有路线A的生态位、LR/TLS或生态位伪bulk形成跨制备稳定表型后，才冻结候选基因集、方向、评分公式和排除规则。', '若空间线M110B-D失败，临床线停止，不重新使用旧F1/TYMP组合代理。'),
        ('P2-A', '单细胞临床投射', '在GSE264586/GSE294109及三队列桥中，仅用基线/预设时间点按患者汇总冻结表型；比较响应或纵向变化，研究分开建模，禁止细胞级 pooled P值。', '患者级方向在至少两个合格队列一致并通过各自FDR；小样本仅报告方向/效应。'),
        ('P2-B', 'bulk ICI响应', '在E-MTAB-3218、CheckMate和JAVELIN中结局盲计算固定程序分数；用预注册协变量的logistic模型检验响应，不优化cutoff。', '至少两个独立队列同方向，主要家族BH-FDR通过；单队列名义结果不升级。'),
        ('P2-C', 'PFS与OS', '在正确CNSR事件编码下，以连续固定分数为主，Cox模型调整年龄/性别/风险分层等可用协变量；队列内估计后再做方向和效应综合。', 'PFS或OS在至少两个队列同方向且FDR通过；PH假设失败则停止该模型。'),
        ('P2-D', '治疗交互', '仅在有随机或明确对照治疗的队列检验固定分数×治疗交互；主效应与交互分开形成FDR家族。', '交互在独立队列重复；否则只能称预后关联，不能称治疗预测。'),
        ('P3', '外部扩展', '仅在固定表型和明确验证角色后考虑PRJNA705464、SCP1288或其他候选；新增数据前先冻结纳入标准、样本单位和主要终点。', '不因方向不理想更换表型或新增队列；元数据不足则降为描述性。'),
    ])
    add_para_callout(vdoc, '临床线推荐输出', '一张跨队列效应森林图、队列级响应/PFS/OS主结果表、治疗交互专表、完整FDR表和失败/不可估计审计。只有跨队列重复后才进入客户主结论。', fill='F4F6F9', accent=DARK_BLUE)

    add_page_break(vdoc)
    add_heading(vdoc, '四、两条路线共用的因果补全模块', 1)
    add_body(vdoc, '计算分析能够深化下游空间结构并筛选临床候选，但不能关闭TYMP→MTA和MTA→SPP1。真正的因果闭环必须由扰动、直接代谢物测量和阻断/救援共同完成。')
    causal = vdoc.add_table(rows=6, cols=4)
    mark_table_header(causal.rows[0])
    set_table_geometry(causal, [1450, 2700, 3260, 1950], indent_dxa=120)
    set_table_borders(causal, color='D8DDE3', size=4)
    for j, text in enumerate(['实验环节', '最低设计', '主要读出', '判定作用']):
        set_cell_shading(causal.cell(0, j), PALE_BLUE)
        add_table_text(causal.cell(0, j), text, bold=True, color=DARK_BLUE, size=8.8, align=WD_ALIGN_PARAGRAPH.CENTER)
    causal_rows = [
        ('TYMP扰动', '敲低/过表达或可验证药理干预；匹配对照；独立生物学重复', 'TYMP RNA/蛋白、细胞状态与活性核查', '确认上游暴露真实发生'),
        ('直接MTA', '相同样本和时间点的靶向LC-MS/MS；保留批次与QC', '胞内/上清5′-MTA绝对或标准化定量', '直接检验TYMP→MTA'),
        ('SPP1读出', '肿瘤细胞条件培养/巨噬细胞体系或共培养；预设时间点', 'SPP1 RNA、蛋白及SPP1-state程序', '直接检验MTA→SPP1'),
        ('通路阻断', 'A2B或STAT3阻断；必要时加入外源MTA救援', 'SPP1与下游血管/免疫功能读出', '区分相容性与机制依赖'),
        ('空间/临床接口', '用计算线冻结的生态位程序，不用实验结果反向挑选签名', '空间定位、固定程序分数与临床端点', '把因果实验与观察队列连接'),
    ]
    for i, row in enumerate(causal_rows, 1):
        for j, text in enumerate(row):
            add_table_text(causal.cell(i, j), text, size=8.15)
            set_cell_margins(causal.cell(i, j), top=55, bottom=55, start=105, end=105)
            if i % 2 == 0:
                set_cell_shading(causal.cell(i, j), 'FAFBFC')
    add_para_callout(vdoc, '样本量原则', '统计样本量按独立生物学重复计算；技术重复、细胞数、spot数和切片内邻域不能代替生物学样本量。', fill='FFF9EB', accent=GOLD)

    add_heading(vdoc, '五、推荐实施顺序：先空间，后临床，因果实验并行', 1)
    sequence = vdoc.add_table(rows=5, cols=4)
    mark_table_header(sequence.rows[0])
    set_table_geometry(sequence, [1100, 2560, 3650, 2050], indent_dxa=120)
    set_table_borders(sequence, color='D8DDE3', size=4)
    for j, text in enumerate(['顺序', '先回答的问题', '工作包', '进入下一步的门槛']):
        set_cell_shading(sequence.cell(0, j), LIGHT_GRAY)
        add_table_text(sequence.cell(0, j), text, bold=True, color=DARK_BLUE, size=8.8, align=WD_ALIGN_PARAGRAPH.CENTER)
    seq_rows = [
        ('1', 'SPP1—内皮是否形成可转移生态位？', '24张切片输入核查、FFPE发现、frozen转移', '跨制备稳定或明确制备依赖'),
        ('2', '生态位是否具有LR、TLS和转录机制深度？', '限定LR、TLS/免疫梯度、切片级伪bulk', '至少一个冻结且可复核的空间表型'),
        ('3', '该固定表型是否具有临床意义？', '单细胞临床投射、bulk响应/PFS/OS、治疗交互', '至少两个独立队列同方向并通过FDR'),
        ('并行', 'TYMP—MTA—SPP1能否形成因果闭环？', '扰动、直接MTA、SPP1 RNA/蛋白、阻断/救援', '直接效应与机制依赖均满足'),
    ]
    for i, row in enumerate(seq_rows, 1):
        for j, text in enumerate(row):
            add_table_text(sequence.cell(i, j), text, size=8.25, align=WD_ALIGN_PARAGRAPH.CENTER if j == 0 else WD_ALIGN_PARAGRAPH.LEFT)
            set_cell_margins(sequence.cell(i, j), top=55, bottom=55, start=105, end=105)
            if i % 2 == 0:
                set_cell_shading(sequence.cell(i, j), 'FAFBFC')

    add_page_break(vdoc)
    add_heading(vdoc, '六、统计、质量控制与报告口径', 1)
    for text in [
        '数据集分开：原始单细胞、免疫治疗单细胞、GSE175540、TJ-RCC、UnitedMet、NCI60及bulk ICI队列保持独立；只在预设综合层连接，不混入同一校正池。',
        '独立单位：单细胞以患者、空间以独立切片、bulk以患者、细胞系以细胞系为统计单位。',
        '小样本：独立n<10只报告方向、效应和一致性；名义显著不得写成FDR显著。',
        '多重检验：基因、通路、LR、空间梯度、响应、生存和治疗交互分别形成预定义FDR家族。',
        '结局盲冻结：生态位、程序基因、评分公式、阈值和排除规则必须在查看临床结局前固定。',
        '阴性保留：上游TYMP→MTA、MTA→SPP1及既有临床阴性必须在报告链条原位保留，下游阳性不能救援上游断点。',
        '失败审计：输入不完整、模型不可识别、PH假设失败、制备不重复和元数据不足均形成正式失败/不可估计记录。',
    ]:
        add_list_item(vdoc, text, bullets, bold_prefix=text.split('：', 1)[0] + '：', after=5)

    add_heading(vdoc, '七、修订后的结论框架', 1)
    add_para_callout(vdoc, '现在可以明确说', 'F1与TYMP在原始单细胞肿瘤状态中相关但存在异质性；SPP1巨噬细胞状态已被操作性定义；SPP1—内皮双向通讯相容性在14名患者中可重复；GSE175540的SPP1—血管空间关联在FFPE和frozen中均获得支持。', fill='EAF4EA', accent='447A4B')
    add_para_callout(vdoc, '现在不能说', '现有数据不支持TYMP驱动MTA、MTA直接诱导SPP1、TYMP普遍极化SPP1巨噬细胞，或完整链条能够预测ICI响应、PFS、OS及治疗获益。', fill='FCEBEC', accent='B94A48')
    add_body(vdoc, '因此，下一阶段不是继续扩张同类相关性，而是沿空间线把已支持的下游界面做深、冻结一个可转移表型；随后再沿临床线进行严格的独立验证。同时由扰动实验补齐中央因果桥。')

    # BUILD_V2_CONTENT

    vdoc.core_properties.title = '项目26210双路线证据与补充分析计划'
    vdoc.core_properties.subject = 'F1/TYMP-MTA-SPP1空间结构与预后治疗两条证据链'
    vdoc.core_properties.author = '项目分析团队'
    vdoc.core_properties.keywords = '26210, F1, TYMP, MTA, SPP1, 空间结构, 预后, 治疗'
    vdoc.save(DOCX_PATH)


build_v2()
print(DOCX_PATH)


def make_seven_feature_workflow(path: Path):
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    chosen = 'Microsoft YaHei' if 'Microsoft YaHei' in fonts else ('SimHei' if 'SimHei' in fonts else 'DejaVu Sans')
    plt.rcParams['font.family'] = chosen
    plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(figsize=(12.4, 5.0))
    ax.set_xlim(0, 12.4); ax.set_ylim(0, 5.0); ax.axis('off')
    boxes = [
        (0.2, 2.55, 2.05, 1.25, '客户图特征注册', '冻结7项方向\n不估计权重'),
        (2.7, 2.55, 2.05, 1.25, '公共数据覆盖核对', 'RC18/RC20完整7项\n缩减队列单独标记'),
        (5.2, 2.55, 2.05, 1.25, '代谢特征迁移检验', '队列内评分\n跨队列方向与FDR'),
        (7.7, 2.55, 2.05, 1.25, 'RNA代理构建', '仅前一步通过\n双向跨队列验证'),
        (10.2, 2.55, 2.0, 1.25, '下游投射', '单细胞→空间→临床\n逐级设置停止门'),
    ]
    for i,(x,y,w,h,title,sub) in enumerate(boxes):
        face = '#E8EEF5' if i < 3 else '#EEF4EA'
        edge = '#2E74B5' if i < 3 else '#447A4B'
        ax.add_patch(plt.Rectangle((x,y),w,h,facecolor=face,edgecolor=edge,linewidth=1.6))
        ax.text(x+w/2,y+h*0.67,title,ha='center',va='center',fontsize=10,fontweight='bold',color='#17365D')
        ax.text(x+w/2,y+h*0.28,sub,ha='center',va='center',fontsize=8.6,color='#3C4650',linespacing=1.25)
        if i < len(boxes)-1:
            ax.annotate('',xy=(boxes[i+1][0]-0.08,y+h/2),xytext=(x+w+0.08,y+h/2),arrowprops=dict(arrowstyle='->',lw=1.7,color='#7E8994'))
    ax.text(0.3,1.55,'核心门槛',fontsize=10.7,fontweight='bold',color='#A97800')
    ax.text(1.55,1.55,'RC18与RC20完整七特征方向一致并达到预设统计门槛，才进入RNA代理构建。',fontsize=9.5,color='#3C4650')
    ax.plot([0.3,12.05],[1.18,1.18],color='#D5DBE2',lw=1)
    ax.text(0.3,0.72,'失败结果同样有效：任何环节未通过即保留“不可迁移/情境依赖”结论，不用后续数据补成阳性。',fontsize=9.5,color='#9B1C1C',fontweight='bold')
    fig.savefig(path,dpi=220,bbox_inches='tight',facecolor='white'); plt.close(fig)


def add_compact_num_definition(doc, num_fmt='bullet', glyph='•'):
    numbering = doc.part.numbering_part.element
    aids = [int(x.get(qn('w:abstractNumId'))) for x in numbering.findall(qn('w:abstractNum'))]
    nids = [int(x.get(qn('w:numId'))) for x in numbering.findall(qn('w:num'))]
    aid, nid = max(aids+[0])+1, max(nids+[0])+1
    abstract=OxmlElement('w:abstractNum'); abstract.set(qn('w:abstractNumId'),str(aid))
    lvl=OxmlElement('w:lvl'); lvl.set(qn('w:ilvl'),'0')
    for tag,val in [('w:start','1'),('w:numFmt',num_fmt),('w:lvlText',glyph if num_fmt=='bullet' else '%1.'),('w:lvlJc','left')]:
        node=OxmlElement(tag); node.set(qn('w:val'),val); lvl.append(node)
    ppr=OxmlElement('w:pPr'); tabs=OxmlElement('w:tabs'); tab=OxmlElement('w:tab'); tab.set(qn('w:val'),'num'); tab.set(qn('w:pos'),'540'); tabs.append(tab); ppr.append(tabs)
    ind=OxmlElement('w:ind'); ind.set(qn('w:left'),'540'); ind.set(qn('w:hanging'),'270'); ppr.append(ind); lvl.append(ppr); abstract.append(lvl); numbering.append(abstract)
    num=OxmlElement('w:num'); num.set(qn('w:numId'),str(nid)); ref=OxmlElement('w:abstractNumId'); ref.set(qn('w:val'),str(aid)); num.append(ref); numbering.append(num)
    return nid


def build_client_plan():
    make_seven_feature_workflow(FIG_PATH)
    pdoc = Document()
    sec = pdoc.sections[0]
    sec.page_width=Inches(8.5); sec.page_height=Inches(11)
    sec.top_margin=Inches(1); sec.bottom_margin=Inches(1); sec.left_margin=Inches(1); sec.right_margin=Inches(1)
    sec.header_distance=Inches(0.492); sec.footer_distance=Inches(0.492)
    normal=pdoc.styles['Normal']; normal.font.name='Calibri'; normal._element.rPr.rFonts.set(qn('w:ascii'),'Calibri'); normal._element.rPr.rFonts.set(qn('w:hAnsi'),'Calibri'); normal._element.rPr.rFonts.set(qn('w:eastAsia'),'Microsoft YaHei'); normal.font.size=Pt(11); normal.font.color.rgb=rgb(BLACK); normal.paragraph_format.space_before=Pt(0); normal.paragraph_format.space_after=Pt(6); normal.paragraph_format.line_spacing=1.25
    for name,size,color,before,after in [('Heading 1',16,BLUE,18,10),('Heading 2',13,BLUE,14,7),('Heading 3',12,DARK_BLUE,10,5)]:
        s=pdoc.styles[name]; s.font.name='Calibri'; s._element.rPr.rFonts.set(qn('w:ascii'),'Calibri'); s._element.rPr.rFonts.set(qn('w:hAnsi'),'Calibri'); s._element.rPr.rFonts.set(qn('w:eastAsia'),'Microsoft YaHei'); s.font.size=Pt(size); s.font.bold=True; s.font.color.rgb=rgb(color); s.paragraph_format.space_before=Pt(before); s.paragraph_format.space_after=Pt(after); s.paragraph_format.keep_with_next=True
    hp=sec.header.paragraphs[0]; hp.alignment=WD_ALIGN_PARAGRAPH.LEFT; hr=hp.add_run('项目26210  |  无权重七特征公共数据迁移分析方案'); set_run_font(hr,size=8.8,color=MUTED,bold=True)
    add_page_number(sec.footer.paragraphs[0])
    bullets=add_compact_num_definition(pdoc,'bullet','•'); decimals=add_compact_num_definition(pdoc,'decimal')

    def cli(text, num_id=bullets, bold_prefix=None, after=4):
        p=pdoc.add_paragraph(); p.paragraph_format.space_before=Pt(0); p.paragraph_format.space_after=Pt(after); p.paragraph_format.line_spacing=1.25
        ppr=p._p.get_or_add_pPr(); numpr=OxmlElement('w:numPr'); ilvl=OxmlElement('w:ilvl'); ilvl.set(qn('w:val'),'0'); nid=OxmlElement('w:numId'); nid.set(qn('w:val'),str(num_id)); numpr.append(ilvl); numpr.append(nid); ppr.append(numpr)
        if bold_prefix and text.startswith(bold_prefix):
            r1=p.add_run(bold_prefix); set_run_font(r1,size=11,color=BLACK,bold=True); r2=p.add_run(text[len(bold_prefix):]); set_run_font(r2,size=11,color=BLACK)
        else:
            r=p.add_run(text); set_run_font(r,size=11,color=BLACK)
        return p

    def body(text, bold_prefix=None, after=6, italic=False):
        p=pdoc.add_paragraph(); p.paragraph_format.space_before=Pt(0); p.paragraph_format.space_after=Pt(after); p.paragraph_format.line_spacing=1.25
        if bold_prefix and text.startswith(bold_prefix):
            r1=p.add_run(bold_prefix); set_run_font(r1,size=11,color=BLACK,bold=True); r2=p.add_run(text[len(bold_prefix):]); set_run_font(r2,size=11,color=BLACK,italic=italic)
        else:
            r=p.add_run(text); set_run_font(r,size=11,color=BLACK,italic=italic)
        return p

    pdoc.add_paragraph().paragraph_format.space_after=Pt(18)
    p=pdoc.add_paragraph(); p.paragraph_format.space_after=Pt(4); r=p.add_run('CLIENT ANALYSIS PLAN'); set_run_font(r,size=9.2,color=GOLD,bold=True)
    p=pdoc.add_paragraph(); p.paragraph_format.space_after=Pt(7); r=p.add_run('项目26210无权重七特征\n公共数据迁移分析方案'); set_run_font(r,size=27,color=NAVY,bold=True)
    p=pdoc.add_paragraph(); p.paragraph_format.space_after=Pt(14); r=p.add_run('从786-O shTYMP代谢扰动特征出发，检验公共数据中的可迁移性与适用边界'); set_run_font(r,size=13.2,color=MID_BLUE,bold=True)
    body('文件用途：客户确认分析范围、步骤、统计口径与停止规则',bold_prefix='文件用途：',after=2)
    body('当前状态：仅完成方案设计；尚未查看七特征联合分析结果',bold_prefix='当前状态：',after=2)
    body('版本日期：V1.0  |  2026年9月1日',bold_prefix='版本日期：',after=10)
    add_para_callout(pdoc,'核心目标','将客户图中已经观察到的shTYMP代谢变化冻结为一个不使用效应权重的七特征方向分数；先在具有配对RNA和代谢组的公共肾癌队列中检验迁移性。公共数据分析用于判断可重复性和适用边界，不替代客户实验，也不修改既有单一MTA阴性结果。',fill='EEF4F8',accent=DARK_BLUE)
    pic=pdoc.add_picture(str(FIG_PATH),width=Inches(6.45)); pic._inline.docPr.set('title','无权重七特征公共数据迁移分析流程'); pic._inline.docPr.set('descr','从客户图特征注册、公共数据覆盖、代谢迁移检验、条件性RNA代理到下游投射的逐级流程，并设置停止门'); pdoc.paragraphs[-1].alignment=WD_ALIGN_PARAGRAPH.CENTER
    cp=pdoc.add_paragraph(); cp.alignment=WD_ALIGN_PARAGRAPH.CENTER; cp.paragraph_format.space_after=Pt(0); rr=cp.add_run('图1  分析步骤、条件门槛与停止规则'); set_run_font(rr,size=8.8,color=MUTED,italic=True)

    add_page_break(pdoc)
    add_heading(pdoc,'一、分析问题与证据边界',1)
    body('本分析不试图利用公共数据库“证明客户实验正确”，而是回答客户图中观察到的代谢变化能否在独立公共样本中表现为一致的分子状态。客户扰动实验、公共观察队列和后续单细胞/空间投射属于不同证据层，分别报告。')
    add_heading(pdoc,'1.1 需要回答的三个问题',2)
    for t in [
        '客户图定义的七项代谢物方向，能否在两个具有完整覆盖的公共肾癌队列中形成一致的样本级综合分数？',
        '该综合分数与TYMP RNA的关联方向是否符合shTYMP图示预期，并在两个队列中重复？',
        '只有上述迁移成立后，能否从配对RNA中建立一个可跨队列转移的转录组代理，用于后续单细胞、空间和临床数据？',
    ]: cli(t,decimals)
    add_heading(pdoc,'1.2 允许与不允许的结论',2)
    boundary=pdoc.add_table(rows=4,cols=2); mark_table_header(boundary.rows[0]); set_table_geometry(boundary,[4680,4680],indent_dxa=120); set_table_borders(boundary,color='D8DDE3',size=4)
    for j,t in enumerate(['允许的结论','不允许的结论']): set_cell_shading(boundary.cell(0,j),PALE_BLUE); add_table_text(boundary.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.5,align=WD_ALIGN_PARAGRAPH.CENTER)
    pairs=[
        ('七特征在公共肾癌样本中方向一致、情境依赖或不支持。','公共数据库复现了786-O shTYMP因果效应。'),
        ('某些代谢物可迁移而另一些不可迁移。','组合分数可以覆盖或推翻既有TYMP-MTA单点阴性。'),
        ('经独立队列验证的RNA代理可用于后续投射。','RNA代理等同于直接代谢物测量或机制闭环。'),
    ]
    for i,row in enumerate(pairs,1):
        for j,t in enumerate(row): add_table_text(boundary.cell(i,j),t,size=9.0); set_cell_margins(boundary.cell(i,j),top=70,bottom=70,start=110,end=110)
    add_para_callout(pdoc,'现有材料边界','目前客户侧仅有图像，因此七项特征只冻结方向，不从箱线图、火山图或VIP图估算权重、效应量或重复数。若后续获得原始矩阵，可另行开展独立的客户数据复算，但不改变本方案预先定义的方向。',fill='FFF9EB',accent=GOLD)

    add_page_break(pdoc)
    add_heading(pdoc,'二、七特征注册与评分定义',1)
    body('七项代谢物及方向在任何联合结果计算前冻结。综合分数统一定义为“shTYMP样代谢分数”；分数越高，表示样本的代谢状态越接近客户图中的shTYMP方向。')
    reg=pdoc.add_table(rows=8,cols=4); mark_table_header(reg.rows[0]); set_table_geometry(reg,[2800,1700,1500,3360],indent_dxa=120); set_table_borders(reg,color='D8DDE3',size=4)
    for j,t in enumerate(['代谢物','shTYMP方向','分数符号','注册说明']): set_cell_shading(reg.cell(0,j),PALE_BLUE); add_table_text(reg.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.1,align=WD_ALIGN_PARAGRAPH.CENTER)
    feature_rows=[
        ('5-methylthioadenosine (5-MTA)','下降','-1','客户靶向验证面板；保留既有单点MTA分析作为独立证据。'),
        ('Thymidine','上升','+1','TYMP经典底物；客户靶向面板中的核心酶学特征。'),
        ('2′-Deoxyuridine','上升','+1','来自客户差异特征图的预设候选；不从图片估算效应量。'),
        ('Methionine','下降','-1','甲硫氨酸循环共响应。'),
        ('S-adenosylmethionine (SAM)','下降','-1','甲硫氨酸循环共响应。'),
        ('Thymine','下降','-1','嘧啶代谢共响应。'),
        ('Spermidine','下降','-1','多胺代谢共响应。'),
    ]
    for i,row in enumerate(feature_rows,1):
        for j,t in enumerate(row): add_table_text(reg.cell(i,j),t,size=8.7,align=WD_ALIGN_PARAGRAPH.CENTER if j in (1,2) else WD_ALIGN_PARAGRAPH.LEFT); set_cell_margins(reg.cell(i,j),top=55,bottom=55,start=105,end=105)
        if i%2==0:
            for c in reg.rows[i].cells: set_cell_shading(c,'FAFBFC')
    add_heading(pdoc,'2.1 无权重评分公式',2)
    add_para_callout(pdoc,'shTYMP样分数','S = [ -Z(5-MTA) + Z(Thymidine) + Z(Deoxyuridine) - Z(Methionine) - Z(SAM) - Z(Thymine) - Z(Spermidine) ] / 7。Z为每个队列内部的秩标准化值；七项权重绝对值均为1，不根据当前数据重新估计。',fill='EEF4F8',accent=DARK_BLUE)
    for t in [
        '期望方向：若公共样本中的TYMP活性状态与客户shTYMP扰动相反，则TYMP RNA与shTYMP样分数预期为负相关。',
        '完整性要求：RC18和RC20只有在同一样本七项均可用时进入主要评分；不插补缺失值。',
        '缩减分数：少于七项的数据集单独命名并降级为敏感性分析，不能称为七特征验证。',
    ]: cli(t)

    add_heading(pdoc,'三、数据源与分析层级',1)
    data=pdoc.add_table(rows=7,cols=5); mark_table_header(data.rows[0]); set_table_geometry(data,[1400,1500,1350,1600,3510],indent_dxa=120); set_table_borders(data,color='D8DDE3',size=4)
    for j,t in enumerate(['数据集','主要样本','特征覆盖','分析层级','用途']): set_cell_shading(data.cell(0,j),LIGHT_GRAY); add_table_text(data.cell(0,j),t,bold=True,color=DARK_BLUE,size=8.8,align=WD_ALIGN_PARAGRAPH.CENTER)
    data_rows=[
        ('RC18','肿瘤n=97','7/7','主要','完整七特征队列1；配对RNA和代谢组。'),
        ('RC20','肿瘤n=52','7/7','主要','完整七特征队列2；独立重复。'),
        ('CPTAC','肿瘤n=50','3/7','敏感性','仅5-MTA、Methionine、Thymine缩减分数。'),
        ('CPTAC_val','肿瘤n=56','3/7','敏感性','独立缩减分数；不升级完整七特征结论。'),
        ('NCI60','全谱系n=58','5/7','外部敏感性','5-MTA、Deoxyuridine、Methionine、SAM、Spermidine。'),
        ('NCI60肾癌','n=8','5/7','探索性','只报告方向、效应和逐细胞系一致性，不做推断检验。'),
    ]
    for i,row in enumerate(data_rows,1):
        for j,t in enumerate(row): add_table_text(data.cell(i,j),t,size=8.25,align=WD_ALIGN_PARAGRAPH.CENTER if j in (1,2,3) else WD_ALIGN_PARAGRAPH.LEFT); set_cell_margins(data.cell(i,j),top=50,bottom=50,start=95,end=95)
        if i%2==0:
            for c in data.rows[i].cells: set_cell_shading(c,'FAFBFC')
    body('UnitedMet来源：https://zenodo.org/records/11286535。四个队列继续使用已验证的样本映射和肿瘤/正常标记；正常样本不进入主要分析。',after=4,italic=True)
    add_para_callout(pdoc,'排除的数据源','CCL180和CCLE缺少可接受的完整MTA及七特征覆盖，不进入本分析。近质量离子、相邻代谢物和名称不明确的峰不作替代。',fill='FFF9EB',accent=GOLD)

    add_page_break(pdoc)
    add_heading(pdoc,'四、立即执行模块：逐步分析内容',1)
    body('本节为第一轮实际工作。所有规则在联合结果计算前固定；执行过程中不根据显著性增删代谢物、改变方向或重新定义主要队列。')

    add_heading(pdoc,'步骤0｜客户特征与输入版本冻结',2)
    body('输入：客户提供的代谢组图、UnitedMet四个配对RNA/代谢文件及样本映射表、NCI60代谢矩阵和已验证TYMP/MTAP RNA映射。',bold_prefix='输入：')
    body('具体分析：记录客户图文件哈希、七项名称和方向；记录所有公共文件的大小、MD5/SHA256、列名、样本数和来源。该步骤只做身份与结构核查，不计算相关性。',bold_prefix='具体分析：')
    body('输出：客户七特征注册表、公共数据文件清单、样本映射审计和分析参数JSON。',bold_prefix='输出：')
    body('进入下一步的门槛：客户图方向已唯一登记；所有主要公共文件哈希匹配既有记录；不存在同名不同物或重复样本ID。',bold_prefix='进入下一步的门槛：')

    add_heading(pdoc,'步骤1｜代谢物身份和跨平台名称映射',2)
    body('输入：各数据集代谢物列名、平台注释及可用HMDB/KEGG标识。',bold_prefix='输入：')
    body('具体分析：为七项代谢物建立规范名、平台原始名、同义名和唯一标识映射；逐队列核对精确身份。只接受明确的5-MTA、Thymidine、2′-Deoxyuridine、Methionine、SAM、Thymine和Spermidine。',bold_prefix='具体分析：')
    for t in ['不使用近质量离子替代5-MTA。','不使用Putrescine或N-acetylputrescine替代Spermidine。','不使用Methionine sulfoxide等衍生物替代Methionine。']: cli(t)
    body('输出：逐数据集七特征覆盖矩阵和排除原因表。',bold_prefix='输出：')
    body('门槛：RC18与RC20必须各自精确覆盖7/7；否则完整七特征主要分析停止，数据集降级为缩减敏感性。',bold_prefix='门槛：')

    add_heading(pdoc,'步骤2｜样本配对、肿瘤筛选与统计单位确认',2)
    body('具体分析：按照现有MasterMapping逐一匹配代谢组ID与RNA ID，检查一一对应、重复ID、缺失样本和肿瘤/正常标签；四个UnitedMet队列分别处理，不合并原始样本。',bold_prefix='具体分析：')
    for t in [
        '主要分析仅保留肿瘤样本：RC18 n=97、RC20 n=52。',
        'CPTAC n=50和CPTAC_val n=56仅用于三特征缩减敏感性。',
        '每个肿瘤样本或细胞系计为一个独立生物学单位；代谢物数量、技术重复数和转录本数量不增加统计n。',
    ]: cli(t)
    body('输出：每队列纳入/排除流程、完整样本表和最终n。',bold_prefix='输出：')
    body('门槛：样本总数、肿瘤/正常数和配对ID须与既有下载验证完全一致；不一致时停止正式计算。',bold_prefix='门槛：')

    add_page_break(pdoc)
    add_heading(pdoc,'步骤3｜队列内标准化与无权重评分',2)
    body('具体分析：每个队列单独进行，不在队列之间做联合标准化。对每项代谢物在合格肿瘤样本中计算平均秩，再转换为均值0、标准差1的秩标准化值。',bold_prefix='具体分析：')
    for t in [
        '按注册符号计算七项等权平均，得到每个样本的shTYMP样分数。',
        '七项必须同时有值；不使用均值、中位数、最近邻或模型插补。',
        'TYMP和MTAP RNA使用已验证的配对表达，转换为log2(TPM+1)；不因结果重新选择表达口径。',
        '评分前检查每项有限性、方差、极端秩和样本顺序；评分后检查均值、标准差和取值范围。',
    ]: cli(t)
    body('输出：逐样本七项秩值、七项有符号贡献和最终shTYMP样分数。',bold_prefix='输出：')
    body('门槛：七项均有非零方差、完整样本数与步骤2一致、评分无NA/Inf；否则停止该队列主要分析。',bold_prefix='门槛：')

    add_heading(pdoc,'步骤4｜RC18与RC20完整七特征主要检验',2)
    body('主要终点：在RC18和RC20中分别计算TYMP RNA与shTYMP样分数的Spearman相关。根据客户图方向，预期rho<0；检验保持双侧。',bold_prefix='主要终点：')
    body('统计证据：报告rho、95%患者重采样Bootstrap区间、双侧P值，并对两个队列的主要终点形成一个两检验BH-FDR家族。',bold_prefix='统计证据：')
    body('组分分析：同一队列内分别计算TYMP与7项代谢物的Spearman相关；每个队列的7项构成独立BH-FDR家族。报告预期方向命中数，不用单项结果修改综合分数。',bold_prefix='组分分析：')
    body('MTAP敏感性：在秩空间中残差化MTAP后计算TYMP与综合分数的偏相关；RC18与RC20两项构成独立敏感性FDR家族。',bold_prefix='MTAP敏感性：')
    body('既有结果保留：原TYMP–MTA单项结果原样并列，不重新命名、删除或由综合分数替代。',bold_prefix='既有结果保留：')
    body('输出：两队列主要结果表、7项组分表、MTAP调整表、样本散点/秩图和贡献热图。',bold_prefix='输出：')

    add_page_break(pdoc)
    add_heading(pdoc,'步骤5｜跨队列综合与正式证据分级',2)
    body('具体分析：将RC18和RC20的主要相关系数进行Fisher z转换，分别报告固定效应和随机效应汇总；同时报告Cochran Q、I²和tau²。由于只有两个完整队列，队列内重复和方向一致性优先于单一汇总P值。',bold_prefix='具体分析：')
    grade=pdoc.add_table(rows=5,cols=3); mark_table_header(grade.rows[0]); set_table_geometry(grade,[2250,4200,2910],indent_dxa=120); set_table_borders(grade,color='D8DDE3',size=4)
    for j,t in enumerate(['证据等级','预设判定条件','报告措辞']): set_cell_shading(grade.cell(0,j),PALE_BLUE); add_table_text(grade.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.0,align=WD_ALIGN_PARAGRAPH.CENTER)
    grade_rows=[
        ('完整迁移支持','RC18与RC20 rho均<0，且两项主要BH q均<0.05。','七特征在两个完整公共肾癌队列中获得支持。'),
        ('方向一致的部分支持','两队列rho均<0，但未同时通过主要FDR；固定效应汇总支持。','方向一致，但证据未达到双队列正式门槛。'),
        ('探索性一致','两队列rho均<0，但队列和汇总检验均未达到门槛。','仅保留方向性候选线索。'),
        ('当前数据不支持','方向相反、主要汇总方向不符，或结果由单个队列驱动。','当前公共数据不支持七特征一般性迁移。'),
    ]
    for i,row in enumerate(grade_rows,1):
        for j,t in enumerate(row): add_table_text(grade.cell(i,j),t,size=8.45); set_cell_margins(grade.cell(i,j),top=55,bottom=55,start=100,end=100)
        if i%2==0:
            for c in grade.rows[i].cells: set_cell_shading(c,'FAFBFC')
    add_para_callout(pdoc,'解释原则','即使综合分数获得支持，也只说明客户扰动方向在公共样本中存在可迁移的代谢状态；不证明这些代谢物由TYMP直接共同调控，也不证明MTA介导SPP1。',fill='FFF9EB',accent=GOLD)

    add_heading(pdoc,'步骤6｜稳健性和影响诊断',2)
    body('所有稳健性分析预先限定，不作为提高显著性的替代模型。',after=4)
    robustness=[
        ('逐特征留一法','依次删除1项，形成7个六特征分数；检查RC18与RC20相关方向是否保持。','判断结论是否被单个代谢物支配。'),
        ('逐样本留一法','每次删除1个样本重算rho，记录范围和方向翻转。','识别高影响样本，不自动删除。'),
        ('置换检验','固定代谢分数，队列内随机置换TYMP标签100,000次。','验证相关性P值不依赖参数近似。'),
        ('原始Z分数敏感性','在同一队列内用常规Z分数替代秩Z，保持方向和等权规则。','仅检查秩变换选择的影响。'),
        ('特征贡献审计','计算每项对总分方差的贡献和相互相关矩阵。','识别冗余，不据此重新赋权。'),
    ]
    rt=pdoc.add_table(rows=len(robustness)+1,cols=3); mark_table_header(rt.rows[0]); set_table_geometry(rt,[2100,4100,3160],indent_dxa=120); set_table_borders(rt,color='D8DDE3',size=4)
    for j,t in enumerate(['分析','执行方法','回答的问题']): set_cell_shading(rt.cell(0,j),LIGHT_GRAY); add_table_text(rt.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.0,align=WD_ALIGN_PARAGRAPH.CENTER)
    for i,row in enumerate(robustness,1):
        for j,t in enumerate(row): add_table_text(rt.cell(i,j),t,size=8.45); set_cell_margins(rt.cell(i,j),top=55,bottom=55,start=100,end=100)

    add_page_break(pdoc)
    add_heading(pdoc,'步骤7｜不完整覆盖队列的缩减敏感性',2)
    body('缩减分析只回答部分特征是否方向相容，不能提升或否定完整七特征主要结论。',after=4)
    add_heading(pdoc,'7.1 CPTAC与CPTAC_val三特征分数',3)
    body('公式：S3 = [ -Z(5-MTA) - Z(Methionine) - Z(Thymine) ] / 3。两个队列分别计算TYMP相关、MTAP调整偏相关和三项组分方向；两个缩减队列的主要终点形成独立BH家族。',bold_prefix='公式：')
    body('报告：明确标注“3/7缩减敏感性”，不得使用“七特征复现”或“独立完整验证”。',bold_prefix='报告：')
    add_heading(pdoc,'7.2 NCI60五特征分数',3)
    body('公式：S5 = [ -Z(5-MTA) + Z(Deoxyuridine) - Z(Methionine) - Z(SAM) - Z(Spermidine) ] / 5。',bold_prefix='公式：')
    for t in [
        '全谱系层：以细胞系为独立单位，报告TYMP相关，并用组织内置换100,000次获得组织校正敏感性P值。',
        '肾癌层：独立n=8，只报告rho、8个细胞系值和leave-one-line-out方向范围；不计算P、FDR或置信区间。',
        '786-0单线：只展示其相对位置，不将一个细胞系解释为shTYMP复现。',
    ]: cli(t)
    body('停止规则：缩减分数无论阳性或阴性，都不能改变RC18/RC20完整七特征的证据等级。',bold_prefix='停止规则：')

    add_page_break(pdoc)
    add_heading(pdoc,'五、并行数据检索与条件执行模块',1)
    body('以下模块不与七特征主要分析同时自动启动。每个模块均有独立进入门槛；上一步失败时，保留失败结论并停止，不更换数据或放宽标准。')

    add_heading(pdoc,'步骤8｜公开扰动数据专项检索',2)
    body('检索目标：寻找与客户实验可比的TYMP遗传敲低/敲除、TYMP药理抑制或明确MTA暴露数据，而不是继续增加普通TYMP表达队列。',bold_prefix='检索目标：')
    body('检索平台：GEO/SRA、ArrayExpress/BioStudies、MetaboLights、Metabolomics Workbench、PRIDE及论文补充数据仓库。',bold_prefix='检索平台：')
    body('纳入条件：必须有可识别对照、独立生物学重复、清楚的细胞/组织背景、可下载原始或处理数据，并至少满足下列之一：测量七项中至少4项代谢物；提供TYMP扰动转录组/蛋白组；提供MTA直接暴露后的细胞响应。',bold_prefix='纳入条件：')
    body('排除条件：只有生存曲线、只有相关性、无对照、技术重复冒充生物学重复、代谢物身份不清或无可下载数据。',bold_prefix='排除条件：')
    body('输出：候选数据注册表，逐项记录模型、扰动、重复数、测量层、可下载性、与客户图的可比程度及是否进入正式分析。',bold_prefix='输出：')
    add_para_callout(pdoc,'当前已知线索','已有会议摘要提示TYMP抑制可导致Thymidine/Deoxyuridine积累并关联铁死亡和免疫治疗，但未发现公开原始数据入口；在获得可验证文件前只列为文献线索，不进入计算验证。',fill='F4F6F9',accent=DARK_BLUE)

    add_heading(pdoc,'步骤9｜将代谢分数转换为可迁移RNA代理（条件模块）',2)
    body('进入门槛：RC18和RC20完整七特征rho均<0，且至少达到“方向一致的部分支持”；若为探索性一致或不支持，本模块不启动。',bold_prefix='进入门槛：')
    body('预测目标：每个样本的七特征shTYMP样分数。预测变量来自配对RNA。',bold_prefix='预测目标：')
    for t in [
        '先冻结候选基因池：七项代谢物直接反应基因及嘧啶、甲硫氨酸循环和多胺代谢通路基因；不根据SPP1、预后或治疗结果选基因。',
        '方向一：RC18训练、RC20外部验证；方向二：RC20训练、RC18外部验证。禁止将两队列混合后随机切分。',
        '训练使用Elastic Net；alpha和lambda仅在训练队列嵌套交叉验证中选择。所有标准化参数由训练队列确定并原样应用到验证队列。',
        '验证报告预测值与真实代谢分数的Spearman rho、R²、校准斜率和置换P值；两个转移方向形成一个BH家族。',
        '只有两个转移方向均为正且均通过FDR，才拟合最终冻结RNA代理；否则停止下游投射。',
    ]: cli(t)
    body('输出：固定基因列表、系数、截距、标准化参数、训练日志、交叉验证轨迹和两个独立验证结果。',bold_prefix='输出：')
    add_para_callout(pdoc,'解释边界','RNA代理代表与七特征相关的转录状态，不等于直接代谢物测量，也不证明TYMP因果驱动七项代谢物。',fill='FFF9EB',accent=GOLD)

    add_page_break(pdoc)
    add_heading(pdoc,'步骤10｜RNA代理向单细胞和空间数据投射（条件模块）',2)
    body('进入门槛：步骤9的双向跨队列验证全部通过，且RNA代理公式完全冻结。',bold_prefix='进入门槛：')
    add_heading(pdoc,'10.1 单细胞肿瘤状态',3)
    for t in [
        '在每个单细胞研究中分别计算恶性/肿瘤样细胞RNA代理分数，不合并不同研究的细胞。',
        '以患者为统计单位，汇总患者内肿瘤细胞中位分数和高分细胞比例。',
        '预设检验4项：SPP1-state占巨噬细胞比例、SPP1/C1QC状态对比、恶性细胞→SPP1通讯特异性、SPP1→恶性细胞通讯特异性；4项构成同一FDR家族。',
        '细胞数量只用于估计患者汇总值，不作为推断性样本量。',
    ]: cli(t)
    add_heading(pdoc,'10.2 空间结构',3)
    for t in [
        '在24张空间切片中计算肿瘤区域RNA代理分数，FFPE和frozen分别分析。',
        '以切片为统计单位检验：肿瘤代理分数与邻近SPP1-state、SPP1—内皮生态位、血管重塑和TLS距离梯度。',
        'spot仅定义空间邻域；不把spot数作为n。两种制备均同方向时才称双制备支持。',
    ]: cli(t)
    body('停止规则：若单细胞和空间结果均不支持，则报告RNA代理未能连接SPP1/空间结构，不进入临床模块。',bold_prefix='停止规则：')

    add_heading(pdoc,'步骤11｜响应、PFS和OS投射（最终条件模块）',2)
    body('进入门槛：RNA代理跨UnitedMet验证通过，并在至少一个独立单细胞或空间数据集中形成预设方向支持。',bold_prefix='进入门槛：')
    for t in [
        '在免疫治疗队列中结局盲计算固定RNA代理，不重新选择基因、权重或cutoff。',
        '响应：患者级连续分数的预注册Logistic模型；响应家族单独FDR。',
        'PFS和OS：采用已校正事件编码的连续分数Cox模型，检查比例风险假设；PFS、OS分别形成FDR家族。',
        '治疗交互：仅在有明确对照治疗的队列检验分数×治疗交互；不得把预后主效应称为治疗预测。',
        '队列分开估计后再做效应综合；不优化分割点，不因单队列方向不理想更换模型。',
    ]: cli(t)
    body('停止规则：临床模块失败不反向否定客户代谢实验，只说明该冻结特征当前不具备稳健临床迁移性。',bold_prefix='停止规则：')

    add_page_break(pdoc)
    add_heading(pdoc,'六、交付内容、验证和决策节点',1)
    add_heading(pdoc,'6.1 第一轮七特征分析交付物',2)
    deliverables=[
        ('特征注册','七项规范名、图示方向、评分符号、来源文件和哈希。'),
        ('数据覆盖审计','各队列精确特征覆盖、样本配对、纳入/排除及缺失原因。'),
        ('样本级评分','RC18/RC20完整七特征分数及七项贡献；缩减队列单独输出。'),
        ('正式统计表','主要相关、组分相关、MTAP敏感性、跨队列综合、FDR和稳健性。'),
        ('结果图','两队列散点/秩图、七项方向热图、leave-one-feature-out稳定性和队列综合图。'),
        ('结论与失败审计','按预设等级给出支持、部分支持、探索性或不支持；保留所有不可估计与失败记录。'),
    ]
    dt=pdoc.add_table(rows=len(deliverables)+1,cols=2); mark_table_header(dt.rows[0]); set_table_geometry(dt,[2100,7260],indent_dxa=120); set_table_borders(dt,color='D8DDE3',size=4)
    for j,t in enumerate(['交付类别','具体内容']): set_cell_shading(dt.cell(0,j),PALE_BLUE); add_table_text(dt.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.2,align=WD_ALIGN_PARAGRAPH.CENTER)
    for i,row in enumerate(deliverables,1):
        for j,t in enumerate(row): add_table_text(dt.cell(i,j),t,size=8.8); set_cell_margins(dt.cell(i,j),top=60,bottom=60,start=110,end=110)
        if i%2==0:
            for c in dt.rows[i].cells: set_cell_shading(c,'FAFBFC')

    add_heading(pdoc,'6.2 独立验证项目',2)
    for t in [
        '重新读取原始文件并复算七项秩标准化值、符号贡献和综合分数。',
        '核对样本ID、肿瘤筛选、完整案例数量和逐行排序，确保无静默丢样本。',
        '独立复算Spearman、Bootstrap区间、BH-FDR、偏相关、置换P值和Meta分析。',
        '检查主要结论是否与原TYMP-MTA单项结果并列而非覆盖。',
        '检查缩减分数在文件名、表头、图题和正文中均明确标记为3/7或5/7。',
        '记录脚本、结果、图和报告的SHA256，生成唯一ID清单。',
    ]: cli(t)

    add_heading(pdoc,'6.3 分阶段决策节点',2)
    decisions=[
        ('节点A','RC18/RC20完整七特征完成','决定是否允许构建RNA代理。'),
        ('节点B','RNA代理双向跨队列验证完成','决定是否允许投射单细胞和空间。'),
        ('节点C','单细胞/空间预设端点完成','决定是否允许进入响应、PFS和OS。'),
        ('任一节点失败','保存阴性/异质/不可估计结果','停止后续链条，不新增同类模型寻找阳性。'),
    ]
    dd=pdoc.add_table(rows=len(decisions)+1,cols=3); mark_table_header(dd.rows[0]); set_table_geometry(dd,[1500,3760,4100],indent_dxa=120); set_table_borders(dd,color='D8DDE3',size=4)
    for j,t in enumerate(['决策点','完成内容','下一步']): set_cell_shading(dd.cell(0,j),LIGHT_GRAY); add_table_text(dd.cell(0,j),t,bold=True,color=DARK_BLUE,size=9.0,align=WD_ALIGN_PARAGRAPH.CENTER)
    for i,row in enumerate(decisions,1):
        for j,t in enumerate(row): add_table_text(dd.cell(i,j),t,size=8.55,align=WD_ALIGN_PARAGRAPH.CENTER if j==0 else WD_ALIGN_PARAGRAPH.LEFT); set_cell_margins(dd.cell(i,j),top=55,bottom=55,start=105,end=105)

    add_page_break(pdoc)
    add_heading(pdoc,'七、客户确认事项与最终执行口径',1)
    add_para_callout(pdoc,'本轮拟执行范围','客户确认后，第一轮只执行步骤0至7：完整七特征主要分析、跨队列综合、稳健性和缩减敏感性。步骤8为并行数据检索；步骤9至11均为条件模块，不在本轮自动启动。',fill='EEF4F8',accent=DARK_BLUE)
    for t in [
        '不要求客户新增湿实验。客户现有图仅用于冻结方向，不作为公共数据的替代样本。',
        '不重新分析或预测客户另外两张未提供的图。',
        '不因公共数据与客户图不一致而删除队列、调整权重或缩短特征列表。',
        '不将完整七特征支持写成TYMP→MTA→SPP1因果闭环；后续每一层仅报告对应数据能支持的结论。',
    ]: cli(t)
    add_heading(pdoc,'计划使用的主要数据来源',2)
    sources=[
        ('客户代谢组图','内部提供；文件SHA256 AABF51942981FBF67CF578259072D4CA9BA84F2EFCFDE2F3855C087C82FFF123','七项方向注册'),
        ('UnitedMet','https://zenodo.org/records/11286535','RC18/RC20完整七特征及配对RNA；CPTAC缩减敏感性'),
        ('NCI60 Metabolon/CellMiner','https://discover.nci.nih.gov/cellminer/','五特征跨细胞系敏感性'),
        ('GEO GSE249082','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249082','MTA处理成纤维细胞候选；仅在专项检索后决定是否使用'),
        ('TYMP抑制ccRCC会议线索','https://registration.uaa2025.org/user/form/view/2112?response=64:53994','方向性文献线索；无公开原始数据前不进入计算'),
    ]
    st=pdoc.add_table(rows=len(sources)+1,cols=3); mark_table_header(st.rows[0]); set_table_geometry(st,[2300,4000,3060],indent_dxa=120); set_table_borders(st,color='D8DDE3',size=4)
    for j,t in enumerate(['来源','地址/身份','本方案中的角色']): set_cell_shading(st.cell(0,j),PALE_BLUE); add_table_text(st.cell(0,j),t,bold=True,color=DARK_BLUE,size=8.9,align=WD_ALIGN_PARAGRAPH.CENTER)
    for i,row in enumerate(sources,1):
        for j,t in enumerate(row): add_table_text(st.cell(i,j),t,size=8.1); set_cell_margins(st.cell(i,j),top=50,bottom=50,start=90,end=90)
    add_para_callout(pdoc,'最终原则','先验证客户七特征能否迁移，再决定是否向下连接SPP1、空间结构和临床结局。分析失败同样形成正式结论；不使用后续数据把前一步失败“补成闭环”。',fill='FFF9EB',accent=GOLD)

    # CLIENT_PLAN_CONTENT
    pdoc.core_properties.title='项目26210无权重七特征公共数据迁移分析方案'
    pdoc.core_properties.subject='客户shTYMP代谢特征注册、公共数据迁移与条件性下游投射'
    pdoc.core_properties.author='项目分析团队'
    pdoc.core_properties.keywords='26210, TYMP, 5-MTA, thymidine, metabolomics, UnitedMet, NCI60'
    pdoc.save(DOCX_PATH)


build_client_plan()
print(DOCX_PATH)
