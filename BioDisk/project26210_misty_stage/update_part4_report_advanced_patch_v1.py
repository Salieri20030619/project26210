#!/usr/bin/env python3
"""Write the complete advanced spatial-analysis patch into the current client DOCX."""

from __future__ import annotations

import hashlib
import json
import struct
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
OUTPUT = STAGE / "project26210_part4_advanced_patch_v1_candidate.docx"
METADATA = STAGE / "project26210_part4_advanced_patch_v1_report_metadata.json"
EXPECTED_SOURCE_SHA256 = "06332668ac51b1513eed799ac802f747be937695275cd90846f9677b4ab1315a"

RBASE = ROOT / "results/part4_spatial"
BAYES_QC = RBASE / "advanced_patch_bayesspace_qc_v1/verification_summary_v1.json"
PROGRAM_QC = RBASE / "advanced_patch_spatial_programs_qc_v1/verification_summary_v1.json"
MISTY_QC = RBASE / "advanced_patch_misty_qc_v1/verification_summary_v1.json"
COMMOT_QC = RBASE / "advanced_patch_commot_cross_sample_v3/commot_v3_independent_qc_summary.json"
MISTY_TABLE = RBASE / "advanced_patch_misty_cross_sample_v3/cross_sample_misty_performance_v3.tsv"
COMMOT_TABLE = RBASE / "advanced_patch_commot_cross_sample_v3/commot_cross_sample_primary_hypothesis_support_v3.tsv"

IMAGES = [
    RBASE / "advanced_patch_bayesspace_glmf1_v1/GLMF1_bayesspace_domains.png",
    RBASE / "advanced_patch_bayesspace_glmf2_v1/GLMF2_bayesspace_domains.png",
    RBASE / "advanced_patch_spatial_programs_cross_sample_v2/cross_sample_module_overlap_v2.png",
    RBASE / "advanced_patch_spatial_programs_cross_sample_v2/cross_sample_gene_rank_concordance_v2.png",
    RBASE / "advanced_patch_misty_cross_sample_v3/cross_sample_misty_performance_gains_v3.png",
    RBASE / "advanced_patch_misty_cross_sample_v3/cross_sample_misty_importance_replication_v3.png",
    RBASE / "advanced_patch_commot_cross_sample_v3/commot_cross_sample_transport_rank_v3.png",
    RBASE / "advanced_patch_commot_cross_sample_v3/commot_cross_sample_hypothesis_proxy_v3.png",
]

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL = "http://schemas.openxmlformats.org/package/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
for prefix, uri in {"w": W, "r": R, "wp": WP, "a": A, "pic": PIC}.items():
    ET.register_namespace(prefix, uri)
ET.register_namespace("", REL)


def q(namespace: str, tag: str) -> str:
    return f"{{{namespace}}}{tag}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def node_text(node: ET.Element) -> str:
    return "".join(value.text or "" for value in node.iter(q(W, "t")))


def paragraph(text: str, style: str | None = None, bold: bool = False, shading: str | None = None) -> ET.Element:
    node = ET.Element(q(W, "p"))
    if style or shading:
        props = ET.SubElement(node, q(W, "pPr"))
        if style:
            ET.SubElement(props, q(W, "pStyle"), {q(W, "val"): style})
        if shading:
            ET.SubElement(props, q(W, "shd"), {q(W, "fill"): shading})
    run = ET.SubElement(node, q(W, "r"))
    if bold:
        rprops = ET.SubElement(run, q(W, "rPr"))
        ET.SubElement(rprops, q(W, "b"))
    ET.SubElement(run, q(W, "t")).text = text
    return node


def page_break() -> ET.Element:
    node = ET.Element(q(W, "p"))
    run = ET.SubElement(node, q(W, "r"))
    ET.SubElement(run, q(W, "br"), {q(W, "type"): "page"})
    return node


def set_paragraph_text(node: ET.Element, text: str, bold: bool = False, shading: str | None = None) -> None:
    props = node.find(q(W, "pPr"))
    for child in list(node):
        if child is not props:
            node.remove(child)
    if shading:
        if props is None:
            props = ET.Element(q(W, "pPr"))
            node.insert(0, props)
        old = props.find(q(W, "shd"))
        if old is not None:
            props.remove(old)
        ET.SubElement(props, q(W, "shd"), {q(W, "fill"): shading})
    run = ET.SubElement(node, q(W, "r"))
    if bold:
        rprops = ET.SubElement(run, q(W, "rPr"))
        ET.SubElement(rprops, q(W, "b"))
    ET.SubElement(run, q(W, "t")).text = text


def direct_paragraph(body: ET.Element, prefix: str, style: str | None = None) -> ET.Element:
    matches = []
    for node in list(body):
        if node.tag != q(W, "p") or not node_text(node).strip().startswith(prefix):
            continue
        if style is not None:
            pstyle = node.find(f"{q(W, 'pPr')}/{q(W, 'pStyle')}")
            if pstyle is None or pstyle.get(q(W, "val")) != style:
                continue
        matches.append(node)
    if len(matches) != 1:
        raise RuntimeError(f"Expected one direct paragraph for {prefix!r}, found {len(matches)}")
    return matches[0]


def any_paragraph(body: ET.Element, prefix: str) -> ET.Element:
    matches = [node for node in body.iter(q(W, "p")) if node_text(node).strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one paragraph for {prefix!r}, found {len(matches)}")
    return matches[0]


def insert_after(body: ET.Element, anchor: ET.Element, nodes: list[ET.Element]) -> None:
    index = list(body).index(anchor) + 1
    for node in nodes:
        body.insert(index, node)
        index += 1


def table(rows: list[list[str]], font_half_points: int = 16) -> ET.Element:
    node = ET.Element(q(W, "tbl"))
    props = ET.SubElement(node, q(W, "tblPr"))
    ET.SubElement(props, q(W, "tblW"), {q(W, "w"): "0", q(W, "type"): "auto"})
    borders = ET.SubElement(props, q(W, "tblBorders"))
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        ET.SubElement(borders, q(W, edge), {
            q(W, "val"): "single", q(W, "sz"): "6", q(W, "space"): "0", q(W, "color"): "808080"
        })
    for row_index, values in enumerate(rows):
        row = ET.SubElement(node, q(W, "tr"))
        row_props = ET.SubElement(row, q(W, "trPr"))
        ET.SubElement(row_props, q(W, "cantSplit"))
        if row_index == 0:
            ET.SubElement(row_props, q(W, "tblHeader"), {q(W, "val"): "true"})
        for value in values:
            cell = ET.SubElement(row, q(W, "tc"))
            cell_props = ET.SubElement(cell, q(W, "tcPr"))
            if row_index == 0:
                ET.SubElement(cell_props, q(W, "shd"), {q(W, "fill"): "D9EAF7"})
            p = ET.SubElement(cell, q(W, "p"))
            run = ET.SubElement(p, q(W, "r"))
            rprops = ET.SubElement(run, q(W, "rPr"))
            ET.SubElement(rprops, q(W, "sz"), {q(W, "val"): str(font_half_points)})
            ET.SubElement(rprops, q(W, "szCs"), {q(W, "val"): str(font_half_points)})
            if row_index == 0:
                ET.SubElement(rprops, q(W, "b"))
            ET.SubElement(run, q(W, "t")).text = str(value)
    return node


def png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()[:24]
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"Not a PNG: {path}")
    return struct.unpack(">II", data[16:24])


def image_paragraph(rid: str, name: str, width_px: int, height_px: int, docpr_id: int, width_emu: int = 6_250_000) -> ET.Element:
    cx = width_emu
    cy = round(cx * height_px / width_px)
    node = ET.Element(q(W, "p"))
    props = ET.SubElement(node, q(W, "pPr"))
    ET.SubElement(props, q(W, "jc"), {q(W, "val"): "center"})
    run = ET.SubElement(node, q(W, "r"))
    drawing = ET.SubElement(run, q(W, "drawing"))
    inline = ET.SubElement(drawing, q(WP, "inline"), {"distT": "0", "distB": "0", "distL": "0", "distR": "0"})
    ET.SubElement(inline, q(WP, "extent"), {"cx": str(cx), "cy": str(cy)})
    ET.SubElement(inline, q(WP, "docPr"), {"id": str(docpr_id), "name": name})
    frame = ET.SubElement(inline, q(WP, "cNvGraphicFramePr"))
    ET.SubElement(frame, q(A, "graphicFrameLocks"), {"noChangeAspect": "1"})
    graphic = ET.SubElement(inline, q(A, "graphic"))
    data = ET.SubElement(graphic, q(A, "graphicData"), {"uri": "http://schemas.openxmlformats.org/drawingml/2006/picture"})
    picture = ET.SubElement(data, q(PIC, "pic"))
    nonvisual = ET.SubElement(picture, q(PIC, "nvPicPr"))
    ET.SubElement(nonvisual, q(PIC, "cNvPr"), {"id": "0", "name": name})
    ET.SubElement(nonvisual, q(PIC, "cNvPicPr"))
    fill = ET.SubElement(picture, q(PIC, "blipFill"))
    ET.SubElement(fill, q(A, "blip"), {q(R, "embed"): rid})
    stretch = ET.SubElement(fill, q(A, "stretch"))
    ET.SubElement(stretch, q(A, "fillRect"))
    shape = ET.SubElement(picture, q(PIC, "spPr"))
    transform = ET.SubElement(shape, q(A, "xfrm"))
    ET.SubElement(transform, q(A, "off"), {"x": "0", "y": "0"})
    ET.SubElement(transform, q(A, "ext"), {"cx": str(cx), "cy": str(cy)})
    geometry = ET.SubElement(shape, q(A, "prstGeom"), {"prst": "rect"})
    ET.SubElement(geometry, q(A, "avLst"))
    return node


def yesno(value: object) -> str:
    return "是" if bool(value) else "否"


def display_target(value: str) -> str:
    mapping = {
        "cell_APOE_C1QC_TAM": "APOE+C1QC+ TAM",
        "cell_B": "B cell",
        "cell_DC": "Dendritic cell",
        "cell_Endothelial": "Endothelial",
        "cell_Epithelial_Normal": "Epithelial Normal",
        "cell_Fibroblast": "Fibroblast",
        "cell_Mast": "Mast cell",
        "cell_Monocyte": "Monocyte",
        "cell_NK": "NK cell",
        "cell_Other_Macrophage": "Other Macrophage",
        "cell_SMC_Pericyte": "SMC/Pericyte",
        "cell_SPP1_TAM": "SPP1+ TAM",
        "cell_T": "T cell",
        "cell_TYMP_High_Malignant": "TYMP-High Malignant",
        "cell_TYMP_Low_Malignant": "TYMP-Low Malignant",
        "cell_other_Malignant": "Other Malignant",
        "cell_TAM_Intermediate": "TAM Intermediate",
    }
    return mapping.get(value, value.replace("program_", "Program: ").replace("_", " "))


required = [SOURCE, BAYES_QC, PROGRAM_QC, MISTY_QC, COMMOT_QC, MISTY_TABLE, COMMOT_TABLE, *IMAGES]
for path in required:
    if not path.is_file():
        raise FileNotFoundError(path)
if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
    raise RuntimeError("Source DOCX hash changed; refusing to generate from an unreviewed source")

bayes_qc = json.loads(BAYES_QC.read_text())
program_qc = json.loads(PROGRAM_QC.read_text())
misty_qc = json.loads(MISTY_QC.read_text())
commot_qc = json.loads(COMMOT_QC.read_text())
assert bayes_qc["checks_passed"] == bayes_qc["checks_total"] == 64
assert bayes_qc["files_manifested"] == 42
assert bayes_qc["manual_visual_qc"]["png_passed"] == bayes_qc["manual_visual_qc"]["png_reviewed"] == 8
assert program_qc["n_pass"] == program_qc["n_checks"] == 115 and program_qc["n_files"] == 80
assert "10/10" in program_qc["manual_figure_review"]
assert misty_qc["checks_passed"] == 82 and misty_qc["checks_failed"] == 0
assert misty_qc["manifest_files"] == 198 and misty_qc["formal_png_manual_review_passed"] == 9
assert commot_qc["checks_passed"] == commot_qc["checks_total"] == 78 and commot_qc["files_hashed"] == 89

misty = pd.read_csv(MISTY_TABLE, sep="\t")
commot = pd.read_csv(COMMOT_TABLE, sep="\t")
assert len(misty) == 20 and len(commot) == 12
assert int(misty["positive_gain_both_sections"].sum()) == 17
assert int(misty["r2_gain_p_lt_0.05_both_sections"].sum()) == 3
assert int(misty["rmse_gain_p_lt_0.05_both_sections"].sum()) == 2
assert int(commot["both_positive_proxy_enrichment"].sum()) == 8
assert int(commot["replicated_proxy_support"].sum()) == 0

misty_rows = [["目标", "ΔR²（G1/G2）", "RMSE改善%（G1/G2）", "两例均正增益", "两例R² P<0.05", "两例RMSE P<0.05"]]
for row in misty.itertuples(index=False):
    misty_rows.append([
        display_target(row.target),
        f"{row._1:.2f} / {row._2:.2f}",
        f"{row._3:.2f} / {row._4:.2f}",
        yesno(row.positive_gain_both_sections),
        yesno(getattr(row, "_18")),
        yesno(getattr(row, "_19")),
    ])

# Avoid pandas itertuples renaming ambiguities for punctuated columns in the COMMOT table.
commot_rows = [["轴（尺度）", "来源代理 → 靶代理", "G1 log2富集/质量排名", "G2 log2富集/质量排名", "两例均正", "严格支持"]]
for _, row in commot.iterrows():
    scale = "110 µm" if row["scale"] == "contact_110um" else "200 µm"
    commot_rows.append([
        f"{row['axis_id']}（{scale}）",
        f"{row['source_proxy']} → {row['target_proxy']}",
        f"{row['offdiagonal_proxy_log2_enrichment_GLMF1']:.3f} / {int(row['offdiagonal_proxy_mass_rank_among_196_GLMF1'])}",
        f"{row['offdiagonal_proxy_log2_enrichment_GLMF2']:.3f} / {int(row['offdiagonal_proxy_mass_rank_among_196_GLMF2'])}",
        yesno(row["both_positive_proxy_enrichment"]),
        yesno(row["replicated_proxy_support"]),
    ])

summary_rows = [
    ["方法", "设计", "核心结果", "证据边界", "独立QC"],
    ["BayesSpace", "两切片独立q调参；q=4", "G1/G2边界spot 105/232", "计算域边界，不是人工病理tumor–NAT边界", "64/64；8/8图"],
    ["SPARK-X + Hotspot", "全基因空间变异筛选与模块化", "两例共同双阳性基因4,242；最强模块重叠37基因", "模块匹配为描述性；基因阳性不等于模块复制", "115/115；10/10图"],
    ["MISTy", "20目标；intra/juxta/para；10折", "两例均正增益17/20；名义显著R² 3/20、RMSE 2/20", "预测关联；P值未多重校正且非患者级重复", "82/82；9/9图"],
    ["COMMOT", "12条先验轴；主尺度110/200 µm", "两例均正富集8/12；严格支持0/12", "代理区室上的最优传输；不证明细胞特异直接通讯", "78/78；9/9图"],
]

with ZipFile(SOURCE) as source_zip:
    document_root = ET.fromstring(source_zip.read("word/document.xml"))
    rels_root = ET.fromstring(source_zip.read("word/_rels/document.xml.rels"))
    source_names = source_zip.namelist()
body = document_root.find(q(W, "body"))
if body is None:
    raise RuntimeError("DOCX document body not found")

# Update front matter and summary without changing the original report identity.
set_paragraph_text(direct_paragraph(body, "正式报告｜"), "正式报告｜基于GLMF1、GLMF2、TJ-RCC 12例、230,618细胞参考对象、三层论文对标验证及四类高级空间分析补丁｜POSTN–APOE增强版", bold=True)
set_paragraph_text(
    direct_paragraph(body, "【增强版新增说明】"),
    "【增强版新增说明】在原第4.11节论文对标验证基础上，新增第4.12节、图4.30–4.37及BayesSpace、SPARK-X/Hotspot、MISTy、COMMOT完整结果与QA审计。新增分析均限定为两张GLMF切片的计算复现或跨切片重复；不得表述为独立患者级重复、直接因果通讯、病理tumor–NAT边界或临床疗效证据。",
    bold=True,
    shading="FFF2CC",
)
set_paragraph_text(
    any_paragraph(body, "GLMF全套分析、"),
    "GLMF全套分析、TJ-RCC 12例扩展、19条定向轴严格验证、跨队列深度敏感性验证、TYMP–SPP1免疫生态、POSTN⁺CAF–APOE⁺TAM论文对标验证，以及BayesSpace、SPARK-X/Hotspot、MISTy和COMMOT高级空间分析补丁均已完成",
)
set_paragraph_text(
    direct_paragraph(body, "技术路线为："),
    "技术路线为：空间数据QC与预处理 → 单细胞参考状态构建 → cell2location spot级丰度映射 → 四类代谢通路GSVA评分 → 样本内独立生态位构建 → 邻域富集、同-spot相关与空间滞后 → 空间分块置换 → BayesSpace计算域与边界 → SPARK-X/Hotspot空间变异基因及模块 → MISTy多视图预测 → COMMOT候选轴最优传输。人工病理tumor–NAT距离仍未获得，因此计算域边界不得替代病理边界。",
)
insert_after(body, direct_paragraph(body, "在原有同-spot"), [paragraph(
    "高级分析补丁显示：BayesSpace在两切片均选择q=4；SPARK-X/Hotspot得到4,242个两例共同双阳性空间基因；MISTy在20个预设目标中有17个两例同时获得正R²及RMSE增益，但更严格的两例名义P<0.05仅R² 3/20、RMSE 2/20；COMMOT 12条主尺度候选轴中8/12两例均为正的非对角代理富集，但没有任何一条同时进入两例传输质量前20%（0/12）。这些结果支持跨切片的描述性重复，不构成患者级显著重复或因果证据。"
)])

# Add image relationships now so section construction can refer to stable IDs.
existing_rids = []
for rel in rels_root.findall(q(REL, "Relationship")):
    rid = rel.get("Id", "")
    if rid.startswith("rId") and rid[3:].isdigit():
        existing_rids.append(int(rid[3:]))
next_rid = max(existing_rids, default=0) + 1
docpr_ids = [int(n.get("id")) for n in document_root.iter(q(WP, "docPr")) if (n.get("id") or "").isdigit()]
next_docpr = max(docpr_ids, default=0) + 1
image_specs = []
for index, image in enumerate(IMAGES, start=1):
    rid = f"rId{next_rid}"
    next_rid += 1
    media_name = f"advanced_patch_v1_{index:02d}_{image.name}"
    ET.SubElement(rels_root, q(REL, "Relationship"), {
        "Id": rid,
        "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
        "Target": f"media/{media_name}",
    })
    width, height = png_size(image)
    image_specs.append((image, rid, media_name, width, height, next_docpr))
    next_docpr += 1

def fig(index: int, caption: str) -> list[ET.Element]:
    path, rid, media_name, width, height, docpr_id = image_specs[index]
    caption_node = paragraph(caption, bold=True)
    props = caption_node.find(q(W, "pPr"))
    if props is None:
        props = ET.Element(q(W, "pPr"))
        caption_node.insert(0, props)
    ET.SubElement(props, q(W, "keepNext"))
    ET.SubElement(props, q(W, "keepLines"))
    return [caption_node, image_paragraph(rid, media_name, width, height, docpr_id)]

section = [
    page_break(),
    paragraph("4.12 【高级分析补丁】BayesSpace、SPARK-X/Hotspot、MISTy与COMMOT", style="Heading2"),
    paragraph("本节在既有cell2location丰度、空间程序和候选轴框架上增加四类互补分析。所有方法均在GLMF1和GLMF2切片内独立运行，再比较跨切片方向或排名；GLMF1与GLMF2不是已确认的独立患者，也不是治疗前后配对，因此本节统一使用“两切片/两例共同”而非“患者级复现”。"),
    table(summary_rows, 15),
    paragraph("4.12.1 BayesSpace计算域与计算边界", style="Heading3"),
    paragraph("BayesSpace对每张切片独立完成q=2–10调参，两例均选择q=4。GLMF1四个域的spot数为507、451、1,920和321，肿瘤样域为1/2；GLMF2为450、898、1,481和1,051，肿瘤样域为2/4。基于相邻计算域得到GLMF1 105个、GLMF2 232个边界spot。签名距离定义为肿瘤样侧≤0、非肿瘤样侧≥0，方向与既有人工病理约定相反；该边界只反映表达驱动的计算分区，不能替代人工标注的tumor–NAT边界。"),
    *fig(0, "图4.30  GLMF1 BayesSpace q=4计算域及域注释。"),
    *fig(1, "图4.31  GLMF2 BayesSpace q=4计算域及域注释。"),
    page_break(),
    paragraph("4.12.2 SPARK-X/Hotspot空间变异基因与模块", style="Heading3"),
    paragraph("GLMF1共检验14,082个基因，SPARK-X经BY校正后13,648个显著，Hotspot正自相关11,269个，两者共同阳性11,136个；Hotspot形成8个模块、纳入427个基因。GLMF2检验13,652个基因，SPARK-X BY显著11,217个，Hotspot正自相关4,824个，共同阳性4,462个；形成10个模块、纳入225个基因。两例共有13,426个共同受检基因，其中4,242个在两例均为双阳性。"),
    paragraph("最强跨样本模块匹配为GLMF1-M1与GLMF2-M1，共享37个基因，Jaccard=0.209，描述性超几何BH P=1.18×10⁻⁶²。POSTN、APOE等目标基因可具有基因层面的空间信号，但未必进入可跨样本匹配的Hotspot模块；因此不能把“基因空间显著”改写成“模块已复现”。"),
    *fig(2, "图4.32  GLMF1与GLMF2 Hotspot模块重叠及描述性匹配。"),
    *fig(3, "图4.33  两切片共同受检基因的空间统计量排名一致性。"),
    page_break(),
    paragraph("4.12.3 MISTy多视图预测增益与重要性重复", style="Heading3"),
    paragraph("MISTy使用20个预设目标：14类cell2location丰度、5个表达程序及37基因共享基质程序。视图为同spot的intra、Delaunay邻接且距离≤110 µm的juxta，以及长度尺度200 µm、近场排除110 µm的Gaussian para；每个模型100棵树、10折交叉验证。与intra基线相比，R²增益为正且RMSE同时改善的目标在GLMF1为18/20、GLMF2为19/20，两例共同为17/20。"),
    paragraph("更严格地要求两例各自名义P<0.05时，R²增益仅3/20，RMSE改善仅2/20。当前mistyR实现使用单侧t检验但未指定paired=TRUE，而且20个目标的P值未作多重校正；交叉验证折或spot也不是患者级生物学重复。因此17/20表示方向性预测增益在两切片共同出现，不能写成17个目标获得患者级显著复现。各视图进入前10%重要性的跨切片重复目标-预测因子对为intra 20、juxta 15、para 13；重要性无符号且不代表因果效应。"),
    table(misty_rows, 14),
    *fig(4, "图4.34  20个MISTy目标在GLMF1/GLMF2的R²与RMSE预测增益。"),
    *fig(5, "图4.35  intra、juxta与para视图的前10%重要性对跨切片重复。"),
    page_break(),
    paragraph("4.12.4 COMMOT候选轴最优传输与代理富集", style="Heading3"),
    paragraph("COMMOT输入为全基因库计数标准化至每spot 10,000后log1p表达。候选配体/受体设置检测门槛max(20 spots, 1% spots)；11条分泌轴以200 µm为主尺度并以500 µm作敏感性分析，APP–CD74接触轴使用110 µm。异源多亚基受体取亚基表达最小值；最优传输参数为ε=0.1、ρ=10、距离权重0.25、最多10,000次迭代。"),
    paragraph("主尺度总传输中，对角同代理区室质量占比很高：GLMF1为0.329–0.958，GLMF2为0.294–0.952，因此主要推断使用排除来源=靶代理区室后的非对角传输。两例主尺度非对角传输质量排名相关Spearman ρ=0.936；500 µm敏感性ρ=0.900。12条预设轴中8/12在两例均呈正的非对角代理富集，但0/12同时进入两例196种代理组合的传输质量前20%；500 µm下为6/11正、0/11严格支持。全4,508种轴×代理组合仅100种满足“两例正富集且两例质量前20%”，主要由SMC/Pericyte与恶性细胞组成驱动，说明总库中存在强区室结构，但未形成预设轴的严格复现证据。"),
    table(commot_rows, 14),
    *fig(6, "图4.36  COMMOT非对角传输质量在两切片间的排名一致性。"),
    *fig(7, "图4.37  12条预设候选轴的非对角代理富集与严格支持判定。"),
    paragraph("4.12.5 QA、版本与证据边界", style="Heading3"),
    paragraph("独立复核结果为：BayesSpace 64/64项通过、42个文件入清单、8/8张终版PNG人工通过；SPARK-X/Hotspot 115/115项通过、80个文件入清单、10/10张终版PNG人工通过；MISTy 82/82项通过、198个文件入清单、9/9张终版PNG人工通过；COMMOT 78/78项通过、89个文件完成哈希、9/9张终版PNG人工通过。SPARK跨样本v1标签重叠图、MISTy被替代图以及COMMOT v2两张总分图的失败记录均保留在审计中，但本报告只使用修订后的终版图。COMMOT v1在末端元数据阶段因变量名遮蔽失败；与正式v2的核心数值差异仅为浮点误差量级（绝对约10⁻¹³、相对约10⁻¹¹），不改变结论。"),
    paragraph("综合解释：四类补丁强化了计算域、空间变异、可预测邻域结构及候选轴传输的描述性证据，同时也给出了明确负结果。两张切片共同出现不等于独立患者复现；名义P值不等于多重校正后显著；预测重要性、空间共现和最优传输均不等于因果、方向性分泌或受体激活。", bold=True, shading="FFF2CC"),
]

heading5 = direct_paragraph(body, "5. 主要结论", style="Heading1")
insert_index = list(body).index(heading5)
for node in section:
    body.insert(insert_index, node)
    insert_index += 1

insert_after(body, direct_paragraph(body, "【论文对标增强结论】"), [paragraph(
    "【高级分析补丁结论】两张GLMF切片均显示可解析的BayesSpace计算域、大量共同空间变异基因、广泛但强度有限的MISTy多视图预测增益，以及COMMOT非对角传输排名的一致性。最稳健的定量描述为：共同双阳性空间基因4,242个；MISTy两例共同正增益17/20，但两例名义显著仅R² 3/20、RMSE 2/20；COMMOT预设轴两例均正8/12，但严格支持0/12。结论应保留“描述性跨切片重复”层级。"
)])
insert_after(body, direct_paragraph(body, "【论文对标增强限制】"), [paragraph(
    "【高级分析补丁限制】BayesSpace边界不是人工病理tumor–NAT边界；SPARK-X/Hotspot的基因与模块证据层级不同；MISTy检验未配对、未跨20目标校正，折和spot不是患者重复；COMMOT结果依赖代理区室、候选轴、距离尺度与最优传输参数，且主尺度存在明显对角质量。GLMF1/GLMF2不能确认独立患者身份，也没有干预、蛋白或临床结局，因此不能从本节推导直接通讯、因果机制、治疗前后变化或患者级显著复现。"
)])

document_bytes = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
rels_bytes = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(SOURCE) as source_zip, ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
    for name in source_names:
        if name == "word/document.xml":
            output_zip.writestr(name, document_bytes)
        elif name == "word/_rels/document.xml.rels":
            output_zip.writestr(name, rels_bytes)
        else:
            output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
    for image, _, media_name, _, _, _ in image_specs:
        output_zip.writestr(f"word/media/{media_name}", image.read_bytes())

with ZipFile(OUTPUT) as check_zip:
    bad = check_zip.testzip()
    if bad is not None:
        raise RuntimeError(f"Corrupt member in staged DOCX: {bad}")
    final_root = ET.fromstring(check_zip.read("word/document.xml"))
    final_text = "\n".join(node_text(p) for p in final_root.iter(q(W, "p")))
    for marker in [
        "4.12 【高级分析补丁】",
        "图4.30", "图4.31", "图4.32", "图4.33", "图4.34", "图4.35", "图4.36", "图4.37",
        "MISTy两例共同正增益17/20", "COMMOT预设轴两例均正8/12", "严格支持0/12",
        "【高级分析补丁结论】", "【高级分析补丁限制】",
    ]:
        if marker not in final_text:
            raise RuntimeError(f"Required report marker missing: {marker}")
    media_names = set(check_zip.namelist())
    for _, _, media_name, _, _, _ in image_specs:
        if f"word/media/{media_name}" not in media_names:
            raise RuntimeError(f"Embedded image missing: {media_name}")

metadata = {
    "status": "staged_and_structurally_verified",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "source": str(SOURCE),
    "source_size": SOURCE.stat().st_size,
    "source_sha256": sha256(SOURCE),
    "staged_output": str(OUTPUT),
    "output_size": OUTPUT.stat().st_size,
    "output_sha256": sha256(OUTPUT),
    "docx_zip_test": "PASS",
    "added_section": "4.12",
    "added_figures": [f"4.{number}" for number in range(30, 38)],
    "added_images": [
        {"source": str(path), "embedded_name": media_name, "size": path.stat().st_size, "sha256": sha256(path)}
        for path, _, media_name, _, _, _ in image_specs
    ],
    "input_qc": {
        "BayesSpace": "64/64 checks; 42 files; 8/8 PNG",
        "SPARK-X_Hotspot": "115/115 checks; 80 files; 10/10 PNG",
        "MISTy": "82/82 checks; 198 files; 9/9 PNG",
        "COMMOT": "78/78 checks; 89 hashed files; 9/9 PNG",
    },
    "key_counts": {
        "misty_targets": 20,
        "misty_positive_gain_both": 17,
        "misty_r2_nominal_both": 3,
        "misty_rmse_nominal_both": 2,
        "commot_primary_axes": 12,
        "commot_both_positive": 8,
        "commot_strict_support": 0,
    },
}
METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(metadata, ensure_ascii=False, indent=2))
