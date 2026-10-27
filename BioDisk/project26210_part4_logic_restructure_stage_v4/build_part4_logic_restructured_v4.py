#!/usr/bin/env python3
from __future__ import annotations

import copy
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.font_manager import fontManager

import sys

UTILS = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_update_stage")
sys.path.insert(0, str(UTILS))
from ooxml_report_utils_v1 import (  # noqa: E402
    A, REL, W, WP, caption_figure, node_text, page_break, paragraph, png_size,
    q, set_paragraph_text, table,
)


SOURCE = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3/"
    "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_"
    "POSTN_APOE论文对标增强版_v3_candidate.docx"
)
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v4")
OUTPUT = STAGE / "project26210_第四部分空间转录组分析报告_主线逻辑重构版_v4_candidate.docx"
FLOWCHART = STAGE / "part4_mainline_framework_v4.png"
METADATA = STAGE / "report_generation_metadata_v4.json"
EXPECTED_SOURCE_SHA256 = "d9ed4172e3a9b675946cdf0fa0b26c4d44154c33e40c72048001c2dee2b19940"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_box(ax, x, y, w, h, title, lines, face, edge):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=1.5, edgecolor=edge, facecolor=face,
    )
    ax.add_patch(box)
    ax.text(x + 0.012, y + h - 0.027, title, fontsize=12.5, fontweight="bold",
            color=edge, va="top")
    ax.text(x + 0.014, y + h - 0.064, "\n".join(f"• {line}" for line in lines),
            fontsize=9.2, color="#243247", va="top", linespacing=1.35)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14,
                                 linewidth=1.5, color="#405A78"))


def build_flowchart() -> None:
    fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc")
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
    fig, ax = plt.subplots(figsize=(16, 9.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    add_box(ax, 0.018, 0.235, 0.165, 0.61, "数据与因果锚点", [
        "GLMF1/GLMF2：发现层",
        "TJ-RCC 12例：患者验证层",
        "单细胞参考：状态与连续谱",
        "客户湿实验：因果锚点",
    ], "#F4F7FB", "#4D6380")

    ax.text(0.215, 0.955, "第一层：建立可信空间地图", fontsize=14, fontweight="bold", color="#1769AA")
    add_box(ax, 0.215, 0.735, 0.225, 0.185, "模块A｜数据基础", [
        "空间QC与有效深度审计", "统一单细胞参考与cell2location", "四类代谢转录程序定位",
    ], "#EAF4FE", "#1769AA")
    add_box(ax, 0.465, 0.735, 0.225, 0.185, "模块B｜组织分区", [
        "分样本独立生态位", "BayesSpace计算域", "计算边界与病理边界严格区分",
    ], "#EAF4FE", "#1769AA")
    add_box(ax, 0.715, 0.735, 0.255, 0.185, "模块C｜空间程序", [
        "SPARK-X/Hotspot无偏筛选", "MISTy多尺度预测结构", "COMMOT候选运输代理",
    ], "#EAF4FE", "#1769AA")

    ax.text(0.215, 0.685, "第二层：GLMF1/GLMF2发现与候选收窄", fontsize=14, fontweight="bold", color="#277B5D")
    add_box(ax, 0.215, 0.485, 0.35, 0.165, "模块D｜主轴发现", [
        "TYMP-High与TAM/Endothelial空间关系", "同spot、邻域滞后与分块置换", "区分绝对共占位与组成特异性",
    ], "#EAF8F1", "#277B5D")
    add_box(ax, 0.595, 0.485, 0.375, 0.165, "模块E｜候选机制筛查", [
        "APOE–TREM2、APP–CD74与MIF轴", "TYMP–Polyamine/MTA转录背景", "保留阳性、异质及阴性结果",
    ], "#EAF8F1", "#277B5D")

    ax.text(0.215, 0.435, "第三层：12名患者外部验证与机制收敛", fontsize=14, fontweight="bold", color="#7B4FA3")
    add_box(ax, 0.215, 0.245, 0.225, 0.155, "模块F｜患者级重复", [
        "12例独立建模", "患者内分块、患者间方向检验", "跨队列深度敏感性",
    ], "#F4EDFA", "#7B4FA3")
    add_box(ax, 0.465, 0.245, 0.225, 0.155, "模块G｜证据收窄", [
        "宽泛空间占位", "19条定向通讯严格检验", "患者级高级方法复核",
    ], "#F4EDFA", "#7B4FA3")
    add_box(ax, 0.715, 0.245, 0.255, 0.155, "模块H｜机制扩展", [
        "TYMP–MTA–SPP1状态", "SPP1–vascular niche与T细胞平衡", "髓系横断面状态排序",
    ], "#F4EDFA", "#7B4FA3")

    ax.text(0.215, 0.198, "第四层：证据分级与可交付结论", fontsize=14, fontweight="bold", color="#B15C16")
    add_box(ax, 0.215, 0.03, 0.755, 0.13, "模块I｜最终判定", [
        "稳定支持：跨患者重复且通过校正；条件支持：样本/尺度/路线敏感；未支持：严格模型或外部队列不通过",
        "空间数据用于人组织定位和外部一致性；直接通讯、代谢中介与因果机制仍由湿实验承担",
    ], "#FFF3E7", "#B15C16")

    arrow(ax, (0.183, 0.55), (0.215, 0.82))
    arrow(ax, (0.58, 0.735), (0.58, 0.65))
    arrow(ax, (0.58, 0.485), (0.58, 0.40))
    arrow(ax, (0.58, 0.245), (0.58, 0.16))
    arrow(ax, (0.44, 0.825), (0.465, 0.825))
    arrow(ax, (0.69, 0.825), (0.715, 0.825))
    arrow(ax, (0.565, 0.565), (0.595, 0.565))
    arrow(ax, (0.44, 0.32), (0.465, 0.32))
    arrow(ax, (0.69, 0.32), (0.715, 0.32))

    fig.savefig(FLOWCHART, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def set_style(node: ET.Element, style: str) -> None:
    props = node.find(q(W, "pPr"))
    if props is None:
        props = ET.Element(q(W, "pPr"))
        node.insert(0, props)
    pstyle = props.find(q(W, "pStyle"))
    if pstyle is None:
        pstyle = ET.SubElement(props, q(W, "pStyle"))
    pstyle.set(q(W, "val"), style)


def clone_range(children: list[ET.Element], start: int, stop: int) -> list[ET.Element]:
    return [copy.deepcopy(node) for node in children[start:stop]]


def retitle(nodes: list[ET.Element], old_prefix: str, new_text: str, style: str | None = None) -> None:
    matches = [node for node in nodes if node.tag == q(W, "p") and node_text(node).strip().startswith(old_prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"Retitle expected one match for {old_prefix!r}, found {len(matches)}")
    set_paragraph_text(matches[0], new_text)
    if style:
        set_style(matches[0], style)


def renumber_source_figures(body: ET.Element) -> dict[str, str]:
    ordered_old = [
        "1", "2", "3A", "3B", "4", "5", "6",
        "30", "31", "32", "33", "34", "35", "36", "37",
        "9", "10", "11", "12", "13", "17", "7", "8", "14", "15", "16",
        "38", "39", "40", "41", "42", "43", "44", "45",
        "18", "19", "20", "21", "22", "23",
        "46", "47", "48", "49", "50", "51", "52", "53",
        "24", "25", "26", "27", "28", "29",
    ]
    new_labels = [
        "2", "3", "4A", "4B", "5", "6", "7",
        *[str(i) for i in range(8, 16)],
        *[str(i) for i in range(16, 55)],
    ]
    if len(ordered_old) != 54 or len(new_labels) != 54:
        raise RuntimeError("Figure map length error")
    mapping = dict(zip(ordered_old, new_labels, strict=True))
    pattern = re.compile(r"((?:图|Figure)\s*)4\.(\d+)([AB]?)")
    for text_node in body.iter(q(W, "t")):
        if not text_node.text:
            continue
        def repl(match: re.Match[str]) -> str:
            old = match.group(2) + match.group(3)
            if old not in mapping:
                return match.group(0)
            return f"{match.group(1)}4.{mapping[old]}"
        text_node.text = pattern.sub(repl, text_node.text)
    return mapping


def add_image_relationship(document_root: ET.Element, rels_root: ET.Element):
    rid_numbers = []
    for rel in rels_root.findall(q(REL, "Relationship")):
        rid = rel.get("Id", "")
        if rid.startswith("rId") and rid[3:].isdigit():
            rid_numbers.append(int(rid[3:]))
    rid = f"rId{max(rid_numbers, default=0) + 1}"
    media_name = "part4_mainline_framework_v4.png"
    ET.SubElement(rels_root, q(REL, "Relationship"), {
        "Id": rid,
        "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
        "Target": f"media/{media_name}",
    })
    docpr_ids = [int(node.get("id")) for node in document_root.iter(q(WP, "docPr"))
                 if (node.get("id") or "").isdigit()]
    width, height = png_size(FLOWCHART)
    return FLOWCHART, rid, media_name, width, height, max(docpr_ids, default=0) + 1


def build_document() -> dict:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("Source report fingerprint changed")
    for path in (OUTPUT, FLOWCHART, METADATA):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite existing output: {path}")
    STAGE.mkdir(parents=True, exist_ok=True)
    build_flowchart()

    with ZipFile(SOURCE) as archive:
        document_root = ET.fromstring(archive.read("word/document.xml"))
        rels_root = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        source_names = archive.namelist()
    body = document_root.find(q(W, "body"))
    if body is None:
        raise RuntimeError("Missing document body")
    children = list(body)
    if len(children) != 338 or children[-1].tag != q(W, "sectPr"):
        raise RuntimeError(f"Unexpected source body structure: {len(children)} blocks")
    section_properties = copy.deepcopy(children[-1])

    cover = clone_range(children, 0, 8)
    set_paragraph_text(cover[1], "主线逻辑重构版最终分析报告")
    set_paragraph_text(cover[3], "项目主线：TYMP相关恶性状态—髓系/血管生态位—患者级外部验证")
    set_paragraph_text(cover[5], "正式报告｜GLMF1/GLMF2发现层｜TJ-RCC 12例验证层｜单细胞髓系状态层｜逻辑重构版 v4")
    set_paragraph_text(cover[6], "【重构原则】保留原v3全部53张分析图与关键数值；仅重排叙事、统一章节和图号，并新增一张总流程图。未新增未经验证的生物学结果。", bold=True, shading="EAF4FE")
    cover_replacements = {
        "单细胞参考构建 → cell2location 映射 → 空间生态位与共定位":
            "空间定位 → GLMF1/GLMF2发现 → TJ-RCC患者验证 → 机制收敛 → 证据分级",
        "2026年8月4日（血管生态位与髓系状态排序更新）":
            "2026年8月5日（主线逻辑重构版）",
    }
    for text_node in cover[4].iter(q(W, "t")):
        if text_node.text in cover_replacements:
            text_node.text = cover_replacements[text_node.text]

    chapter2_mapping = clone_range(children, 38, 48)
    retitle(chapter2_mapping, "2. cell2location", "2.1 统一参考映射与代谢转录程序定位", "Heading2")
    chapter2_niches = clone_range(children, 48, 59)
    retitle(chapter2_niches, "3. 空间生态位", "2.2 分样本空间生态位与组织区室", "Heading2")
    retitle(chapter2_niches, "3.1", "2.2.1 GLMF1与GLMF2独立构建", "Heading3")

    discovery_core = clone_range(children, 60, 82)
    retitle(discovery_core, "4.1 ", "3.1 TYMP-High与SPP1+ TAM：同位、邻近与组成对照", "Heading2")
    retitle(discovery_core, "4.2 ", "3.2 APOE+C1QC+ TAM：区室分离信号", "Heading2")
    retitle(discovery_core, "4.3 ", "3.3 分块置换收窄候选空间关系", "Heading2")
    retitle(discovery_core, "4.4 ", "3.4 MIF–CD74/CXCR4候选轴及髓系广泛性", "Heading2")

    # The final blank block in each source section already carries a page break.
    # Exclude it because the restructured document inserts chapter page breaks explicitly.
    glmf_advanced = clone_range(children, 180, 218)
    retitle(glmf_advanced, "4.12 ", "3.5 多方法深化：计算域、空间程序、多尺度预测与运输代理", "Heading2")
    for old, new in [
        ("4.12.1", "3.5.1 BayesSpace计算域与计算边界"),
        ("4.12.2", "3.5.2 SPARK-X/Hotspot空间变异基因与模块"),
        ("4.12.3", "3.5.3 MISTy多视图预测增益与重要性重复"),
        ("4.12.4", "3.5.4 COMMOT候选轴最优传输与代理富集"),
        ("4.12.5", "3.5.5 QA、版本与证据边界"),
    ]:
        retitle(glmf_advanced, old, new, "Heading3")

    tj_c2l = clone_range(children, 90, 111)
    retitle(tj_c2l, "4.6 ", "4.1 TJ-RCC 12例统一映射与队列空间图谱", "Heading2")
    tj_depth = clone_range(children, 125, 135)
    retitle(tj_depth, "4.8 ", "4.2 跨队列有效深度敏感性", "Heading2")
    tj_p0 = clone_range(children, 82, 90)
    retitle(tj_p0, "4.5 ", "4.3 基因层面P0外部重复验证", "Heading2")
    tj_strict = clone_range(children, 111, 125)
    retitle(tj_strict, "4.7 ", "4.4 从空间邻近收窄到19条定向通讯", "Heading2")
    for node in tj_strict:
        for text_node in node.iter(q(W, "t")):
            if text_node.text:
                text_node.text = text_node.text.replace("第4.5–4.6节", "第4.3–4.4节")
    tj_advanced = clone_range(children, 219, 254)
    retitle(tj_advanced, "4.13 ", "4.5 12例患者级高级方法复核", "Heading2")
    for old, new in [
        ("4.13.1", "4.5.1 BayesSpace计算域及程序分布"),
        ("4.13.2", "4.5.2 SPARK-X/Hotspot空间程序重复"),
        ("4.13.3", "4.5.3 MISTy多尺度生态预测"),
        ("4.13.4", "4.5.4 靶向COMMOT与可估计性边界"),
        ("4.13.5", "4.5.5 完成性、QA与患者图"),
    ]:
        retitle(tj_advanced, old, new, "Heading3")

    mech_immune = clone_range(children, 135, 156)
    retitle(mech_immune, "4.9 ", "5.1 TYMP–SPP1免疫生态与Polyamine–MTA背景", "Heading2")
    mech_state = clone_range(children, 255, 276)
    retitle(mech_state, "4.14 ", "5.2 TYMP–MTA-context–SPP1巨噬细胞状态与下游表型", "Heading2")
    retitle(mech_state, "4.14.1", "5.2.1 MTA转录代理中介与下游表型", "Heading3")
    retitle(mech_state, "4.14.2", "5.2.2 QA与机制证据边界", "Heading3")
    mech_vascular = clone_range(children, 276, 297)
    retitle(mech_vascular, "4.15 ", "5.3 血管生态位与髓系横断面状态排序", "Heading2")
    retitle(mech_vascular, "4.15.1", "5.3.1 TYMP–SPP1血管生态位的患者级空间关系", "Heading3")
    retitle(mech_vascular, "4.15.2", "5.3.2 髓系横断面状态排序", "Heading3")
    retitle(mech_vascular, "4.15.3", "5.3.3 QA与证据边界", "Heading3")
    mech_paper = clone_range(children, 159, 177)
    retitle(mech_paper, "4.11 ", "5.4 平行机制对标：TGFβ1–POSTN–APOE–CD8轴", "Heading2")
    mech_boundary = clone_range(children, 156, 158)
    retitle(mech_boundary, "4.10 ", "5.5 尚未完成：病理tumor–NAT边界", "Heading2")

    for child in list(body):
        body.remove(child)
    for node in cover:
        body.append(node)

    body.append(paragraph("执行摘要", style="Heading1"))
    body.append(paragraph(
        "本报告围绕一个主问题展开：TYMP相关恶性细胞状态是否对应可重复的髓系与血管空间生态，并能否在独立患者中收敛为可信的人组织证据。分析不再按工具或补丁时间排列，而按‘发现—外部验证—机制收敛—证据分级’组织。"
    ))
    body.append(table([
        ["证据层", "主要对象", "本层回答的问题", "当前结论"],
        ["空间发现", "GLMF1/GLMF2", "候选关系在哪里、是否跨切片同向", "形成TYMP–TAM/血管、APOE–TREM2和MIF等候选；部分关系对尺度或组成敏感"],
        ["患者验证", "TJ-RCC 12名独立患者", "候选能否跨患者重复", "APOE–TREM2空间占位、TYMP-High→Endothelial及两条血管生态位关系较稳定"],
        ["机制收敛", "12例空间数据+单细胞参考+客户湿实验", "相关性是否支持TYMP–MTA–SPP1机制", "空间层仅提供趋势/异质性；直接因果仍由湿实验锚定"],
        ["证据边界", "严格通讯、代谢中介、论文闭环", "能否升级为直接机制", "19条定向通讯0/19通过；MTA中介与完整POSTN–APOE–CD8闭环未获支持"],
    ], 15))
    body.append(paragraph(
        "一句话结论：当前最可靠的人组织空间证据不是一条已经闭合的配体–受体因果链，而是‘TYMP相关恶性区域—血管/髓系生态位—T细胞平衡’的患者级关联框架；其中若干具体环节获得支持，另一些环节表现为异质、路线敏感或明确未支持。",
        bold=True, shading="FFF2CC",
    ))
    body.append(paragraph("图FLOW 第四部分空间分析主线：从空间发现到患者验证，再到机制边界。", bold=True))
    body.append(paragraph("[FLOWCHART_IMAGE_PLACEHOLDER]"))
    body.append(paragraph("目录与阅读路径", style="Heading1"))
    for line in [
        "1. 研究问题与证据策略",
        "2. 数据基础与空间图谱",
        "3. GLMF1/GLMF2发现层：候选关系与方法深化",
        "4. TJ-RCC 12名患者验证层：从空间重复到严格收窄",
        "5. 机制收敛层：TYMP–MTA–SPP1、血管生态位与平行对标",
        "6. 证据分级、主要结论与下一步",
    ]:
        body.append(paragraph(line, style="Normal"))

    body.append(page_break())
    body.append(paragraph("1. 研究问题与证据策略", style="Heading1"))
    body.append(paragraph(
        "第三部分提出TYMP相关恶性状态及髓系互作假设；第四部分的职责是把该假设放回真实组织空间，并判断哪些环节能够跨切片、跨患者重复。报告依次回答：目标状态能否可靠定位；候选关系是同spot共占位、邻域相邻还是组成变化；发现层结果能否在12名独立患者中重复；空间证据能否支持具体通讯、代谢中介或因果机制。"
    ))
    body.append(table([
        ["层级", "统计单位", "核心方法", "允许解释"],
        ["定位层", "spot/切片", "QC、cell2location、GSVA、生态位", "组织内位置和异质性"],
        ["发现层", "GLMF1/GLMF2分别分析", "同spot、spatial lag、分块置换、BayesSpace、SPARK/Hotspot、MISTy、COMMOT", "候选空间关系与跨切片方向"],
        ["验证层", "12名患者", "患者内分块、患者间符号检验、BH校正", "患者级重复或异质性"],
        ["机制层", "患者+单细胞样本/供者+湿实验", "状态评分、横断面DPT、对标检验", "空间相关性；因果由干预实验承担"],
    ], 15))
    body.append(paragraph(
        "证据判定规则：患者级、校正后重复优先于spot级显著；空间占位不等于定向通讯；转录代谢代理不等于代谢物或通量；横断面DPT不等于谱系追踪。",
        bold=True, shading="EAF4FE",
    ))
    for node in clone_range(children, 28, 38):
        body.append(node)

    body.append(page_break())
    body.append(paragraph("2. 数据基础与空间图谱", style="Heading1"))
    body.append(paragraph(
        "本章只回答‘分析对象是否可比、状态能否定位、组织是否存在稳定区室’。在此基础上，后续章节才讨论具体候选轴。GLMF1/GLMF2均使用相同14类参考状态；TJ-RCC 12例沿用统一参考签名，但所有患者分别建模和汇总。"
    ))
    for node in chapter2_mapping + chapter2_niches:
        body.append(node)

    body.append(page_break())
    body.append(paragraph("3. GLMF1/GLMF2发现层：候选关系与方法深化", style="Heading1"))
    body.append(paragraph(
        "发现层先用可解释的同spot、邻域和分块统计筛选候选，再用计算域、无偏空间程序、多尺度预测和最优传输代理判断其稳健性。两张切片共同出现只说明跨切片方向一致，不等同于患者级复现。"
    ))
    for node in discovery_core + glmf_advanced:
        body.append(node)

    body.append(page_break())
    body.append(paragraph("4. TJ-RCC 12名患者验证层：从空间重复到严格收窄", style="Heading1"))
    body.append(paragraph(
        "验证层以患者为生物学重复，先建立12例统一空间图谱并审计有效深度，再从不依赖反卷积的基因层P0关系进入cell2location状态关系、19条定向通讯和四类高级方法。验证标准逐步提高，因此‘宽泛空间关系重复、严格通讯未通过’并不矛盾。"
    ))
    for node in tj_c2l + tj_depth + tj_p0 + tj_strict + tj_advanced:
        body.append(node)

    body.append(page_break())
    body.append(paragraph("5. 机制收敛层：TYMP–MTA–SPP1、血管生态位与平行对标", style="Heading1"))
    body.append(paragraph(
        "本章围绕客户湿实验主轴进行人组织空间验证：先区分TYMP-High与SPP1+ TAM的绝对共占位和组成特异性，再检验巨噬细胞极化、MTA转录代理、血管生态位与T细胞平衡，最后用单细胞横断面状态排序和POSTN–APOE论文轴作正交对照。"
    ))
    for node in mech_immune + mech_state + mech_vascular + mech_paper + mech_boundary:
        body.append(node)

    body.append(page_break())
    body.append(paragraph("6. 证据分级、主要结论与下一步", style="Heading1"))
    body.append(paragraph("6.1 证据矩阵", style="Heading2"))
    body.append(table([
        ["证据等级", "结论", "关键依据", "报告口径"],
        ["稳定支持", "APOE–TREM2空间占位/检出概率具有跨患者重复", "P0与cell2location两条路线均11/12正向并通过患者级FDR", "空间占位或微环境共现；不是直接诱导"],
        ["稳定支持", "TYMP-High恶性状态与邻近Endothelial跨患者正向", "4×4为12/12、5×5为11/12，方向FDR通过", "血管相关空间生态；不是已证实VEGF通讯"],
        ["稳定支持", "SPP1–vascular niche与T-cell balance负向；TYMP–MTA core与SPP1–vascular niche正向", "两尺度均12/12方向支持，全局FDR=0.00220", "患者级空间关联；不是免疫排斥或治疗耐药因果"],
        ["条件支持", "TYMP-High与SPP1+ TAM/极化状态", "GLMF部分样本/尺度正向；TJ-RCC总体异质；参考轻依赖路线FDR=0.154", "趋势性人组织支持，注明样本、尺度和路线敏感"],
        ["条件支持", "MIF–CD74/CXCR4与POSTN–APOE上游环节", "GLMF局部重复或部分样本/患者支持；外部完整链不重复", "候选环节，不拼接为闭环"],
        ["未支持", "19条具体定向通讯", "严格0/19；宽松队列P<0.05亦0/19", "不得升级为空间通讯机制"],
        ["未支持", "MTA转录代理中介与完整POSTN–APOE–CD8闭环", "间接效应接近0；反馈和CD8排斥未复现或方向相反", "需要代谢物、蛋白和干预验证"],
        ["未完成", "病理tumor–NAT有符号距离", "缺少人工病理多边形", "不展示伪边界，不给出梯度结论"],
    ], 14))
    body.append(paragraph("6.2 主线结论", style="Heading2"))
    body.append(paragraph(
        "第一，空间定位和队列扩展均已完成：GLMF1/GLMF2用于发现和方法深化，TJ-RCC 12名独立患者用于患者级复核，单细胞参考用于髓系横断面状态描述。第二，跨患者最稳定的不是单一配体–受体链，而是TYMP相关恶性状态、血管区室、APOE–TREM2占位及SPP1–vascular niche/T细胞平衡之间的空间生态联系。第三，TYMP–5′-MTA–SPP1主轴在人组织中获得趋势性和背景性支持，但SPP1特异性、MTA中介及促恶性因果没有被空间数据独立证明。"
    ))
    body.append(paragraph(
        "因此，第四部分对项目主机制的贡献应定位为‘真实人组织中的空间锚定、患者异质性评估和可转化验证优先级排序’，而不是重复客户湿实验已经承担的因果证明。",
        bold=True, shading="EAF8F1",
    ))
    body.append(paragraph("6.3 当前限制", style="Heading2"))
    for text_value in [
        "GLMF1/GLMF2不能解释为治疗前后；两切片共同方向不等于患者级重复。",
        "Visium spot混合多类细胞；分块置换可降低但不能完全消除空间自相关和组织区室混杂。",
        "GSVA与Polyamine–MTA context均为转录代理，不是代谢通量或5′-MTA实测。",
        "BayesSpace边界不是病理边界；MISTy重要性和COMMOT最优传输不是因果效应。",
        "当前无治疗响应、生存、纵向取样、RNA velocity或谱系追踪，不能推导进展、耐药或谱系转换。",
    ]:
        body.append(paragraph("• " + text_value))
    body.append(paragraph("6.4 建议的下一步验证顺序", style="Heading2"))
    body.append(table([
        ["优先级", "下一步", "目的"],
        ["P1", "用多重IF/IHC或高分辨空间技术同时定位TYMP肿瘤细胞、SPP1/C1QC巨噬细胞、Endothelial与T细胞", "验证细胞来源、邻近方向与混合spot偏差"],
        ["P1", "在既有湿实验体系中加入5′-MTA定量、MTAP/TYMP干预与SPP1/C1QC状态读出", "直接检验代谢中介和巨噬细胞状态特异性"],
        ["P2", "围绕SPP1–vascular niche/T-cell balance进行病例分层和独立队列复核", "验证当前最稳定的新患者级空间关系"],
        ["P2", "完成GLMF1/GLMF2人工病理tumor–NAT标注", "建立真实边界梯度，替代计算型边界"],
        ["P3", "若获得治疗或生存标签，再开展疗效/预后关联", "避免用横断面空间关系替代临床终点"],
    ], 15))
    body.append(paragraph(
        "最终总述：本次重构后，报告从‘按分析补丁罗列’转为‘发现—验证—机制—边界’证据链。全部原始图表仍保留，但阳性、条件性和阴性结果被置于同一逻辑框架中，便于客户沿主线阅读并判断哪些结论可以进入正文、哪些仅适合作为补充或下一步假设。",
        bold=True, shading="FFF2CC",
    ))
    body.append(section_properties)

    figure_mapping = renumber_source_figures(body)
    flow_spec = add_image_relationship(document_root, rels_root)
    placeholder = [node for node in list(body) if node.tag == q(W, "p") and node_text(node) == "[FLOWCHART_IMAGE_PLACEHOLDER]"]
    if len(placeholder) != 1:
        raise RuntimeError("Flowchart placeholder missing")
    index = list(body).index(placeholder[0])
    body.remove(placeholder[0])
    flow_nodes = caption_figure(flow_spec, "图4.1  第四部分空间分析主线：从空间发现到患者验证，再到机制边界。", width_emu=6_250_000)
    # caption_figure returns caption + image; an explicit caption already precedes the placeholder.
    existing_caption = [node for node in list(body) if node.tag == q(W, "p") and node_text(node).startswith("图FLOW")]
    if len(existing_caption) != 1:
        raise RuntimeError("Flowchart caption placeholder missing")
    set_paragraph_text(existing_caption[0], "图4.1  第四部分空间分析主线：从空间发现到患者验证，再到机制边界。", bold=True)
    body.insert(index, flow_nodes[1])

    document_bytes = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
    rels_bytes = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)

    shanghai = timezone(timedelta(hours=8))
    now = datetime.now(shanghai)
    core_bytes = None
    with ZipFile(SOURCE) as source_zip:
        core_root = ET.fromstring(source_zip.read("docProps/core.xml"))
        dc = "http://purl.org/dc/elements/1.1/"
        dcterms = "http://purl.org/dc/terms/"
        title_node = core_root.find(q(dc, "title"))
        if title_node is not None:
            title_node.text = "project26210 第四部分空间转录组分析报告（主线逻辑重构版v4）"
        modified = core_root.find(q(dcterms, "modified"))
        if modified is not None:
            modified.text = now.isoformat(timespec="seconds")
        core_bytes = ET.tostring(core_root, encoding="utf-8", xml_declaration=True)

        with ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
            for name in source_names:
                if name == "word/document.xml":
                    output_zip.writestr(name, document_bytes)
                elif name == "word/_rels/document.xml.rels":
                    output_zip.writestr(name, rels_bytes)
                elif name == "docProps/core.xml":
                    output_zip.writestr(name, core_bytes)
                else:
                    output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
            output_zip.writestr(f"word/media/{flow_spec[2]}", FLOWCHART.read_bytes())

    with ZipFile(OUTPUT) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise RuntimeError(f"ZIP CRC failure: {bad}")

    metadata = {
        "status": "generated_pending_render_qc",
        "generated_at": now.isoformat(timespec="seconds"),
        "source": str(SOURCE),
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
        "output_size": OUTPUT.stat().st_size,
        "flowchart": str(FLOWCHART),
        "flowchart_sha256": sha256(FLOWCHART),
        "structure_change": "chronological patch stack -> discovery / patient validation / mechanism convergence / evidence grading",
        "original_numbered_figures_retained": 53,
        "original_captioned_figures_retained": 54,
        "new_flowchart_figures": 1,
        "figure_mapping": figure_mapping,
        "content_policy": "Reordered and synthesized only; no unverified biological result added.",
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return metadata


if __name__ == "__main__":
    result = build_document()
    print(json.dumps(result, ensure_ascii=False, indent=2))
