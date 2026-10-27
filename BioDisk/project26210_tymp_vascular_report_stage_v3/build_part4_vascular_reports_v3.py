#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd


UTILS = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_update_stage")
sys.path.insert(0, str(UTILS))
from ooxml_report_utils_v1 import (  # noqa: E402
    A, REL, W, WP, any_paragraph, caption_figure, direct_paragraph, insert_after,
    node_text, page_break, paragraph, png_size, q, set_paragraph_text, sha256,
    table, table_row,
)


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3")
FULL_SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
CURATED_SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
FULL_OUTPUT = STAGE / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3_candidate.docx"
CURATED_OUTPUT = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3_candidate.docx"
METADATA = STAGE / "report_generation_metadata_v3.json"

EXPECTED_FULL_SHA = "4c830304b52b7e7ad9c4e40902cc81558e5b50690e4089ac4ee52ffd0b4db69b"
EXPECTED_CURATED_SHA = "55b95f7a15b3f33ba9e3d10f65bb88a76d9104882bbc26585508790f5bc55f39"

VASCULAR = ROOT / "results/part4_spatial/tymp_vascular_niche_v1/tables/tjrcc_cohort_direction_summary_v1.tsv"
TRAJECTORY = ROOT / "results/part4_spatial/tymp_myeloid_state_ordering_v1/tables/cohort_pseudotime_direction_summary_v1.tsv"
SYNTHESIS = ROOT / "results/part4_spatial/tymp_vascular_trajectory_synthesis_v1"
FIGURES = [
    SYNTHESIS / "figures/vascular_niche_patient_effects_v1.png",
    SYNTHESIS / "figures/vascular_niche_cohort_summary_v1.png",
    SYNTHESIS / "figures/myeloid_cross_sectional_state_overview_v2.png",
    SYNTHESIS / "figures/myeloid_state_patient_directions_v1.png",
]
NUMERIC_QC = SYNTHESIS / "qc/independent_verification_summary_v2.json"
FIGURE_QC = SYNTHESIS / "qc/joint_figure_verification_summary_v1.json"
MANUAL_QC = SYNTHESIS / "qc/joint_figure_manual_review_v1.tsv"


RELATION_LABELS = {
    "SPP1_TAM_to_vascular_compartment": "SPP1 TAM vs vascular compartment",
    "SPP1_vascular_niche_to_T_cell_balance": "SPP1–vascular niche vs T-cell balance",
    "SPP1_vascular_niche_to_cytotoxicity": "SPP1–vascular niche vs cytotoxicity",
    "SPP1_vascular_niche_to_hypoxia_angiogenesis": "SPP1–vascular niche vs hypoxia/angiogenesis",
    "SPP1_vascular_niche_to_invasion_EMT": "SPP1–vascular niche vs invasion/EMT",
    "TYMP_MTA_core_to_SPP1_monocyte_balance": "TYMP–MTA core vs SPP1/monocyte balance",
    "TYMP_MTA_core_to_SPP1_vascular_niche": "TYMP–MTA core vs SPP1–vascular niche",
}
ROUTE_LABELS = {
    "route_independent": "route-independent",
    "reference_light": "reference-light",
    "reference_dependent": "reference-dependent",
}
SCORE_LABELS = {
    "Monocyte_program": "Monocyte program",
    "SPP1_state_no_SPP1": "SPP1-state（排除SPP1）",
    "C1QC_state": "C1QC-state",
    "Vascular_remodeling": "Vascular-remodeling",
    "Immunosuppression": "Immunosuppression",
    "MTA_context": "MTA-context",
    "SPP1_vs_C1QC": "SPP1-state minus C1QC-state",
}


def validate_sources() -> tuple[pd.DataFrame, pd.DataFrame]:
    required = [FULL_SOURCE, CURATED_SOURCE, VASCULAR, TRAJECTORY, NUMERIC_QC, FIGURE_QC, MANUAL_QC, *FIGURES]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(FULL_SOURCE) != EXPECTED_FULL_SHA:
        raise RuntimeError("Full v2 report fingerprint changed")
    if sha256(CURATED_SOURCE) != EXPECTED_CURATED_SHA:
        raise RuntimeError("Curated v2 report fingerprint changed")
    if FULL_OUTPUT.exists() or CURATED_OUTPUT.exists() or METADATA.exists():
        raise RuntimeError("Refusing to overwrite v3 staging outputs")
    numeric_qc = json.loads(NUMERIC_QC.read_text())
    figure_qc = json.loads(FIGURE_QC.read_text())
    manual = pd.read_csv(MANUAL_QC, sep="\t")
    if numeric_qc.get("status") != "passed" or numeric_qc.get("checks_passed") != 35 or numeric_qc.get("checks_failed") != 0:
        raise RuntimeError("Numeric QC is not 35/35 PASS")
    if figure_qc.get("status") != "passed" or figure_qc.get("figures_passed") != 4 or figure_qc.get("figures_failed") != 0:
        raise RuntimeError("Figure QC is not 4/4 PASS")
    if len(manual) != 4 or not manual["status"].eq("PASS").all():
        raise RuntimeError("Manual figure review is not 4/4 PASS")

    vascular = pd.read_csv(VASCULAR, sep="\t")
    trajectory = pd.read_csv(TRAJECTORY, sep="\t")
    if len(vascular) != 18 or set(vascular["grid"]) != {"4x4", "5x5"}:
        raise RuntimeError("Unexpected vascular cohort table")
    if len(trajectory) != 7 or set(trajectory["n_samples"]) != {63}:
        raise RuntimeError("Unexpected trajectory cohort table")

    key_t = vascular[(vascular["relation"] == "SPP1_vascular_niche_to_T_cell_balance") &
                     (vascular["route"] == "route_independent")].sort_values("grid")
    key_v = vascular[(vascular["relation"] == "TYMP_MTA_core_to_SPP1_vascular_niche") &
                     (vascular["route"] == "reference_dependent")].sort_values("grid")
    if not (len(key_t) == len(key_v) == 2 and (key_t["direction_supported_patients"] == 12).all()
            and (key_v["direction_supported_patients"] == 12).all()
            and np.allclose(key_t["sign_test_fdr_bh_global"], 0.002197265625)
            and np.allclose(key_v["sign_test_fdr_bh_global"], 0.002197265625)):
        raise RuntimeError("Key vascular results changed")
    t = trajectory.set_index("score")
    expected = {
        "C1QC_state": (57, 0.47482326630318866, 5.738490487915832e-11),
        "SPP1_state_no_SPP1": (41, 0.09120838747808518, 0.031604961742560056),
        "SPP1_vs_C1QC": (11, -0.27607333776555015, 3.904831277027714e-7),
    }
    for score, (positive, median, fdr) in expected.items():
        row = t.loc[score]
        if int(row.positive_samples) != positive or not np.isclose(row.median_rho, median) or not np.isclose(row.sign_test_fdr_bh, fdr):
            raise RuntimeError(f"Trajectory result changed: {score}")
    return vascular, trajectory


def relationship_image_specs(source: Path, images: list[Path], prefix: str):
    with ZipFile(source) as archive:
        document_root = ET.fromstring(archive.read("word/document.xml"))
        rels_root = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        names = archive.namelist()
    existing_rids = [int(x.get("Id")[3:]) for x in rels_root.findall(q(REL, "Relationship"))
                     if (x.get("Id") or "").startswith("rId") and x.get("Id")[3:].isdigit()]
    next_rid = max(existing_rids, default=0) + 1
    docpr_ids = [int(node.get("id")) for node in document_root.iter(q(WP, "docPr"))
                 if (node.get("id") or "").isdigit()]
    next_docpr = max(docpr_ids, default=0) + 1
    specs = []
    for index, image in enumerate(images, start=1):
        rid = f"rId{next_rid}"
        next_rid += 1
        media_name = f"{prefix}_{index:02d}_{image.name}"
        ET.SubElement(rels_root, q(REL, "Relationship"), {
            "Id": rid,
            "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
            "Target": f"media/{media_name}",
        })
        width, height = png_size(image)
        specs.append((image, rid, media_name, width, height, next_docpr))
        next_docpr += 1
    return document_root, rels_root, names, specs


def write_docx(source: Path, output: Path, document_root: ET.Element, rels_root: ET.Element,
               source_names: list[str], specs) -> None:
    document_bytes = ET.tostring(document_root, encoding="utf-8", xml_declaration=True)
    rels_bytes = ET.tostring(rels_root, encoding="utf-8", xml_declaration=True)
    with ZipFile(source) as source_zip, ZipFile(output, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
        for name in source_names:
            if name == "word/document.xml":
                output_zip.writestr(name, document_bytes)
            elif name == "word/_rels/document.xml.rels":
                output_zip.writestr(name, rels_bytes)
            else:
                output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
        for image, _, media_name, _, _, _ in specs:
            output_zip.writestr(f"word/media/{media_name}", image.read_bytes())
    with ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise RuntimeError(f"ZIP CRC failure: {output}")
        for _, _, media_name, _, _, _ in specs:
            if f"word/media/{media_name}" not in archive.namelist():
                raise RuntimeError(f"Missing embedded image: {media_name}")


def build_tables(vascular: pd.DataFrame, trajectory: pd.DataFrame):
    vascular_5 = vascular[vascular["grid"] == "5x5"].copy()
    vascular_rows = [["空间关系", "路线", "预期方向", "方向支持患者", "中位rho", "全局BH FDR"]]
    for row in vascular_5.itertuples(index=False):
        vascular_rows.append([
            RELATION_LABELS[row.relation], ROUTE_LABELS[row.route], row.expected_direction,
            f"{int(row.direction_supported_patients)}/{int(row.n_patients)}",
            f"{row.median_rho:.3f}", f"{row.sign_test_fdr_bh_global:.4g}",
        ])
    trajectory_rows = [["状态评分", "正向样本", "中位rho", "方向P", "BH FDR"]]
    for row in trajectory.itertuples(index=False):
        trajectory_rows.append([
            SCORE_LABELS[row.score], f"{int(row.positive_samples)}/{int(row.n_samples)}",
            f"{row.median_rho:.3f}", f"{row.sign_test_p:.4g}", f"{row.sign_test_fdr_bh:.4g}",
        ])
    curated_vascular = [["患者级空间关系", "尺度", "方向支持患者", "中位rho", "BH FDR"]]
    selected = vascular[
        ((vascular["relation"] == "SPP1_vascular_niche_to_T_cell_balance") & (vascular["route"] == "route_independent")) |
        ((vascular["relation"] == "TYMP_MTA_core_to_SPP1_vascular_niche") & (vascular["route"] == "reference_dependent"))
    ].sort_values(["relation", "grid"])
    for row in selected.itertuples(index=False):
        label = ("SPP1–vascular niche vs T-cell balance（预期负向）" if row.relation == "SPP1_vascular_niche_to_T_cell_balance"
                 else "TYMP–MTA core vs SPP1–vascular niche（reference-dependent）")
        curated_vascular.append([
            label, row.grid, f"{int(row.direction_supported_patients)}/{int(row.n_patients)}",
            f"{row.median_rho:.3f}", f"{row.sign_test_fdr_bh_global:.4g}",
        ])
    curated_trajectory = [["横断面DPT状态", "正向样本", "中位rho", "BH FDR", "解释"]]
    for score in ["Monocyte_program", "SPP1_state_no_SPP1", "C1QC_state", "SPP1_vs_C1QC"]:
        row = trajectory[trajectory["score"] == score].iloc[0]
        interpretation = {
            "Monocyte_program": "随DPT下降",
            "SPP1_state_no_SPP1": "弱正向、存在异质性",
            "C1QC_state": "最稳定上升状态",
            "SPP1_vs_C1QC": "多数样本向C1QC侧偏移",
        }[score]
        curated_trajectory.append([
            SCORE_LABELS[score], f"{int(row.positive_samples)}/{int(row.n_samples)}",
            f"{row.median_rho:.3f}", f"{row.sign_test_fdr_bh:.4g}", interpretation,
        ])
    return vascular_rows, trajectory_rows, curated_vascular, curated_trajectory


def build_full(vascular_rows, trajectory_rows) -> dict:
    document_root, rels_root, names, specs = relationship_image_specs(FULL_SOURCE, FIGURES, "vascular_trajectory_v3")
    body = document_root.find(q(W, "body"))
    if body is None:
        raise RuntimeError("Full report body missing")
    full_text = "\n".join(node_text(p) for p in document_root.iter(q(W, "p")))
    if "4.14 【TYMP机制状态分析】" not in full_text or "4.15 【血管生态位与髓系横断面状态排序扩展】" in full_text:
        raise RuntimeError("Unexpected full report section state")

    set_paragraph_text(
        direct_paragraph(body, "正式报告｜"),
        "正式报告｜GLMF1/GLMF2与TJ-RCC 12例全套空间分析｜TYMP血管生态位与髓系状态排序增强整合版 v3",
        bold=True,
    )
    set_paragraph_text(
        direct_paragraph(body, "【本次整合新增说明】"),
        "【本次整合新增说明】在既有4.12–4.14节基础上，新增4.15节TYMP–SPP1血管生态位与髓系横断面状态排序。新增分析覆盖GLMF1/GLMF2共7,079 spots、TJ-RCC 12例34,609 spots及单细胞参考中的53,794个髓系细胞；患者级空间关联和横断面DPT均保留非因果、非谱系追踪边界。",
        bold=True,
        shading="FFF2CC",
    )
    date_paragraphs = [p for p in body.iter(q(W, "p")) if "2026年8月3日（12例高级分析与TYMP状态更新）" in node_text(p)]
    if len(date_paragraphs) != 1:
        raise RuntimeError(f"Expected one full report date, found {len(date_paragraphs)}")
    for text_node in date_paragraphs[0].iter(q(W, "t")):
        if text_node.text and "2026年8月3日" in text_node.text:
            text_node.text = text_node.text.replace("2026年8月3日（12例高级分析与TYMP状态更新）", "2026年8月4日（血管生态位与髓系状态排序更新）")

    section = [
        page_break(),
        paragraph("4.15 【血管生态位与髓系横断面状态排序扩展】", style="Heading2"),
        paragraph("本节补充检验TYMP–SPP1方向是否伴随血管区室、T细胞平衡、缺氧/血管生成及侵袭程序，并以单细胞参考对象的髓系细胞进行横断面状态排序。空间部分覆盖GLMF1、GLMF2和TJ-RCC 12例；患者级汇总以12名TJ-RCC患者为生物学重复。DPT部分覆盖53,794个Monocyte、Other Macrophage、APOE+C1QC+ TAM和SPP1+ TAM细胞，63个样本达到样本内相关检验门槛。"),
        paragraph("4.15.1 TYMP–SPP1血管生态位的患者级空间关系", style="Heading3"),
        paragraph("在5×5主尺度的9条预设关系中，两条获得全局BH FDR<0.05。route-independent SPP1–vascular niche与T-cell balance在4×4和5×5均为12/12患者预期负向，中位rho分别为−0.456和−0.383；reference-dependent TYMP–MTA core与SPP1–vascular niche在两尺度均为12/12正向，中位rho为0.485和0.545；四项全局FDR均为0.00220。SPP1 TAM与血管区室、细胞毒性、缺氧/血管生成、侵袭/EMT及SPP1/monocyte balance未形成统一FDR支持。因此结果支持两条稳定空间关系，不支持把全部预设环节串联成完整机制链。"),
        table(vascular_rows, 14),
        *caption_figure(specs[0], "图4.50  TJ-RCC 12例九条关系–路线组合的5×5患者级空间分块效应。每点为一名患者，黑线为患者中位数；箭头仅表示预设检验方向，不表示因果作用。", width_emu=6_100_000),
        page_break(),
        *caption_figure(specs[1], "图4.51  TYMP–SPP1血管生态位预设关系的患者级方向汇总。只有SPP1–vascular niche vs T-cell balance和reference-dependent TYMP–MTA core vs SPP1–vascular niche达到全局FDR<0.05。", width_emu=6_100_000),
        paragraph("4.15.2 髓系横断面状态排序", style="Heading3"),
        paragraph("以单核细胞中心根、scVI 30近邻图、diffusion map和DPT构建横断面状态排序。C1QC-state在57/63样本中随DPT上升，中位rho=0.475、方向FDR=5.74×10⁻¹¹；不含SPP1基因本身的SPP1-state为41/63正向，中位rho=0.091、FDR=0.0316。SPP1-state minus C1QC-state仅11/63正向，中位rho=−0.276、FDR=3.90×10⁻⁷，说明排序后段整体更偏C1QC状态，而不是形成普遍的SPP1相对优势。"),
        table(trajectory_rows, 14),
        page_break(),
        *caption_figure(specs[2], "图4.52  53,794个髓系细胞的横断面状态排序。图示DPT、SPP1-state、C1QC-state、两状态差值及血管重塑评分；DPT不是RNA velocity或谱系示踪。", width_emu=6_100_000),
        page_break(),
        *caption_figure(specs[3], "图4.53  63个合格样本中七类状态评分与DPT的样本内Spearman方向。每点为一个样本，黑线为中位数；方向取决于所选单核细胞根。", width_emu=6_100_000),
        paragraph("4.15.3 QA与证据边界", style="Heading3"),
        paragraph("新增数值结果通过35/35独立检查，4张正式图通过4/4文件检查和4/4人工图检。血管生态位关系属于空间关联，不能证明免疫排斥、治疗耐药或直接通讯；当前数据无治疗响应、生存和纵向标签。DPT为横断面状态连续谱，不能证明单核细胞向SPP1+ TAM的谱系转化，也不能替代RNA velocity、谱系示踪或干预实验。"),
    ]
    heading5 = direct_paragraph(body, "5. 主要结论", style="Heading1")
    index = list(body).index(heading5)
    for node in section:
        body.insert(index, node)
        index += 1

    insert_after(body, direct_paragraph(body, "【TYMP机制状态结论】"), [
        paragraph("【血管生态位与髓系状态排序扩展结论】12例患者级证据稳定支持SPP1–vascular niche与T-cell balance负向关联，以及reference-dependent TYMP–MTA core与SPP1–vascular niche正向关联；两条关系在4×4和5×5均为12/12方向支持、全局FDR=0.00220。其余预设血管、缺氧、侵袭和SPP1/monocyte路线未形成统一支持。单细胞DPT中C1QC-state的上升比SPP1-state更稳定，故仅支持异质髓系状态连续谱，不支持普遍的单核细胞→SPP1+ TAM转化。"),
    ])
    insert_after(body, direct_paragraph(body, "【TYMP状态分析限制】"), [
        paragraph("【血管生态位与DPT限制】SPP1–vascular niche为三个细胞丰度分位数组成的复合指标，T-cell balance和表达程序虽经质量及空间多项式校正，仍不能排除残余组织区室混杂。DPT根与邻居图来自同一横断面参考对象，没有spliced/unspliced层、纵向取样、克隆追踪或疗效标签，因此不得用于谱系、进展或治疗耐药因果结论。"),
    ])

    status_tables = [node for node in body.iter(q(W, "tbl")) if "TJ-RCC 12例cell2location及全套下游分析" in node_text(node)]
    if len(status_tables) != 1:
        raise RuntimeError(f"Expected one status table, found {len(status_tables)}")
    status_tables[0].append(table_row([
        "TYMP血管生态位与髓系横断面状态排序", "已完成；35/35数值检查、4/4文件检查与4/4图检通过", "已纳入第4.15节和第5章"
    ], font_half_points=16))

    set_paragraph_text(
        any_paragraph(body, "最终总述"),
        "最终总述　第四部分现已形成四层证据：GLMF1/GLMF2用于空间发现和方法深化，TJ-RCC 12名独立患者用于全套高级方法及血管生态位关系重复，单细胞参考用于髓系横断面状态排序，客户湿实验作为TYMP–5′-MTA–SPP1机制的因果锚点。新增结果稳定支持SPP1–vascular niche与T-cell balance负向关系，以及reference-dependent TYMP–MTA core与SPP1–vascular niche正向关系；但缺氧、侵袭和SPP1/monocyte等其余环节未形成完整链式验证。DPT后段更稳定地偏向C1QC-state而非普遍SPP1相对优势。故报告支持人体ccRCC中的空间相容性与髓系状态异质性，不能把空间相关或横断面排序改写为已证明的直接通讯、谱系转化、治疗耐药或因果进展。",
        bold=True,
    )

    write_docx(FULL_SOURCE, FULL_OUTPUT, document_root, rels_root, names, specs)
    return {
        "source": str(FULL_SOURCE), "source_sha256": sha256(FULL_SOURCE),
        "candidate": str(FULL_OUTPUT), "candidate_sha256": sha256(FULL_OUTPUT),
        "candidate_size": FULL_OUTPUT.stat().st_size, "added_section": "4.15",
        "added_figures": ["4.50", "4.51", "4.52", "4.53"], "added_tables": 2,
    }


def replace_text_node(root: ET.Element, old: str, new: str) -> int:
    count = 0
    for node in root.iter(q(W, "t")):
        if node.text and old in node.text:
            node.text = node.text.replace(old, new)
            count += 1
    return count


def build_curated(curated_vascular, curated_trajectory) -> dict:
    curated_images = [FIGURES[1], FIGURES[3]]
    document_root, rels_root, names, specs = relationship_image_specs(CURATED_SOURCE, curated_images, "curated_v3")
    body = document_root.find(q(W, "body"))
    if body is None:
        raise RuntimeError("Curated report body missing")
    curated_text = "\n".join(node_text(p) for p in document_root.iter(q(W, "p")))
    if "3. 12名独立患者中的方向重复" not in curated_text or "3.2 血管生态位与髓系横断面状态排序" in curated_text:
        raise RuntimeError("Unexpected curated report section state")

    set_paragraph_text(direct_paragraph(body, "空间转录组图谱增强精选验证报告"), "空间转录组血管生态位与髓系状态增强精选验证报告", bold=True)
    if replace_text_node(document_root, "客户方向空间图谱增强精选版", "客户方向血管生态位与髓系状态增强精选版") != 1:
        raise RuntimeError("Curated version label replacement failed")
    if replace_text_node(document_root, " v2", " v3") != 1:
        raise RuntimeError("Curated v2 to v3 replacement failed")
    if replace_text_node(document_root, "报告日期：2026年8月3日", "报告日期：2026年8月4日") != 1:
        raise RuntimeError("Curated date replacement failed")

    exec_anchor = direct_paragraph(body, "在TJ-RCC 12名独立患者中")
    insert_after(body, exec_anchor, [
        paragraph("新增血管生态位检验显示，两条关系在4×4和5×5均获得患者级稳定支持：SPP1–vascular niche与T-cell balance为12/12预期负向，reference-dependent TYMP–MTA core与SPP1–vascular niche为12/12正向，四项全局BH FDR均为0.00220。其余血管、缺氧、侵袭和SPP1/monocyte路线未形成统一FDR支持。"),
        paragraph("单细胞参考的53,794个髓系细胞DPT排序中，C1QC-state为57/63样本正向，而排除SPP1本身的SPP1-state仅41/63正向，SPP1-state minus C1QC-state仅11/63正向。该结果支持异质髓系状态连续谱，不证明单核细胞向SPP1+ TAM的谱系转化。"),
    ])

    section = [
        page_break(),
        paragraph("3.2 血管生态位与髓系横断面状态排序", style="Heading2"),
        paragraph("患者级空间分析进一步将TYMP–SPP1方向放入血管区室和T细胞平衡背景中；单细胞参考DPT用于判断SPP1样与C1QC样状态在横断面连续谱中的相对方向。仅保留达到患者级全局FDR支持的空间关系及直接限制机制解释的DPT结果。"),
        table(curated_vascular, 16),
        *caption_figure(specs[0], "图3  TJ-RCC 12例TYMP–SPP1血管生态位关系的患者级方向汇总。两条关系在两种尺度均为12/12方向支持、全局FDR=0.00220；其余预设路线未获统一支持。", width_emu=5_950_000),
        paragraph("图3解释：SPP1–vascular niche为SPP1+ TAM、Endothelial和SMC/Pericyte丰度分位数组成的复合空间指标。与T-cell balance的负向关系支持空间区室中的T细胞相对不足，但无治疗响应或生存标签，不能称为免疫治疗耐药。reference-dependent TYMP–MTA core关系依赖既有细胞状态参考，不能作为独立因果验证。"),
        page_break(),
        table(curated_trajectory, 16),
        *caption_figure(specs[1], "图4  63个合格样本中髓系状态评分与横断面DPT的方向。C1QC-state上升最稳定；SPP1-state呈弱正向，而SPP1相对C1QC多数为负向。", width_emu=5_950_000),
        paragraph("图4解释：DPT使用单核细胞根、scVI邻居图和diffusion map。它不使用spliced/unspliced信息，也没有纵向或谱系示踪，因此只能描述横断面状态排序。图中方向不能写成单核细胞必然向SPP1+ TAM分化。"),
    ]
    heading4 = direct_paragraph(body, "4. 面向客户方向的结论", style="Heading1")
    index = list(body).index(heading4)
    for node in section:
        body.insert(index, node)
        index += 1

    insert_after(body, direct_paragraph(body, "4.1 可以由本轮空间分析支持的表述", style="Heading2"), [
        paragraph("新增证据支持两条患者级空间关系：SPP1–vascular niche与T-cell balance稳定负向，reference-dependent TYMP–MTA core与SPP1–vascular niche稳定正向；两者均跨4×4/5×5尺度重复。横断面DPT显示C1QC-state比SPP1-state更稳定上升，提示髓系状态重塑具有分支与异质性。"),
    ])
    insert_after(body, direct_paragraph(body, "4.2 证据边界", style="Heading2"), [
        paragraph("血管生态位空间关系不能单独证明直接细胞通讯、T细胞功能性排斥或治疗耐药；DPT不能证明单核细胞向SPP1+ TAM的谱系转化。当前没有RNA velocity、谱系示踪、治疗应答、生存或纵向进展标签。"),
    ])
    set_paragraph_text(
        direct_paragraph(body, "结合客户既有湿实验"),
        "结合客户既有湿实验，本空间转录组分析支持：TYMP相关恶性区域与MTA相关转录背景、SPP1⁺巨噬细胞及血管生态位在人体ccRCC组织中具有空间相容性；SPP1–vascular niche与T-cell balance负向关系及reference-dependent TYMP–MTA core与SPP1–vascular niche正向关系在12名患者和两种尺度中稳定重复。横断面DPT同时显示髓系状态并非普遍向SPP1相对优势推进，而更稳定地伴随C1QC-state上升。该结果为湿实验机制提供人组织空间相关性和患者异质性补充，但不替代5′-MTA实测、直接通讯、谱系示踪或干预因果证据。",
    )
    set_paragraph_text(
        direct_paragraph(body, "本报告只呈现"),
        "本报告呈现一张8-panel组织空间生态位图、一张6-panel患者空间图谱、两张血管生态位/髓系状态患者级汇总图及其直接支撑结果。全部12例均进入患者效应分析；完整9条血管关系、7类DPT评分、阴性结果、方法敏感性和全部审计文件保留在项目正式结果目录和完整Part 4报告中。精选版不替代完整技术报告。",
    )

    write_docx(CURATED_SOURCE, CURATED_OUTPUT, document_root, rels_root, names, specs)
    return {
        "source": str(CURATED_SOURCE), "source_sha256": sha256(CURATED_SOURCE),
        "candidate": str(CURATED_OUTPUT), "candidate_sha256": sha256(CURATED_OUTPUT),
        "candidate_size": CURATED_OUTPUT.stat().st_size, "added_section": "3.2",
        "added_figures": ["3", "4"], "added_tables": 2,
    }


def structural_audit(path: Path, expected_tables: int, expected_drawings: int, expected_media: int, markers: list[str]) -> dict:
    with ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise RuntimeError(f"CRC failure: {path}")
        root = ET.fromstring(archive.read("word/document.xml"))
        text = "\n".join(node_text(p) for p in root.iter(q(W, "p")))
        counts = {
            "paragraphs": len(list(root.iter(q(W, "p")))),
            "tables": len(list(root.iter(q(W, "tbl")))),
            "drawings": len(list(root.iter(q(W, "drawing")))),
            "media": len([name for name in archive.namelist() if name.startswith("word/media/")]),
        }
        if counts["tables"] != expected_tables or counts["drawings"] != expected_drawings or counts["media"] != expected_media:
            raise RuntimeError(f"Unexpected structure for {path.name}: {counts}")
        for marker in markers:
            if marker not in text:
                raise RuntimeError(f"Missing marker in {path.name}: {marker}")
        if "•" in text:
            raise RuntimeError(f"Forbidden bullet character in {path.name}")
        return counts


def main() -> None:
    STAGE.mkdir(parents=True, exist_ok=True)
    vascular, trajectory = validate_sources()
    tables = build_tables(vascular, trajectory)
    full = build_full(tables[0], tables[1])
    curated = build_curated(tables[2], tables[3])
    full["structure"] = structural_audit(
        FULL_OUTPUT, 31, 57, 57,
        ["4.15 【血管生态位与髓系横断面状态排序扩展】", "图4.50", "图4.53", "35/35", "12/12", "谱系转化"],
    )
    curated["structure"] = structural_audit(
        CURATED_OUTPUT, 7, 4, 4,
        ["3.2 血管生态位与髓系横断面状态排序", "图3", "图4", "0.00220", "57/63", "谱系转化"],
    )
    metadata = {
        "status": "staged_and_structurally_verified",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "full_report": full,
        "curated_report": curated,
        "source_tables": [
            {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)}
            for path in [VASCULAR, TRAJECTORY]
        ],
        "embedded_images": [
            {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in FIGURES
        ],
        "source_qc": {"numeric": "35/35", "figure_file": "4/4", "manual": "4/4"},
        "inference_boundary": "Spatial associations are not direct communication or treatment resistance; DPT is cross-sectional ordering, not lineage tracing.",
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
