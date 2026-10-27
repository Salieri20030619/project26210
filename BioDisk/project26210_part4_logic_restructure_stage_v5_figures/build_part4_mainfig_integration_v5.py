#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


UTILS = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_update_stage")
sys.path.insert(0, str(UTILS))
from ooxml_report_utils_v1 import A, REL, W, WP, caption_figure, node_text, png_size, q  # noqa: E402


SOURCE = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v4/"
    "project26210_第四部分空间转录组分析报告_主线逻辑重构版_v4_candidate.docx"
)
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v5_figures")
OUTPUT = STAGE / "project26210_第四部分空间转录组分析报告_主线逻辑重构版_三张主图整合_v5_candidate.docx"
METADATA = STAGE / "report_generation_metadata_v5.json"
FIGROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_three_figures_final_v4")
FIGURES = [
    {
        "number": "4.1",
        "path": FIGROOT / "Figure1_TYMP_MTA_SPP1_evidence_narrowing_v4.png",
        "media": "part4_main_figure_4_1_v4.png",
        "anchor": "主评分以GPNMB、LGALS3、FABP5、CD9和LPL构成",
        "reference": "整体证据链见图4.1。",
        "caption": (
            "图4.1  TYMP–MTA–SPP1巨噬细胞假说的空间证据逐级收窄。"
            "A，发现、患者验证与严格审计设计；B，GLMF2双变量空间结构；"
            "C，GLMF1/GLMF2跨尺度发现效应；D，12例患者中无SPP1状态与SPP1相对C1QC效应；"
            "E，从空间相容性到MTA中介和定向通讯的证据分级。空间数据支持人组织背景，"
            "但不独立证明5′-MTA中介或直接通讯。"
        ),
    },
    {
        "number": "4.2",
        "path": FIGROOT / "Figure2_SPP1_vascular_niche_Tcell_balance_v4.png",
        "media": "part4_main_figure_4_2_v4.png",
        "anchor": "在5×5主尺度的9条预设关系中，两条获得全局BH FDR<0.05",
        "reference": "队列空间关系及尺度稳定性见图4.2。",
        "caption": (
            "图4.2  患者重复的SPP1–vascular niche与T-cell balance空间关系。"
            "A，变量构建；B，代表性患者Z43_T双变量空间图；C，两个分块尺度的患者效应；"
            "D，九条预设关系–路线组合的队列审计；E，计算型肿瘤界面两侧的评分分布；"
            "F，患者级方向检验支持的关联网络。上述为关联结构，不证明直接通讯、免疫排斥或治疗耐药。"
        ),
    },
    {
        "number": "4.3",
        "path": FIGROOT / "Figure3_myeloid_state_continuum_v4.png",
        "media": "part4_main_figure_4_3_v4.png",
        "anchor": "以单核细胞中心根、scVI 30近邻图、diffusion map和DPT构建横断面状态排序",
        "reference": "髓系状态连续谱的综合结果见图4.3。",
        "caption": (
            "图4.3  髓系细胞形成以C1QC状态增强为主的横断面连续谱，而非已证明的SPP1谱系。"
            "A，四类髓系状态UMAP；B，各状态的DPT分布；C，样本内标准化后的状态评分轨迹；"
            "D，共享近邻图连通性；E，63个合格样本的方向一致性；F，观察结果与不可外推结论。"
            "DPT不替代RNA velocity、纵向取样或谱系示踪。"
        ),
    },
]
EXPECTED_SOURCE_SHA256 = "980f6f54a1d6e0af851ed568b8d248b7ba022b6db2b477ae04bf2772df4d5083"
EXPECTED_FIGURE_SHA256 = {
    "4.1": "cbf75c1996fc79d02d116807b4bfac8be760fc6c09a80d951007f6d2e2dbf8cd",
    "4.2": "b320e12041593f40ab92646e49c0bb2457d8708199bf32d479f1c1ae39ee1811",
    "4.3": "f709bd890227ffebe0c91e7cb390913a192e46c5aa17ae7998d4840d71756950",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def replace_all_text(root: ET.Element, old: str, new: str) -> int:
    count = 0
    for text_node in root.iter(q(W, "t")):
        if text_node.text and old in text_node.text:
            count += text_node.text.count(old)
            text_node.text = text_node.text.replace(old, new)
    return count


def relabel_existing_figures(root: ET.Element) -> dict[str, int]:
    chinese = re.compile(r"(?<!补充)图\s*4\.(\d+)([A-Z]?)")
    english = re.compile(r"(?<!Supplementary )Figure\s*4\.(\d+)([A-Z]?)", re.IGNORECASE)
    chinese_count = 0
    english_count = 0
    for text_node in root.iter(q(W, "t")):
        if not text_node.text:
            continue
        text_node.text, n_cn = chinese.subn(lambda m: f"补充图 S4.{m.group(1)}{m.group(2)}", text_node.text)
        text_node.text, n_en = english.subn(lambda m: f"Supplementary Figure S4.{m.group(1)}{m.group(2)}", text_node.text)
        chinese_count += n_cn
        english_count += n_en
    if chinese_count < 54:
        raise RuntimeError(f"Too few existing figure labels relabeled: {chinese_count}")
    return {"chinese_replacements": chinese_count, "english_replacements": english_count}


def append_reference(anchor: ET.Element, reference: str) -> None:
    texts = list(anchor.iter(q(W, "t")))
    if not texts:
        raise RuntimeError(f"Anchor has no text: {node_text(anchor)}")
    if reference not in node_text(anchor):
        texts[-1].text = (texts[-1].text or "") + reference


def find_direct_paragraph(body: ET.Element, prefix: str) -> ET.Element:
    matches = [node for node in list(body)
               if node.tag == q(W, "p") and node_text(node).strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one direct paragraph for {prefix!r}, found {len(matches)}")
    return matches[0]


def add_image_spec(document_root: ET.Element, rels_root: ET.Element, figure: dict, offset: int):
    rid_numbers = []
    for rel in rels_root.findall(q(REL, "Relationship")):
        rid = rel.get("Id", "")
        if rid.startswith("rId") and rid[3:].isdigit():
            rid_numbers.append(int(rid[3:]))
    rid = f"rId{max(rid_numbers, default=0) + 1}"
    ET.SubElement(rels_root, q(REL, "Relationship"), {
        "Id": rid,
        "Type": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
        "Target": f"media/{figure['media']}",
    })
    docpr_ids = [int(node.get("id")) for node in document_root.iter(q(WP, "docPr"))
                 if (node.get("id") or "").isdigit()]
    width, height = png_size(figure["path"])
    return figure["path"], rid, figure["media"], width, height, max(docpr_ids, default=0) + 1 + offset


def insert_after(body: ET.Element, anchor: ET.Element, nodes: list[ET.Element]) -> None:
    index = list(body).index(anchor) + 1
    for node in nodes:
        body.insert(index, node)
        index += 1


def main() -> None:
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("Source report fingerprint changed")
    for figure in FIGURES:
        if sha256(figure["path"]) != EXPECTED_FIGURE_SHA256[figure["number"]]:
            raise RuntimeError(f"Figure fingerprint changed: {figure['path']}")
    STAGE.mkdir(parents=True, exist_ok=True)
    for path in (OUTPUT, METADATA):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite existing output: {path}")

    with ZipFile(SOURCE) as archive:
        source_names = archive.namelist()
        document_root = ET.fromstring(archive.read("word/document.xml"))
        rels_root = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
    body = document_root.find(q(W, "body"))
    if body is None:
        raise RuntimeError("Missing document body")

    relabel_counts = relabel_existing_figures(document_root)
    replace_all_text(document_root, "主线逻辑重构版最终分析报告", "主线逻辑重构版最终分析报告（主图整合版）")
    replace_all_text(document_root, "逻辑重构版 v4", "逻辑重构版 v5｜三张主图整合")
    replace_all_text(
        document_root,
        "【重构原则】保留原v3全部53张分析图与关键数值；仅重排叙事、统一章节和图号，并新增一张总流程图。未新增未经验证的生物学结果。",
        "【主图规则】正文新增3张主题组合主图；原报告全部分析图保持原位置并统一标注为补充图S4.x。未改变原图数据或生物学结论。",
    )
    replace_all_text(
        document_root,
        "全部原始图表仍保留",
        "全部原始图表仍在原位置保留并标注为补充图",
    )

    inserted = []
    media_payloads = []
    for offset, figure in enumerate(FIGURES):
        anchor = find_direct_paragraph(body, figure["anchor"])
        append_reference(anchor, figure["reference"])
        spec = add_image_spec(document_root, rels_root, figure, offset)
        nodes = caption_figure(spec, figure["caption"], width_emu=6_250_000)
        insert_after(body, anchor, nodes)
        media_payloads.append((figure["media"], figure["path"]))
        inserted.append({"number": figure["number"], "anchor": figure["anchor"], "media": figure["media"]})

    full_text = "\n".join(node_text(node) for node in body.iter(q(W, "p")))
    main_caption_count = sum(full_text.count(f"图{number}  ") for number in ("4.1", "4.2", "4.3"))
    if main_caption_count != 3:
        raise RuntimeError(f"Expected 3 main captions, found {main_caption_count}")
    supplement_captions = [node_text(node).strip() for node in list(body)
                           if node.tag == q(W, "p") and node_text(node).strip().startswith("补充图 S4.")]
    if len(supplement_captions) != 55:
        raise RuntimeError(f"Expected 55 existing captions marked supplementary, found {len(supplement_captions)}")

    shanghai = timezone(timedelta(hours=8))
    now = datetime.now(shanghai)
    with ZipFile(SOURCE) as source_zip:
        core_root = ET.fromstring(source_zip.read("docProps/core.xml"))
        dc = "http://purl.org/dc/elements/1.1/"
        dcterms = "http://purl.org/dc/terms/"
        title_node = core_root.find(q(dc, "title"))
        if title_node is not None:
            title_node.text = "project26210 第四部分空间转录组分析报告（主线逻辑重构版·三张主图整合v5）"
        modified = core_root.find(q(dcterms, "modified"))
        if modified is not None:
            modified.text = now.isoformat(timespec="seconds")
        replacements = {
            "word/document.xml": ET.tostring(document_root, encoding="utf-8", xml_declaration=True),
            "word/_rels/document.xml.rels": ET.tostring(rels_root, encoding="utf-8", xml_declaration=True),
            "docProps/core.xml": ET.tostring(core_root, encoding="utf-8", xml_declaration=True),
        }
        with ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
            for name in source_names:
                if name in replacements:
                    output_zip.writestr(name, replacements[name])
                else:
                    output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))
            for media_name, media_path in media_payloads:
                output_zip.writestr(f"word/media/{media_name}", media_path.read_bytes())

    with ZipFile(OUTPUT) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("DOCX ZIP CRC failure")

    metadata = {
        "status": "generated_pending_render_qc",
        "generated_at": now.isoformat(timespec="seconds"),
        "source": str(SOURCE),
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
        "output_size": OUTPUT.stat().st_size,
        "main_figures": inserted,
        "main_figure_count": 3,
        "existing_caption_count_marked_supplementary": len(supplement_captions),
        "existing_figure_positions": "unchanged; labels and in-text references only",
        "relabel_counts": relabel_counts,
        "figure_hashes": EXPECTED_FIGURE_SHA256,
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
