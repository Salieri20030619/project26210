#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


SOURCE = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part4_logic_restructure_stage_v5_figures/"
    "project26210_第四部分空间转录组分析报告_主线逻辑重构版_三张主图整合_v5_candidate.docx"
)
STAGE = SOURCE.parent
OUTPUT = STAGE / "project26210_第四部分空间转录组分析报告_仅三张主图精简版_v5_candidate.docx"
METADATA = STAGE / "compact_report_generation_metadata_v5.json"
EXPECTED_SOURCE_SHA256 = "fc371d99c36d0f68b892d524bb594577115e940648518794c465d29ddea1a64a"
MAIN_MEDIA = {
    "media/part4_main_figure_4_1_v4.png",
    "media/part4_main_figure_4_2_v4.png",
    "media/part4_main_figure_4_3_v4.png",
}
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"
for prefix, uri in {
    "w": W[1:-1],
    "r": R[1:-1],
    "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
    "a": A[1:-1],
    "pic": "http://schemas.openxmlformats.org/drawingml/2006/picture",
}.items():
    ET.register_namespace(prefix, uri)
ET.register_namespace("", REL[1:-1])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def node_text(node: ET.Element) -> str:
    return "".join(text.text or "" for text in node.iter(W + "t"))


def replace_all_text(root: ET.Element, old: str, new: str) -> int:
    count = 0
    for text_node in root.iter(W + "t"):
        if text_node.text and old in text_node.text:
            count += text_node.text.count(old)
            text_node.text = text_node.text.replace(old, new)
    return count


def clear_paragraph_content(node: ET.Element) -> int:
    if node.tag != W + "p":
        raise RuntimeError(f"Expected paragraph block, got {node.tag}")
    removed_page_break_before = 0
    props = node.find(W + "pPr")
    if props is not None:
        for page_break_before in list(props.findall(W + "pageBreakBefore")):
            props.remove(page_break_before)
            removed_page_break_before += 1
    for child in list(node):
        if child.tag != W + "pPr":
            node.remove(child)
    return removed_page_break_before


def has_explicit_page_break(node: ET.Element) -> bool:
    return any(br.get(W + "type") == "page" for br in node.iter(W + "br"))


def main() -> None:
    if sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise RuntimeError("Integrated report fingerprint changed")
    for path in (OUTPUT, METADATA):
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite existing output: {path}")

    with ZipFile(SOURCE) as archive:
        source_names = archive.namelist()
        document_root = ET.fromstring(archive.read("word/document.xml"))
        rels_root = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
    body = document_root.find(W + "body")
    if body is None:
        raise RuntimeError("Missing document body")

    rel_targets = {rel.get("Id"): rel.get("Target") for rel in rels_root.findall(REL + "Relationship")}
    removed_supplement_captions = 0
    removed_supplement_drawings = 0
    removed_subpanel_labels = 0
    cleared_page_break_before = 0
    subpanel_labels = {
        "A. GLMF1：spot 级 QC 指标分布",
        "B. GLMF2：spot 级 QC 指标分布",
        "C. GLMF1：QC 指标的组织空间分布",
        "D. GLMF2：QC 指标的组织空间分布",
    }
    for node in list(body):
        text = node_text(node).strip()
        if node.tag == W + "p" and text.startswith("补充图 S4."):
            cleared_page_break_before += clear_paragraph_content(node)
            removed_supplement_captions += 1
            continue
        if node.tag == W + "p" and text in subpanel_labels:
            cleared_page_break_before += clear_paragraph_content(node)
            removed_subpanel_labels += 1
            continue
        blips = list(node.iter(A + "blip"))
        if blips:
            targets = {rel_targets.get(blip.get(R + "embed"), "") for blip in blips}
            if not targets.issubset(MAIN_MEDIA):
                cleared_page_break_before += clear_paragraph_content(node)
                removed_supplement_drawings += len(blips)

    redundant_explicit_page_breaks = 0
    meaningful_since_break = True
    for node in list(body):
        if has_explicit_page_break(node):
            if not meaningful_since_break:
                clear_paragraph_content(node)
                redundant_explicit_page_breaks += 1
            else:
                meaningful_since_break = False
            continue
        if node_text(node).strip() or node.tag == W + "tbl" or list(node.iter(W + "drawing")):
            meaningful_since_break = True

    if removed_supplement_captions != 55:
        raise RuntimeError(f"Expected 55 supplementary captions removed, got {removed_supplement_captions}")
    if removed_supplement_drawings != 58:
        raise RuntimeError(f"Expected 58 supplementary drawings removed, got {removed_supplement_drawings}")
    if removed_subpanel_labels != 4:
        raise RuntimeError(f"Expected 4 composite subpanel labels removed, got {removed_subpanel_labels}")

    replace_all_text(document_root, "主线逻辑重构版最终分析报告（主图整合版）",
                     "主线逻辑重构版最终分析报告（仅三张主图精简版）")
    replace_all_text(document_root, "逻辑重构版 v5｜三张主图整合",
                     "逻辑重构版 v5｜仅三张主图精简版")
    replace_all_text(
        document_root,
        "【主图规则】正文新增3张主题组合主图；原报告全部分析图保持原位置并统一标注为补充图S4.x。未改变原图数据或生物学结论。",
        "【精简版规则】正文仅保留图4.1–图4.3三张主题组合主图；原始单项分析图未在本版本嵌入，完整图集见主图整合完整版。未改变数据、表格或生物学结论。",
    )
    replace_all_text(
        document_root,
        "正文展示8张队列图，逐病例图保留在正式结果目录作为补充审计材料。",
        "主图整合完整版展示8张队列图；本精简版仅保留3张组合主图，逐病例图仍在正式结果目录作为审计材料。",
    )
    replace_all_text(
        document_root,
        "全部原始图表仍在原位置保留并标注为补充图",
        "本精简版仅嵌入3张组合主图；完整原始图集保留在主图整合完整版",
    )

    used_image_rids = {blip.get(R + "embed") for blip in document_root.iter(A + "blip")}
    removed_media_targets = set()
    for rel in list(rels_root.findall(REL + "Relationship")):
        if rel.get("Type", "").endswith("/image") and rel.get("Id") not in used_image_rids:
            removed_media_targets.add("word/" + rel.get("Target", ""))
            rels_root.remove(rel)
    remaining_targets = {rel.get("Target") for rel in rels_root.findall(REL + "Relationship")
                         if rel.get("Type", "").endswith("/image")}
    if remaining_targets != MAIN_MEDIA:
        raise RuntimeError(f"Unexpected remaining image relationships: {remaining_targets}")
    if len(list(document_root.iter(W + "drawing"))) != 3:
        raise RuntimeError("Compact report must contain exactly three drawings")

    shanghai = timezone(timedelta(hours=8))
    now = datetime.now(shanghai)
    with ZipFile(SOURCE) as source_zip:
        core_root = ET.fromstring(source_zip.read("docProps/core.xml"))
        dc = "http://purl.org/dc/elements/1.1/"
        dcterms = "http://purl.org/dc/terms/"
        title_node = core_root.find("{" + dc + "}title")
        if title_node is not None:
            title_node.text = "project26210 第四部分空间转录组分析报告（仅三张主图精简版v5）"
        modified = core_root.find("{" + dcterms + "}modified")
        if modified is not None:
            modified.text = now.isoformat(timespec="seconds")
        replacements = {
            "word/document.xml": ET.tostring(document_root, encoding="utf-8", xml_declaration=True),
            "word/_rels/document.xml.rels": ET.tostring(rels_root, encoding="utf-8", xml_declaration=True),
            "docProps/core.xml": ET.tostring(core_root, encoding="utf-8", xml_declaration=True),
        }
        with ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
            for name in source_names:
                if name in removed_media_targets:
                    continue
                if name in replacements:
                    output_zip.writestr(name, replacements[name])
                else:
                    output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))

    with ZipFile(OUTPUT) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("DOCX ZIP CRC failure")
        remaining_media = {name for name in archive.namelist() if name.startswith("word/media/")}
    expected_media = {"word/" + target for target in MAIN_MEDIA}
    if remaining_media != expected_media:
        raise RuntimeError(f"Unexpected compact media payload: {remaining_media}")

    metadata = {
        "status": "generated_pending_render_qc",
        "generated_at": now.isoformat(timespec="seconds"),
        "source": str(SOURCE),
        "source_sha256": EXPECTED_SOURCE_SHA256,
        "output": str(OUTPUT),
        "output_sha256": sha256(OUTPUT),
        "output_size": OUTPUT.stat().st_size,
        "main_figures_retained": 3,
        "supplementary_captions_removed": removed_supplement_captions,
        "supplementary_drawings_removed": removed_supplement_drawings,
        "composite_subpanel_labels_removed": removed_subpanel_labels,
        "cleared_page_break_before_from_removed_figures": cleared_page_break_before,
        "redundant_explicit_page_breaks_removed": redundant_explicit_page_breaks,
        "visible_drawings": 3,
        "embedded_media_assets": len(remaining_media),
        "content_scope": "full logical narrative and tables retained; only three composite main figures embedded",
    }
    METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
