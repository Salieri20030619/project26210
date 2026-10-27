#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_curated_report_stage_v2")
MARKDOWN = STAGE / "project26210_tymp_curated_report_v2.md"
LUA = STAGE / "pagebreak.lua"
RAW = STAGE / "project26210_tymp_curated_report_v2_raw.docx"
OUTPUT = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2_candidate.docx"
METADATA = STAGE / "curated_report_generation_metadata_v2.json"
FIG1 = STAGE / "Figure1_TYMP_MTA_SPP1_spatial_niche_8panel_v2.png"
FIG2 = STAGE / "Figure2_TYMP_SPP1_patient_spatial_atlas_6panel_v2.png"

IMMUNE = ROOT / "results/part4_spatial/tymp_immune_ecology_v2/tables/original_tymp_immune_spatial_tests.tsv"
SIGNATURE = ROOT / "results/part4_spatial/tymp_immune_ecology_v2/tables/original_tymp_signature_spatial_tests.tsv"
COHORT_STATE = ROOT / "results/part4_spatial/tymp_spp1_state_v1/tables/cohort_state_direction_summary_v1.tsv"
COHORT_DOWN = ROOT / "results/part4_spatial/tymp_spp1_state_v1/tables/cohort_interface_downstream_summary_v1.tsv"
SELECTED = STAGE / "selected_cohort_evidence_v2.tsv"
PANEL_MANIFEST = STAGE / "panel_source_manifest_v2.tsv"

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
CP = "http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
DC = "http://purl.org/dc/elements/1.1/"
DCTERMS = "http://purl.org/dc/terms/"
XSI = "http://www.w3.org/2001/XMLSchema-instance"
for prefix, uri in {
    "w": W,
    "r": R,
    "wp": WP,
    "a": A,
    "pic": PIC,
    "cp": CP,
    "dc": DC,
    "dcterms": DCTERMS,
    "xsi": XSI,
}.items():
    ET.register_namespace(prefix, uri)


def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text_of(node: ET.Element) -> str:
    return "".join(item.text or "" for item in node.iter(q(W, "t"))).strip()


def child(parent: ET.Element, tag: str) -> ET.Element:
    item = parent.find(q(W, tag))
    if item is None:
        item = ET.SubElement(parent, q(W, tag))
    return item


def paragraph_properties(paragraph: ET.Element) -> ET.Element:
    props = paragraph.find(q(W, "pPr"))
    if props is None:
        props = ET.Element(q(W, "pPr"))
        paragraph.insert(0, props)
    return props


def run_properties(run: ET.Element) -> ET.Element:
    props = run.find(q(W, "rPr"))
    if props is None:
        props = ET.Element(q(W, "rPr"))
        run.insert(0, props)
    return props


def replace_property(parent: ET.Element, tag: str, attributes: dict[str, str]) -> ET.Element:
    old = parent.find(q(W, tag))
    if old is not None:
        parent.remove(old)
    return ET.SubElement(parent, q(W, tag), {q(W, key): value for key, value in attributes.items()})


def format_run(run: ET.Element, *, size: int, bold: bool = False, color: str = "222222") -> None:
    props = run_properties(run)
    replace_property(props, "rFonts", {"ascii": "Arial", "hAnsi": "Arial", "eastAsia": "Microsoft YaHei", "cs": "Arial"})
    replace_property(props, "sz", {"val": str(size)})
    replace_property(props, "szCs", {"val": str(size)})
    replace_property(props, "color", {"val": color})
    existing_bold = props.find(q(W, "b"))
    if bold and existing_bold is None:
        ET.SubElement(props, q(W, "b"))
    elif not bold and existing_bold is not None:
        props.remove(existing_bold)


def format_paragraph(
    paragraph: ET.Element,
    *,
    size: int = 21,
    bold: bool = False,
    color: str = "222222",
    center: bool = False,
    after: int = 90,
    line: int = 300,
    keep_next: bool = False,
    shading: str | None = None,
) -> None:
    props = paragraph_properties(paragraph)
    replace_property(props, "spacing", {"after": str(after), "line": str(line), "lineRule": "auto"})
    if center:
        replace_property(props, "jc", {"val": "center"})
    if keep_next and props.find(q(W, "keepNext")) is None:
        ET.SubElement(props, q(W, "keepNext"))
    if shading:
        replace_property(props, "shd", {"fill": shading})
    for run in paragraph.findall(q(W, "r")):
        format_run(run, size=size, bold=bold, color=color)


for required in [MARKDOWN, LUA, FIG1, FIG2, IMMUNE, SIGNATURE, COHORT_STATE, COHORT_DOWN, SELECTED, PANEL_MANIFEST]:
    if not required.is_file():
        raise FileNotFoundError(required)
if RAW.exists() or OUTPUT.exists():
    raise FileExistsError("v2 raw or candidate already exists; choose a new version")

subprocess.run(
    [
        "pandoc",
        str(MARKDOWN),
        "--from=markdown+raw_tex",
        "--to=docx",
        f"--lua-filter={LUA}",
        "--resource-path",
        str(STAGE),
        "-o",
        str(RAW),
    ],
    check=True,
)

with ZipFile(RAW) as archive:
    names = archive.namelist()
    if archive.testzip() is not None:
        raise RuntimeError("Raw DOCX CRC failure")
    document = ET.fromstring(archive.read("word/document.xml"))
    core = ET.fromstring(archive.read("docProps/core.xml"))

body = document.find(q(W, "body"))
if body is None:
    raise RuntimeError("Missing document body")

# A4 portrait with 0.70-inch side margins, leaving 6.87 inches for the two master figures.
section = body.find(q(W, "sectPr"))
if section is None:
    raise RuntimeError("Missing section properties")
replace_property(section, "pgSz", {"w": "11907", "h": "16840"})
replace_property(section, "pgMar", {"top": "900", "right": "1008", "bottom": "900", "left": "1008", "header": "360", "footer": "360", "gutter": "0"})

cover_formats = {
    "TYMP–5′-MTA–SPP1⁺巨噬细胞生态位机制": dict(size=38, bold=True, color="17365D", center=True, after=180, line=320),
    "空间转录组精选验证报告": dict(size=30, bold=True, color="1F4E78", center=True, after=180, line=320),
    "围绕客户湿实验方向的人组织空间相关性与患者重复证据": dict(size=22, bold=True, color="404040", center=True, after=180, line=300),
}

for paragraph in body.iter(q(W, "p")):
    value = text_of(paragraph)
    props = paragraph_properties(paragraph)
    style = props.find(q(W, "pStyle"))
    style_id = style.get(q(W, "val"), "") if style is not None else ""
    in_table = any(ancestor.tag == q(W, "tc") for ancestor in [])
    if value in cover_formats:
        format_paragraph(paragraph, **cover_formats[value])
    elif value.startswith(("数据集：", "报告版本：", "报告日期：", "报告状态：")):
        format_paragraph(paragraph, size=20, center=True, after=40, line=280)
    elif value.startswith("证据定位："):
        format_paragraph(paragraph, size=20, bold=True, color="7F6000", center=True, after=100, line=280, shading="FFF2CC")
    elif style_id in {"Title", "title"}:
        format_paragraph(paragraph, size=38, bold=True, color="17365D", center=True, after=180, line=320)
    elif style_id in {"Subtitle", "subtitle"}:
        format_paragraph(paragraph, size=28, bold=True, color="1F4E78", center=True, after=160, line=300)
    elif style_id in {"Heading1", "heading 1", "1"} or value in {"执行摘要"} or value.startswith(("1. ", "2. ", "3. ", "4. ", "5. ")):
        format_paragraph(paragraph, size=30, bold=True, color="1F4E78", after=130, line=300, keep_next=True)
    elif style_id in {"Heading2", "heading 2", "2"} or value.startswith(("1.1 ", "1.2 ", "2.1 ", "3.1 ", "4.1 ", "4.2 ", "4.3 ")):
        format_paragraph(paragraph, size=24, bold=True, color="2F75B5", after=100, line=290, keep_next=True)
    elif style_id.lower().replace(" ", "") in {"imagecaption", "caption"} or value.startswith(("图1", "图2")):
        format_paragraph(paragraph, size=18, bold=True, color="333333", center=True, after=80, line=260, keep_next=True)
    elif paragraph.find(q(W, "r") + "/" + q(W, "drawing")) is not None or paragraph.find(".//" + q(WP, "inline")) is not None:
        format_paragraph(paragraph, size=18, center=True, after=80, line=240)
    else:
        format_paragraph(paragraph, size=21, bold=False, color="222222", after=90, line=300)

# Enlarge both composite figures to 6.80 inches and keep their aspect ratios.
target_width = int(6.80 * 914400)
for inline in body.iter(q(WP, "inline")):
    extent = inline.find(q(WP, "extent"))
    if extent is None:
        continue
    old_width = int(extent.get("cx", "1"))
    old_height = int(extent.get("cy", "1"))
    target_height = round(old_height * target_width / old_width)
    extent.set("cx", str(target_width))
    extent.set("cy", str(target_height))
    transform_extent = inline.find(".//" + q(A, "xfrm") + "/" + q(A, "ext"))
    if transform_extent is not None:
        transform_extent.set("cx", str(target_width))
        transform_extent.set("cy", str(target_height))

tables = list(body.iter(q(W, "tbl")))
for table in tables:
    props = table.find(q(W, "tblPr"))
    if props is None:
        props = ET.Element(q(W, "tblPr"))
        table.insert(0, props)
    replace_property(props, "tblW", {"w": "5000", "type": "pct"})
    replace_property(props, "tblLayout", {"type": "autofit"})
    old_borders = props.find(q(W, "tblBorders"))
    if old_borders is not None:
        props.remove(old_borders)
    borders = ET.SubElement(props, q(W, "tblBorders"))
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        ET.SubElement(borders, q(W, edge), {q(W, "val"): "single", q(W, "sz"): "5", q(W, "space"): "0", q(W, "color"): "A6A6A6"})
    rows = table.findall(q(W, "tr"))
    for row_index, row in enumerate(rows):
        row_props = row.find(q(W, "trPr"))
        if row_props is None:
            row_props = ET.Element(q(W, "trPr"))
            row.insert(0, row_props)
        if row_props.find(q(W, "cantSplit")) is None:
            ET.SubElement(row_props, q(W, "cantSplit"))
        if row_index == 0 and row_props.find(q(W, "tblHeader")) is None:
            ET.SubElement(row_props, q(W, "tblHeader"))
        for cell in row.findall(q(W, "tc")):
            cell_props = cell.find(q(W, "tcPr"))
            if cell_props is None:
                cell_props = ET.Element(q(W, "tcPr"))
                cell.insert(0, cell_props)
            replace_property(cell_props, "vAlign", {"val": "center"})
            if row_index == 0:
                replace_property(cell_props, "shd", {"fill": "D9EAF7"})
            for paragraph in cell.findall(q(W, "p")):
                format_paragraph(paragraph, size=18, bold=row_index == 0, color="17365D" if row_index == 0 else "222222", after=35, line=250)

for tag, value in {
    "title": "TYMP–5′-MTA–SPP1⁺巨噬细胞生态位机制空间转录组精选验证报告",
    "subject": "客户湿实验方向的人组织空间相关性与患者重复证据",
}.items():
    node = core.find(q(DC, tag))
    if node is None:
        node = ET.SubElement(core, q(DC, tag))
    node.text = value

document_bytes = ET.tostring(document, encoding="utf-8", xml_declaration=True)
core_bytes = ET.tostring(core, encoding="utf-8", xml_declaration=True)
with ZipFile(RAW) as source, ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as destination:
    for name in names:
        if name == "word/document.xml":
            destination.writestr(name, document_bytes)
        elif name == "docProps/core.xml":
            destination.writestr(name, core_bytes)
        else:
            destination.writestr(source.getinfo(name), source.read(name))

with ZipFile(OUTPUT) as archive:
    if archive.testzip() is not None:
        raise RuntimeError("Final DOCX CRC failure")
    final_document = ET.fromstring(archive.read("word/document.xml"))
    report_text = "".join(item.text or "" for item in final_document.iter(q(W, "t")))
    drawings = list(final_document.iter(q(W, "drawing")))
    final_tables = list(final_document.iter(q(W, "tbl")))
    media = [name for name in archive.namelist() if name.startswith("word/media/")]
    if len(drawings) != 2 or len(final_tables) != 5 or len(media) != 2:
        raise RuntimeError(f"Unexpected structure: drawings={len(drawings)}, tables={len(final_tables)}, media={len(media)}")
    for marker in ["执行摘要", "图1", "图2", "推荐客户表述", "10/12", "0.154", "0.193"]:
        if marker not in report_text:
            raise RuntimeError(f"Missing report marker: {marker}")

metadata = {
    "status": "staged",
    "generator": str(Path(__file__)),
    "report": str(OUTPUT),
    "report_size": OUTPUT.stat().st_size,
    "report_sha256": sha256(OUTPUT),
    "figures": [
        {"path": str(FIG1), "sha256": sha256(FIG1), "panels": 8},
        {"path": str(FIG2), "sha256": sha256(FIG2), "panels": 6},
    ],
    "source_tables": [
        {"path": str(path), "sha256": sha256(path)}
        for path in [IMMUNE, SIGNATURE, COHORT_STATE, COHORT_DOWN, SELECTED, PANEL_MANIFEST]
    ],
    "docx_structure": {"tables": len(final_tables), "drawings": len(drawings), "media": len(media)},
    "layout": {"page": "A4 portrait", "side_margin_in": 0.70, "figure_width_in": 6.80},
    "scope": "Positive, mechanism-aligned spatial evidence only; explicit inference boundary retained.",
}
METADATA.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(metadata, ensure_ascii=False, indent=2))
