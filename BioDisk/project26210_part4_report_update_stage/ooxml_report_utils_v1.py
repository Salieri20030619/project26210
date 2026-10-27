#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

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


def table_row(values: list[str], header: bool = False, font_half_points: int = 16) -> ET.Element:
    row = ET.Element(q(W, "tr"))
    row_props = ET.SubElement(row, q(W, "trPr"))
    ET.SubElement(row_props, q(W, "cantSplit"))
    if header:
        ET.SubElement(row_props, q(W, "tblHeader"), {q(W, "val"): "true"})
    for value in values:
        cell = ET.SubElement(row, q(W, "tc"))
        cell_props = ET.SubElement(cell, q(W, "tcPr"))
        if header:
            ET.SubElement(cell_props, q(W, "shd"), {q(W, "fill"): "D9EAF7"})
        p = ET.SubElement(cell, q(W, "p"))
        run = ET.SubElement(p, q(W, "r"))
        rprops = ET.SubElement(run, q(W, "rPr"))
        ET.SubElement(rprops, q(W, "sz"), {q(W, "val"): str(font_half_points)})
        ET.SubElement(rprops, q(W, "szCs"), {q(W, "val"): str(font_half_points)})
        if header:
            ET.SubElement(rprops, q(W, "b"))
        ET.SubElement(run, q(W, "t")).text = str(value)
    return row


def table(rows: list[list[str]], font_half_points: int = 16) -> ET.Element:
    node = ET.Element(q(W, "tbl"))
    props = ET.SubElement(node, q(W, "tblPr"))
    ET.SubElement(props, q(W, "tblW"), {q(W, "w"): "0", q(W, "type"): "auto"})
    borders = ET.SubElement(props, q(W, "tblBorders"))
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        ET.SubElement(borders, q(W, edge), {q(W, "val"): "single", q(W, "sz"): "6", q(W, "space"): "0", q(W, "color"): "808080"})
    for row_index, values in enumerate(rows):
        node.append(table_row(values, header=row_index == 0, font_half_points=font_half_points))
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


def caption_figure(spec: tuple[Path, str, str, int, int, int], caption: str, width_emu: int = 6_250_000) -> list[ET.Element]:
    _, rid, media_name, width, height, docpr_id = spec
    caption_node = paragraph(caption, bold=True)
    props = caption_node.find(q(W, "pPr"))
    if props is None:
        props = ET.Element(q(W, "pPr"))
        caption_node.insert(0, props)
    ET.SubElement(props, q(W, "keepNext"))
    ET.SubElement(props, q(W, "keepLines"))
    return [caption_node, image_paragraph(rid, media_name, width, height, docpr_id, width_emu=width_emu)]
