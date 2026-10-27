#!/usr/bin/env python3
"""Correct the report date after the advanced-patch content was finalized."""

from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


SOURCE = Path("/mnt/BioDisk/ludingyu/project26210/docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
OUTPUT = STAGE / "project26210_part4_advanced_patch_v1_datefix_candidate.docx"
AUDIT = STAGE / "project26210_part4_advanced_patch_v1_datefix_audit.json"
EXPECTED_SOURCE = "646b1abf7d7580e9770e60021bc26634bace6055f2b969e0cbb0b6c13919ec0d"
OLD = "2026年7月31日（本版更新）"
NEW = "2026年8月2日（高级分析补丁更新）"
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
ET.register_namespace("w", W)


def q(tag: str) -> str:
    return f"{{{W}}}{tag}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text(node: ET.Element) -> str:
    return "".join(x.text or "" for x in node.iter(q("t")))


if sha256(SOURCE) != EXPECTED_SOURCE:
    raise RuntimeError("Official report changed before the date correction")
with ZipFile(SOURCE) as source_zip:
    names = source_zip.namelist()
    root = ET.fromstring(source_zip.read("word/document.xml"))
    before_text = [text(x) for x in root.iter(q("p"))]
    matches = [x for x in root.iter(q("p")) if text(x).strip() == OLD]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one old report-date paragraph, found {len(matches)}")
    node = matches[0]
    props = node.find(q("pPr"))
    for child in list(node):
        if child is not props:
            node.remove(child)
    run = ET.SubElement(node, q("r"))
    ET.SubElement(run, q("t")).text = NEW
    after_text = [text(x) for x in root.iter(q("p"))]
    differences = [(i, a, b) for i, (a, b) in enumerate(zip(before_text, after_text)) if a != b]
    if differences != [(before_text.index(OLD), OLD, NEW)]:
        raise RuntimeError(f"Unexpected paragraph changes: {differences}")
    document_bytes = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    with ZipFile(OUTPUT, "w", ZIP_DEFLATED, allowZip64=True) as output_zip:
        for name in names:
            if name == "word/document.xml":
                output_zip.writestr(name, document_bytes)
            else:
                output_zip.writestr(source_zip.getinfo(name), source_zip.read(name))

with ZipFile(SOURCE) as before_zip, ZipFile(OUTPUT) as after_zip:
    if after_zip.testzip() is not None:
        raise RuntimeError("Date-corrected candidate failed ZIP test")
    if before_zip.namelist() != after_zip.namelist():
        raise RuntimeError("DOCX member list changed during date correction")
    changed_members = []
    for name in before_zip.namelist():
        if hashlib.sha256(before_zip.read(name)).digest() != hashlib.sha256(after_zip.read(name)).digest():
            changed_members.append(name)
    if changed_members != ["word/document.xml"]:
        raise RuntimeError(f"Unexpected changed DOCX members: {changed_members}")
    final_root = ET.fromstring(after_zip.read("word/document.xml"))
    final_text = "\n".join(text(x) for x in final_root.iter(q("p")))
    if OLD in final_text or final_text.count(NEW) != 1:
        raise RuntimeError("Final report-date text check failed")
    for marker in ["4.12 【高级分析补丁】", "图4.37", "【高级分析补丁结论】", "【高级分析补丁限制】"]:
        if marker not in final_text:
            raise RuntimeError(f"Existing advanced-patch marker missing after date correction: {marker}")

audit = {
    "status": "PASS",
    "source": str(SOURCE),
    "source_size": SOURCE.stat().st_size,
    "source_sha256": sha256(SOURCE),
    "candidate": str(OUTPUT),
    "candidate_size": OUTPUT.stat().st_size,
    "candidate_sha256": sha256(OUTPUT),
    "old_text": OLD,
    "new_text": NEW,
    "changed_members": ["word/document.xml"],
    "changed_paragraphs": 1,
    "zip_test": "PASS",
}
AUDIT.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(audit, ensure_ascii=False, indent=2))
