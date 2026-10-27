#!/usr/bin/env python3
"""Independently verify the staged Part 4 advanced-patch DOCX."""

from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_report_stage")
DOCX = STAGE / "project26210_part4_advanced_patch_v1_candidate.docx"
META = STAGE / "project26210_part4_advanced_patch_v1_report_metadata.json"
OUT_JSON = STAGE / "project26210_part4_advanced_patch_v1_verification.json"
OUT_TSV = STAGE / "project26210_part4_advanced_patch_v1_verification.tsv"
MISTY_PATH = ROOT / "results/part4_spatial/advanced_patch_misty_cross_sample_v3/cross_sample_misty_performance_v3.tsv"
COMMOT_PATH = ROOT / "results/part4_spatial/advanced_patch_commot_cross_sample_v3/commot_cross_sample_primary_hypothesis_support_v3.tsv"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
REL = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text(node: ET.Element) -> str:
    return "".join(x.text or "" for x in node.iter(W + "t"))


def table_rows(node: ET.Element) -> list[list[str]]:
    return [[text(cell) for cell in row.findall(W + "tc")] for row in node.findall(W + "tr")]


checks: list[dict[str, object]] = []


def check(name: str, passed: bool, observed: object, expected: object) -> None:
    checks.append({"check": name, "passed": bool(passed), "observed": str(observed), "expected": str(expected)})


for path in (DOCX, META, MISTY_PATH, COMMOT_PATH):
    check(f"file_exists::{path.name}", path.is_file(), path.is_file(), True)
if not all(bool(x["passed"]) for x in checks):
    raise SystemExit("Required input missing")

meta = json.loads(META.read_text())
check("docx_sha256_matches_metadata", sha256(DOCX) == meta["output_sha256"], sha256(DOCX), meta["output_sha256"])

with ZipFile(DOCX) as archive:
    bad = archive.testzip()
    names = archive.namelist()
    document = ET.fromstring(archive.read("word/document.xml"))
    rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
    check("zip_crc", bad is None, bad, None)
    check("zip_entries_unique", len(names) == len(set(names)), len(names), len(set(names)))
    check("zip_entry_count", len(names) == 60, len(names), 60)
    check("media_count", sum(x.startswith("word/media/") for x in names) == 41, sum(x.startswith("word/media/") for x in names), 41)
    check("document_xml_parse", True, "parsed", "parsed")
    check("relationships_xml_parse", True, "parsed", "parsed")

    rel_map = {x.get("Id"): x.get("Target") for x in rels.findall(REL + "Relationship")}
    embeds = [x.get(R + "embed") for x in document.iter("{http://schemas.openxmlformats.org/drawingml/2006/main}blip")]
    unresolved = [rid for rid in embeds if rid not in rel_map]
    missing_targets = []
    for rid in embeds:
        target = rel_map.get(rid, "")
        if target.startswith("media/") and f"word/{target}" not in names:
            missing_targets.append((rid, target))
    check("all_image_relationships_resolved", not unresolved, unresolved, [])
    check("all_embedded_media_present", not missing_targets, missing_targets, [])

    metadata_images = {x["embedded_name"]: x for x in meta["added_images"]}
    new_media = sorted(x for x in names if x.startswith("word/media/advanced_patch_v1_"))
    check("advanced_patch_media_count", len(new_media) == 8, len(new_media), 8)
    for member in new_media:
        name = Path(member).name
        observed = hashlib.sha256(archive.read(member)).hexdigest()
        expected = metadata_images[name]["sha256"]
        check(f"embedded_hash::{name}", observed == expected, observed, expected)

body = document.find(W + "body")
assert body is not None
body_nodes = list(body)
body_text = [text(x).strip() for x in body_nodes]
all_text = "\n".join(text(x) for x in document.iter(W + "p"))
check("table_count", len(list(document.iter(W + "tbl"))) == 23, len(list(document.iter(W + "tbl"))), 23)
check("drawing_count", len(list(document.iter(W + "drawing"))) == 41, len(list(document.iter(W + "drawing"))), 41)

advanced_positions = [i for i, value in enumerate(body_text) if value.startswith("4.12 【高级分析补丁】")]
heading5_positions = []
for i, node in enumerate(body_nodes):
    if node.tag != W + "p" or not body_text[i].startswith("5. 主要结论"):
        continue
    style = node.find(f"{W}pPr/{W}pStyle")
    if style is not None and style.get(W + "val") == "Heading1":
        heading5_positions.append(i)
check("advanced_section_unique", len(advanced_positions) == 1, advanced_positions, "one position")
check("chapter5_heading_unique", len(heading5_positions) == 1, heading5_positions, "one position")
check("advanced_before_chapter5", advanced_positions[0] < heading5_positions[0], (advanced_positions[0], heading5_positions[0]), "advanced < chapter5")

required_markers = [
    "BayesSpace在两切片均选择q=4",
    "共同双阳性空间基因4,242个",
    "MISTy两例共同正增益17/20",
    "两例名义显著仅R² 3/20、RMSE 2/20",
    "COMMOT预设轴两例均正8/12",
    "严格支持0/12",
    "不能写成17个目标获得患者级显著复现",
    "0/12同时进入两例196种代理组合的传输质量前20%",
    "【高级分析补丁结论】",
    "【高级分析补丁限制】",
]
for marker in required_markers:
    check(f"marker::{marker[:28]}", marker in all_text, marker in all_text, True)
for number in range(30, 38):
    marker = f"图4.{number}  "
    check(f"figure_caption::{marker}", all_text.count(marker) == 1, all_text.count(marker), 1)

tables = [table_rows(x) for x in document.iter(W + "tbl")]
misty_tables = [x for x in tables if x and x[0] and x[0][0] == "目标" and "ΔR²（G1/G2）" in x[0]]
commot_tables = [x for x in tables if x and x[0] and x[0][0] == "轴（尺度）"]
check("misty_table_unique", len(misty_tables) == 1, len(misty_tables), 1)
check("commot_table_unique", len(commot_tables) == 1, len(commot_tables), 1)
if misty_tables:
    rows = misty_tables[0]
    check("misty_table_rows", len(rows) == 21, len(rows), 21)
    check("misty_both_positive_count", sum(x[3] == "是" for x in rows[1:]) == 17, sum(x[3] == "是" for x in rows[1:]), 17)
    check("misty_r2_nominal_both_count", sum(x[4] == "是" for x in rows[1:]) == 3, sum(x[4] == "是" for x in rows[1:]), 3)
    check("misty_rmse_nominal_both_count", sum(x[5] == "是" for x in rows[1:]) == 2, sum(x[5] == "是" for x in rows[1:]), 2)
if commot_tables:
    rows = commot_tables[0]
    check("commot_table_rows", len(rows) == 13, len(rows), 13)
    check("commot_both_positive_count", sum(x[4] == "是" for x in rows[1:]) == 8, sum(x[4] == "是" for x in rows[1:]), 8)
    check("commot_strict_count", sum(x[5] == "是" for x in rows[1:]) == 0, sum(x[5] == "是" for x in rows[1:]), 0)

check("misty_source_rows", len(pd.read_csv(MISTY_PATH, sep="\t")) == 20, len(pd.read_csv(MISTY_PATH, sep="\t")), 20)
check("commot_source_rows", len(pd.read_csv(COMMOT_PATH, sep="\t")) == 12, len(pd.read_csv(COMMOT_PATH, sep="\t")), 12)

failed = [x for x in checks if not x["passed"]]
summary = {
    "status": "PASS" if not failed else "FAIL",
    "checks_passed": len(checks) - len(failed),
    "checks_total": len(checks),
    "checks_failed": len(failed),
    "docx": str(DOCX),
    "docx_size": DOCX.stat().st_size,
    "docx_sha256": sha256(DOCX),
    "failed_checks": failed,
}
OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
pd.DataFrame(checks).to_csv(OUT_TSV, sep="\t", index=False)
print(json.dumps(summary, ensure_ascii=False, indent=2))
if failed:
    raise SystemExit(1)
