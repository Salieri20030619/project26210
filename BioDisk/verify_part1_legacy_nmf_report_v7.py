#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree as ET


DEFAULT_REPORT = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/"
    "05_nmf_result_report_v1/NMF_client_results_report_v7.docx"
)
REPORT = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_REPORT
OUT_DIR = Path("/mnt/BioDisk/wangmingyuan/part1_legacy_nmf_report_v7_qc_staging")
QC_TSV = OUT_DIR / "NMF_client_results_report_v7_verification.tsv"
METADATA_JSON = OUT_DIR / "NMF_client_results_report_v7_verification_metadata.json"

SOURCE_IMAGES = [
    Path("/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/06_nmf_customer_requested_redraw_v1/figures/Figure_NMF_customer_requested_main_v1.png"),
    Path("/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/08_nmf_customer_tree_style_v3_1/figures/Figure_NMF_clone_state_trees_customer_style_v3_1.png"),
    Path("/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/06_nmf_customer_requested_redraw_v1/figures/C_metaprogram_expression_heatmap_reference_style_customer_v1.png"),
    Path("/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/06_nmf_customer_requested_redraw_v1/figures/E_normalized_gene_contribution_pale_blue_customer_v1.png"),
    Path("/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/04_nmf_sensitivity_v1/figures/NMF_sensitivity_overview_v1.png"),
]

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PR = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"w": W, "a": A, "r": R, "pr": PR}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def add(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append(
            {
                "check": name,
                "pass": "true" if passed else "false",
                "observed": str(observed),
                "expected": str(expected),
            }
        )

    add("report_exists", REPORT.is_file(), REPORT.is_file(), True)
    add("report_size_gt_1mb", REPORT.stat().st_size > 1_000_000, REPORT.stat().st_size, ">1000000")
    source_missing = [str(path) for path in SOURCE_IMAGES if not path.is_file()]
    add("all_five_source_images_exist", not source_missing, len(SOURCE_IMAGES) - len(source_missing), 5)

    with zipfile.ZipFile(REPORT) as archive:
        names = archive.namelist()
        add("zip_crc", archive.testzip() is None, archive.testzip(), None)
        required_parts = {
            "[Content_Types].xml",
            "word/document.xml",
            "word/_rels/document.xml.rels",
        }
        add("required_ooxml_parts", required_parts.issubset(names), len(required_parts & set(names)), len(required_parts))

        xml_names = [name for name in names if name.endswith((".xml", ".rels"))]
        xml_failures = []
        for name in xml_names:
            try:
                ET.fromstring(archive.read(name))
            except ET.ParseError:
                xml_failures.append(name)
        add("all_xml_parse", not xml_failures, len(xml_names) - len(xml_failures), len(xml_names))

        document = ET.fromstring(archive.read("word/document.xml"))
        rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        text = "".join(node.text or "" for node in document.findall(".//w:t", NS))
        table_count = len(document.findall(".//w:tbl", NS))
        drawing_count = len(document.findall(".//w:drawing", NS))
        blips = document.findall(".//a:blip", NS)
        embed_ids = [node.attrib.get(f"{{{R}}}embed", "") for node in blips]
        page_breaks = [
            node for node in document.findall(".//w:br", NS)
            if node.attrib.get(f"{{{W}}}type") == "page"
        ]

        rel_map = {
            rel.attrib["Id"]: rel.attrib.get("Target", "")
            for rel in rels.findall("pr:Relationship", NS)
        }
        image_targets = [rel_map.get(rid, "") for rid in embed_ids]
        normalized_media = [
            "word/" + target.lstrip("/") if not target.startswith("/") else target.lstrip("/")
            for target in image_targets
        ]
        missing_targets = [name for name in normalized_media if name not in names]
        embedded_hashes = {
            sha256_bytes(archive.read(name)) for name in normalized_media if name in names
        }
        source_hashes = {sha256_file(path) for path in SOURCE_IMAGES if path.is_file()}

        add("table_count", table_count == 4, table_count, 4)
        add("drawing_count", drawing_count == 5, drawing_count, 5)
        add("image_reference_count", len(embed_ids) == 5, len(embed_ids), 5)
        add("unique_image_relationships", len(set(embed_ids)) == 5, len(set(embed_ids)), 5)
        add("all_image_relationships_resolve", not missing_targets, len(normalized_media) - len(missing_targets), len(normalized_media))
        add("embedded_images_match_sources", embedded_hashes == source_hashes, len(embedded_hashes & source_hashes), 5)
        add("explicit_page_breaks", len(page_breaks) >= 2, len(page_breaks), ">=2")
        add("no_solid_bullet_character", "•" not in text, text.count("•"), 0)

        required_text = [
            "恶性上皮细胞 cNMF 分析结果说明书",
            "29,338 个 CopyKAT aneuploid 恶性上皮细胞",
            "Factor-TYMP 候选程序",
            "rank 5 组件轮廓系数",
            "MP5 主导细胞的最大单研究占比",
            "当前数据尚不支持的结论",
            "尚不能将 MP5 定义为已经验证的跨队列稳定 Factor-TYMP 亚型",
            "聚类继续使用 5 维 cNMF 平均使用度、Euclidean 距离和 average linkage",
            "不代表系统发育时间",
            "每个可见叶节点和分叉节点均使用字母编号",
            "同一分叉处的入边与两条出边使用互不相同的颜色",
            "不可将图上的几何枝长直接当作原始 Euclidean 距离",
            "不是系统发育树",
            "不混用补充细胞前仅覆盖部分细胞的 metaVIPER 结果",
            "所有场景中，候选主导细胞的最大单研究占比均为 89.2%–100%",
            "去掉 GSE207493 后",
            "细胞数 | 6,944 | 29,338",
        ]
        for phrase in required_text:
            compact_phrase = phrase.replace(" | ", "")
            present = phrase in text or compact_phrase in text
            add("required_text::" + phrase[:42], present, present, True)

    all_passed = all(row["pass"] == "true" for row in checks)
    with QC_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "pass", "observed", "expected"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)

    metadata = {
        "status": "PASS" if all_passed else "FAIL",
        "checked_at": datetime.now().astimezone().isoformat(),
        "report": str(REPORT),
        "report_bytes": REPORT.stat().st_size,
        "report_sha256": sha256_file(REPORT),
        "checks_passed": sum(row["pass"] == "true" for row in checks),
        "checks_total": len(checks),
        "source_images": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in SOURCE_IMAGES if path.is_file()
        ],
        "qc_tsv": str(QC_TSV),
    }
    METADATA_JSON.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"NMF_REPORT_V7_{metadata['status']} checks={metadata['checks_passed']}/{metadata['checks_total']}")
    print(f"report_sha256={metadata['report_sha256']}")
    if not all_passed:
        for row in checks:
            if row["pass"] != "true":
                print("FAIL", row["check"], row["observed"], row["expected"], sep="\t", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
