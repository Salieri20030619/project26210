#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_report_stage_v3")
FULL_SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v2.docx"
CURATED_SOURCE = ROOT / "docs/client_reports/part4_spatial/project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2.docx"
FULL = STAGE / "project26210_第四部分空间转录组阶段性报告_客户样式综合图版_POSTN_APOE论文对标增强版_v3_candidate.docx"
CURATED = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v3_candidate.docx"
METADATA = STAGE / "report_generation_metadata_v3.json"
OUT = STAGE / "independent_verification_v3.tsv"
SUMMARY = STAGE / "independent_verification_summary_v3.json"
VASCULAR = ROOT / "results/part4_spatial/tymp_vascular_niche_v1/tables/tjrcc_cohort_direction_summary_v1.tsv"
TRAJECTORY = ROOT / "results/part4_spatial/tymp_myeloid_state_ordering_v1/tables/cohort_pseudotime_direction_summary_v1.tsv"
FIG_DIR = ROOT / "results/part4_spatial/tymp_vascular_trajectory_synthesis_v1/figures"
FIGURES = [
    FIG_DIR / "vascular_niche_patient_effects_v1.png",
    FIG_DIR / "vascular_niche_cohort_summary_v1.png",
    FIG_DIR / "myeloid_cross_sectional_state_overview_v2.png",
    FIG_DIR / "myeloid_state_patient_directions_v1.png",
]

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def q(namespace: str, tag: str) -> str:
    return f"{{{namespace}}}{tag}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def node_text(node: ET.Element) -> str:
    return "".join(item.text or "" for item in node.iter(q(W, "t")))


def archive_info(path: Path) -> dict:
    with ZipFile(path) as archive:
        bad = archive.testzip()
        root = ET.fromstring(archive.read("word/document.xml"))
        rels = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
        media = {name: archive.read(name) for name in archive.namelist() if name.startswith("word/media/")}
        text = "\n".join(node_text(p) for p in root.iter(q(W, "p")))
        image_rel = {}
        for rel in rels.findall(q(REL, "Relationship")):
            if (rel.get("Type") or "").endswith("/image"):
                image_rel[rel.get("Id")] = "word/" + rel.get("Target").lstrip("/")
        embeds = [node.get(q(R, "embed")) for node in root.iter() if node.tag.endswith("}blip") and node.get(q(R, "embed"))]
        tables = list(root.iter(q(W, "tbl")))
        return {
            "bad": bad, "root": root, "rels": rels, "media": media, "text": text,
            "image_rel": image_rel, "embeds": embeds, "tables": tables,
            "drawings": list(root.iter(q(W, "drawing"))),
        }


def main() -> None:
    if OUT.exists() or SUMMARY.exists():
        raise SystemExit("Refusing to overwrite independent verification")
    gate = subprocess.run([str(ROOT / "check_memory_gate.sh")], capture_output=True, text=True, check=False)
    gate_line = (gate.stdout or gate.stderr).strip()
    if gate.returncode != 0:
        raise RuntimeError(gate_line)
    metadata = json.loads(METADATA.read_text())
    vascular = pd.read_csv(VASCULAR, sep="\t")
    trajectory = pd.read_csv(TRAJECTORY, sep="\t")
    checks = []

    def check(report: str, item: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({"report": report, "check": item, "status": "PASS" if passed else "FAIL",
                       "observed": str(observed), "expected": str(expected)})

    check("sources", "full v2 fingerprint", sha256(FULL_SOURCE) == metadata["full_report"]["source_sha256"], sha256(FULL_SOURCE), metadata["full_report"]["source_sha256"])
    check("sources", "curated v2 fingerprint", sha256(CURATED_SOURCE) == metadata["curated_report"]["source_sha256"], sha256(CURATED_SOURCE), metadata["curated_report"]["source_sha256"])
    check("candidates", "full candidate fingerprint", sha256(FULL) == metadata["full_report"]["candidate_sha256"], sha256(FULL), metadata["full_report"]["candidate_sha256"])
    check("candidates", "curated candidate fingerprint", sha256(CURATED) == metadata["curated_report"]["candidate_sha256"], sha256(CURATED), metadata["curated_report"]["candidate_sha256"])

    source_full = archive_info(FULL_SOURCE)
    source_curated = archive_info(CURATED_SOURCE)
    reports = {
        "full": (FULL, archive_info(FULL), source_full, 31, 57, 57),
        "curated": (CURATED, archive_info(CURATED), source_curated, 7, 4, 4),
    }
    for label, (path, info, source, tables, drawings, media) in reports.items():
        check(label, "ZIP CRC", info["bad"] is None, info["bad"], None)
        check(label, "table count", len(info["tables"]) == tables, len(info["tables"]), tables)
        check(label, "drawing count", len(info["drawings"]) == drawings, len(info["drawings"]), drawings)
        check(label, "media count", len(info["media"]) == media, len(info["media"]), media)
        missing_targets = [rid for rid in info["embeds"] if rid not in info["image_rel"] or info["image_rel"][rid] not in info["media"]]
        check(label, "image relationship targets", not missing_targets, missing_targets, [])
        check(label, "no forbidden bullet character", "•" not in info["text"], "present" if "•" in info["text"] else "absent", "absent")
        bordered = 0
        for table_node in info["tables"]:
            if table_node.find(f"{q(W, 'tblPr')}/{q(W, 'tblBorders')}") is not None:
                bordered += 1
        check(label, "all tables bordered", bordered == tables, bordered, tables)
        source_media = source["media"]
        unchanged = sum(name in info["media"] and info["media"][name] == payload for name, payload in source_media.items())
        check(label, "pre-existing media byte-identical", unchanged == len(source_media), unchanged, len(source_media))

    full_markers = {
        "version title": "增强整合版 v3",
        "date": "2026年8月4日（血管生态位与髓系状态排序更新）",
        "section": "4.15 【血管生态位与髓系横断面状态排序扩展】",
        "figure 4.50": "图4.50",
        "figure 4.53": "图4.53",
        "numeric QC": "35/35独立检查",
        "T-cell result": "中位rho分别为−0.456和−0.383",
        "vascular result": "中位rho为0.485和0.545",
        "C1QC result": "57/63样本中随DPT上升",
        "SPP1-vs-C1QC result": "仅11/63正向",
        "lineage boundary": "不能证明单核细胞向SPP1+ TAM的谱系转化",
        "chapter5 conclusion": "【血管生态位与髓系状态排序扩展结论】",
        "chapter5 limitation": "【血管生态位与DPT限制】",
        "status row": "已纳入第4.15节和第5章",
    }
    for item, marker in full_markers.items():
        check("full", item, marker in reports["full"][1]["text"], "found" if marker in reports["full"][1]["text"] else "missing", "found")

    curated_markers = {
        "version title": "血管生态位与髓系状态增强精选验证报告",
        "version": "血管生态位与髓系状态增强精选版 v3",
        "date": "报告日期：2026年8月4日",
        "section": "3.2 血管生态位与髓系横断面状态排序",
        "figure 3": "图3  TJ-RCC 12例TYMP–SPP1血管生态位关系",
        "figure 4": "图4  63个合格样本中髓系状态评分",
        "global FDR": "全局BH FDR均为0.00220",
        "C1QC result": "C1QC-state为57/63样本正向",
        "lineage boundary": "不证明单核细胞向SPP1+ TAM的谱系转化",
        "four figures": "两张血管生态位/髓系状态患者级汇总图",
    }
    for item, marker in curated_markers.items():
        check("curated", item, marker in reports["curated"][1]["text"], "found" if marker in reports["curated"][1]["text"] else "missing", "found")

    # New embedded image bytes must exactly match the formal source PNGs.
    full_new = {name: payload for name, payload in reports["full"][1]["media"].items() if name not in source_full["media"]}
    curated_new = {name: payload for name, payload in reports["curated"][1]["media"].items() if name not in source_curated["media"]}
    formal_hashes = {sha256(path) for path in FIGURES}
    full_hashes = {hashlib.sha256(payload).hexdigest() for payload in full_new.values()}
    curated_expected = {sha256(FIGURES[1]), sha256(FIGURES[3])}
    curated_hashes = {hashlib.sha256(payload).hexdigest() for payload in curated_new.values()}
    check("full", "four formal images embedded byte-identically", full_hashes == formal_hashes and len(full_new) == 4, sorted(full_hashes), sorted(formal_hashes))
    check("curated", "two selected images embedded byte-identically", curated_hashes == curated_expected and len(curated_new) == 2, sorted(curated_hashes), sorted(curated_expected))

    # Recompute the key numbers represented in both reports from formal tables.
    key = vascular[
        (((vascular.relation == "SPP1_vascular_niche_to_T_cell_balance") & (vascular.route == "route_independent")) |
         ((vascular.relation == "TYMP_MTA_core_to_SPP1_vascular_niche") & (vascular.route == "reference_dependent")))
    ]
    key_ok = len(key) == 4 and (key.direction_supported_patients == 12).all() and np.allclose(key.sign_test_fdr_bh_global, 0.002197265625)
    check("sources", "key vascular recomputation", key_ok, f"rows={len(key)}; supported={key.direction_supported_patients.tolist()}; fdr={key.sign_test_fdr_bh_global.tolist()}", "4 rows; all 12/12; all 0.002197265625")
    t = trajectory.set_index("score")
    t_ok = int(t.loc["C1QC_state", "positive_samples"]) == 57 and int(t.loc["SPP1_state_no_SPP1", "positive_samples"]) == 41 and int(t.loc["SPP1_vs_C1QC", "positive_samples"]) == 11
    check("sources", "trajectory direction recomputation", t_ok, {score: int(t.loc[score, "positive_samples"]) for score in ["C1QC_state", "SPP1_state_no_SPP1", "SPP1_vs_C1QC"]}, {"C1QC_state": 57, "SPP1_state_no_SPP1": 41, "SPP1_vs_C1QC": 11})

    # Pandoc must independently reverse-parse both DOCX candidates.
    for label, path in [("full", FULL), ("curated", CURATED)]:
        result = subprocess.run(["pandoc", str(path), "-t", "plain"], capture_output=True, text=True, check=False)
        parsed_ok = result.returncode == 0 and "谱系转化" in result.stdout and "0.00220" in result.stdout
        check(label, "Pandoc reverse parse", parsed_ok, f"exit={result.returncode}; chars={len(result.stdout)}", "exit=0; required markers")

    frame = pd.DataFrame(checks)
    frame.to_csv(OUT, sep="\t", index=False)
    passed = int(frame.status.eq("PASS").sum())
    failed = int(frame.status.eq("FAIL").sum())
    summary = {
        "status": "passed" if failed == 0 else "failed",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "memory_gate": gate_line,
        "checks_passed": passed,
        "checks_failed": failed,
        "full_sha256": sha256(FULL),
        "curated_sha256": sha256(CURATED),
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(frame.to_string(index=False))
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
