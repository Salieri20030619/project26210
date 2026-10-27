#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
STAGE = Path("/mnt/BioDisk/wangmingyuan/project26210_tymp_curated_report_stage_v2")
REPORT = STAGE / "project26210_TYMP_5MTA_SPP1机制空间验证精选版_v2_candidate.docx"
META = STAGE / "curated_report_generation_metadata_v2.json"
FIG_META = STAGE / "figure_build_metadata_v2.json"
PANEL = STAGE / "panel_source_manifest_v2.tsv"
SELECTED = STAGE / "selected_cohort_evidence_v2.tsv"
EFFECT_MATRIX = STAGE / "selected_patient_effect_matrix_v2.tsv"
REPRESENTATIVE = STAGE / "representative_patient_selection_v2.tsv"
IMMUNE = ROOT / "results/part4_spatial/tymp_immune_ecology_v2/tables/original_tymp_immune_spatial_tests.tsv"
SIGNATURE = ROOT / "results/part4_spatial/tymp_immune_ecology_v2/tables/original_tymp_signature_spatial_tests.tsv"
OUT = STAGE / "curated_report_independent_verification_v2.tsv"
SUMMARY = STAGE / "curated_report_independent_verification_summary_v2.json"

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"

def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

checks: list[tuple[str, bool, str]] = []
def check(name: str, condition: bool, detail: str) -> None:
    checks.append((name, bool(condition), detail))

for path in [REPORT, META, FIG_META, PANEL, SELECTED, EFFECT_MATRIX, REPRESENTATIVE, IMMUNE, SIGNATURE]:
    check(f"exists_{path.name}", path.is_file(), str(path))

meta = json.loads(META.read_text(encoding="utf-8"))
fig_meta = json.loads(FIG_META.read_text(encoding="utf-8"))
check("report_hash_matches_metadata", sha256(REPORT) == meta["report_sha256"], sha256(REPORT))
check("report_size_matches_metadata", REPORT.stat().st_size == meta["report_size"], str(REPORT.stat().st_size))
check("figure1_eight_panels", fig_meta["figure_1"]["panels"] == 8, str(fig_meta["figure_1"]["panels"]))
check("figure2_six_panels", fig_meta["figure_2"]["panels"] == 6, str(fig_meta["figure_2"]["panels"]))

with ZipFile(REPORT) as archive:
    names = archive.namelist()
    check("zip_crc", archive.testzip() is None, str(archive.testzip()))
    check("zip_unique_entries", len(names) == len(set(names)), str(len(names)))
    media = [name for name in names if name.startswith("word/media/")]
    check("embedded_media_count", len(media) == 2, str(len(media)))
    embedded_hashes = {sha256_bytes(archive.read(name)) for name in media}
    expected_hashes = {fig_meta["figure_1"]["sha256"], fig_meta["figure_2"]["sha256"]}
    check("embedded_figure_hashes", embedded_hashes == expected_hashes, str(embedded_hashes))
    root = ET.fromstring(archive.read("word/document.xml"))
    document_xml = archive.read("word/document.xml").decode("utf-8")
    text = "".join(t.text or "" for t in root.iter(q(W, "t")))
    tables = list(root.iter(q(W, "tbl")))
    drawings = list(root.iter(q(W, "drawing")))
    check("table_count", len(tables) == 5, str(len(tables)))
    check("drawing_count", len(drawings) == 2, str(len(drawings)))
    all_cells = list(root.iter(q(W, "tc")))
    empty_cells = [cell for cell in all_cells if not "".join(t.text or "" for t in cell.iter(q(W, "t"))).strip()]
    check("all_table_cells_have_text", len(empty_cells) == 0, f"empty={len(empty_cells)};total={len(all_cells)}")
    bordered = sum(table.find(f"{q(W, 'tblPr')}/{q(W, 'tblBorders')}") is not None for table in tables)
    check("explicit_table_borders", bordered == 5, f"{bordered}/5")
    section = root.find(f".//{q(W, 'sectPr')}")
    page_size = section.find(q(W, "pgSz")) if section is not None else None
    check("a4_page_size", page_size is not None and page_size.get(q(W, "w")) == "11907" and page_size.get(q(W, "h")) == "16840", str(page_size.attrib if page_size is not None else None))
    extents = list(root.iter("{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}extent"))
    check("two_full_width_figures", len(extents) == 2 and all(abs(int(item.get("cx")) / 914400 - 6.8) < 0.001 for item in extents), str([item.attrib for item in extents]))
    check("no_unicode_bullet", "•" not in text, "forbidden bullet absent")
    check("no_list_numbering", "<w:numPr" not in document_xml, "numPr absent")
    for marker in [
        "TYMP–5′-MTA–SPP1⁺巨噬细胞生态位机制", "执行摘要", "1. 分析范围与证据规则",
        "2. 发现样本中的机制链空间相容性", "3. 12名独立患者中的方向重复",
        "4. 面向客户方向的结论", "推荐客户表述", "5. 精选版交付说明",
        "0.0162", "0.0423", "10/12", "0.154", "9/12", "0.305", "0.193",
    ]:
        check(f"marker_{marker[:20]}", marker in text, marker)
    check("figure1_caption_unique", text.count("图1 GLMF1/GLMF2中") == 1, str(text.count("图1 GLMF1/GLMF2中")))
    check("figure2_caption_unique", text.count("图2 TJ-RCC 12例中") == 1, str(text.count("图2 TJ-RCC 12例中")))

panel = pd.read_csv(PANEL, sep="\t")
check("panel_manifest_14", len(panel) == 14, str(len(panel)))
for figure, expected_letters in [(1, set("ABCDEFGH")), (2, set("ABCDEF"))]:
    letters = set(panel.loc[panel.figure == figure, "panel"])
    check(f"figure{figure}_panel_letters", letters == expected_letters, "".join(sorted(letters)))
check("panel_sources_exist", all(Path(path.split(";", 1)[0]).is_file() for path in panel.source), "all first source paths")

effect_matrix = pd.read_csv(EFFECT_MATRIX, sep="\t", index_col=0)
check("effect_matrix_12_by_7", effect_matrix.shape == (12, 7), str(effect_matrix.shape))
representative = pd.read_csv(REPRESENTATIVE, sep="\t")
check("representative_is_x49", len(representative) == 12 and representative.iloc[0]["sample"] == "X49_T", representative.head(1).to_json(orient="records"))

immune = pd.read_csv(IMMUNE, sep="\t")
signature = pd.read_csv(SIGNATURE, sep="\t")
def one(frame: pd.DataFrame, **filters) -> pd.Series:
    mask = pd.Series(True, index=frame.index)
    for key, value in filters.items(): mask &= frame[key].eq(value)
    rows = frame[mask]
    if len(rows) != 1: raise RuntimeError(f"Source recomputation failed for {filters}: {len(rows)} rows")
    return rows.iloc[0]

for sample, expected_rho, expected_fdr in [("GLMF1", 0.7352941176470588, 0.0162), ("GLMF2", 0.6470588235294118, 0.0423)]:
    row = one(signature, sample=sample, source="TYMP-High Malignant", target="Polyamine-MTA context", scale="spatial_block_permutation", grid="4x4")
    check(f"{sample}_tymp_mta", abs(row.rho - expected_rho) < 1e-12 and abs(row.fdr_bh - expected_fdr) < 5e-5, f"rho={row.rho};fdr={row.fdr_bh}")
for sample, expected_fdr in [("GLMF1", 0.0072), ("GLMF2", 0.0054)]:
    row = one(signature, sample=sample, source="SPP1+ TAM", target="Polyamine-MTA context", scale="spatial_block_permutation", grid="5x5")
    check(f"{sample}_spp1_mta", abs(row.fdr_bh - expected_fdr) < 1e-12, f"rho={row.rho};fdr={row.fdr_bh}")
for sample, expected_fdr in [("GLMF1", 0.0086666666666666), ("GLMF2", 0.0011428571428571)]:
    row = one(immune, sample=sample, metric="absolute_q05", source="TYMP-High Malignant", target="SPP1+ TAM", scale="spatial_block_permutation", grid="5x5")
    check(f"{sample}_tymp_spp1_absolute", abs(row.fdr_bh - expected_fdr) < 1e-12, f"rho={row.rho};fdr={row.fdr_bh}")

selected = pd.read_csv(SELECTED, sep="\t")
check("selected_six_rows", len(selected) == 6, str(len(selected)))
state = selected[selected.endpoint == "SPP1-state (without SPP1)"]
direct = selected[selected.endpoint == "Direct SPP1 expression (ref.-dependent)"]
emt = selected[selected.endpoint == "Invasion-EMT context"]
check("state_10_of_12_both", len(state) == 2 and (state.positive == 10).all() and (state.bh_fdr.round(6) == 0.154297).all(), state.to_json(orient="records"))
check("direct_spp1_9_of_12_both", len(direct) == 2 and (direct.positive == 9).all(), direct.to_json(orient="records"))
check("emt_10_of_12_both", len(emt) == 2 and (emt.positive == 10).all() and emt.bh_fdr.between(0.1928, 0.1930).all(), emt.to_json(orient="records"))

passed = sum(ok for _, ok, _ in checks)
OUT.write_text("check\tstatus\tdetail\n" + "".join(f"{name}\t{'PASS' if ok else 'FAIL'}\t{detail}\n" for name, ok, detail in checks), encoding="utf-8")
summary = {"status": "passed" if passed == len(checks) else "failed", "passed": passed, "total": len(checks), "failed": [name for name, ok, _ in checks if not ok], "report": str(REPORT), "report_sha256": sha256(REPORT), "verification": str(OUT)}
SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
if passed != len(checks): raise SystemExit(1)
