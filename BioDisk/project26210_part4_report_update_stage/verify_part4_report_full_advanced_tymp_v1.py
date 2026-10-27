#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import pandas as pd

DOCX = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/project26210_part4_full_advanced_tymp_v2_candidate.docx")
META = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/project26210_part4_full_advanced_tymp_v2_metadata.json")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/independent_verification_v2.tsv")
SUMMARY = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_report_full_update_stage/independent_verification_summary_v2.json")
ADV = Path("/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2/cohort_synthesis_v1")
STATE = Path("/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_spp1_state_v1")

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006/relationships"

def q(ns: str, tag: str) -> str:
    return f"{{{ns}}}{tag}"

def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

checks: list[tuple[str, bool, str]] = []
def check(name: str, condition: bool, detail: str) -> None:
    checks.append((name, bool(condition), detail))

meta = json.loads(META.read_text(encoding="utf-8"))
check("candidate_exists", DOCX.is_file(), str(DOCX))
check("candidate_hash_matches_metadata", sha256_file(DOCX) == meta["candidate_sha256"], sha256_file(DOCX))
check("candidate_size_matches_metadata", DOCX.stat().st_size == meta["candidate_size"], str(DOCX.stat().st_size))

with ZipFile(DOCX) as archive:
    names = archive.namelist()
    check("zip_crc", archive.testzip() is None, str(archive.testzip()))
    check("zip_names_unique", len(names) == len(set(names)), f"{len(names)} entries")
    check("zip_entry_count", len(names) == 72, str(len(names)))
    media = [n for n in names if n.startswith("word/media/")]
    check("media_count", len(media) == 53, str(len(media)))
    doc_root = ET.fromstring(archive.read("word/document.xml"))
    rel_root = ET.fromstring(archive.read("word/_rels/document.xml.rels"))
    doc_xml = archive.read("word/document.xml").decode("utf-8")
    visible_text = "".join(t.text or "" for t in doc_root.iter(q(W, "t")))
    tables = list(doc_root.iter(q(W, "tbl")))
    drawings = list(doc_root.iter(q(W, "drawing")))
    check("table_count", len(tables) == 29, str(len(tables)))
    check("drawing_count", len(drawings) == 53, str(len(drawings)))
    bordered = sum(t.find(f"{q(W, 'tblPr')}/{q(W, 'tblBorders')}") is not None for t in tables)
    check("all_tables_have_explicit_borders", bordered == len(tables), f"{bordered}/{len(tables)}")
    rels = {r.get("Id"): r.get("Target") for r in rel_root.findall(q(PKG, "Relationship"))}
    embeds = [b.get(q(R, "embed")) for b in doc_root.iter(q(A, "blip"))]
    unresolved = [rid for rid in embeds if rid not in rels]
    missing_targets = [rels[rid] for rid in embeds if rid in rels and f"word/{rels[rid]}" not in names]
    check("drawing_relationships_resolved", not unresolved and not missing_targets, f"unresolved={unresolved};missing={missing_targets}")
    embedded_meta_ok = True
    for image in meta["embedded_images"]:
        member = f"word/media/{image['embedded_name']}"
        embedded_meta_ok &= member in names and sha256_bytes(archive.read(member)) == image["sha256"]
    check("new_image_hashes", embedded_meta_ok, "12/12 metadata hashes")

check("no_unicode_bullet", "•" not in visible_text, "forbidden bullet absent")
check("no_numbering_properties", "<w:numPr" not in doc_xml, "numPr absent")
check("report_date", visible_text.count("2026年8月3日（12例高级分析与TYMP状态更新）") == 1, "new date exactly once")
heading_413 = "4.13 【TJ-RCC 12例全套高级分析】"
heading_414 = "4.14 【TYMP机制状态分析】"
check("section_4_13_unique", visible_text.count(heading_413) == 1, "4.13 exactly once")
check("section_4_14_unique", visible_text.count(heading_414) == 1, "4.14 exactly once")
check("section_order", visible_text.index(heading_413) < visible_text.index(heading_414) < visible_text.rindex("5. 主要结论"), "4.13 < 4.14 < chapter 5 body heading")
for number in range(38, 50):
    check(f"figure_{number}_unique", visible_text.count(f"图4.{number}") == 1, f"图4.{number} exactly once")
for marker in [
    "【12例全套高级分析结论】", "【TYMP机制状态结论】", "【12例高级分析限制】", "【TYMP状态分析限制】",
    "客户湿实验应作为因果锚点", "不能独立证明5′-MTA中介", "345个基因获12/12", "18/20 MISTy目标",
    "参考轻依赖的无SPP1极化评分", "5×5 FDR=0.00781", "BH FDR=0.775–0.826", "侵袭/EMT背景",
]:
    check(f"marker_{marker[:18]}", marker in visible_text, marker)

# Recompute the principal cohort claims from formal source tables.
bayes = pd.read_csv(ADV / "bayesspace_patient_summary_v1.tsv", sep="\t")
effects = pd.read_csv(ADV / "bayesspace_patient_program_effects_v1.tsv", sep="\t")
genes = pd.read_csv(ADV / "spatial_program_gene_replication_v1.tsv", sep="\t")
candidates = pd.read_csv(ADV / "candidate_program_gene_replication_v1.tsv", sep="\t")
misty = pd.read_csv(ADV / "misty_patient_performance_summary_v1.tsv", sep="\t")
commot = pd.read_csv(ADV / "commot_patient_hypothesis_proxy_summary_v1.tsv", sep="\t")
check("advanced_patient_count", len(bayes) == 12 and bayes["sample"].nunique() == 12, f"rows={len(bayes)}")
check("advanced_spot_count", int(bayes["n_spots"].sum()) == 34609, str(int(bayes["n_spots"].sum())))
check("bayesspace_q_range", bayes["selected_q"].between(4, 6).all(), f"{bayes.selected_q.min()}-{bayes.selected_q.max()}")
te = effects[effects.feature == "Tumor-epithelial"]
check("tumor_program_12_of_12", len(te) == 12 and int((te.tumor_minus_non_tumor_median > 0).sum()) == 12, f"median={te.tumor_minus_non_tumor_median.median():.6f}")
check("strict_gene_345", int((genes.n_both_strict == 12).sum()) == 345, str(int((genes.n_both_strict == 12).sum())))
expected_programs = {"Tumor-epithelial": 2, "POSTN-ECM": 7, "APOE-TAM": 4, "TYMP-context": 5, "CD8-effector": 0}
for program, count in expected_programs.items():
    got = int(((candidates.program == program) & (candidates.n_both_strict == 12)).sum())
    check(f"candidate_{program}", got == count, str(got))
gain = misty[misty.measure == "gain.R2"]
check("misty_gain_fdr_18_of_20", len(gain) == 20 and int((gain.wilcoxon_bh_fdr_within_measure < 0.05).sum()) == 18, str(int((gain.wilcoxon_bh_fdr_within_measure < 0.05).sum())))
apoe = commot[(commot.axis_id == "APOE_TREM2_TYROBP")]
check("commot_apoe_two_scales", len(apoe) == 2 and (apoe.positive_enrichment_n == 12).all(), f"rows={len(apoe)}")
check("commot_itgb3_not_estimable", len(commot) == 19 and not commot.axis_id.str.contains("ITGB3", na=False).any(), "19 estimable rows; ITGB3 axes absent")

state = pd.read_csv(STATE / "tables/cohort_state_direction_summary_v1.tsv", sep="\t")
med = pd.read_csv(STATE / "tables/cohort_exploratory_mediation_summary_v1.tsv", sep="\t")
down = pd.read_csv(STATE / "tables/cohort_interface_downstream_summary_v1.tsv", sep="\t")
primary = state[(state.adjustment == "composition_qc") & (state.subset == "macrophage_enriched")]
ref_light = primary[(primary.route == "reference_light") & (primary.outcome == "SPP1_state_no_SPP1")]
check("state_ref_light_10_of_12", len(ref_light) == 2 and (ref_light.positive_patients == 10).all(), f"rows={len(ref_light)}")
check("state_ref_light_fdr", len(ref_light) == 2 and (ref_light.sign_test_fdr_bh.round(6) == 0.154297).all(), str(ref_light.sign_test_fdr_bh.tolist()))
contrast = primary[(primary.route == "reference_light") & (primary.outcome == "SPP1_vs_C1QC_contrast") & (primary.grid == "5x5")]
check("state_contrast_5x5_negative", len(contrast) == 1 and int(contrast.iloc[0].positive_patients) == 0 and abs(float(contrast.iloc[0].sign_test_fdr_bh) - 0.0078125) < 1e-10, contrast.to_json(orient="records"))
check("mediation_indirect_unsupported", len(med) == 4 and med.indirect_a_times_b_positive_patients.between(4, 5).all() and med.indirect_a_times_b_sign_fdr_bh.between(0.775, 0.827).all(), f"rows={len(med)}")
emt = down[down.outcome == "Invasion-EMT context"]
check("downstream_emt_trend", len(emt) == 2 and (emt.positive_patients == 10).all() and emt.sign_test_fdr_bh.between(0.1928, 0.1930).all(), f"rows={len(emt)}")

passed = sum(ok for _, ok, _ in checks)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("check\tstatus\tdetail\n" + "".join(f"{n}\t{'PASS' if ok else 'FAIL'}\t{d}\n" for n, ok, d in checks), encoding="utf-8")
summary = {"status": "passed" if passed == len(checks) else "failed", "passed": passed, "total": len(checks), "failed": [n for n, ok, _ in checks if not ok], "docx": str(DOCX), "docx_sha256": sha256_file(DOCX), "verification_tsv": str(OUT)}
SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
if passed != len(checks):
    raise SystemExit(1)
