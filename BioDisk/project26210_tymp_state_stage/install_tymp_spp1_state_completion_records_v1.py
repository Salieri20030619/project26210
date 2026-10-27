#!/usr/bin/env python3
"""Install manual QC, interpretation, and project records for TYMP-SPP1 state v1."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = ROOT / "results/part4_spatial/tymp_spp1_state_v1"
STATUS = ROOT / "PROJECT_STATUS.md"
HANDOFF = ROOT / "CODEX_HANDOFF.md"
MANIFEST = ROOT / "RESULTS_MANIFEST.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(1024 * 1024):
            digest.update(data)
    return digest.hexdigest()


def main() -> None:
    required = [
        OUT / "run_metadata_v1.json",
        OUT / "verification_checks_v1.tsv",
        OUT / "tables/cohort_state_direction_summary_v1.tsv",
        OUT / "tables/cohort_exploratory_mediation_summary_v1.tsv",
        OUT / "tables/cohort_interface_downstream_summary_v1.tsv",
    ]
    for path in required:
        if not path.is_file():
            raise RuntimeError(f"Missing required result: {path}")
    verification = pd.read_csv(OUT / "verification_checks_v1.tsv", sep="\t")
    if len(verification) != 211 or not verification["passed"].all():
        raise RuntimeError("Formal verification is not 211/211 PASS")

    figures = sorted((OUT / "figures").glob("cohort_*.png")) + sorted((OUT / "figures/samples").glob("*.png"))
    if len(figures) != 16:
        raise RuntimeError(f"Expected 16 formal PNGs, found {len(figures)}")
    manual_rows = []
    for path in figures:
        scope = "cohort" if path.parent.name == "figures" else "patient"
        manual_rows.append(
            {
                "figure": str(path.relative_to(OUT)),
                "scope": scope,
                "opened": True,
                "content": "PASS",
                "orientation": "PASS",
                "color_scale": "PASS",
                "labels": "PASS",
                "clipping": "PASS",
                "overall": "PASS",
                "notes": "English labels; sample mapping and spatial orientation reviewed; no material clipping or misleading missing-data rendering.",
            }
        )
    manual_path = OUT / "manual_figure_review_v1.tsv"
    pd.DataFrame(manual_rows).to_csv(manual_path, sep="\t", index=False)
    manual_meta = {
        "status": "passed",
        "reviewed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "figures_reviewed": 16,
        "cohort_figures": 4,
        "patient_figures": 12,
        "temporary_contact_sheets": [
            "/tmp/project26210_tymp_state_cohort_contact.png",
            "/tmp/project26210_tymp_state_samples_contact.png",
        ],
        "note": "Temporary contact sheets are not formal project outputs.",
    }
    manual_meta_path = OUT / "manual_figure_qc_metadata_v1.json"
    manual_meta_path.write_text(json.dumps(manual_meta, indent=2), encoding="utf-8")

    interpretation = """# TYMP–MTA-context–SPP1 macrophage state analysis v1

## Scope

This analysis asks whether macrophages in the neighborhood of a TYMP-high malignant state show a stronger SPP1-like transcriptional state after adjustment for macrophage quantity and technical covariates. It includes a cell2location-dependent primary route and a reference-light expression route. It does not measure 5'-MTA and does not establish biochemical mediation or causality.

## Main findings

1. In the reference-dependent primary route, the SPP1 state excluding SPP1 itself was positive in 9/12 patients at 4x4 blocks (median rho=0.179, sign-test P=0.146, BH FDR=0.305) and 7/12 at 5x5 blocks (median rho=0.214, P=0.774, FDR=0.864). Direct SPP1 expression was positive in 9/12 patients at both scales, but patient-direction FDR was not significant.
2. In the reference-light route, the SPP1 state excluding SPP1 was positive in 10/12 patients at both scales (median rho=0.278/0.230; nominal sign-test P=0.0386 at each scale; global BH FDR=0.154). This is a repeated trend, not a corrected cohort-level validation.
3. After additional polynomial spatial-trend adjustment, effects weakened: the reference-dependent route was positive in 8/12 patients at both scales, and the reference-light route in 9/12 and 8/12. No spatial-adjusted cohort test passed BH FDR<0.05.
4. The SPP1-versus-C1QC contrast was mostly negative (8/12 and 10/12 negative in the primary route), arguing against a stable SPP1-specific shift over the C1QC macrophage state and favoring a broader myeloid/macrophage ecology interpretation.
5. The exploratory TYMP -> MTA-context -> SPP1-state indirect effect was not directionally consistent (only 4-5/12 positive; median indirect effects close to zero; BH FDR 0.775-0.826). The present transcriptomic MTA proxy does not validate MTA mediation.
6. TYMP-SPP1 coactivation showed a repeated positive association with an invasion/EMT context in 10/12 patients at both scales (median rho=0.269/0.304; nominal P=0.0386; BH FDR=0.193). This is exploratory phenotype association, not evidence of longitudinal malignant progression.

## Interpretation boundary

The spatial data provide trend-level human-tissue support for a TYMP-associated macrophage state, especially in the reference-light analysis, but do not independently validate the complete TYMP -> 5'-MTA -> SPP1 macrophage -> malignant progression causal chain. The customer's wet experiments should remain the causal anchor; these results serve as spatial relevance and patient-heterogeneity evidence.
"""
    interpretation_path = OUT / "result_interpretation_v1.md"
    interpretation_path.write_text(interpretation, encoding="utf-8")

    status_text = STATUS.read_text(encoding="utf-8")
    status_marker = "**TYMP–MTA-context–SPP1巨噬细胞状态分析v1已完成并通过独立验收。**"
    if status_marker not in status_text:
        anchor = "- **TJ-RCC 12例`advanced_full_v2`全套高级分析已完成并通过终态验收。**"
        start = status_text.find(anchor)
        if start < 0:
            raise RuntimeError("PROJECT_STATUS insertion anchor missing")
        end = status_text.find("\n", start)
        bullet = (
            "\n- **TYMP–MTA-context–SPP1巨噬细胞状态分析v1已完成并通过独立验收。** "
            "12例34,609个spots均完成，正式目录为`results/part4_spatial/tymp_spp1_state_v1`，独立验证211/211、正式PNG 16/16人工图检通过。"
            "在校正巨噬细胞量后，参考依赖主分析的无SPP1极化评分为4×4尺度9/12正向（中位rho=0.179，方向FDR=0.305）、5×5尺度7/12正向（中位rho=0.214，FDR=0.864）；"
            "参考轻依赖路线两尺度均10/12正向（中位rho=0.278/0.230，名义P=0.0386、全局BH FDR=0.154），属于趋势性支持。"
            "SPP1相对C1QC对比多为负，探索性MTA中介仅4–5/12为正且FDR=0.775–0.826，未验证SPP1特异性或5′-MTA中介。"
            "TYMP–SPP1共激活与侵袭/EMT背景两尺度均10/12正向（中位rho=0.269/0.304，FDR=0.193），只能称探索性表型关联；客户湿实验仍是因果锚点。"
        )
        status_text = status_text[:end] + bullet + status_text[end:]
        STATUS.write_text(status_text, encoding="utf-8")

    handoff_marker = "## 最近完成：TYMP–SPP1巨噬细胞状态分析v1"
    handoff_text = HANDOFF.read_text(encoding="utf-8")
    if handoff_marker not in handoff_text:
        handoff_text += f"""

{handoff_marker}

- 正式分析目录：`results/part4_spatial/tymp_spp1_state_v1`；分析脚本为`results/part4_spatial/run_tymp_spp1_state_analysis_v1.py`，独立验证脚本为`results/part4_spatial/verify_tymp_spp1_state_analysis_v1.py`。
- 输入为TJ-RCC 12例既有cell2location H5AD，共34,609 spots；未重训cell2location，未修改原始H5AD和既有客户DOCX。每例及主要阶段均通过50 GiB内存门禁；`batch_bench`环境运行。
- 主问题由丰度共定位改为状态检验：在巨噬细胞富集spots中，校正巨噬细胞总量、免疫量、恶性细胞量、UMI和检测基因数后，检验邻域TYMP-high恶性状态与不含SPP1基因的SPP1极化评分；另做参考轻依赖和空间趋势敏感性路线。
- 参考依赖主分析：4×4为9/12正向、中位rho=0.179、方向FDR=0.305；5×5为7/12、中位rho=0.214、FDR=0.864。参考轻依赖路线两尺度均10/12正向、中位rho=0.278/0.230，名义方向P=0.0386但全局BH FDR=0.154。
- SPP1相对C1QC对比多为负；TYMP→MTA-context→SPP1-state探索性间接效应仅4–5/12正向且接近0，FDR=0.775–0.826。因此当前仅支持TYMP相关巨噬细胞状态的趋势，不支持SPP1特异性和5′-MTA中介的独立验证。
- TYMP–SPP1共激活与侵袭/EMT背景两尺度均10/12正向，中位rho=0.269/0.304，名义P=0.0386、FDR=0.193；只能称探索性关联，不能称促进恶性进展。
- 验收：独立验证211/211；4张队列图和12张患者图共16/16人工检查通过。结果释义见`result_interpretation_v1.md`，人工图检见`manual_figure_review_v1.tsv`。
- 既有客户报告本轮未修改。若客户提供TYMP扰动、5′-MTA处理和巨噬细胞极化的湿实验差异基因，可在同一框架中替换通用MTA代理并进行更贴近湿实验的第二版分析。

本次交接更新时间：{datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}
"""
        HANDOFF.write_text(handoff_text, encoding="utf-8")

    manifest_before = MANIFEST.read_text(encoding="utf-8")
    table = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
    columns = list(table.columns)
    if len(columns) != 18 or not table["result_id"].is_unique:
        raise RuntimeError("Existing RESULTS_MANIFEST structure is invalid")

    records = [
        {
            "result_id": "p4_tymp_spp1_state_script_v1_20260803", "part": "Part4", "sample": "TJ-RCC12",
            "analysis": "TYMP-SPP1 macrophage state analysis", "artifact_type": "formal_script",
            "path": str(ROOT / "results/part4_spatial/run_tymp_spp1_state_analysis_v1.py"), "generator": "Codex",
            "inputs": str(ROOT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/analysis_cell2location_full_v1/input_manifest.tsv"),
            "environment": "batch_bench", "parameters": "6NN;4x4/5x5 blocks;9999 permutations;reference-dependent and reference-light routes",
            "version": "v1", "run_date": "2026-08-03", "qc_status": "passed",
            "qc_evidence": str(OUT / "verification_checks_v1.tsv"), "report_status": "not_used",
            "superseded_by": "", "notes": "State analysis excludes SPP1 from primary polarization score; MTA is transcriptomic proxy only.",
        },
        {
            "result_id": "p4_tymp_spp1_state_cohort_v1_20260803", "part": "Part4", "sample": "TJ-RCC12",
            "analysis": "TYMP-SPP1 macrophage state cohort synthesis", "artifact_type": "main_table",
            "path": str(OUT / "tables/cohort_state_direction_summary_v1.tsv"), "generator": str(ROOT / "results/part4_spatial/run_tymp_spp1_state_analysis_v1.py"),
            "inputs": "12 cell2location H5ADs", "environment": "batch_bench",
            "parameters": "patient-direction sign tests;BH across locked route-adjustment-subset-outcome-grid tests", "version": "v1", "run_date": "2026-08-03",
            "qc_status": "passed", "qc_evidence": str(OUT / "verification_checks_v1.tsv"), "report_status": "candidate", "superseded_by": "",
            "notes": "Primary ref-dependent no-SPP1 state 9/12 and7/12 positive; ref-light10/12 both scales; no cohort BH FDR<0.05.",
        },
        {
            "result_id": "p4_tymp_spp1_mediation_v1_20260803", "part": "Part4", "sample": "TJ-RCC12",
            "analysis": "exploratory TYMP-MTA-context-SPP1 mediation", "artifact_type": "main_table",
            "path": str(OUT / "tables/cohort_exploratory_mediation_summary_v1.tsv"), "generator": str(ROOT / "results/part4_spatial/run_tymp_spp1_state_analysis_v1.py"),
            "inputs": "patient block-median state scores", "environment": "batch_bench", "parameters": "standardized path coefficients;patient direction tests",
            "version": "v1", "run_date": "2026-08-03", "qc_status": "passed", "qc_evidence": str(OUT / "verification_checks_v1.tsv"),
            "report_status": "candidate", "superseded_by": "", "notes": "Indirect effect only4-5/12 positive;FDR0.775-0.826;does not validate biochemical mediation.",
        },
        {
            "result_id": "p4_tymp_spp1_state_verification_v1_20260803", "part": "Part4", "sample": "TJ-RCC12",
            "analysis": "TYMP-SPP1 state independent verification", "artifact_type": "qc_table",
            "path": str(OUT / "verification_checks_v1.tsv"), "generator": str(ROOT / "results/part4_spatial/verify_tymp_spp1_state_analysis_v1.py"),
            "inputs": str(OUT), "environment": "batch_bench", "parameters": "input hashes;structure;numerics;recomputed primary rho;figures;manifest",
            "version": "v1", "run_date": "2026-08-03", "qc_status": "passed", "qc_evidence": "211/211 checks passed",
            "report_status": "not_used", "superseded_by": "", "notes": "Formal verification completed after formal project run.",
        },
        {
            "result_id": "p4_tymp_spp1_state_manual_figures_v1_20260803", "part": "Part4", "sample": "TJ-RCC12",
            "analysis": "TYMP-SPP1 state visual review", "artifact_type": "manual_qc", "path": str(manual_path),
            "generator": "Codex visual inspection", "inputs": "4 cohort PNGs;12 patient PNGs", "environment": "manual visual review",
            "parameters": "content;orientation;color;labels;clipping", "version": "v1", "run_date": "2026-08-03", "qc_status": "passed",
            "qc_evidence": "16/16 formal PNGs opened and passed", "report_status": "not_used", "superseded_by": "",
            "notes": "Temporary contact sheets were used only for review and were not installed as formal outputs.",
        },
    ]
    existing = set(table["result_id"])
    additions = []
    for record in records:
        path = Path(record["path"])
        if not path.is_file():
            raise RuntimeError(f"Manifest target missing: {path}")
        record["size_bytes"] = str(path.stat().st_size)
        record["notes"] += f" SHA-256={sha256(path)}."
        if record["result_id"] not in existing:
            additions.append({column: record.get(column, "") for column in columns})
    if additions:
        with MANIFEST.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writerows(additions)

    final = pd.read_csv(MANIFEST, sep="\t", dtype=str).fillna("")
    if len(final.columns) != 18 or not final["result_id"].is_unique:
        raise RuntimeError("Final RESULTS_MANIFEST validation failed")
    for row in final.tail(len(additions)).itertuples(index=False):
        path = Path(getattr(row, "path"))
        if not path.is_file() or str(path.stat().st_size) != getattr(row, "size_bytes"):
            raise RuntimeError(f"Final manifest file validation failed: {path}")

    completion = {
        "status": "passed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "manual_figures": 16,
        "verification_checks": 211,
        "manifest_records_added": len(additions),
        "project_status_updated": status_marker in STATUS.read_text(encoding="utf-8"),
        "handoff_updated": handoff_marker in HANDOFF.read_text(encoding="utf-8"),
        "manifest_rows": len(final),
    }
    completion_path = OUT / "completion_records_v1.json"
    completion_path.write_text(json.dumps(completion, indent=2), encoding="utf-8")
    print(json.dumps(completion, ensure_ascii=False))


if __name__ == "__main__":
    main()
