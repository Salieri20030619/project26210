#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs_v1"
FINAL = ROOT / "results_v2"
PART = ROOT / "results_v2.part"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (INPUT / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    temporary = path.with_name(path.name + ".part")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def f(value: str) -> float:
    return float(value)


def fmt(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}g}"


def top_unique_tf(rows: list[dict[str, str]], endpoint: str, limit: int = 15) -> list[str]:
    selected = [row for row in rows if row["endpoint"] == endpoint]
    selected.sort(key=lambda row: f(row["mean_neglog10_q"]), reverse=True)
    seen: set[str] = set()
    result: list[str] = []
    for row in selected:
        tf = row["TF"]
        if tf not in seen:
            seen.add(tf)
            result.append(tf)
        if len(result) == limit:
            break
    return result


if FINAL.exists() or PART.exists():
    raise SystemExit("Refusing to overwrite an existing Phase65 result directory")
PART.mkdir(parents=True)

phase63 = read_tsv("IMMUNE_EVIDENCE_REGISTER_v1.tsv")
phase63_validation = read_tsv("PHASE63_VALIDATION_v1.tsv")
large = read_tsv("large_cohort_convergence_grading_v2.tsv")
cycle_gate = read_tsv("gse207493_cell_cycle_data_gate_v1.tsv")
cycle = read_tsv("gse207493_cell_cycle_adjusted_effect_summary_v1.tsv")
met_counts = read_tsv("gse207493_metabolic_result_counts_v1.tsv")
met_contrast = read_tsv("gse207493_metabolic_contrast_summary_retained_v1.tsv")
met_corr = read_tsv("gse207493_metabolic_correlation_summary_retained_v1.tsv")
celltype = read_tsv("gse207493_celltype_feature_summary_v1.tsv")
composition = read_tsv("gse207493_author_composition_metrics_v1.tsv")
tf_retained = read_tsv("de_novo_tf_integrated_retained_v1.tsv")
tf_cross = read_tsv("de_novo_tf_cross_method_FDR_v1.tsv")

frozen_tf_lists: dict[str, list[str]] = {}
for line in (INPUT / "gse207493_de_novo_tf_regulon_record_v1.md").read_text(encoding="utf-8").splitlines():
    for endpoint in ("tumor_f1", "myeloid_dc_tymp", "treg_tymp"):
        prefix = f"- {endpoint}: "
        if line.startswith(prefix):
            frozen_tf_lists[endpoint] = [value.strip() for value in line[len(prefix):].split(",")]

with (INPUT / "de_novo_integrated_status_v1.json").open(encoding="utf-8") as handle:
    tf_status = json.load(handle)
with (INPUT / "de_novo_scenic_status_v1.json").open(encoding="utf-8") as handle:
    scenic_status = json.load(handle)
with (INPUT / "de_novo_collectri_viper_status_v1.json").open(encoding="utf-8") as handle:
    viper_status = json.load(handle)

large_by_domain = {row["domain"]: row for row in large}
cycle_by_contrast = {row["contrast"]: row for row in cycle}
cell_by_type = {row["cell_type"]: row for row in celltype}
met_by_contrast = {row["contrast"]: row for row in met_counts}

tf_endpoint_summary: list[dict[str, object]] = []
for endpoint in ["tumor_f1", "myeloid_dc_tymp", "treg_tymp"]:
    cross = [row for row in tf_cross if row["endpoint"] == endpoint]
    retained = [row for row in tf_retained if row["endpoint"] == endpoint]
    n_values = [int(float(row[column])) for row in cross for column in ("n_SCENIC", "n_VIPER")]
    tf_endpoint_summary.append(
        {
            "endpoint": endpoint,
            "cross_method_fdr_rows": len(cross),
            "unique_cross_method_tfs": len({row["TF"] for row in cross}),
            "retained_rows": len(retained),
            "eligible_sample_n_min": min(n_values),
            "eligible_sample_n_max": max(n_values),
            "positive_same_direction_rows": sum(f(row["effect_SCENIC"]) > 0 and f(row["effect_VIPER"]) > 0 for row in cross),
            "negative_same_direction_rows": sum(f(row["effect_SCENIC"]) < 0 and f(row["effect_VIPER"]) < 0 for row in cross),
            "top_15_unique_tfs": ",".join(frozen_tf_lists[endpoint]),
        }
    )

met_summary: list[dict[str, object]] = []
for row in met_counts:
    contrast = row["contrast"]
    crows = [x for x in met_contrast if x["contrast"] == contrast and f(x["adjusted_q_value"]) < 0.05]
    rrows = [x for x in met_corr if x["contrast"] == contrast and f(x["adjusted_q_value"]) < 0.05]
    crows.sort(key=lambda x: abs(f(x["median_cycle_adjusted_difference"])), reverse=True)
    rrows.sort(key=lambda x: abs(f(x["median_cycle_adjusted_spearman_rho"])), reverse=True)
    met_summary.append(
        {
            **row,
            "top_5_adjusted_contrast_pathways": ",".join(x["pathway"] for x in crows[:5]),
            "top_5_continuous_pathways": ",".join(x["pathway"] for x in rrows[:5]),
        }
    )


def evidence(
    evidence_id: str,
    domain: str,
    finding: str,
    unit: str,
    n: str,
    effect: str,
    statistic: str,
    global_grade: str,
    external_grade: str,
    interpretation: str,
    limitation: str,
    source: str,
) -> dict[str, str]:
    return locals()


evidence_rows: list[dict[str, str]] = [
    evidence(
        "I01", "direct_clinical", "F1 as response/PFS/OS predictor", "patient", "three external clinical cohorts",
        "No robust validation; treatment-interaction directions are inconsistent",
        large_by_domain["F1_direct_clinical"]["E_MTAB_3218"] + "; " + large_by_domain["F1_direct_clinical"]["CheckMate_009_010_025"] + "; " + large_by_domain["F1_direct_clinical"]["JAVELIN_Renal101"],
        "D", large_by_domain["F1_direct_clinical"]["grade"],
        "F1 cannot be reported as a validated stand-alone immunotherapy predictor.",
        "Cohorts, endpoints and treatment interactions are heterogeneous.", "large_cohort_convergence_grading_v2.tsv",
    ),
    evidence(
        "I02", "direct_clinical", "TYMP as response/PFS/OS predictor", "patient", "three external clinical cohorts",
        "No nominal/FDR direct validation in the tested core models",
        large_by_domain["TYMP_direct_clinical"]["E_MTAB_3218"] + "; " + large_by_domain["TYMP_direct_clinical"]["CheckMate_009_010_025"] + "; " + large_by_domain["TYMP_direct_clinical"]["JAVELIN_Renal101"],
        "D", large_by_domain["TYMP_direct_clinical"]["grade"],
        "TYMP cannot be reported as a validated stand-alone immunotherapy predictor.",
        "Bulk clinical cohorts do not provide cell-resolved source information.", "large_cohort_convergence_grading_v2.tsv",
    ),
    evidence(
        "I03", "cell_context", "Bulk F1 association with Treg-rich immune state", "patient/sample and cell-resolved", "726 JAVELIN RNA samples plus 19 GSE207493 samples",
        "Bulk association replicates, but F1 is tumor-high and Treg-low at cell resolution",
        large_by_domain["F1_Treg_context"]["JAVELIN_Renal101"] + "; " + large_by_domain["F1_Treg_context"]["GSE207493_scRNA"],
        "A", large_by_domain["F1_Treg_context"]["grade"],
        "The repeated bulk F1-Treg relation is contextual/compositional, not evidence of Treg-intrinsic F1.",
        "Cell-resolved validation is one large ccRCC cohort.", "large_cohort_convergence_grading_v2.tsv;gse207493_cell_cycle_adjusted_effect_summary_v1.tsv",
    ),
    evidence(
        "I04", "cell_source", "TYMP expression in Treg and myeloid/DC", "patient and cell-resolved", "17-19 paired samples; 19 total samples",
        "TYMP is enriched in Treg versus other T cells, while myeloid/DC are the strongest overall sources",
        f"Treg adjusted difference={fmt(f(cycle_by_contrast['TYMP_Treg_vs_otherT']['median_cycle_adjusted_mean_difference']))}, q={fmt(f(cycle_by_contrast['TYMP_Treg_vs_otherT']['adjusted_q_value']))}; myeloid/DC adjusted difference={fmt(f(cycle_by_contrast['TYMP_MyeloidDC_vs_other']['median_cycle_adjusted_mean_difference']))}, q={fmt(f(cycle_by_contrast['TYMP_MyeloidDC_vs_other']['adjusted_q_value']))}",
        "A", large_by_domain["TYMP_Treg_biology"]["grade"],
        "TYMP has replicated Treg-associated biology but Treg is not its dominant tissue source.",
        "Expression/source evidence is not a clinical-response validation.", "large_cohort_convergence_grading_v2.tsv;gse207493_cell_cycle_adjusted_effect_summary_v1.tsv",
    ),
    evidence(
        "I05", "cell_source", "F1 tumor-cell localization", "patient-level paired cell contrast", "19 samples; 102,723 cells audited",
        "F1 remains higher in tumor than non-tumor cells after cycle adjustment",
        f"adjusted median difference={fmt(f(cycle_by_contrast['F1_Tumor_vs_nonTumor']['median_cycle_adjusted_mean_difference']))}, q={fmt(f(cycle_by_contrast['F1_Tumor_vs_nonTumor']['adjusted_q_value']))}",
        "B", "CELL_RESOLVED_SINGLE_COHORT", "F1 is primarily a tumor-cell transcriptional program.",
        "One cell-resolved validation cohort; not a response endpoint.", "gse207493_cell_cycle_adjusted_effect_summary_v1.tsv",
    ),
    evidence(
        "I06", "confounding_audit", "Cell-cycle robustness of F1/TYMP source contrasts", "patient-level paired residual contrast", "17-19 samples",
        "All four prespecified source contrasts retain direction and FDR after S/G2M adjustment",
        "4/4 adjusted q<0.05; 102,723/102,723 cells assigned cycle scores", "B", "CONFOUNDING_AUDIT_PASS",
        "Cell cycle is not the driver of the principal F1/TYMP cell-source findings.",
        "Residualization addresses measured S/G2M programs, not every possible confounder.", "gse207493_cell_cycle_data_gate_v1.tsv;gse207493_cell_cycle_adjusted_effect_summary_v1.tsv",
    ),
    evidence(
        "I07", "metabolism", "F1-high tumor metabolic state", "patient-specific within-lineage contrast", "19 samples",
        "Oxidative respiration, hypoxia/ROS handling and glycolysis are enriched after cycle and overlap removal",
        f"{met_by_contrast['F1_high_vs_low_Tumor']['adjusted_q_lt_0.05']}/31 contrast pathways and {met_by_contrast['F1_high_vs_low_Tumor']['continuous_FDR']}/31 continuous pathways at FDR<0.05",
        "C1", "MECHANISTIC_SINGLE_COHORT_FDR", "F1 marks a metabolically heterogeneous oxidative/hypoxic tumor program.",
        "Single cohort and expression-derived pathway scores; no matched metabolomics.", "gse207493_metabolic_*_v1.tsv",
    ),
    evidence(
        "I08", "metabolism", "TYMP-high myeloid/DC metabolic state", "patient-specific within-lineage contrast", "19 samples",
        "Pyrimidine salvage/catabolism plus glycolytic, redox and hypoxia programs remain after TYMP removal",
        f"{met_by_contrast['TYMP_high_vs_low_MyeloidDC']['adjusted_q_lt_0.05']}/31 contrast pathways and {met_by_contrast['TYMP_high_vs_low_MyeloidDC']['continuous_FDR']}/31 continuous pathways at FDR<0.05",
        "C1", "MECHANISTIC_SINGLE_COHORT_FDR", "TYMP marks a broad pyrimidine-remodeling myeloid/DC state.",
        "Single cohort and expression-derived pathway scores.", "gse207493_metabolic_*_v1.tsv",
    ),
    evidence(
        "I09", "metabolism", "TYMP-high Treg metabolic state", "patient-specific within-lineage contrast", "16 grouped; 17 continuous samples",
        "Hypoxia/redox/energy and nucleoside transport are enriched; broad pyrimidine salvage/catabolism is not",
        f"{met_by_contrast['TYMP_high_vs_low_Treg']['adjusted_q_lt_0.05']}/31 contrast pathways and {met_by_contrast['TYMP_high_vs_low_Treg']['continuous_FDR']}/31 continuous pathways at FDR<0.05",
        "C1", "MECHANISTIC_SINGLE_COHORT_FDR", "Treg TYMP denotes adaptation/transport rather than the full myeloid pyrimidine program.",
        "Smaller Treg subset and single cohort.", "gse207493_metabolic_*_v1.tsv",
    ),
]

tf_labels = {
    "tumor_f1": ("I10", "F1-high tumor TF program", "Stress/AP-1/CEBP regulatory state"),
    "myeloid_dc_tymp": ("I11", "TYMP-high myeloid/DC TF program", "IRF/SPI1/NF-kB-centered regulatory state"),
    "treg_tymp": ("I12", "TYMP-high Treg TF program", "GATA3/NF-kB/STAT/TFEB regulatory state"),
}
for summary in tf_endpoint_summary:
    eid, finding, interpretation = tf_labels[str(summary["endpoint"])]
    evidence_rows.append(
        evidence(
            eid, "tf_regulon", finding, "patient-level SCENIC and CollecTRI-VIPER", f"n={summary['eligible_sample_n_min']}-{summary['eligible_sample_n_max']}",
            f"{summary['cross_method_fdr_rows']} same-direction cross-method FDR rows; top TFs: {summary['top_15_unique_tfs']}",
            f"SCENIC and CollecTRI-VIPER q<0.05 with concordant direction; {summary['unique_cross_method_tfs']} unique TFs",
            "C1", "A_CROSS_METHOD_FDR_INTERNAL", interpretation,
            "Two methods use the same cohort and expression matrix; method convergence is not independent-cohort replication.",
            "de_novo_tf_cross_method_FDR_v1.tsv",
        )
    )

evidence_rows.extend(
    [
        evidence(
            "I13", "clinical_context", "Broad immune-state association with outcome", "patient", "three external clinical cohorts",
            "Outcome association is present in selected settings but endpoint/model dependent",
            large_by_domain["Broad_immune_state"]["E_MTAB_3218"] + "; " + large_by_domain["Broad_immune_state"]["CheckMate_009_010_025"] + "; " + large_by_domain["Broad_immune_state"]["JAVELIN_Renal101"],
            "B", large_by_domain["Broad_immune_state"]["grade"], "Immune state is clinically contextual rather than a universal predictor.",
            "Different assays, endpoints and treatment settings.", "large_cohort_convergence_grading_v2.tsv",
        ),
        evidence(
            "I14", "clinical_context", "Randomized angiogenesis/metabolism benchmark", "randomized treatment interaction", "JAVELIN Renal 101 RNA subset",
            "Angiogenesis, lipid/organic-acid metabolism and signaling interactions pass FDR",
            large_by_domain["Angiogenesis_metabolism"]["JAVELIN_Renal101"], "A", large_by_domain["Angiogenesis_metabolism"]["grade"],
            "The analysis reproduces a randomized treatment-context benchmark, but this is not F1/TYMP-specific evidence.",
            "Benchmark supports context validity, not the target axis directly.", "large_cohort_convergence_grading_v2.tsv",
        ),
        evidence(
            "I15", "foundational_small_cohorts", "Phase60-63 patient-level immune analyses", "patient or paired patient", "10 evidence families; mostly small n",
            "One B-grade F1-TYMP/Treg association, contextual pathway remodeling, and no validated direct predictor",
            f"Phase63 evidence rows={len(phase63)}; validation checks={sum(row['status']=='PASS' for row in phase63_validation)} PASS",
            "B", "FOUNDATIONAL_CONVERGENCE", "Earlier small-cohort findings remain supportive/contextual and are not discarded.",
            "Several models have only 3-4 pairs or three no-effect patients.", "IMMUNE_EVIDENCE_REGISTER_v1.tsv",
        ),
    ]
)

axis_rows = [
    {"axis": "F1", "dominant_cell_context": "Tumor", "cell_source_conclusion": "Tumor-high; Treg-low", "metabolic_program": "Oxidative respiration; hypoxia/ROS; glycolysis", "tf_program": tf_endpoint_summary[0]["top_15_unique_tfs"], "clinical_status": "Not validated; inconsistent exploratory interactions", "overall_interpretation": "Tumor-state marker, not established ICB predictor"},
    {"axis": "TYMP", "dominant_cell_context": "Myeloid/DC", "cell_source_conclusion": "Strongest overall source; Treg enriched versus other T cells", "metabolic_program": "Myeloid pyrimidine remodeling; Treg hypoxia/redox/transport", "tf_program": tf_endpoint_summary[1]["top_15_unique_tfs"] + " | Treg: " + tf_endpoint_summary[2]["top_15_unique_tfs"], "clinical_status": "No direct validation", "overall_interpretation": "Lineage-dependent metabolic/immune state marker"},
    {"axis": "Broad immune context", "dominant_cell_context": "Multiple immune lineages", "cell_source_conclusion": "Composition and state are context dependent", "metabolic_program": "Treatment setting prioritizes angiogenesis/metabolism in randomized data", "tf_program": "No single universal TF program", "clinical_status": "Supported only in selected endpoints/settings", "overall_interpretation": "Contextual modifier rather than universal biomarker"},
]

evidence_fields = list(evidence_rows[0])
write_tsv(PART / "IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v2.tsv", evidence_rows, evidence_fields)
write_tsv(PART / "IMMUNE_AXIS_CONVERGENCE_MATRIX_v1.tsv", axis_rows, list(axis_rows[0]))
write_tsv(PART / "TF_CROSS_METHOD_ENDPOINT_SUMMARY_v1.tsv", tf_endpoint_summary, list(tf_endpoint_summary[0]))
write_tsv(PART / "METABOLIC_ENDPOINT_SUMMARY_v1.tsv", met_summary, list(met_summary[0]))

input_provenance = []
for path in sorted(INPUT.iterdir()):
    if path.is_file():
        input_provenance.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
write_tsv(PART / "INPUT_PROVENANCE_v1.tsv", input_provenance, ["file", "size_bytes", "sha256"])

grade_counts = Counter(row["global_grade"] for row in evidence_rows)
tf_counts = {row["endpoint"]: row["cross_method_fdr_rows"] for row in tf_endpoint_summary}
met_counts_map = {row["contrast"]: row["adjusted_q_lt_0.05"] for row in met_summary}
summary = f"""# Part 1 immune series: mechanistic convergence v2

## Outcome

The immune series now supports a strong cell-state interpretation but still does not validate F1 or TYMP as a stand-alone immunotherapy-response biomarker. This distinction is the central Phase65 conclusion.

## Evidence base

- Foundational phases 60-63: {len(phase63)} patient-level evidence families covering response, longitudinal change, composition, state, pseudobulk differential expression and Hallmark pathways.
- Large-cohort Phase64: 19 GSE207493 samples and 102,723 QC-passing cells; author-composition Pearson r={f(composition[1]['value']):.4f} across nine broad classes.
- Cell-cycle audit: 4/4 prespecified F1/TYMP source contrasts retained direction and FDR after within-patient S/G2M adjustment.
- De novo regulation: {tf_status['n_integrated_rows']:,} integrated rows, {tf_status['n_retained_rows']:,} retained rows and {tf_status['n_grade_A_cross_method_FDR']} same-direction SCENIC/CollecTRI-VIPER FDR rows; old Part-1 networks were not used as the primary analysis.

## Converged biology

1. F1 is primarily a tumor-cell program. The repeated bulk F1-Treg association is real but is best interpreted as tissue composition/shared microenvironment rather than Treg-intrinsic F1.
2. TYMP is strongest overall in myeloid and dendritic cells, while Treg still shows reproducible enrichment relative to other T cells.
3. F1-high tumor cells couple oxidative respiration, hypoxia/ROS handling and glycolysis. {met_counts_map['F1_high_vs_low_Tumor']}/31 prespecified metabolic pathways pass the adjusted high/low FDR gate.
4. TYMP-high myeloid/DC cells show a broad pyrimidine-remodeling state ({met_counts_map['TYMP_high_vs_low_MyeloidDC']}/31 FDR pathways). TYMP-high Treg instead shows hypoxia/redox/energy adaptation and nucleoside transport ({met_counts_map['TYMP_high_vs_low_Treg']}/31 FDR pathways), without a complete TYMP-independent pyrimidine program.
5. Regulatory convergence is lineage-specific: tumor F1 has {tf_counts['tumor_f1']} cross-method FDR rows, myeloid/DC TYMP has {tf_counts['myeloid_dc_tymp']}, and Treg TYMP has {tf_counts['treg_tymp']}. These are within-cohort method-convergent results, not independent-cohort replication.

## Clinical boundary

- F1 direct clinical validation remains exploratory and directionally inconsistent across external cohorts.
- TYMP has no supported direct response/PFS/OS validation in the tested external models.
- Broad immune-state associations are treatment- and endpoint-dependent.
- Randomized JAVELIN evidence supports angiogenesis/metabolic treatment context, but it is not direct evidence for F1/TYMP.

## Global evidence grades

- A: {grade_counts['A']} rows (replicated/cell-resolved context or randomized benchmark).
- B: {grade_counts['B']} rows (limited or contextual support).
- C1: {grade_counts['C1']} rows (single-cohort FDR-supported mechanisms, including cross-method TF convergence).
- D: {grade_counts['D']} rows (currently unsupported as direct clinical predictors).

## Reporting language

Supported: “F1 and TYMP define lineage-dependent tumor and immune metabolic/regulatory states.”

Not supported: “F1 or TYMP is a validated stand-alone immunotherapy-response biomarker,” “the TFs prove causal regulation,” or “the metabolic scores constitute matched metabolomic validation.”

Structured evidence is in `IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v2.tsv`; axis-level synthesis is in `IMMUNE_AXIS_CONVERGENCE_MATRIX_v1.tsv`.
"""
(PART / "PART1_IMMUNE_SERIES_CONVERGENCE_v2.md").write_text(summary, encoding="utf-8")

summary_cn = f"""# 第一部分免疫系列：机制证据收敛（v2）

## 一、总结果

免疫系列目前已经形成较强的细胞状态与机制收敛证据，但仍未验证 F1 或 TYMP 是独立的免疫治疗疗效标志物。机制证据增强与临床预测结论不升级，是 Phase65 最重要的边界。

## 二、证据基础

- Phase60–63：{len(phase63)} 个患者级证据家族，覆盖疗效、纵向变化、组成、状态、pseudobulk 差异表达和 Hallmark 通路。
- Phase64 大队列：GSE207493 共 19 个样本、102,723 个质控合格细胞；九类广义细胞组成与作者结果 Pearson r={f(composition[1]['value']):.4f}。
- 细胞周期审计：4/4 个预设 F1/TYMP 细胞来源对比在患者内 S/G2M 校正后均保持方向并通过 FDR。
- de novo 调控分析：{tf_status['n_integrated_rows']:,} 条整合结果、{tf_status['n_retained_rows']:,} 条保留结果、{tf_status['n_grade_A_cross_method_FDR']} 条 SCENIC/CollecTRI‑VIPER 同方向双方法 FDR 结果；旧 Part 1 调控网络未作为主分析输入。

## 三、收敛生物学结论

1. F1 主要是肿瘤细胞程序。多队列 bulk F1–Treg 关联可以复现，但细胞分辨率结果支持其主要反映组织组成或共享微环境，而不是 Treg 内源性 F1。
2. TYMP 的总体最强来源是髓系和树突细胞；与此同时，Treg 相对其他 T 细胞仍存在可复现的 TYMP 富集。
3. F1-high 肿瘤细胞同时呈现氧化呼吸、缺氧/ROS 处理和糖酵解程序；{met_counts_map['F1_high_vs_low_Tumor']}/31 条预设代谢通路通过校正后的高低组 FDR 门槛。
4. TYMP-high 髓系/DC 呈现广泛嘧啶重塑（{met_counts_map['TYMP_high_vs_low_MyeloidDC']}/31 条 FDR 通路）。TYMP-high Treg 主要表现为缺氧、氧化还原、能量适应和核苷转运（{met_counts_map['TYMP_high_vs_low_Treg']}/31 条 FDR 通路），未形成完整的、独立于 TYMP 锚定基因的嘧啶程序。
5. 调控收敛具有谱系特异性：Tumor F1 有 {tf_counts['tumor_f1']} 条双方法 FDR 结果，Myeloid/DC TYMP 有 {tf_counts['myeloid_dc_tymp']} 条，Treg TYMP 有 {tf_counts['treg_tymp']} 条。这是同一队列内的方法收敛，不是独立队列复现。

## 四、临床解释边界

- F1 的直接临床验证仍属探索性，且外部队列治疗交互方向不一致。
- TYMP 在已检验的外部模型中没有获得直接的 response/PFS/OS 支持。
- 广义免疫状态与结局的关系依赖治疗背景、终点和统计模型。
- 随机化 JAVELIN 证据支持血管生成/代谢治疗背景，但不是 F1/TYMP 的直接证据。

## 五、全局证据等级

- A：{grade_counts['A']} 项，代表复现并获得细胞分辨率支持的情境结论或随机化基准。
- B：{grade_counts['B']} 项，代表有限或情境性支持。
- C1：{grade_counts['C1']} 项，代表单队列 FDR 支持的机制结果，包括双方法 TF 收敛。
- D：{grade_counts['D']} 项，代表当前不支持作为直接临床预测指标。

## 六、推荐写法

可以写：“F1 与 TYMP 界定了谱系依赖的肿瘤及免疫代谢/调控状态。”

不能写：“F1 或 TYMP 已经是独立、经过验证的免疫治疗疗效标志物”“TF 结果证明了因果调控”或“表达评分等同于配对代谢组验证”。

结构化证据见 `IMMUNE_MECHANISTIC_EVIDENCE_REGISTER_v2.tsv`，轴级收敛见 `IMMUNE_AXIS_CONVERGENCE_MATRIX_v1.tsv`。
"""
(PART / "PART1_IMMUNE_SERIES_CONVERGENCE_v2.md").write_text(summary_cn, encoding="utf-8")

readme = """# Phase65: Part 1 immune mechanistic convergence

This stage consolidates validated Phase60-63 patient-level immune evidence with Phase64 large-cohort cell-source, cell-cycle, metabolism and de novo TF/regulon results.

No raw counts, H5AD, normalization, integration, clustering, annotation or model training is performed. The stage operates only on compact, versioned result tables and preserves the distinction between mechanistic convergence and clinical prediction.
"""
(PART / "README.md").write_text(readme, encoding="utf-8")

supersession = """# Phase65 results v2 supersession record

`results_v2` supersedes the unsynchronized local `results_v1` draft. The numerical evidence, row counts and grades are unchanged. Version 2 restores the representative TF lists directly from the frozen Phase64 record instead of re-ranking tied cross-method rows.
"""
(PART / "SUPERSESSION_RECORD_v2.md").write_text(supersession, encoding="utf-8")

metadata = {
    "status": "BUILD_COMPLETE_VALIDATION_PENDING",
    "created_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
    "phase63_evidence_rows": len(phase63),
    "phase65_evidence_rows": len(evidence_rows),
    "grade_counts": dict(sorted(grade_counts.items())),
    "tf_cross_method_rows": len(tf_cross),
    "metabolic_contrasts": len(met_summary),
    "result_version": "v2",
    "raw_data_read": False,
    "model_training_performed": False,
}
(PART / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

os.replace(PART, FINAL)
print(json.dumps(metadata, ensure_ascii=False, indent=2))
