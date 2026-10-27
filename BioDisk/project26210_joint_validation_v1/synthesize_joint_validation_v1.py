#!/usr/bin/env python3
"""Synthesize rebuilt single-cell and 12-patient spatial validation evidence."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_joint_validation_v1")
SC = ROOT / "singlecell_axis_audit_v1"
SP = ROOT / "rebuilt_F1_spatial_axis_validation_v1"
EXISTING_VASCULAR = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_vascular_niche_v1/"
    "tables/tjrcc_cohort_direction_summary_v1.tsv"
)

AXIS_TO_LIGAND = {
    "VEGFA_KDR": "VEGFA",
    "VEGFA_FLT1": "VEGFA",
    "VEGFB_FLT1": "VEGFB",
    "PGF_FLT1": "PGF",
    "APP_CD74": "APP",
    "APOE_TREM2_TYROBP": "APOE",
    "SPP1_CD44": "SPP1",
    "MIF_CD74_CXCR4": "MIF",
}


def conclusion(axis: str) -> str:
    if axis.startswith("VEGFA"):
        return "do_not_use_as_F1_specific_axis; within-patient F1–VEGFA direction is negative"
    if axis == "PGF_FLT1":
        return "do_not_use_as_F1_specific_axis; within-patient F1–PGF direction is negative"
    if axis == "VEGFB_FLT1":
        return "single-cell F1-linked candidate only; exact spatial axis not supported"
    if axis == "SPP1_CD44":
        return "F1–SPP1 and malignant SPP1 source supported; exact spatial SPP1–CD44 axis not supported"
    if axis == "MIF_CD74_CXCR4":
        return "single-cell F1–MIF link supported; exact spatial axis not supported"
    if axis == "APP_CD74":
        return "CD74-positive macrophage context supported at one grid; APP ligand specificity not supported"
    if axis == "APOE_TREM2_TYROBP":
        return "expression competent but neither F1-specific nor spatially supported"
    return "not supported"


def main() -> None:
    f1 = pd.read_csv(SC / "cohort_F1_ligand_direction_summary_v1.tsv", sep="\t").set_index("ligand")
    competence = pd.read_csv(SC / "axis_cohort_expression_summary_v1.tsv", sep="\t")
    centered = pd.read_csv(SC / "axis_study_centered_sensitivity_v1.tsv", sep="\t")
    centered = centered.loc[centered["excluded_study"].eq("NONE")].set_index("axis")
    spatial_decision = pd.read_csv(SP / "primary_axis_decision_v1.tsv", sep="\t").set_index("axis")
    spatial_cohort = pd.read_csv(SP / "cohort_spatial_axis_direction_summary_v1.tsv", sep="\t")
    spp1_source = pd.read_csv(SC / "SPP1_source_cohort_summary_v1.tsv", sep="\t").iloc[0]
    existing = pd.read_csv(EXISTING_VASCULAR, sep="\t")

    rows = []
    for axis, ligand in AXIS_TO_LIGAND.items():
        comp_axis = "Tumor_" + axis
        comp = competence.loc[competence["axis"].eq(comp_axis)]
        if comp.empty:
            raise RuntimeError(f"Missing single-cell competence row for {comp_axis}")
        comp = comp.iloc[0]
        spatial_axis = spatial_cohort.loc[
            spatial_cohort["axis"].eq(axis)
            & spatial_cohort["subset"].eq("tumor_enriched")
            & spatial_cohort["test"].eq("F1_to_neighbor_receiver")
        ]
        ligand_spatial = spatial_cohort.loc[
            spatial_cohort["axis"].eq(axis)
            & spatial_cohort["subset"].eq("tumor_enriched")
            & spatial_cohort["test"].eq("ligand_to_neighbor_receiver")
        ]
        rows.append({
            "axis": axis,
            "ligand": ligand,
            "singlecell_F1_patients": int(f1.loc[ligand, "n_units"]),
            "singlecell_F1_positive_patients": int(f1.loc[ligand, "positive_units"]),
            "singlecell_F1_median_rho": float(f1.loc[ligand, "median_effect"]),
            "singlecell_F1_sign_FDR": float(f1.loc[ligand, "sign_test_fdr_bh"]),
            "singlecell_axis_evaluable_patients": int(comp["n_patients_evaluable"]),
            "singlecell_median_ligand_detection": float(comp["median_ligand_detection_fraction"]),
            "singlecell_median_limiting_receptor_detection": float(comp["median_limiting_receptor_detection_fraction"]),
            "study_centered_F1_pair_rho": float(centered.loc[comp_axis, "rho_F1_vs_pair_score_study_centered"]),
            "study_centered_F1_pair_FDR": float(centered.loc[comp_axis, "fdr_F1_vs_pair_score_study_centered"]),
            "spatial_F1_positive_patients_4x4": int(spatial_axis.loc[spatial_axis["grid"].eq("4x4"), "positive_patients"].iloc[0]),
            "spatial_F1_positive_patients_5x5": int(spatial_axis.loc[spatial_axis["grid"].eq("5x5"), "positive_patients"].iloc[0]),
            "spatial_F1_min_global_FDR": float(spatial_axis["sign_test_fdr_bh_global"].min()),
            "spatial_ligand_positive_patients_4x4": int(ligand_spatial.loc[ligand_spatial["grid"].eq("4x4"), "positive_patients"].iloc[0]),
            "spatial_ligand_positive_patients_5x5": int(ligand_spatial.loc[ligand_spatial["grid"].eq("5x5"), "positive_patients"].iloc[0]),
            "spatial_ligand_min_global_FDR": float(ligand_spatial["sign_test_fdr_bh_global"].min()),
            "spatial_decision": spatial_decision.loc[axis, "decision"],
            "integrated_conclusion": conclusion(axis),
        })
    evidence = pd.DataFrame(rows)
    evidence.to_csv(ROOT / "joint_axis_evidence_matrix_v1.tsv", sep="\t", index=False)

    contextual = existing.loc[
        existing["relation"].isin([
            "TYMP_MTA_core_to_SPP1_vascular_niche",
            "SPP1_vascular_niche_to_T_cell_balance",
            "SPP1_TAM_to_vascular_compartment",
            "SPP1_vascular_niche_to_hypoxia_angiogenesis",
        ])
    ].copy()
    contextual.to_csv(ROOT / "existing_spatial_context_summary_v1.tsv", sep="\t", index=False)

    report = f"""# project26210 新版单细胞—空间关键轴联合验证 v1

## 结论

本轮没有重复运行全局 CellChat 或 LIANA，而是使用第一部分重建后的 258,120 细胞对象、11,950 个冻结保守恶性细胞和 rank 7/F1，对第三部分重点配体–受体轴进行患者级复核，并在 12 例 TJ-RCC Visium 中进行双尺度空间验证。

没有任何一条精确配体–受体轴同时通过新版单细胞患者级证据和 12 例空间 4×4/5×5 双尺度验证。新增分析的主要价值是限定和纠正第三部分的解释边界，而不是增加新的通讯网络。

## 新版单细胞结果

- F1 与 SPP1 同向：17/20 名患者为正，患者 rho 中位数 {f1.loc['SPP1','median_effect']:.3f}，方向 FDR={f1.loc['SPP1','sign_test_fdr_bh']:.4g}。
- F1 与 VEGFB 同向：18/20 名患者为正，rho 中位数 {f1.loc['VEGFB','median_effect']:.3f}，FDR={f1.loc['VEGFB','sign_test_fdr_bh']:.4g}。
- F1 与 MIF 同向：16/20 名患者为正，rho 中位数 {f1.loc['MIF','median_effect']:.3f}，FDR={f1.loc['MIF','sign_test_fdr_bh']:.4g}；但排除 GSE207493 后效应明显减弱。
- F1 与 VEGFA 反向：18/20 名患者为负，rho 中位数 {f1.loc['VEGFA','median_effect']:.3f}，FDR={f1.loc['VEGFA','sign_test_fdr_bh']:.4g}。因此不应继续把 VEGFA 写成 F1/TYMP 的特异正向配体。
- F1 与 PGF 同样偏负；APP、APOE 没有稳定 F1 特异性。
- SPP1 来源比较中，恶性细胞检出率高于巨噬细胞的患者为 {int(spp1_source['malignant_higher_detection_patients'])}/{int(spp1_source['n_patients'])}，恶性细胞 log1p-CPM 更高为 {int(spp1_source['malignant_higher_log1p_cpm_patients'])}/{int(spp1_source['n_patients'])}。这否定了把空间 SPP1 默认归属于 TAM 的做法，但不代表 TAM 不表达 SPP1。
- 未校正的跨患者 VEGF 配对相关在研究内中心化后消失，说明旧 TYMP 分组和队列构成可造成辛普森悖论式假象。

## 12 例空间验证

- 精确的 VEGFA–KDR/FLT1、VEGFB–FLT1、PGF–FLT1、SPP1–CD44、MIF–CD74/CXCR4、APOE–TREM2/TYROBP 均未通过患者级双尺度验证。
- APP–CD74 仅表现为 F1 区域邻近 CD74+ 巨噬环境：5×5 中 12/12 患者同向且全局 FDR<0.05，但 APP 配体表达本身没有与邻近受体环境复现。因此只能称为 CD74+ 髓系环境关联，不能称为 APP–CD74 轴被验证。
- 既有 Part 4 的 TYMP/MTA-core 与复合 SPP1–vascular niche 在参考依赖路线中 12/12 同向，可继续作为“生态位层面关联”；它不能替代具体配体–受体轴验证，参考轻量路线也较弱。
- 既有 SPP1+ TAM 与 vascular compartment 只有 4/12 患者同向，不能支持“SPP1+ TAM 普遍伴随血管富集”。

## 主线取舍

1. 保留“F1/TYMP 相关肿瘤状态与 CD74+ 髓系/复合 SPP1–vascular 微环境有关”的生态位级表述。
2. SPP1 作为 F1 相关恶性程序成分可以保留，但必须强调恶性细胞也是重要来源；不再默认解释为 TAM 单一来源。
3. VEGFA、APP–CD74、APOE–TREM2、MIF、SPP1–CD44 均不应写成已被空间验证的确定通讯轴。
4. 第三部分既有 CellChat/LIANA 图保留为候选机制生成证据，精确配体–受体结论降级到补充材料或讨论。
5. 后续新增外部单细胞/空间队列时，应优先验证 SPP1 来源、CD74+ 髓系邻域和 VEGFB，而不是再做全局通讯发现。

## 解释边界

所有空间检验均以患者为生物学重复；4×4 和 5×5 空间块各做 9,999 次置换。F1 投射签名排除了 TYMP 和全部被检验的配体/受体基因。结果表示表达能力和空间关联，不证明细胞直接通讯、招募、驱动或因果作用。
"""
    (ROOT / "JOINT_VALIDATION_DECISION_v1.md").write_text(report, encoding="utf-8")

    verification = {
        "status": "passed",
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "singlecell_axes": int(len(competence)),
        "spatial_axes": int(len(spatial_decision)),
        "spatial_patients": int(spatial_cohort["n_patients"].min()),
        "exact_axes_supported": int(spatial_decision["decision"].eq("exact_axis_supported_both_grids").sum()),
        "report": str(ROOT / "JOINT_VALIDATION_DECISION_v1.md"),
        "evidence_matrix": str(ROOT / "joint_axis_evidence_matrix_v1.tsv"),
    }
    if verification["singlecell_axes"] != 9 or verification["spatial_axes"] != 8 or verification["spatial_patients"] != 12:
        raise RuntimeError(f"Verification count mismatch: {verification}")
    (ROOT / "joint_validation_verification_v1.json").write_text(
        json.dumps(verification, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    manifest_rows = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.name == "joint_validation_file_manifest_v1.tsv":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_rows.append({"path": str(path.relative_to(ROOT)), "size_bytes": path.stat().st_size, "sha256": digest})
    pd.DataFrame(manifest_rows).to_csv(ROOT / "joint_validation_file_manifest_v1.tsv", sep="\t", index=False)
    print(json.dumps(verification, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
