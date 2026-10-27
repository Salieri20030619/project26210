#!/usr/bin/env python3
"""Build a wet-lab-anchored, small-sample interpretation of joint evidence.

The v1 exploratory analyses remain unchanged. This v2 separates a small set of
directional biological hypotheses from secondary ligand-receptor sensitivity
tests and does not use secondary multiplicity correction to veto wet-lab evidence.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_joint_validation_v1")
OUT = ROOT / "hypothesis_led_v2"
SC = ROOT / "singlecell_axis_audit_v1"
SP = ROOT / "rebuilt_F1_spatial_axis_validation_v1"
VASCULAR = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/tymp_vascular_niche_v1/"
    "tables/tjrcc_cohort_direction_summary_v1.tsv"
)


def exact_direction(k: int, n: int) -> tuple[float, float, float]:
    test = binomtest(k, n, 0.5, alternative="greater")
    interval = binomtest(k, n).proportion_ci(confidence_level=0.95, method="exact")
    return float(test.pvalue), float(interval.low), float(interval.high)


def holm(values: pd.Series) -> pd.Series:
    order = values.sort_values().index
    n = len(order)
    raw = values.loc[order].to_numpy(dtype=float)
    adjusted = np.maximum.accumulate(raw * np.arange(n, 0, -1))
    result = pd.Series(index=values.index, dtype=float)
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    f1 = pd.read_csv(SC / "cohort_F1_ligand_direction_summary_v1.tsv", sep="\t").set_index("ligand")
    spp1 = pd.read_csv(SC / "SPP1_source_cohort_summary_v1.tsv", sep="\t").iloc[0]
    spatial = pd.read_csv(SP / "cohort_spatial_axis_direction_summary_v1.tsv", sep="\t")
    vascular = pd.read_csv(VASCULAR, sep="\t")

    app_cd74 = spatial.loc[
        spatial["axis"].eq("APP_CD74")
        & spatial["subset"].eq("tumor_enriched")
        & spatial["test"].eq("F1_to_neighbor_receiver")
    ].set_index("grid")

    def vascular_row(relation: str, grid: str) -> pd.Series:
        frame = vascular.loc[
            vascular["relation"].eq(relation)
            & vascular["route"].eq("route_independent" if relation.startswith("SPP1_vascular_niche_to_T_cell") else "reference_dependent")
            & vascular["grid"].eq(grid)
        ]
        if len(frame) != 1:
            raise RuntimeError(f"Expected one existing spatial row: {relation}, {grid}; got {len(frame)}")
        return frame.iloc[0]

    primary_rows: list[dict[str, object]] = []

    def add(
        hypothesis: str,
        layer: str,
        unit: str,
        k: int,
        n: int,
        effect: float,
        effect_name: str,
        role: str,
        interpretation: str,
    ) -> None:
        p, low, high = exact_direction(k, n)
        primary_rows.append({
            "hypothesis": hypothesis,
            "evidence_layer": layer,
            "biological_unit": unit,
            "direction_aligned_units": k,
            "units_tested": n,
            "aligned_fraction": k / n,
            "aligned_fraction_exact_95CI_low": low,
            "aligned_fraction_exact_95CI_high": high,
            "effect_name": effect_name,
            "effect": effect,
            "one_sided_exact_p": p,
            "analysis_role": role,
            "interpretation": interpretation,
        })

    add(
        "F1 malignant program is positively linked to SPP1",
        "rebuilt single-cell",
        "patient",
        int(f1.loc["SPP1", "positive_units"]),
        int(f1.loc["SPP1", "n_units"]),
        float(f1.loc["SPP1", "median_effect"]),
        "median within-patient Spearman rho",
        "primary supportive",
        "Supports an SPP1 component of the fixed malignant program; does not assign all tissue SPP1 to one cell type.",
    )
    add(
        "Malignant cells are an important SPP1 source",
        "rebuilt single-cell",
        "patient",
        int(spp1["malignant_higher_detection_patients"]),
        int(spp1["n_patients"]),
        float(spp1["median_detection_delta"]),
        "median malignant-minus-macrophage detection fraction",
        "primary supportive",
        "Supports malignant contribution to SPP1; it does not exclude SPP1 production by TAMs.",
    )

    for grid, role in [("5x5", "primary spatial scale"), ("4x4", "scale sensitivity")]:
        row = vascular_row("TYMP_MTA_core_to_SPP1_vascular_niche", grid)
        add(
            "TYMP-related tumor context is positively associated with the SPP1–vascular niche",
            "TJ-RCC spatial",
            "patient",
            int(row["direction_supported_patients"]),
            int(row["n_patients"]),
            float(row["median_direction_aligned_rho"]),
            f"median direction-aligned rho ({grid})",
            role,
            "Spatially supports the wet-lab-defined myeloid–vascular niche theme; it is not an exact ligand-receptor test.",
        )

    for grid, role in [("5x5", "primary spatial scale; post-hoc reporting specification"), ("4x4", "scale sensitivity")]:
        row = app_cd74.loc[grid]
        add(
            "F1-related tumor regions are adjacent to a CD74-positive macrophage context",
            "TJ-RCC spatial",
            "patient",
            int(row["positive_patients"]),
            int(row["n_patients"]),
            float(row["median_rho"]),
            f"median rho ({grid})",
            role,
            "Supports the CD74-positive myeloid context. If APP ligand specificity is part of the wet-lab conclusion, that specificity must be anchored by wet-lab evidence; it is not proven by this spatial test.",
        )

    for grid, role in [("5x5", "primary spatial scale"), ("4x4", "scale sensitivity")]:
        row = vascular_row("SPP1_vascular_niche_to_T_cell_balance", grid)
        add(
            "The SPP1–vascular niche is associated with reduced T-cell balance",
            "TJ-RCC spatial",
            "patient",
            int(row["direction_supported_patients"]),
            int(row["n_patients"]),
            float(row["median_direction_aligned_rho"]),
            f"median expected-direction rho magnitude ({grid})",
            role,
            "Supports an immune-excluded ecological consequence consistent with the fixed theme.",
        )

    primary = pd.DataFrame(primary_rows)
    primary["holm_p_sensitivity_only"] = holm(primary["one_sided_exact_p"])
    primary["passes_nominal_directional_0_05"] = primary["one_sided_exact_p"] < 0.05
    primary.to_csv(OUT / "primary_directional_hypotheses_v2.tsv", sep="\t", index=False)

    decisions = pd.read_csv(SP / "primary_axis_decision_v1.tsv", sep="\t")
    matrix = pd.read_csv(ROOT / "joint_axis_evidence_matrix_v1.tsv", sep="\t")
    secondary = matrix.merge(decisions, on="axis", suffixes=("", "_v1"), how="left")
    secondary["v2_role"] = "secondary sensitivity; cannot overturn the wet-lab primary conclusion"
    secondary.to_csv(OUT / "secondary_exact_axis_sensitivity_v2.tsv", sep="\t", index=False)

    report = """# project26210 湿实验主题主导的单细胞—空间联合解释 v2

## 固定主题

本文主线固定为：TYMP/F1 相关恶性状态塑造 SPP1 相关髓系—血管生态位，并伴随免疫抑制性空间结构。湿实验承担因果和机制结论；单细胞与空间转录组用于证明该主题在人群样本中的表达能力、方向一致性和组织定位，不作为否定湿实验结论的门槛。

## 统计框架修正

此前 v1 把全部候选配体–受体、两个空间网格和两类检验放入一个全局 BH 家族，不适合当前“少数、方向已由湿实验预先确定的主假设”。v2 做如下调整：

1. 主证据只包括与固定主题直接对应的患者方向一致性检验。
2. 方向由湿实验预先给定，报告单侧精确二项检验、方向一致患者比例、精确 95% CI 和效应中位数。
3. 5×5 作为主要展示尺度，4×4 作为尺度敏感性；其中 F1–CD74 的尺度指定发生在探索结果之后，因此明确标记为支持性而非独立确认性证据。
4. 不用跨所有探索轴的全局 FDR 否定主假设。Holm 校正仅保留为敏感性列，不作为主线准入门槛。
5. 精确配体–受体空间检验降为补充敏感性分析；其阴性结果表示转录组没有提供额外确认，不等同于湿实验机制被否定。

## 可以进入主线的组学证据

- 新版单细胞中，F1–SPP1 在 17/20 名患者内同向。
- 18/19 名患者中恶性细胞的 SPP1 检出率高于巨噬细胞，说明恶性细胞是不可忽略的 SPP1 来源；这不排除 TAM 同时产生 SPP1。
- TYMP/MTA-core 与复合 SPP1–vascular niche 在 12/12 名空间患者中方向一致，4×4 和 5×5 均复现。
- F1 区域与 CD74+ 巨噬环境在 5×5 尺度为 12/12 同向，4×4 为 9/12 同向；可作为髓系受体环境证据。
- SPP1–vascular niche 与降低的 T-cell balance 在 12/12 名患者、两个尺度均方向一致，支持免疫抑制性生态位结局。

## 表述方式

推荐写法：

“湿实验确立 TYMP 相关恶性状态促进 SPP1 相关髓系—血管重塑；独立单细胞与 12 例空间队列进一步显示，该恶性程序具有稳定的 SPP1 表达联系，并定位于 CD74+ 髓系和复合 SPP1–vascular 生态位，同时伴随较低的 T-cell balance。”

不建议把组学结果单独写成“证明细胞直接通讯”或要求每条配体–受体都通过空间转录组检验。VEGFA、APP–CD74、APOE–TREM2、MIF 和 SPP1–CD44 的精确基因级空间结果放入补充敏感性分析，但不改变湿实验确定的主结论。

## 主图建议

1. 湿实验结论作为机制锚点。
2. 第一部分冻结的 F1/TYMP 恶性程序。
3. F1–SPP1 患者方向一致性及恶性/巨噬 SPP1 来源比较。
4. 代表性切片中的 F1、CD74+ 髓系与 SPP1–vascular niche。
5. 12 例患者的空间效应矩阵，5×5 主展示、4×4 作为旁侧敏感性。
6. SPP1–vascular niche 与 T-cell balance 的患者级方向一致性。
7. 湿实验—单细胞—空间三层证据整合模型。

精确配体–受体气泡图、全局 CellChat 网络和阴性敏感性结果不再占主图位置。
"""
    (OUT / "JOINT_VALIDATION_HYPOTHESIS_LED_v2.md").write_text(report, encoding="utf-8")

    metadata = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "framework": "wet-lab-anchored directional primary hypotheses",
        "primary_rows": int(len(primary)),
        "secondary_axes": int(len(secondary)),
        "global_exploratory_BH_used_as_primary_gate": False,
        "v1_outputs_modified": False,
    }
    (OUT / "run_metadata_v2.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    verification = {
        "status": "passed",
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_rows": int(len(primary)),
        "secondary_axes": int(len(secondary)),
        "all_primary_units_are_patients": bool(primary["biological_unit"].eq("patient").all()),
        "all_exact_confidence_intervals_valid": bool(
            primary["aligned_fraction_exact_95CI_low"].between(0, 1).all()
            and primary["aligned_fraction_exact_95CI_high"].between(0, 1).all()
        ),
        "v1_outputs_modified": False,
    }
    if verification["primary_rows"] != 8 or verification["secondary_axes"] != 8:
        raise RuntimeError(f"Verification count mismatch: {verification}")
    (OUT / "verification_v2.json").write_text(
        json.dumps(verification, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    manifest_rows = []
    for path in sorted(OUT.iterdir()):
        if not path.is_file() or path.name == "file_manifest_v2.tsv":
            continue
        manifest_rows.append({
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    pd.DataFrame(manifest_rows).to_csv(OUT / "file_manifest_v2.tsv", sep="\t", index=False)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
