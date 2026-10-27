from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"D:\project26210_external_joint_v1")
PHASE = Path(__file__).resolve().parent
OUT = PHASE / "analysis_local_v1"
RNG_SEED = 2621058
N_RESAMPLES = 200
GENES = ["TYMP", "ODC1", "TK1", "SRM", "PCYT2"]

P57 = ROOT / "phase57_part1_by_cohort_signal_audit_local_v1" / "analysis_local_v1"
CELL_FILE = (
    ROOT
    / "phase31_part1_m17_enzyme_matrix_local_v1"
    / "results_local_v1"
    / "m17_cell_expression_v1.tsv.gz"
)


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(pd.to_numeric(values, errors="coerce"))
    p = values[valid].astype(float)
    if p.empty:
        return result
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result.loc[p.index.to_numpy()[order]] = adjusted
    return result


def wilcoxon_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    values = values[values != 0]
    if len(values) == 0:
        return 1.0
    return float(stats.wilcoxon(values, alternative="two-sided", method="auto").pvalue)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan, np.nan
    result = stats.spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def patient_jackknife(effects: pd.DataFrame, cohort_q: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    q_lookup = cohort_q.set_index(["signal_id", "study"])[
        "wilcoxon_fdr_bh_within_study"
    ].to_dict()
    for (signal_id, study), group in effects.groupby(["signal_id", "study"], sort=True):
        values = group["effect"].to_numpy(float)
        full_median = float(np.median(values))
        full_sign = np.sign(full_median)
        jackknife: list[tuple[str, float, float]] = []
        if len(group) >= 3:
            for idx, row in group.iterrows():
                kept = group.drop(index=idx)["effect"].to_numpy(float)
                jk_median = float(np.median(kept))
                jackknife.append(
                    (
                        str(row["curated_patient_id"]),
                        jk_median,
                        jk_median - full_median,
                    )
                )
        if jackknife:
            jk_values = np.array([x[1] for x in jackknife])
            shifts = np.array([x[2] for x in jackknife])
            max_index = int(np.argmax(np.abs(shifts)))
            sign_stability = float(np.mean(np.sign(jk_values) == full_sign))
            influential_patient = jackknife[max_index][0]
            max_shift = float(shifts[max_index])
            jk_min = float(np.min(jk_values))
            jk_max = float(np.max(jk_values))
        else:
            sign_stability = np.nan
            influential_patient = ""
            max_shift = np.nan
            jk_min = np.nan
            jk_max = np.nan
        q = float(q_lookup.get((signal_id, study), np.nan))
        if len(group) < 3:
            status = "DIRECTION_ONLY_LT3_PATIENTS"
        elif sign_stability < 1:
            status = "PATIENT_OUTLIER_SENSITIVE"
        elif np.isfinite(q) and q < 0.05:
            status = "INTERNALLY_STABLE_FDR_SUPPORTED"
        else:
            status = "DIRECTIONALLY_STABLE_NOT_FDR_SUPPORTED"
        rows.append(
            {
                "signal_id": signal_id,
                "signal_label": group.iloc[0]["signal_label"],
                "study": study,
                "n_patients": len(group),
                "n_positive": int(np.sum(values > 0)),
                "n_negative": int(np.sum(values < 0)),
                "full_median_effect": full_median,
                "full_wilcoxon_p": wilcoxon_p(values),
                "within_study_bh_q_18signals": q,
                "jackknife_median_min": jk_min,
                "jackknife_median_max": jk_max,
                "jackknife_sign_stability_fraction": sign_stability,
                "most_influential_patient": influential_patient,
                "median_shift_when_removed": max_shift,
                "status": status,
            }
        )
    return pd.DataFrame(rows)


def cell_count_dependency(effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (signal_id, study), group in effects.groupby(["signal_id", "study"], sort=True):
        if len(group) < 3:
            continue
        x = np.log10(group["n_cells"].to_numpy(float))
        y = group["effect"].to_numpy(float)
        rho, p = safe_spearman(x, y)
        order = np.argsort(x)
        split = len(group) // 2
        low = y[order[:split]]
        high = y[order[-split:]]
        mw_p = (
            float(stats.mannwhitneyu(low, high, alternative="two-sided").pvalue)
            if len(low) >= 2 and len(high) >= 2
            else np.nan
        )
        if len(group) < 5:
            risk = "SMALL_N_NOT_INTERPRETABLE"
        elif np.isfinite(rho) and abs(rho) >= 0.7:
            risk = "POSSIBLE_CELL_COUNT_DEPENDENCY"
        else:
            risk = "NO_LARGE_CELL_COUNT_ASSOCIATION"
        rows.append(
            {
                "signal_id": signal_id,
                "signal_label": group.iloc[0]["signal_label"],
                "study": study,
                "n_patients": len(group),
                "spearman_effect_vs_log10_cells": rho,
                "spearman_p": p,
                "low_cell_patients_n": len(low),
                "high_cell_patients_n": len(high),
                "low_cell_median_effect": float(np.median(low)),
                "high_cell_median_effect": float(np.median(high)),
                "high_minus_low_median": float(np.median(high) - np.median(low)),
                "mannwhitney_p": mw_p,
                "dependency_flag": risk,
            }
        )
    out = pd.DataFrame(rows)
    p_for_bh = out["spearman_p"].where(out["n_patients"] >= 5, np.nan)
    out["spearman_bh_q_across_signal_study"] = bh_adjust(p_for_bh)
    return out


def patient_rho(data: pd.DataFrame, gene: str) -> float:
    x = data["F1"].to_numpy(float)
    y = data[f"{gene}_log1p_cp10k"].to_numpy(float)
    if len(data) < 20 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan
    return float(stats.spearmanr(x, y).statistic)


def balanced_resampling(cells: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows: list[dict] = []
    for gene in GENES:
        for study, study_data in cells.groupby("study", sort=True):
            patient_groups = []
            original_effects = []
            for patient, group in study_data.groupby("curated_patient_id", sort=True):
                rho = patient_rho(group, gene)
                if np.isfinite(rho):
                    patient_groups.append((patient, group))
                    original_effects.append(rho)
            if not patient_groups:
                continue
            resampled_medians = []
            for _ in range(N_RESAMPLES):
                patient_effects = []
                for _, group in patient_groups:
                    n_take = min(100, len(group))
                    selected = rng.choice(len(group), size=n_take, replace=False)
                    rho = patient_rho(group.iloc[selected], gene)
                    if np.isfinite(rho):
                        patient_effects.append(rho)
                if patient_effects:
                    resampled_medians.append(float(np.median(patient_effects)))
            original_median = float(np.median(original_effects))
            resampled = np.asarray(resampled_medians)
            rows.append(
                {
                    "gene": gene,
                    "signal_id": f"F1__{gene}_RNA",
                    "study": study,
                    "n_eligible_patients": len(patient_groups),
                    "original_patient_equal_median": original_median,
                    "resampling_cap_cells_per_patient": 100,
                    "n_resamples": len(resampled),
                    "resampled_median_effect": float(np.median(resampled)),
                    "resampled_ci_low": float(np.quantile(resampled, 0.025)),
                    "resampled_ci_high": float(np.quantile(resampled, 0.975)),
                    "same_sign_fraction": float(
                        np.mean(np.sign(resampled) == np.sign(original_median))
                    ),
                }
            )
    return pd.DataFrame(rows)


def fixed_count_resampling(cells: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows: list[dict] = []
    for gene in GENES:
        for study, study_data in cells.groupby("study", sort=True):
            if study != "GSE207493":
                continue
            for fixed_n in [50, 100, 200]:
                groups = []
                for patient, group in study_data.groupby("curated_patient_id", sort=True):
                    if len(group) >= fixed_n and np.isfinite(patient_rho(group, gene)):
                        groups.append((patient, group))
                if len(groups) < 2:
                    continue
                medians = []
                for _ in range(N_RESAMPLES):
                    effects = []
                    for _, group in groups:
                        selected = rng.choice(len(group), size=fixed_n, replace=False)
                        rho = patient_rho(group.iloc[selected], gene)
                        if np.isfinite(rho):
                            effects.append(rho)
                    medians.append(float(np.median(effects)))
                medians = np.asarray(medians)
                rows.append(
                    {
                        "gene": gene,
                        "signal_id": f"F1__{gene}_RNA",
                        "study": study,
                        "fixed_cells_per_patient": fixed_n,
                        "n_patients": len(groups),
                        "n_resamples": len(medians),
                        "median_effect": float(np.median(medians)),
                        "ci_low": float(np.quantile(medians, 0.025)),
                        "ci_high": float(np.quantile(medians, 0.975)),
                        "positive_resample_fraction": float(np.mean(medians > 0)),
                        "negative_resample_fraction": float(np.mean(medians < 0)),
                    }
                )
    return pd.DataFrame(rows)


def gse207493_patient_detail(cells: pd.DataFrame) -> pd.DataFrame:
    data = cells[cells["study"].eq("GSE207493")]
    rows: list[dict] = []
    for patient, group in data.groupby("curated_patient_id", sort=True):
        row = {
            "study": "GSE207493",
            "curated_patient_id": patient,
            "n_cells": len(group),
            "median_F1": float(group["F1"].median()),
            "F1_nonzero_fraction": float(np.mean(group["F1"] > 0)),
            "F1_dominant_fraction": float(
                np.mean(group["dominant_factor_rank7_v1"].eq("F1"))
            ),
        }
        for gene in GENES:
            row[f"{gene}_spearman_rho"] = patient_rho(group, gene)
            row[f"{gene}_detection_fraction"] = float(
                np.mean(group[f"{gene}_raw_count"] > 0)
            )
        rows.append(row)
    return pd.DataFrame(rows)


def gse207493_summary(
    detail: pd.DataFrame,
    jackknife: pd.DataFrame,
    dependency: pd.DataFrame,
    balanced: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []
    for gene in GENES:
        signal_id = f"F1__{gene}_RNA"
        effects = detail[f"{gene}_spearman_rho"].dropna().to_numpy(float)
        jk = jackknife[
            jackknife["signal_id"].eq(signal_id)
            & jackknife["study"].eq("GSE207493")
        ]
        dep = dependency[
            dependency["signal_id"].eq(signal_id)
            & dependency["study"].eq("GSE207493")
        ]
        bal = balanced[
            balanced["signal_id"].eq(signal_id)
            & balanced["study"].eq("GSE207493")
        ]
        rows.append(
            {
                "gene": gene,
                "n_patients": len(effects),
                "n_positive": int(np.sum(effects > 0)),
                "n_negative": int(np.sum(effects < 0)),
                "median_patient_rho": float(np.median(effects)),
                "wilcoxon_p": wilcoxon_p(effects),
                "jackknife_median_min": float(jk.iloc[0]["jackknife_median_min"]),
                "jackknife_median_max": float(jk.iloc[0]["jackknife_median_max"]),
                "jackknife_sign_stability_fraction": float(
                    jk.iloc[0]["jackknife_sign_stability_fraction"]
                ),
                "effect_vs_log_cells_rho": float(
                    dep.iloc[0]["spearman_effect_vs_log10_cells"]
                ),
                "effect_vs_log_cells_p": float(dep.iloc[0]["spearman_p"]),
                "balanced100_median_rho": float(bal.iloc[0]["resampled_median_effect"]),
                "balanced100_ci_low": float(bal.iloc[0]["resampled_ci_low"]),
                "balanced100_ci_high": float(bal.iloc[0]["resampled_ci_high"]),
                "balanced100_same_sign_fraction": float(bal.iloc[0]["same_sign_fraction"]),
            }
        )
    out = pd.DataFrame(rows)
    out["wilcoxon_bh_q_5genes"] = bh_adjust(out["wilcoxon_p"])
    return out


def fmt(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def write_summary(
    gse: pd.DataFrame,
    jackknife: pd.DataFrame,
    dependency: pd.DataFrame,
) -> None:
    tymp = gse[gse["gene"].eq("TYMP")].iloc[0]
    stable = int((jackknife["status"] == "INTERNALLY_STABLE_FDR_SUPPORTED").sum())
    sensitive = int((jackknife["status"] == "PATIENT_OUTLIER_SENSITIVE").sum())
    count_flags = dependency[dependency["dependency_flag"].eq("POSSIBLE_CELL_COUNT_DEPENDENCY")]
    count_flags_fdr = count_flags[count_flags["spearman_bh_q_across_signal_study"] < 0.05]
    lines = [
        "# 第一部分额外单队列深挖 v1",
        "",
        "## 分析内容",
        "",
        "本阶段在Phase 57逐队列审计基础上增加三类单队列内部检查：患者留一法、患者效应与细胞数关系，以及对TYMP/ODC1/TK1/SRM/PCYT2五个RNA端点的每患者最多100细胞平衡重抽样。所有推断仍以患者为单位；细胞重抽样只用于检查单个大细胞量患者是否造成效应。",
        "",
        "## GSE207493核心结果",
        "",
        f"- F1–TYMP在10名患者中10/10正向，中位rho={fmt(tymp['median_patient_rho'])}，五基因内BH q={fmt(tymp['wilcoxon_bh_q_5genes'], 4)}。",
        f"- 逐一删除任一患者后，中位rho范围为{fmt(tymp['jackknife_median_min'])}至{fmt(tymp['jackknife_median_max'])}，方向稳定率={fmt(tymp['jackknife_sign_stability_fraction'])}。",
        f"- 将每名患者最多平衡抽取100个细胞、重复{N_RESAMPLES}次后，中位rho={fmt(tymp['balanced100_median_rho'])}，95%重抽样区间{fmt(tymp['balanced100_ci_low'])}至{fmt(tymp['balanced100_ci_high'])}，同方向比例={fmt(tymp['balanced100_same_sign_fraction'])}。",
        f"- 患者效应与log10细胞数的Spearman rho={fmt(tymp['effect_vs_log_cells_rho'])}，P={fmt(tymp['effect_vs_log_cells_p'], 4)}；未见效应由大细胞量患者单调拉高的证据。",
        "",
        "因此，GSE207493内部的F1–TYMP RNA关系不是由单一患者或少数超大细胞量患者造成；它是该队列内较稳定的共同方向。跨队列外推限制仍然存在，因为其他队列严格恶性患者数较少，且F1主导状态的细胞组成高度偏向GSE207493。",
        "",
        "## 全部单队列稳定性概览",
        "",
        f"- 在88个“信号–队列”组合中，{stable}个同时通过队列内18信号BH校正并具有100%患者留一方向稳定性。",
        f"- {sensitive}个组合在删除某名患者后中位方向改变，被标记为患者异常值敏感。",
        f"- 效应与细胞数的探索性筛查产生{len(count_flags)}个|rho|≥0.7的提示，其中{len(count_flags_fdr)}个通过全表BH校正；这些只作为潜在精度/组成依赖提示，不作为生物学结论。",
        "",
        "## 解释",
        "",
        "新增分析支持把问题拆成两层：GSE207493内部的核心F1–TYMP关系稳定；但严格TYMP功能轴及部分调控端点是否能跨队列推广，仍取决于其他队列能否提供足够的严格恶性患者，而不是继续增加GSE207493细胞数。",
        "",
    ]
    (PHASE / "SINGLE_COHORT_DEEP_DIVE_SUMMARY_ZH_v1.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    effects = pd.read_csv(P57 / "patient_effects_long_v1.tsv", sep="\t")
    cohort_q = pd.read_csv(P57 / "per_cohort_effects_v1.tsv", sep="\t")
    usecols = [
        "study",
        "curated_patient_id",
        "F1",
        "dominant_factor_rank7_v1",
        *[f"{gene}_raw_count" for gene in GENES],
        *[f"{gene}_log1p_cp10k" for gene in GENES],
    ]
    cells = pd.read_csv(CELL_FILE, sep="\t", usecols=usecols)

    jackknife = patient_jackknife(effects, cohort_q)
    dependency = cell_count_dependency(effects)
    balanced = balanced_resampling(cells, rng)
    fixed = fixed_count_resampling(cells, rng)
    detail = gse207493_patient_detail(cells)
    gse_summary = gse207493_summary(detail, jackknife, dependency, balanced)

    jackknife.to_csv(OUT / "within_cohort_patient_jackknife_v1.tsv", sep="\t", index=False)
    dependency.to_csv(OUT / "effect_vs_patient_cell_count_v1.tsv", sep="\t", index=False)
    balanced.to_csv(OUT / "balanced_cell_resampling_v1.tsv", sep="\t", index=False)
    fixed.to_csv(OUT / "fixed_cell_count_resampling_v1.tsv", sep="\t", index=False)
    detail.to_csv(OUT / "gse207493_patient_rna_detail_v1.tsv", sep="\t", index=False)
    gse_summary.to_csv(OUT / "gse207493_internal_summary_v1.tsv", sep="\t", index=False)
    write_summary(gse_summary, jackknife, dependency)

    checks = [
        ("phase57_signal_count", effects["signal_id"].nunique() == 18),
        ("phase57_patient_effect_rows", len(effects) == 355),
        ("cell_count", len(cells) == 11950),
        ("five_studies", cells["study"].nunique() == 5),
        ("gse207493_twelve_strict_patients", len(detail) == 12),
        ("gse207493_ten_tymp_eligible_patients", int(detail["TYMP_spearman_rho"].notna().sum()) == 10),
        ("gse207493_five_gene_summary", len(gse_summary) == 5),
        ("tymp_all_gse207493_positive", int(gse_summary.loc[gse_summary.gene.eq("TYMP"), "n_positive"].iloc[0]) == 10),
        ("tymp_jackknife_sign_stable", float(gse_summary.loc[gse_summary.gene.eq("TYMP"), "jackknife_sign_stability_fraction"].iloc[0]) == 1.0),
        ("tymp_balanced_sign_stable", float(gse_summary.loc[gse_summary.gene.eq("TYMP"), "balanced100_same_sign_fraction"].iloc[0]) == 1.0),
        ("finite_balanced_results", bool(np.isfinite(balanced["resampled_median_effect"]).all())),
    ]
    validation = pd.DataFrame(checks, columns=["check", "pass"])
    validation.to_csv(PHASE / "validation_checks_v1.tsv", sep="\t", index=False)
    if not validation["pass"].all():
        raise RuntimeError("Phase 58 validation failed")

    files = [*OUT.glob("*"), PHASE / "SINGLE_COHORT_DEEP_DIVE_SUMMARY_ZH_v1.md", PHASE / "validation_checks_v1.tsv"]
    manifest = []
    for path in sorted(files):
        manifest.append(
            {
                "file": str(path.relative_to(PHASE)),
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    pd.DataFrame(manifest).to_csv(PHASE / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    print(
        f"PASS: {len(jackknife)} signal-study jackknife rows; "
        f"{len(balanced)} balanced-resampling rows; {len(detail)} GSE207493 patients"
    )


if __name__ == "__main__":
    main()
