#!/usr/bin/env python3
"""Patient-aware analysis of lineage-specific de novo SCENIC regulons."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
INPUT = PROJECT / "inputs_v1/GSE207493/tf_regulon_de_novo_v1/scenic_inputs_v1"
ROOT = PROJECT / "results_v1/GSE207493/tf_regulon_de_novo_v1/scenic_de_novo_v1"
LINEAGES = {
    "tumor": ("tumor_f1", "f1_score"),
    "myeloid_dc": ("myeloid_dc_tymp", "TYMP_log1p"),
    "treg": ("treg_tymp", "TYMP_log1p"),
}
TOPIC = {
    "HIF1A", "EPAS1", "NFE2L2", "SREBF1", "SREBF2", "MYC", "PPARA", "PPARG",
    "STAT1", "STAT2", "STAT3", "IRF1", "IRF3", "IRF5", "IRF7", "IRF9",
    "RELA", "RELB", "NFKB1", "NFKB2", "JUN", "JUNB", "JUND", "FOS", "FOSL1",
    "FOXO1", "FOXO3", "TFEB", "MITF", "ATF3", "ATF4", "XBP1", "CREB1",
    "ETS1", "BATF", "BATF3", "NFATC1", "NFATC2", "TOX", "FOXP3", "IKZF2",
    "E2F1", "E2F8", "CEBPB", "CEBPD", "SPI1", "HNF1B", "HNF4A",
}


def bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    mask = values.notna() & np.isfinite(values)
    p = values.loc[mask].to_numpy(float)
    if not len(p): return out
    order = np.argsort(p); adjusted = p[order] * len(p) / np.arange(1, len(p)+1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted); restored[order] = np.minimum(adjusted, 1)
    out.loc[mask] = restored
    return out


def residualize(values: np.ndarray, s: np.ndarray, g: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(s)), s, g])
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def wilcox(values: pd.Series) -> float:
    x = values.dropna().to_numpy(float)
    if len(x) < 3 or np.allclose(x, 0): return np.nan
    return float(stats.wilcoxon(x, method="approx").pvalue)


def main() -> None:
    rows = []
    for lineage, (endpoint, exposure_column) in LINEAGES.items():
        auc_path = ROOT / lineage / "auc_grnboost2_all_cells_v1.csv"
        meta = pd.read_csv(INPUT / lineage / "cell_metadata_v1.tsv.gz", sep="\t",
                           dtype={"cell_barcode": str}).set_index("cell_barcode", drop=False)
        auc = pd.read_csv(auc_path).rename(columns=lambda x: "cell_barcode" if x == "Cell" else x)
        if "cell_barcode" not in auc.columns:
            auc = auc.rename(columns={auc.columns[0]: "cell_barcode"})
        auc = auc.set_index("cell_barcode").loc[meta.index]
        regulons = auc.columns.astype(str).tolist(); matrix = auc.to_numpy(float)
        for sample, one_meta in meta.groupby("sample", sort=True):
            if len(one_meta) < 40 or one_meta[exposure_column].std() == 0: continue
            positions = meta.index.get_indexer(one_meta.index)
            s = one_meta.S_score.to_numpy(float); g = one_meta.G2M_score.to_numpy(float)
            exposure = residualize(one_meta[[exposure_column]].to_numpy(float), s, g).ravel()
            activity = residualize(matrix[positions], s, g)
            xrank = stats.rankdata(exposure)
            yrank = np.apply_along_axis(stats.rankdata, 0, activity)
            continuous = np.corrcoef(xrank, yrank, rowvar=False)[0, 1:]
            low = exposure <= np.quantile(exposure, .25); high = exposure >= np.quantile(exposure, .75)
            contrast = activity[high].mean(0) - activity[low].mean(0)
            for metric, effect in [("cycle_adjusted_spearman", continuous),
                                   ("cycle_adjusted_high_minus_low", contrast)]:
                rows.append(pd.DataFrame({
                    "lineage": lineage, "sample": sample, "endpoint": endpoint, "metric": metric,
                    "regulon": regulons, "effect": effect, "n_cells": len(one_meta),
                    "n_low": int(low.sum()), "n_high": int(high.sum()),
                }))
    effects = pd.concat(rows, ignore_index=True)
    effects.to_csv(ROOT / "de_novo_scenic_effects_by_sample_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    keys = ["lineage", "endpoint", "metric", "regulon"]
    summary = effects.groupby(keys, observed=True).effect.agg(
        n_eligible_samples="size", median_effect="median", mean_effect="mean",
        positive_fraction=lambda x: float(np.mean(x > 0)), wilcoxon_p_two_sided=wilcox,
    ).reset_index()
    summary["tf"] = summary.regulon.str.replace(r"\([+-]\)$", "", regex=True)
    summary["BH_FDR"] = np.nan
    for _, idx in summary[summary.n_eligible_samples >= 10].groupby(
        ["endpoint", "metric"], observed=True
    ).groups.items():
        summary.loc[idx, "BH_FDR"] = bh(summary.loc[idx, "wilcoxon_p_two_sided"])
    summary["topic_related"] = summary.tf.isin(TOPIC)
    summary["nominal_rank"] = summary.groupby(
        ["endpoint", "metric"], observed=True
    ).wilcoxon_p_two_sided.rank(method="min")
    summary["retained"] = (
        (summary.BH_FDR < .05)
        | (summary.topic_related & (summary.wilcoxon_p_two_sided < .05))
        | (summary.nominal_rank <= 15)
    )
    summary["evidence_tier"] = np.select(
        [summary.BH_FDR < .05,
         summary.topic_related & (summary.wilcoxon_p_two_sided < .05),
         summary.nominal_rank <= 15],
        ["FDR_significant", "topic_nominal_p_lt_0.05", "top15_nominal"],
        default="not_retained",
    )
    summary["abs_median_effect"] = summary.median_effect.abs()
    summary = summary.sort_values(
        ["endpoint", "metric", "wilcoxon_p_two_sided", "abs_median_effect"],
        ascending=[True, True, True, False], na_position="last",
    )
    summary.to_csv(ROOT / "de_novo_scenic_summary_full_v1.tsv", sep="\t", index=False)
    summary[summary.retained].to_csv(ROOT / "de_novo_scenic_summary_retained_v1.tsv", sep="\t", index=False)
    summary[summary.topic_related].to_csv(ROOT / "de_novo_scenic_topic_TFs_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_DE_NOVO_SCENIC_ANALYSIS_V1", "n_cells": 39912,
        "n_regulons": {lineage: int(summary[summary.lineage == lineage].regulon.nunique()) for lineage in LINEAGES},
        "n_summary_rows": len(summary), "n_retained_rows": int(summary.retained.sum()),
        "old_part1_regulons_used": False,
    }
    (ROOT / "analysis_status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    (ROOT / "ANALYSIS_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
