#!/usr/bin/env python3
"""Patient-aware analysis of fixed Part-1 SCENIC regulons in GSE207493."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
INPUT = PROJECT / "inputs_v1/GSE207493/tf_regulon_v1"
OUT = PROJECT / "results_v1/GSE207493/tf_regulon_v1/scenic_fixed_v1"
INTERMEDIATE = OUT / "intermediate_v1"
METHODS = {
    "GENIE3": INTERMEDIATE / "auc_genie3_fixed_external_v1.csv",
    "GRNBoost2": INTERMEDIATE / "auc_grnboost2_fixed_external_v1.csv",
}
ENDPOINTS = {
    "tumor_f1": ({"Tumor"}, "f1_score"),
    "myeloid_dc_tymp": ({"Myeloid", "Dendritic"}, "TYMP_log1p"),
    "treg_tymp": ({"Treg"}, "TYMP_log1p"),
}
TOPIC_TFS = {
    "ILF2", "HOXD3", "HNF4A", "HIF1A", "EPAS1", "NFE2L2", "SREBF1", "SREBF2",
    "MYC", "PPARA", "PPARG", "STAT1", "STAT3", "IRF1", "IRF5", "RELA", "RELB",
    "NFKB1", "NFKB2", "FOXO1", "FOXO3", "TFEB", "MITF", "ATF4", "XBP1", "CREB1",
    "SP1", "ETS1", "BATF", "NFATC1", "NFATC2", "TOX", "FOXP3", "IKZF2", "E2F1", "E2F8",
}


def bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.notna() & np.isfinite(values)
    p = values.loc[finite].to_numpy(float)
    if not len(p):
        return out
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out.loc[finite] = restored
    return out


def residualize(values: np.ndarray, s_score: np.ndarray, g_score: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(s_score)), s_score, g_score])
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    return values - design @ coefficients


def safe_wilcox(values: pd.Series) -> float:
    x = values.dropna().to_numpy(float)
    if len(x) < 3 or np.allclose(x, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(x, alternative="two-sided", method="approx").pvalue)
    except ValueError:
        return np.nan


def analyze_method(method: str, path: Path, meta: pd.DataFrame) -> pd.DataFrame:
    auc = pd.read_csv(path)
    cell_column = auc.columns[0]
    auc = auc.rename(columns={cell_column: "cell_barcode"}).set_index("cell_barcode")
    if auc.index.has_duplicates or set(auc.index) != set(meta.index):
        raise RuntimeError(f"{method}: AUCell cell identifiers mismatch")
    auc = auc.loc[meta.index]
    regulons = auc.columns.astype(str).tolist()
    matrix = auc.to_numpy(np.float64)
    rows: list[pd.DataFrame] = []
    for endpoint, (cell_types, exposure_column) in ENDPOINTS.items():
        for sample, sample_meta in meta.groupby("sample", sort=True):
            use = sample_meta["cell_type"].isin(cell_types).to_numpy()
            if use.sum() < 40:
                continue
            one_meta = sample_meta.loc[use]
            positions = meta.index.get_indexer(one_meta.index)
            exposure = one_meta[exposure_column].to_numpy(float)
            if np.std(exposure) == 0:
                continue
            s_score = one_meta["S_score"].to_numpy(float)
            g_score = one_meta["G2M_score"].to_numpy(float)
            exposure_resid = residualize(exposure[:, None], s_score, g_score).ravel()
            activity_resid = residualize(matrix[positions, :], s_score, g_score)
            exposure_rank = stats.rankdata(exposure_resid)
            activity_rank = np.apply_along_axis(stats.rankdata, 0, activity_resid)
            continuous = np.corrcoef(exposure_rank, activity_rank, rowvar=False)[0, 1:]
            low = exposure_resid <= np.quantile(exposure_resid, 0.25)
            high = exposure_resid >= np.quantile(exposure_resid, 0.75)
            contrast = activity_resid[high].mean(axis=0) - activity_resid[low].mean(axis=0)
            for metric, effect in [
                ("cycle_adjusted_spearman", continuous),
                ("cycle_adjusted_high_minus_low", contrast),
            ]:
                rows.append(pd.DataFrame({
                    "method": method, "sample": sample, "endpoint": endpoint, "metric": metric,
                    "regulon": regulons, "effect": effect, "n_cells": int(use.sum()),
                    "n_low": int(low.sum()), "n_high": int(high.sum()),
                }))
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "RUN_COMPLETE").exists():
        print("GSE207493_FIXED_SCENIC_V1_REUSED")
        return
    missing = [str(path) for path in METHODS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(missing)
    meta = pd.read_csv(INPUT / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"cell_barcode": str})
    meta = meta.set_index("cell_barcode", drop=False)
    all_effects = []
    for method, path in METHODS.items():
        print(f"[SCENIC analysis] {method}", flush=True)
        all_effects.append(analyze_method(method, path, meta))
    effects = pd.concat(all_effects, ignore_index=True)
    effects.to_csv(OUT / "scenic_effects_by_sample_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    group_columns = ["method", "endpoint", "metric", "regulon"]
    summary = effects.groupby(group_columns, observed=True)["effect"].agg(
        n_eligible_samples="size", median_effect="median", mean_effect="mean",
        positive_fraction=lambda x: float(np.mean(x > 0)), wilcoxon_p_two_sided=safe_wilcox,
    ).reset_index()
    summary["tf"] = summary["regulon"].str.replace(r"\([+-]\)$", "", regex=True)
    summary["BH_FDR"] = np.nan
    for _, index in summary[summary.n_eligible_samples >= 10].groupby(
        ["method", "endpoint", "metric"], observed=True
    ).groups.items():
        summary.loc[index, "BH_FDR"] = bh(summary.loc[index, "wilcoxon_p_two_sided"])
    summary["topic_related"] = summary["tf"].isin(TOPIC_TFS)
    summary["nominal_rank"] = summary.groupby(
        ["method", "endpoint", "metric"], observed=True
    )["wilcoxon_p_two_sided"].rank(method="min")
    summary["retained"] = (
        (summary.BH_FDR < 0.05)
        | (summary.topic_related & (summary.wilcoxon_p_two_sided < 0.05))
        | (summary.nominal_rank <= 15)
    )
    summary["evidence_tier"] = np.select(
        [summary.BH_FDR < 0.05,
         summary.topic_related & (summary.wilcoxon_p_two_sided < 0.05),
         summary.nominal_rank <= 15],
        ["FDR_significant", "topic_nominal_p_lt_0.05", "top15_nominal"],
        default="not_retained",
    )
    summary = summary.sort_values(
        ["method", "endpoint", "metric", "wilcoxon_p_two_sided", "median_effect"],
        ascending=[True, True, True, True, False], na_position="last",
    )
    summary.to_csv(OUT / "scenic_regulon_summary_full_v1.tsv", sep="\t", index=False)
    summary[summary.retained].to_csv(OUT / "scenic_regulon_summary_retained_v1.tsv", sep="\t", index=False)
    summary[summary.tf.isin(["ILF2", "HOXD3", "HNF4A", "HIF1A"])].to_csv(
        OUT / "scenic_old_part1_driver_replication_v1.tsv", sep="\t", index=False
    )

    pivot = summary.pivot_table(
        index=["endpoint", "metric", "tf"], columns="method",
        values=["median_effect", "wilcoxon_p_two_sided", "BH_FDR", "retained"], aggfunc="first",
    )
    pivot.columns = [f"{value}__{method}" for value, method in pivot.columns]
    consensus = pivot.reset_index()
    effect_columns = [column for column in consensus if column.startswith("median_effect__")]
    consensus["same_direction"] = np.sign(consensus[effect_columns[0]]) == np.sign(consensus[effect_columns[1]])
    consensus["n_methods_fdr"] = sum(
        (consensus[column] < 0.05).astype(int) for column in consensus if column.startswith("BH_FDR__")
    )
    consensus["n_methods_retained"] = sum(
        consensus[column].fillna(False).astype(int) for column in consensus if column.startswith("retained__")
    )
    consensus["consensus_tier"] = np.select(
        [consensus.same_direction & (consensus.n_methods_fdr == 2),
         consensus.same_direction & (consensus.n_methods_retained == 2),
         consensus.n_methods_retained >= 1],
        ["two_method_FDR", "two_method_retained", "single_method_retained"],
        default="not_retained",
    )
    consensus.to_csv(OUT / "scenic_cross_method_consensus_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_FIXED_SCENIC_EXTERNAL_VALIDATION",
        "n_cells": len(meta), "n_methods": len(METHODS), "n_summary_rows": len(summary),
        "n_retained": int(summary.retained.sum()),
        "network_policy": "reuse original Part-1 GENIE3 and GRNBoost2 regulons; AUCell only",
        "cycle_adjustment": "residualize activity and exposure on S_score and G2M_score within sample",
    }
    (OUT / "status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    (OUT / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
