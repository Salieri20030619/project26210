#!/usr/bin/env python3
"""Finalize patient-level de novo CollecTRI-VIPER evidence."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
ROOT = PROJECT / "results_v1/GSE207493/tf_regulon_de_novo_v1/collectri_viper_v1"
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
    if not len(p):
        return out
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted); restored[order] = np.minimum(adjusted, 1)
    out.loc[mask] = restored
    return out


def wilcox(values: pd.Series) -> float:
    x = values.dropna().to_numpy(float)
    if len(x) < 3 or np.allclose(x, 0):
        return np.nan
    return float(stats.wilcoxon(x, method="approx").pvalue)


def main() -> None:
    effects = pd.read_csv(ROOT / "collectri_viper_effects_by_sample_v1.tsv.gz", sep="\t")
    keys = ["endpoint", "metric", "source"]
    summary = effects.groupby(keys, observed=True)["effect"].agg(
        n_eligible_samples="size", median_effect="median", mean_effect="mean",
        positive_fraction=lambda x: float(np.mean(x > 0)), wilcoxon_p_two_sided=wilcox,
    ).reset_index()
    summary["BH_FDR"] = np.nan
    for _, idx in summary[summary.n_eligible_samples >= 10].groupby(
        ["endpoint", "metric"], observed=True
    ).groups.items():
        summary.loc[idx, "BH_FDR"] = bh(summary.loc[idx, "wilcoxon_p_two_sided"])
    summary["topic_related"] = summary.source.isin(TOPIC)
    summary["nominal_rank"] = summary.groupby(
        ["endpoint", "metric"], observed=True
    ).wilcoxon_p_two_sided.rank(method="min")
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
    summary["abs_median_effect"] = summary.median_effect.abs()
    summary = summary.sort_values(
        ["endpoint", "metric", "wilcoxon_p_two_sided", "abs_median_effect"],
        ascending=[True, True, True, False], na_position="last",
    )
    summary.to_csv(ROOT / "collectri_viper_summary_full_v1.tsv", sep="\t", index=False)
    summary[summary.retained].to_csv(ROOT / "collectri_viper_summary_retained_v1.tsv", sep="\t", index=False)
    summary[summary.topic_related].to_csv(ROOT / "collectri_viper_topic_TFs_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_COLLECTRI_VIPER_DE_NOVO_V1", "n_samples": int(effects["sample"].nunique()),
        "n_sources": int(effects.source.nunique()), "n_summary_rows": len(summary),
        "n_retained_rows": int(summary.retained.sum()), "old_part1_network_used": False,
    }
    (ROOT / "status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    (ROOT / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
