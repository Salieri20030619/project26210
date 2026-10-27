#!/usr/bin/env python3
"""Finalize patient-level metaVIPER tables from completed per-sample effects."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
OUT = PROJECT / "results_v1/GSE207493/tf_regulon_v1/metaviper_v1"
TF_FILE = Path("/mnt/BioDisk/ludingyu/project26210/references/pyscenic/hg38/allTFs_hg38.txt")
TOPIC = {
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


def safe_wilcox(values: pd.Series) -> float:
    x = values.dropna().to_numpy(float)
    if len(x) < 3 or np.allclose(x, 0):
        return np.nan
    try:
        return float(stats.wilcoxon(x, alternative="two-sided", method="approx").pvalue)
    except ValueError:
        return np.nan


def main() -> None:
    effects = pd.read_csv(OUT / "metaviper_effects_by_sample_v1.tsv.gz", sep="\t",
                          dtype={"entrez_id": str})
    if effects["sample"].nunique() != 19 or effects["entrez_id"].nunique() != 5341:
        raise RuntimeError("Incomplete metaVIPER per-sample effects")
    keys = ["endpoint", "metric", "entrez_id", "symbol"]
    summary = effects.groupby(keys, observed=True, dropna=False)["effect"].agg(
        n_eligible_samples="size", median_effect="median", mean_effect="mean",
        positive_fraction=lambda x: float(np.mean(x > 0)), wilcoxon_p_two_sided=safe_wilcox,
    ).reset_index()
    summary["BH_FDR"] = np.nan
    for _, index in summary[summary.n_eligible_samples >= 10].groupby(
        ["endpoint", "metric"], observed=True
    ).groups.items():
        summary.loc[index, "BH_FDR"] = bh(summary.loc[index, "wilcoxon_p_two_sided"])
    tf_symbols = {line.strip() for line in TF_FILE.read_text().splitlines() if line.strip()}
    summary["is_TF"] = summary.symbol.isin(tf_symbols)
    summary["topic_related"] = summary.symbol.isin(TOPIC)
    summary["nominal_rank_all"] = summary.groupby(
        ["endpoint", "metric"], observed=True
    )["wilcoxon_p_two_sided"].rank(method="min")
    summary["nominal_rank_TF"] = summary["wilcoxon_p_two_sided"].where(summary.is_TF).groupby(
        [summary.endpoint, summary.metric], observed=True
    ).rank(method="min")
    summary["retained"] = summary.is_TF & (
        (summary.BH_FDR < 0.05)
        | (summary.topic_related & (summary.wilcoxon_p_two_sided < 0.05))
        | (summary.nominal_rank_TF <= 15)
    )
    summary["evidence_tier"] = np.select(
        [summary.BH_FDR < 0.05,
         summary.topic_related & (summary.wilcoxon_p_two_sided < 0.05),
         summary.is_TF & (summary.nominal_rank_TF <= 15)],
        ["FDR_significant", "topic_nominal_p_lt_0.05", "top15_nominal_TF"],
        default="not_retained",
    )
    summary["abs_median_effect"] = summary.median_effect.abs()
    summary = summary.sort_values(
        ["endpoint", "metric", "wilcoxon_p_two_sided", "abs_median_effect"],
        ascending=[True, True, True, False], na_position="last",
    )
    summary.to_csv(OUT / "metaviper_regulator_summary_full_v1.tsv", sep="\t", index=False)
    summary[summary.retained].to_csv(OUT / "metaviper_tf_summary_retained_v1.tsv", sep="\t", index=False)
    summary[summary.topic_related].to_csv(OUT / "metaviper_topic_regulators_v1.tsv", sep="\t", index=False)
    old = summary[summary.symbol.isin({"ILF2", "HOXD3", "HNF4A", "HIF1A"})].copy()
    old["prespecified_role"] = np.where(
        old.symbol == "HIF1A", "contextual", "old_part1_SCENIC_driver"
    )
    old.to_csv(OUT / "metaviper_old_part1_driver_replication_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_METAVIPER_EXTERNAL_VALIDATION",
        "n_cells": 39912, "n_samples": 19, "n_regulators": 5341,
        "n_summary_rows": len(summary), "n_retained_TF_rows": int(summary.retained.sum()),
        "network": "aracne.networks KIRC plus KIRP; identical to rebuilt Part 1",
        "cycle_adjustment": "within-sample residualization on S_score and G2M_score",
        "finalizer_note": "R completed all 19 checkpoints/effects; Python finalized after R setorder(abs()) syntax error",
    }
    (OUT / "status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    (OUT / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
