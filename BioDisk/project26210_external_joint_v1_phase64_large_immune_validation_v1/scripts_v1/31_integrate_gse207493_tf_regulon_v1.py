#!/usr/bin/env python3
"""Integrate fixed SCENIC and metaVIPER evidence and write a reproducible record."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
ROOT = PROJECT / "results_v1/GSE207493/tf_regulon_v1"
SCENIC = ROOT / "scenic_fixed_v1"
VIPER = ROOT / "metaviper_v1"
LOGS = PROJECT / "logs_v1"
OUT = ROOT / "integrated_v1"
METHOD_LOGS = {
    "GENIE3": LOGS / "gse207493_scenic_genie3_aucell_v1.log",
    "GRNBoost2": LOGS / "gse207493_scenic_grnboost2_aucell_v1.log",
}
OLD = {"ILF2", "HOXD3", "HNF4A", "HIF1A"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_metric(metric: str) -> str:
    return "continuous" if metric in {"cycle_adjusted_spearman", "cycle_adjusted_continuous"} else "high_minus_low"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scenic = pd.read_csv(SCENIC / "scenic_regulon_summary_full_v1.tsv", sep="\t")
    viper = pd.read_csv(VIPER / "metaviper_regulator_summary_full_v1.tsv", sep="\t",
                        dtype={"entrez_id": str})
    warning_rows = []
    low_coverage: dict[str, set[str]] = {}
    pattern = re.compile(r"genes in ([^()]+)\(\+\) are present")
    for method, log in METHOD_LOGS.items():
        tfs = []
        for line in log.read_text(errors="replace").splitlines():
            match = pattern.search(line)
            if match:
                tfs.append(match.group(1))
        low_coverage[method] = set(tfs)
        warning_rows.extend({"method": method, "tf": tf, "coverage_warning_lt_80pct": True} for tf in sorted(set(tfs)))
    coverage = pd.DataFrame(warning_rows)
    coverage.to_csv(OUT / "scenic_low_coverage_regulons_v1.tsv", sep="\t", index=False)

    scenic_long = scenic.rename(columns={
        "tf": "symbol", "median_effect": "effect", "wilcoxon_p_two_sided": "p_value",
        "BH_FDR": "q_value",
    })[["method", "endpoint", "metric", "symbol", "effect", "p_value", "q_value",
        "n_eligible_samples", "retained", "evidence_tier", "topic_related"]].copy()
    scenic_long["canonical_metric"] = scenic_long.metric.map(canonical_metric)
    scenic_long["low_coverage"] = [
        symbol in low_coverage[method] for method, symbol in zip(scenic_long.method, scenic_long.symbol)
    ]
    viper_long = viper[viper.is_TF].rename(columns={
        "median_effect": "effect", "wilcoxon_p_two_sided": "p_value", "BH_FDR": "q_value",
    })[["endpoint", "metric", "symbol", "effect", "p_value", "q_value",
        "n_eligible_samples", "retained", "evidence_tier", "topic_related"]].copy()
    viper_long.insert(0, "method", "metaVIPER")
    viper_long["canonical_metric"] = viper_long.metric.map(canonical_metric)
    viper_long["low_coverage"] = False
    long = pd.concat([scenic_long, viper_long], ignore_index=True)
    long["fdr_significant"] = long.q_value < 0.05
    long["direction"] = np.sign(long.effect)
    long.to_csv(OUT / "tf_method_evidence_long_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    rows = []
    keys = ["endpoint", "canonical_metric", "symbol"]
    for key, group in long.groupby(keys, observed=True, dropna=False):
        record = dict(zip(keys, key))
        for method in ["GENIE3", "GRNBoost2", "metaVIPER"]:
            one = group[group.method == method]
            if len(one):
                one = one.sort_values(["p_value", "q_value"], na_position="last").iloc[0]
                for column in ["effect", "p_value", "q_value", "n_eligible_samples", "retained",
                               "evidence_tier", "low_coverage", "topic_related"]:
                    record[f"{column}__{method}"] = one[column]
        rows.append(record)
    combined = pd.DataFrame(rows)
    for method in ["GENIE3", "GRNBoost2", "metaVIPER"]:
        for column, default in [("retained", False), ("low_coverage", False)]:
            name = f"{column}__{method}"
            if name not in combined:
                combined[name] = default
            combined[name] = combined[name].fillna(default).astype(bool)
        for column in ["effect", "p_value", "q_value"]:
            name = f"{column}__{method}"
            if name not in combined:
                combined[name] = np.nan

    combined["n_methods_FDR"] = sum((combined[f"q_value__{m}"] < 0.05).astype(int)
                                      for m in ["GENIE3", "GRNBoost2", "metaVIPER"])
    combined["n_methods_retained"] = sum(combined[f"retained__{m}"].astype(int)
                                           for m in ["GENIE3", "GRNBoost2", "metaVIPER"])
    viper_sign = np.sign(combined["effect__metaVIPER"])
    scenic_fdr_same = []
    scenic_retained_same = []
    for index in combined.index:
        n_fdr = 0; n_retained = 0
        for method in ["GENIE3", "GRNBoost2"]:
            if combined.at[index, f"low_coverage__{method}"]:
                continue
            same = (
                np.isfinite(combined.at[index, f"effect__{method}"])
                and np.isfinite(combined.at[index, "effect__metaVIPER"])
                and np.sign(combined.at[index, f"effect__{method}"]) == viper_sign[index]
                and viper_sign[index] != 0
            )
            if same and combined.at[index, f"q_value__{method}"] < 0.05:
                n_fdr += 1
            if same and combined.at[index, f"retained__{method}"]:
                n_retained += 1
        scenic_fdr_same.append(n_fdr)
        scenic_retained_same.append(n_retained)
    combined["n_reliable_SCENIC_FDR_same_as_VIPER"] = scenic_fdr_same
    combined["n_reliable_SCENIC_retained_same_as_VIPER"] = scenic_retained_same
    viper_fdr = combined["q_value__metaVIPER"] < 0.05
    genie_fdr_reliable = (combined["q_value__GENIE3"] < 0.05) & ~combined["low_coverage__GENIE3"]
    grn_fdr_reliable = (combined["q_value__GRNBoost2"] < 0.05) & ~combined["low_coverage__GRNBoost2"]
    scenic_same = (
        np.sign(combined["effect__GENIE3"]) == np.sign(combined["effect__GRNBoost2"])
    ) & (np.sign(combined["effect__GENIE3"]) != 0)
    combined["evidence_grade"] = np.select(
        [viper_fdr & (combined.n_reliable_SCENIC_FDR_same_as_VIPER == 2),
         viper_fdr & (combined.n_reliable_SCENIC_FDR_same_as_VIPER >= 1),
         genie_fdr_reliable & grn_fdr_reliable & scenic_same,
         viper_fdr & (combined.n_reliable_SCENIC_retained_same_as_VIPER >= 1),
         combined.n_methods_retained >= 1],
        ["A_three_method_FDR", "A_cross_platform_FDR", "B_two_SCENIC_networks_FDR",
         "B_VIPER_FDR_plus_SCENIC_retained", "C_single_platform_or_direction_conflict"],
        default="not_retained",
    )
    q_columns = [f"q_value__{m}" for m in ["GENIE3", "GRNBoost2", "metaVIPER"]]
    clipped = combined[q_columns].clip(lower=1e-300)
    combined["mean_neglog10_q"] = -np.log10(clipped).mean(axis=1, skipna=True)
    grade_order = {
        "A_three_method_FDR": 1, "A_cross_platform_FDR": 2,
        "B_two_SCENIC_networks_FDR": 3, "B_VIPER_FDR_plus_SCENIC_retained": 4,
        "C_single_platform_or_direction_conflict": 5, "not_retained": 6,
    }
    combined["grade_order"] = combined.evidence_grade.map(grade_order)
    combined = combined.sort_values(
        ["endpoint", "canonical_metric", "grade_order", "mean_neglog10_q", "symbol"],
        ascending=[True, True, True, False, True], na_position="last",
    )
    combined.to_csv(OUT / "tf_three_method_integrated_full_v1.tsv", sep="\t", index=False)
    retained = combined[combined.evidence_grade != "not_retained"].copy()
    retained.to_csv(OUT / "tf_three_method_integrated_retained_v1.tsv", sep="\t", index=False)
    combined[combined.symbol.isin(OLD)].to_csv(OUT / "old_part1_driver_integrated_replication_v1.tsv", sep="\t", index=False)

    count_rows = [
        {"metric": "target_cells", "value": 39912},
        {"metric": "samples", "value": 19},
        {"metric": "GENIE3_regulons_scored", "value": int(scenic[scenic.method == "GENIE3"].regulon.nunique())},
        {"metric": "GRNBoost2_regulons_scored", "value": int(scenic[scenic.method == "GRNBoost2"].regulon.nunique())},
        {"metric": "GENIE3_low_coverage_warning", "value": len(low_coverage["GENIE3"])},
        {"metric": "GRNBoost2_low_coverage_warning", "value": len(low_coverage["GRNBoost2"])},
        {"metric": "metaVIPER_regulators", "value": 5341},
        {"metric": "integrated_A_three_method_rows", "value": int((combined.evidence_grade == "A_three_method_FDR").sum())},
        {"metric": "integrated_A_cross_platform_rows", "value": int((combined.evidence_grade == "A_cross_platform_FDR").sum())},
    ]
    counts = pd.DataFrame(count_rows)
    counts.to_csv(OUT / "tf_regulon_result_counts_v1.tsv", sep="\t", index=False)

    top_lines = []
    for endpoint in ["tumor_f1", "myeloid_dc_tymp", "treg_tymp"]:
        top = retained[(retained.endpoint == endpoint) & retained.evidence_grade.str.startswith("A_")].head(10)
        values = ", ".join(f"{row.symbol} ({row.evidence_grade})" for row in top.itertuples())
        top_lines.append(f"- {endpoint}: {values if values else 'no grade-A TF'}")
    old_table = combined[combined.symbol.isin(OLD)][
        ["endpoint", "canonical_metric", "symbol", "effect__GENIE3", "q_value__GENIE3",
         "low_coverage__GENIE3", "effect__GRNBoost2", "q_value__GRNBoost2",
         "low_coverage__GRNBoost2", "effect__metaVIPER", "q_value__metaVIPER", "evidence_grade"]
    ]
    record = f"""# GSE207493 transcription-factor/regulon validation record v1

## Decision

Both SCENIC and metaVIPER were retained because they answer complementary questions. SCENIC tests the
original Part-1 data-derived GENIE3 and GRNBoost2 regulons by fixed-network AUCell in an external cohort;
metaVIPER tests regulator/protein activity with the independent KIRC/KIRP ARACNe interactomes. No network
was relearned in GSE207493.

## Frozen scope and statistics

- 39,912 target cells from 19 tumors: Tumor, Myeloid/Dendritic and Treg.
- Within each sample and lineage, both regulator activity and F1/TYMP exposure were residualized on S and
  G2M scores. Continuous association and upper-minus-lower quartile contrasts were then summarized with
  sample as the biological unit.
- At >=10 eligible samples BH FDR is used. Topic TFs with nominal p<0.05 and the top 15 nominal TFs per
  endpoint/method are retained as lower tiers.
- Fixed SCENIC scored all 299 GENIE3 and 314 GRNBoost2 regulons. Low-target-coverage warnings affected
  {len(low_coverage['GENIE3'])} and {len(low_coverage['GRNBoost2'])}, respectively; affected regulons cannot
  contribute to grade-A cross-platform evidence.
- metaVIPER inferred 5,341 regulators in every sample using the same KIRC/KIRP networks as rebuilt Part 1.

## Main interpretation

- ILF2 is the clearest old-Part-1 replication: positive activity with tumor F1 and with TYMP-high myeloid/DC
  and Treg states, supported by both SCENIC networks and metaVIPER in the principal endpoints.
- HOXD3 and HNF4A are not clean SCENIC external replications because their fixed regulons had poor/zero
  usable activity in this dataset. metaVIPER does estimate them, often with negative TYMP associations;
  those are method-specific findings, not confirmation of the old positive SCENIC signal.
- HIF1A is contextual rather than an old primary SCENIC driver. metaVIPER supports positive HIF1A activity
  in TYMP-high myeloid/DC but not Treg, and a negative association with tumor F1.

## Top cross-platform results (first 10 grade-A rows per endpoint)

{chr(10).join(top_lines)}

## Files

- `tf_three_method_integrated_full_v1.tsv`: complete three-method table.
- `tf_three_method_integrated_retained_v1.tsv`: FDR/topic nominal/top-15 retained evidence with grades.
- `old_part1_driver_integrated_replication_v1.tsv`: ILF2, HOXD3, HNF4A and HIF1A audit.
- `scenic_low_coverage_regulons_v1.tsv`: regulons downgraded for <80% target coverage.
"""
    (OUT / "gse207493_tf_regulon_validation_record_v1.md").write_text(record)
    status = {
        "status": "PASS_TF_REGULON_INTEGRATED_V1", "n_cells": 39912, "n_samples": 19,
        "n_integrated_rows": len(combined), "n_retained_rows": len(retained),
        "n_grade_A_three_method": int((combined.evidence_grade == "A_three_method_FDR").sum()),
        "n_grade_A_cross_platform": int((combined.evidence_grade == "A_cross_platform_FDR").sum()),
    }
    (OUT / "status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2))
    print("\n".join(top_lines))


if __name__ == "__main__":
    main()
