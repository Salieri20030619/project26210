#!/usr/bin/env python3
"""Integrate primary de novo SCENIC with independent CollecTRI-VIPER evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
ROOT = PROJECT / "results_v1/GSE207493/tf_regulon_de_novo_v1"
SCENIC = ROOT / "scenic_de_novo_v1/de_novo_scenic_summary_full_v1.tsv"
VIPER = ROOT / "collectri_viper_v1/collectri_viper_summary_full_v1.tsv"
OUT = ROOT / "integrated_primary_v1"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def metric_name(value: str) -> str:
    return "continuous" if value in {"cycle_adjusted_spearman", "cycle_adjusted_continuous"} else "high_minus_low"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scenic = pd.read_csv(SCENIC, sep="\t")
    viper = pd.read_csv(VIPER, sep="\t")
    scenic["canonical_metric"] = scenic.metric.map(metric_name)
    viper["canonical_metric"] = viper.metric.map(metric_name)
    scenic_one = scenic.rename(columns={
        "tf": "TF", "median_effect": "effect_SCENIC", "wilcoxon_p_two_sided": "p_SCENIC",
        "BH_FDR": "q_SCENIC", "retained": "retained_SCENIC", "evidence_tier": "tier_SCENIC",
        "n_eligible_samples": "n_SCENIC",
    })[["endpoint", "canonical_metric", "TF", "effect_SCENIC", "p_SCENIC", "q_SCENIC",
        "n_SCENIC", "retained_SCENIC", "tier_SCENIC", "topic_related"]]
    viper_one = viper.rename(columns={
        "source": "TF", "median_effect": "effect_VIPER", "wilcoxon_p_two_sided": "p_VIPER",
        "BH_FDR": "q_VIPER", "retained": "retained_VIPER", "evidence_tier": "tier_VIPER",
        "n_eligible_samples": "n_VIPER",
    })[["endpoint", "canonical_metric", "TF", "effect_VIPER", "p_VIPER", "q_VIPER",
        "n_VIPER", "retained_VIPER", "tier_VIPER", "topic_related"]].rename(
            columns={"topic_related": "topic_related_VIPER"}
        )
    scenic_one = scenic_one.rename(columns={"topic_related": "topic_related_SCENIC"})
    merged = scenic_one.merge(viper_one, on=["endpoint", "canonical_metric", "TF"], how="outer")
    for column in ["retained_SCENIC", "retained_VIPER", "topic_related_SCENIC", "topic_related_VIPER"]:
        merged[column] = merged[column].fillna(False).astype(bool)
    same_direction = (
        np.isfinite(merged.effect_SCENIC) & np.isfinite(merged.effect_VIPER)
        & (np.sign(merged.effect_SCENIC) == np.sign(merged.effect_VIPER))
        & (np.sign(merged.effect_SCENIC) != 0)
    )
    merged["same_direction"] = same_direction
    scenic_fdr = merged.q_SCENIC < .05
    viper_fdr = merged.q_VIPER < .05
    merged["evidence_grade"] = np.select(
        [same_direction & scenic_fdr & viper_fdr,
         same_direction & ((scenic_fdr & merged.retained_VIPER) | (viper_fdr & merged.retained_SCENIC)),
         same_direction & merged.retained_SCENIC & merged.retained_VIPER,
         scenic_fdr | viper_fdr,
         merged.retained_SCENIC | merged.retained_VIPER],
        ["A_cross_method_FDR", "B_one_FDR_one_retained", "B_both_retained",
         "C_single_method_FDR_or_conflict", "D_single_method_retained"],
        default="not_retained",
    )
    merged["mean_neglog10_q"] = -np.log10(merged[["q_SCENIC", "q_VIPER"]].clip(lower=1e-300)).mean(axis=1)
    order = {
        "A_cross_method_FDR": 1, "B_one_FDR_one_retained": 2, "B_both_retained": 3,
        "C_single_method_FDR_or_conflict": 4, "D_single_method_retained": 5, "not_retained": 6,
    }
    merged["grade_order"] = merged.evidence_grade.map(order)
    merged = merged.sort_values(
        ["endpoint", "canonical_metric", "grade_order", "mean_neglog10_q", "TF"],
        ascending=[True, True, True, False, True], na_position="last",
    )
    merged.to_csv(OUT / "de_novo_tf_integrated_full_v1.tsv", sep="\t", index=False)
    retained = merged[merged.evidence_grade != "not_retained"].copy()
    retained.to_csv(OUT / "de_novo_tf_integrated_retained_v1.tsv", sep="\t", index=False)
    grade_a = merged[merged.evidence_grade == "A_cross_method_FDR"].copy()
    grade_a.to_csv(OUT / "de_novo_tf_cross_method_FDR_v1.tsv", sep="\t", index=False)

    top_lines = []
    for endpoint in ["tumor_f1", "myeloid_dc_tymp", "treg_tymp"]:
        one = grade_a[grade_a.endpoint == endpoint].drop_duplicates("TF").head(15)
        values = ", ".join(one.TF.tolist()) if len(one) else "none"
        top_lines.append(f"- {endpoint}: {values}")
    record = f"""# GSE207493 de novo TF/regulon primary analysis v1

## Primary-analysis reset

The prior fixed Part-1 regulon/KIRC-KIRP analysis is superseded as a primary analysis and retained only as
a sensitivity audit. This analysis was restarted without old Part-1 TF candidates or networks.

## Methods

- SCENIC networks were relearned separately from raw counts in Tumor, Myeloid/DC and Treg using GRNBoost2,
  sample-balanced cells, hg38 cisTarget motif pruning and lineage-wide AUCell.
- VIPER used the independently frozen public CollecTRI human network (43,159 edges, 1,186 sources), not
  the old KIRC/KIRP interactomes. CollecTRI was fetched through the required server proxy.
- Both methods used within-sample S/G2M adjustment. Patient/sample is the biological unit.
- FDR-significant results are primary; theme-related nominal p<0.05 and top-15 nominal results are lower tiers.

## Cross-method grade-A TFs (first 15 unique TFs per endpoint)

{chr(10).join(top_lines)}

Grade A requires SCENIC and CollecTRI-VIPER to have the same direction and BH FDR below 0.05.
"""
    (OUT / "gse207493_de_novo_tf_regulon_record_v1.md").write_text(record)
    status = {
        "status": "PASS_DE_NOVO_TF_REGULON_PRIMARY_V1", "n_cells": 39912, "n_samples": 19,
        "n_integrated_rows": len(merged), "n_retained_rows": len(retained),
        "n_grade_A_cross_method_FDR": len(grade_a), "old_part1_primary_used": False,
    }
    (OUT / "status_v1.json").write_text(json.dumps(status, indent=2) + "\n")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text("PASS\n")
    print(json.dumps(status, indent=2)); print("\n".join(top_lines))


if __name__ == "__main__":
    main()
