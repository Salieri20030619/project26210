#!/usr/bin/env python3
"""Validate malignant evidence tiers and minimum-cell thresholds at patient level."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
TIER_DIR = HERE / "posthoc_malignancy_evidence_tier_programs_v3"
EXPANDED_RUN = HERE / "runs" / "20260812_expanded_malignancy_v2"
OUT = HERE / "posthoc_tier_threshold_validation_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
METHODS = [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]
THRESHOLDS = [10, 20, 30, 50]
RELATIONS = [
    ("H1_F1_to_SPP1", "f1_top100_weighted", "spp1_state_no_spp1", "myeloid", 1),
    ("H1_TYMP_to_SPP1", "tymp_expression", "spp1_state_no_spp1", "myeloid", 1),
    ("H2_F1_to_APOE_C1QC", "f1_top100_weighted", "apoe_tam_expression_proxy", "myeloid", 1),
    ("H2_F1_to_C1QC", "f1_top100_weighted", "c1qc_state", "myeloid", 1),
    ("H2_F1_to_HLAII", "f1_top100_weighted", "hla_ii_supporting", "myeloid", 0),
    ("H3_F1_to_hypoxia", "f1_top100_weighted", "hypoxia_angiogenesis", "malignant", 1),
    ("H3_TYMP_to_hypoxia", "tymp_expression", "hypoxia_angiogenesis", "malignant", 1),
    ("H5_F1_to_MTA", "f1_top100_weighted", "polyamine_mta_context", "malignant", 1),
    ("H5_TYMP_to_MTA", "tymp_expression", "polyamine_mta_context", "malignant", 1),
]


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(values)
    if not valid.any():
        return result
    p = values.loc[valid].to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    result.loc[valid] = adjusted
    return result


def cross_sectional_nonmalignant(dataset: str) -> pd.DataFrame:
    sample = pd.read_csv(EXPANDED_RUN / dataset / "tables" / "sample_program_summary_v1.tsv", sep="\t")
    sample = sample[sample["analysis_compartment"] != "malignant_provisional"].copy()
    if dataset == "GSE178481":
        sample = sample[~sample["sample_id"].astype(str).str.contains("Normal", case=False, regex=False)]
    if dataset == "GSE264586":
        sample = sample[sample["timepoint"] == "Baseline"]
    rows = []
    keys = ["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"]
    for group_keys, group in sample.groupby(keys, sort=True, dropna=False):
        weights = group["n_cells"].to_numpy(float)
        rows.append(dict(
            zip(keys, group_keys), n_cells=int(weights.sum()),
            mean_primary=float(np.average(group["mean_primary"], weights=weights)),
            mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=weights)),
        ))
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir()
    malignant = pd.read_csv(TIER_DIR / "cross_sectional_patient_evidence_tier_programs_v1.tsv", sep="\t")
    audit = pd.read_csv(TIER_DIR / "cell_evidence_method_and_annotation_audit_v1.tsv", sep="\t")

    composition = (
        audit.groupby(["dataset", "evidence_tier", "evidence_detail", "phase7_annotation_confidence"], as_index=False)["n_cells"].sum()
    )
    tier_totals = composition.groupby(["dataset", "evidence_tier"])["n_cells"].transform("sum")
    composition["within_dataset_tier_fraction"] = composition["n_cells"] / tier_totals
    atomic_tsv(composition, OUT / "evidence_method_annotation_composition_v1.tsv")

    rows = []
    for dataset_index, dataset in enumerate(DATASETS):
        nonmal = cross_sectional_nonmalignant(dataset)
        mal_dataset = malignant[malignant["dataset"] == dataset]
        for threshold in THRESHOLDS:
            nwide = nonmal[nonmal["n_cells"] >= threshold].copy()
            nwide["key"] = nwide["signature_id"] + "__" + nwide["analysis_compartment"]
            for method_index, (method, value) in enumerate(METHODS):
                nw = nwide.pivot(index=["dataset", "patient_key", "patient_id"], columns="key", values=value).reset_index()
                for tier in ["high_confidence", "supportive", "expanded"]:
                    mt = mal_dataset[(mal_dataset["evidence_tier"] == tier) & (mal_dataset["n_cells"] >= threshold)]
                    mw = mt.pivot(index=["dataset", "patient_key", "patient_id"], columns="signature_id", values=value).reset_index()
                    wide = mw.merge(nw, on=["dataset", "patient_key", "patient_id"], how="outer")
                    for relation_index, (relation, x_sig, y_sig, y_comp, expected) in enumerate(RELATIONS):
                        y_col = y_sig if y_comp == "malignant" else f"{y_sig}__{y_comp}"
                        columns = [c for c in [x_sig, y_col] if c in wide.columns]
                        sub = wide[["dataset", "patient_key", "patient_id", *columns]].dropna() if len(columns) == 2 else pd.DataFrame()
                        n = len(sub)
                        if n >= 5 and sub[x_sig].nunique() > 1 and sub[y_col].nunique() > 1:
                            x = sub[x_sig].to_numpy(float)
                            y = sub[y_col].to_numpy(float)
                            test = stats.spearmanr(x, y)
                            rho, p = float(test.statistic), float(test.pvalue)
                            loo = []
                            for i in range(n):
                                xi, yi = np.delete(x, i), np.delete(y, i)
                                if np.unique(xi).size > 1 and np.unique(yi).size > 1:
                                    loo.append(float(stats.spearmanr(xi, yi).statistic))
                            loo_same = float(np.mean(np.sign(loo) == np.sign(rho))) if loo and rho != 0 else np.nan
                        else:
                            rho = p = loo_same = np.nan
                        rows.append({
                            "dataset": dataset, "method": method, "evidence_tier": tier,
                            "minimum_cells": threshold, "relation": relation,
                            "n_patients": n, "spearman_rho": rho, "p_value": p,
                            "loo_sign_stability_fraction": loo_same,
                            "expected_direction": expected,
                        })
    results = pd.DataFrame(rows)
    results["p_fdr_bh_uniform_9"] = results.groupby(
        ["dataset", "method", "evidence_tier", "minimum_cells"], group_keys=False
    )["p_value"].apply(bh)
    atomic_tsv(results, OUT / "tier_threshold_malignant_relations_v1.tsv")

    key_names = {"H2_F1_to_C1QC", "H2_F1_to_APOE_C1QC", "H3_F1_to_hypoxia", "H1_TYMP_to_SPP1"}
    key = results[results["relation"].isin(key_names)].sort_values(
        ["relation", "dataset", "method", "minimum_cells", "evidence_tier"]
    )
    atomic_tsv(key, OUT / "key_relations_tier_threshold_v1.tsv")

    stability_rows = []
    for keys, group in results[results["spearman_rho"].notna()].groupby(["dataset", "method", "evidence_tier", "relation"]):
        signs = np.sign(group["spearman_rho"].to_numpy(float))
        stability_rows.append({
            "dataset": keys[0], "method": keys[1], "evidence_tier": keys[2], "relation": keys[3],
            "thresholds_tested": len(group), "minimum_n_patients": int(group["n_patients"].min()),
            "maximum_n_patients": int(group["n_patients"].max()),
            "median_rho": float(group["spearman_rho"].median()),
            "minimum_rho": float(group["spearman_rho"].min()),
            "maximum_rho": float(group["spearman_rho"].max()),
            "threshold_sign_stability": float(np.mean(signs == signs[0])),
        })
    stability = pd.DataFrame(stability_rows)
    atomic_tsv(stability, OUT / "threshold_stability_summary_v1.tsv")

    support = audit[audit["evidence_tier"] == "supportive"]
    supportive_copykat = int(support.loc[support["evidence_detail"].str.contains("copykat_only"), "n_cells"].sum())
    supportive_infercnv = int(support.loc[support["evidence_detail"].str.contains("infercnv_only"), "n_cells"].sum())
    expanded20 = results[(results["evidence_tier"] == "expanded") & (results["minimum_cells"] == 20)]
    status = {
        "schema_version": "phase14_tier_threshold_validation_v1",
        "status": "PASS_POSTHOC_TIER_THRESHOLD_VALIDATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "supportive_cells": int(support["n_cells"].sum()),
        "supportive_copykat_only_cells": supportive_copykat,
        "supportive_infercnv_only_cells": supportive_infercnv,
        "supportive_copykat_only_fraction": supportive_copykat / int(support["n_cells"].sum()),
        "thresholds": THRESHOLDS,
        "relation_rows": int(len(results)),
        "expanded_threshold20_tested": int(expanded20["spearman_rho"].notna().sum()),
        "expanded_threshold20_nominal": int((expanded20["p_value"] < 0.05).sum()),
        "expanded_threshold20_uniform9_fdr": int((expanded20["p_fdr_bh_uniform_9"] < 0.05).sum()),
        "formal_strict_v6_replaced": False,
    }
    atomic_json(status, OUT / "TIER_THRESHOLD_VALIDATION_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
