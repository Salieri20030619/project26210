#!/usr/bin/env python3
"""Post-hoc diagnostic decomposition of the Phase 14 v6 zero-FDR result.

This analysis is descriptive/exploratory and never replaces the frozen within-dataset
confirmatory analysis. It quantifies testability, CNV-gate attrition, v3-to-v6 changes,
and cross-study correlation consistency.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
V3 = HERE / "patient_analysis_v3"
V6 = HERE / "patient_analysis_final_cnv_v6"
RUN3 = HERE / "runs" / "20260812_082000_provisional_v3"
RUN6 = HERE / "runs" / "20260812_final_cnv_gated_v6"
OUT = HERE / "posthoc_zero_fdr_diagnostic_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def meta_correlations(group: pd.DataFrame) -> dict:
    use = group.loc[group["spearman_rho"].notna() & (group["n_patients"] >= 5)].copy()
    k = len(use)
    if k < 2:
        return {"studies": k}
    rho = np.clip(use["spearman_rho"].to_numpy(float), -0.999999, 0.999999)
    n = use["n_patients"].to_numpy(float)
    yi = np.arctanh(rho)
    vi = 1.0 / (n - 3.0)
    wi = 1.0 / vi
    fixed_y = float(np.sum(wi * yi) / np.sum(wi))
    fixed_se = float(math.sqrt(1.0 / np.sum(wi)))
    fixed_z = fixed_y / fixed_se
    fixed_p = float(2.0 * stats.norm.sf(abs(fixed_z)))
    q = float(np.sum(wi * (yi - fixed_y) ** 2))
    q_df = k - 1
    q_p = float(stats.chi2.sf(q, q_df))
    c = float(np.sum(wi) - np.sum(wi**2) / np.sum(wi))
    tau2 = float(max(0.0, (q - q_df) / c)) if c > 0 else 0.0
    wr = 1.0 / (vi + tau2)
    random_y = float(np.sum(wr * yi) / np.sum(wr))
    random_se = float(math.sqrt(1.0 / np.sum(wr)))
    random_z = random_y / random_se
    random_p = float(2.0 * stats.norm.sf(abs(random_z)))
    i2 = float(max(0.0, (q - q_df) / q) * 100.0) if q > 0 else 0.0
    return {
        "studies": k,
        "patients_total": int(n.sum()),
        "rho_min": float(rho.min()),
        "rho_max": float(rho.max()),
        "positive_studies": int((rho > 0).sum()),
        "negative_studies": int((rho < 0).sum()),
        "fixed_rho": float(np.tanh(fixed_y)),
        "fixed_ci_low": float(np.tanh(fixed_y - 1.96 * fixed_se)),
        "fixed_ci_high": float(np.tanh(fixed_y + 1.96 * fixed_se)),
        "fixed_p": fixed_p,
        "random_rho": float(np.tanh(random_y)),
        "random_ci_low": float(np.tanh(random_y - 1.96 * random_se)),
        "random_ci_high": float(np.tanh(random_y + 1.96 * random_se)),
        "random_p": random_p,
        "q": q,
        "q_p": q_p,
        "i2_percent": i2,
        "tau2": tau2,
    }


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)

    assoc3_path = V3 / "within_dataset_patient_associations_v3.tsv"
    assoc6_path = V6 / "within_dataset_patient_associations_v6.tsv"
    assoc3 = pd.read_csv(assoc3_path, sep="\t")
    assoc6 = pd.read_csv(assoc6_path, sep="\t")

    tested6 = assoc6[assoc6["status"].str.startswith("TESTED")].copy()
    testability_rows = []
    for (dataset, method), group in assoc6.groupby(["dataset", "method"], sort=True):
        tested = group[group["status"].str.startswith("TESTED")]
        testability_rows.append(
            {
                "dataset": dataset,
                "method": method,
                "hypotheses_total": len(group),
                "hypotheses_tested": len(tested),
                "hypotheses_not_tested": len(group) - len(tested),
                "bh_family_valid_p": int(tested["p_value"].notna().sum()),
                "min_n_patients_tested": tested["n_patients"].min() if len(tested) else np.nan,
                "max_n_patients_tested": tested["n_patients"].max() if len(tested) else np.nan,
                "min_nominal_p": tested["p_value"].min() if len(tested) else np.nan,
                "min_bh_q": tested["p_fdr_bh_within_dataset_method"].min() if len(tested) else np.nan,
                "nominal_p_lt_0_05": int((tested["p_value"] < 0.05).sum()),
                "bh_q_lt_0_05": int((tested["p_fdr_bh_within_dataset_method"] < 0.05).sum()),
            }
        )
    testability = pd.DataFrame(testability_rows)
    atomic_tsv(testability, OUT / "testability_by_dataset_method_v1.tsv")

    top = tested6.sort_values(["p_value", "dataset", "method", "relation"]).copy()
    top["bootstrap_ci_excludes_zero"] = (top["bootstrap_ci_low"] > 0) | (top["bootstrap_ci_high"] < 0)
    top["nominal_p_lt_0_05"] = top["p_value"] < 0.05
    atomic_tsv(top, OUT / "all_tested_associations_ranked_v1.tsv")

    keys = ["dataset", "method", "relation"]
    compare = assoc3.merge(assoc6, on=keys, how="outer", suffixes=("_v3", "_v6"), validate="one_to_one")
    compare["delta_n_v6_minus_v3"] = compare["n_patients_v6"] - compare["n_patients_v3"]
    compare["delta_rho_v6_minus_v3"] = compare["spearman_rho_v6"] - compare["spearman_rho_v3"]
    compare["lost_testability_after_cnv"] = compare["status_v3"].str.startswith("TESTED") & ~compare["status_v6"].str.startswith("TESTED")
    atomic_tsv(compare, OUT / "v3_vs_v6_association_diagnostics_v1.tsv")

    retention_rows = []
    for dataset in DATASETS:
        p3 = pd.read_csv(RUN3 / dataset / "tables" / "patient_program_summary_v1.tsv", sep="\t")
        p6 = pd.read_csv(RUN6 / dataset / "tables" / "patient_program_summary_v1.tsv", sep="\t")
        f3 = p3[(p3["analysis_compartment"] == "malignant_provisional") & (p3["signature_id"] == "f1_top100_weighted")]
        f6 = p6[(p6["analysis_compartment"] == "malignant_provisional") & (p6["signature_id"] == "f1_top100_weighted")]
        left = f3[["patient_key", "patient_id", "n_cells"]].rename(columns={"n_cells": "provisional_malignant_cells"})
        right = f6[["patient_key", "patient_id", "n_cells"]].rename(columns={"n_cells": "final_cnv_malignant_cells"})
        merged = left.merge(right, on=["patient_key", "patient_id"], how="outer")
        merged[["provisional_malignant_cells", "final_cnv_malignant_cells"]] = merged[["provisional_malignant_cells", "final_cnv_malignant_cells"]].fillna(0).astype(int)
        merged.insert(0, "dataset", dataset)
        merged["retention_fraction"] = np.where(
            merged["provisional_malignant_cells"] > 0,
            merged["final_cnv_malignant_cells"] / merged["provisional_malignant_cells"],
            np.nan,
        )
        merged["provisional_ge20"] = merged["provisional_malignant_cells"] >= 20
        merged["final_ge20"] = merged["final_cnv_malignant_cells"] >= 20
        retention_rows.append(merged)
    retention = pd.concat(retention_rows, ignore_index=True).sort_values(["dataset", "patient_key"])
    atomic_tsv(retention, OUT / "malignant_patient_cell_retention_v1.tsv")

    retention_summary = (
        retention.groupby("dataset", as_index=False)
        .agg(
            provisional_patients_ge20=("provisional_ge20", "sum"),
            final_patients_ge20=("final_ge20", "sum"),
            provisional_cells=("provisional_malignant_cells", "sum"),
            final_cells=("final_cnv_malignant_cells", "sum"),
            median_patient_retention=("retention_fraction", "median"),
        )
    )
    retention_summary["pooled_cell_retention"] = retention_summary["final_cells"] / retention_summary["provisional_cells"]
    atomic_tsv(retention_summary, OUT / "malignant_retention_summary_v1.tsv")

    meta_rows = []
    for (method, relation), group in assoc6.groupby(["method", "relation"], sort=True):
        result = meta_correlations(group)
        result.update({"method": method, "relation": relation})
        meta_rows.append(result)
    meta = pd.DataFrame(meta_rows)
    meta["fixed_p_fdr_bh_within_method"] = meta.groupby("method", group_keys=False)["fixed_p"].apply(bh)
    meta["random_p_fdr_bh_within_method"] = meta.groupby("method", group_keys=False)["random_p"].apply(bh)
    meta["analysis_status"] = "POSTHOC_EXPLORATORY_NOT_CONFIRMATORY"
    atomic_tsv(meta, OUT / "exploratory_cross_study_correlation_meta_v1.tsv")

    manifest_rows = []
    for path in [assoc3_path, assoc6_path, V6 / "PATIENT_ANALYSIS_STATUS_v6.json", V6 / "INDEPENDENT_VALIDATION_STATUS_v6.tsv"]:
        manifest_rows.append({"input": str(path.relative_to(HERE.parent)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), OUT / "INPUT_MANIFEST_v1.tsv")

    nominal = tested6[tested6["p_value"] < 0.05]
    lost = compare[compare["lost_testability_after_cnv"]]
    meta_random_sig = meta[meta["random_p_fdr_bh_within_method"] < 0.05]
    status = {
        "schema_version": "phase14_zero_fdr_diagnostic_v1",
        "status": "PASS_POSTHOC_DIAGNOSTIC",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "confirmatory_result_unchanged": "zero v6 within-dataset FDR-significant associations",
        "tested_associations": int(len(tested6)),
        "nominal_p_lt_0_05": int(len(nominal)),
        "bootstrap_ci_excludes_zero": int(((tested6["bootstrap_ci_low"] > 0) | (tested6["bootstrap_ci_high"] < 0)).sum()),
        "lost_testability_after_final_cnv_gate": int(len(lost)),
        "final_malignant_cells": int(retention["final_cnv_malignant_cells"].sum()),
        "provisional_malignant_cells": int(retention["provisional_malignant_cells"].sum()),
        "exploratory_random_effect_meta_fdr_lt_0_05": int(len(meta_random_sig)),
        "warning": "Cross-study meta-analysis was added after inspecting the zero-FDR result and is exploratory, not confirmatory.",
    }
    atomic_json(status, OUT / "DIAGNOSTIC_STATUS_v1.json")


if __name__ == "__main__":
    main()

