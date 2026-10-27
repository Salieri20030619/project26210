#!/usr/bin/env python3
"""Post-hoc malignant-relation sensitivity across phenotype, any-method, and strict CNV tiers."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
ANY = HERE / "posthoc_any_method_positive_inputs_v2"
V3 = HERE / "patient_analysis_v3" / "within_dataset_patient_associations_v3.tsv"
V6 = HERE / "patient_analysis_final_cnv_v6" / "within_dataset_patient_associations_v6.tsv"
OUT = HERE / "posthoc_malignancy_gate_tiers_v3"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
METHODS = ["primary_gene_z_mean", "sensitivity_mean_log_normalized"]
MIN_CELLS = 20
MIN_PATIENTS = 5
N_BOOT = 9999
SEED = 20260812
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


def bootstrap_ci(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(N_BOOT, len(x)))
    xr = stats.rankdata(x[idx], axis=1)
    yr = stats.rankdata(y[idx], axis=1)
    xc = xr - xr.mean(axis=1, keepdims=True)
    yc = yr - yr.mean(axis=1, keepdims=True)
    den = np.sqrt(np.sum(xc**2, axis=1) * np.sum(yc**2, axis=1))
    values = np.divide(np.sum(xc * yc, axis=1), den, out=np.full(N_BOOT, np.nan), where=den > 0)
    values = values[np.isfinite(values)]
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)), len(values) / N_BOOT


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    phenotype = pd.read_csv(V3, sep="\t")
    strict = pd.read_csv(V6, sep="\t")
    malignant_names = {x[0] for x in RELATIONS}
    reused = []
    for tier, frame in [("phenotype_candidate", phenotype), ("strict_conservative", strict)]:
        sub = frame[frame["relation"].isin(malignant_names)].copy()
        sub.insert(0, "gate_tier", tier)
        reused.append(sub)

    any_rows = []
    for dataset_index, dataset in enumerate(DATASETS):
        malignant = pd.read_csv(ANY / dataset / "malignant_patient_program_summary_any_method_v1.tsv", sep="\t")
        nonmal = pd.read_csv(HERE / "runs" / "20260812_final_cnv_gated_v6" / dataset / "tables" / "patient_program_summary_v1.tsv", sep="\t")
        nonmal = nonmal[nonmal["analysis_compartment"] != "malignant_provisional"]
        for method_index, method in enumerate(METHODS):
            value = "mean_primary" if method == "primary_gene_z_mean" else "mean_sensitivity"
            mw = malignant[malignant["n_cells"] >= MIN_CELLS].pivot(index=["dataset", "patient_key", "patient_id"], columns="signature_id", values=value).reset_index()
            nw = nonmal[nonmal["n_cells"] >= MIN_CELLS].pivot(index=["dataset", "patient_key", "patient_id"], columns=["signature_id", "analysis_compartment"], values=value)
            nw.columns = [f"{a}__{b}" for a, b in nw.columns]
            wide = mw.merge(nw.reset_index(), on=["dataset", "patient_key", "patient_id"], how="outer")
            for relation_index, (relation, x_sig, y_sig, y_comp, expected) in enumerate(RELATIONS):
                y_col = y_sig if y_comp == "malignant" else f"{y_sig}__{y_comp}"
                columns = [c for c in [x_sig, y_col] if c in wide.columns]
                sub = wide[["dataset", "patient_key", "patient_id", *columns]].dropna()
                n = len(sub) if len(columns) == 2 else 0
                if n >= MIN_PATIENTS and sub[x_sig].nunique() > 1 and sub[y_col].nunique() > 1:
                    x = sub[x_sig].to_numpy(float)
                    y = sub[y_col].to_numpy(float)
                    test = stats.spearmanr(x, y)
                    rho = float(test.statistic)
                    p = float(test.pvalue)
                    low, high, valid = bootstrap_ci(x, y, SEED + dataset_index * 1000 + method_index * 100 + relation_index)
                    loo = [float(stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic) for i in range(n)]
                    status = "TESTED_POSTHOC_ANY_METHOD_POSITIVE"
                else:
                    rho = p = low = high = valid = np.nan
                    loo = []
                    status = "NOT_TESTED_LT5_ELIGIBLE_PATIENTS"
                any_rows.append({
                    "gate_tier": "any_method_positive",
                    "dataset": dataset,
                    "method": method,
                    "relation": relation,
                    "n_patients": n,
                    "spearman_rho": rho,
                    "bootstrap_ci_low": low,
                    "bootstrap_ci_high": high,
                    "bootstrap_valid_fraction": valid,
                    "p_value": p,
                    "direction_matches_expected": bool(np.sign(rho) == expected) if expected and np.isfinite(rho) else np.nan,
                    "loo_min_rho": min(loo) if loo else np.nan,
                    "loo_max_rho": max(loo) if loo else np.nan,
                    "loo_sign_stability_fraction": float(np.mean(np.sign(loo) == np.sign(rho))) if loo and rho != 0 else np.nan,
                    "status": status,
                })
    any_frame = pd.DataFrame(any_rows)
    any_frame["p_fdr_bh_within_dataset_method"] = any_frame.groupby(["dataset", "method"], group_keys=False)["p_value"].apply(bh)

    common = ["gate_tier", "dataset", "method", "relation", "n_patients", "spearman_rho", "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_valid_fraction", "p_value", "direction_matches_expected", "loo_min_rho", "loo_max_rho", "loo_sign_stability_fraction", "status", "p_fdr_bh_within_dataset_method"]
    combined = pd.concat([reused[0][common], any_frame[common], reused[1][common]], ignore_index=True)
    combined = combined.sort_values(["dataset", "method", "relation", "gate_tier"])
    combined["p_fdr_bh_within_dataset_method_uniform_9_malignant"] = combined.groupby(
        ["gate_tier", "dataset", "method"], group_keys=False
    )["p_value"].apply(bh)
    atomic_tsv(combined, OUT / "malignant_relations_by_gate_tier_v1.tsv")

    summary = (
        combined.groupby("gate_tier", as_index=False)
        .agg(
            relation_rows=("relation", "size"),
            rows_tested=("spearman_rho", lambda x: int(x.notna().sum())),
            nominal_p_lt_0_05=("p_value", lambda x: int((x < 0.05).sum())),
            fdr_q_lt_0_05=("p_fdr_bh_within_dataset_method", lambda x: int((x < 0.05).sum())),
            median_tested_n=("n_patients", lambda x: float(x[combined.loc[x.index, "spearman_rho"].notna()].median()) if combined.loc[x.index, "spearman_rho"].notna().any() else np.nan),
        )
    )
    atomic_tsv(summary, OUT / "gate_tier_summary_v1.tsv")
    status = {
        "schema_version": "phase14_malignancy_gate_tiers_posthoc_v3",
        "status": "PASS_POSTHOC_SENSITIVITY",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "warning": "Any-method-positive tier was defined after inspecting final results; it is exploratory and cannot replace strict conservative v6.",
        "tiers": ["phenotype_candidate", "any_method_positive", "strict_conservative"],
    }
    atomic_json(status, OUT / "GATE_TIER_SENSITIVITY_STATUS_v1.json")


if __name__ == "__main__":
    main()
