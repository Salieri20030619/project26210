#!/usr/bin/env python3
"""Patient-level, within-dataset analysis of frozen Phase 13 programs."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUN = Path(os.environ.get("PROJECT26210_ANALYSIS_RUN", str(HERE / "runs" / "20260812_082000_provisional_v3"))).resolve()
OUT = Path(os.environ.get("PROJECT26210_ANALYSIS_OUT", str(HERE / "patient_analysis_v3"))).resolve()
OUTPUT_VERSION = os.environ.get("PROJECT26210_ANALYSIS_VERSION", "v3")
FINAL_MALIGNANCY = os.environ.get("PROJECT26210_FINAL_MALIGNANCY", "0") == "1"
PASS_STATUS = os.environ.get("PROJECT26210_PASS_STATUS", "PASS" if FINAL_MALIGNANCY else "PASS_PROVISIONAL_MALIGNANCY")
MALIGNANT_TEST_STATUS = os.environ.get("PROJECT26210_MALIGNANT_TEST_STATUS", "TESTED" if FINAL_MALIGNANCY else "TESTED_PROVISIONAL_MALIGNANCY")
MALIGNANCY_GATE_LABEL = os.environ.get(
    "PROJECT26210_MALIGNANCY_GATE_LABEL",
    "strict_conservative_rebuild_priority" if FINAL_MALIGNANCY else "phenotype_candidate_provisional",
)
MALIGNANCY_LIMITATION = os.environ.get(
    "PROJECT26210_MALIGNANCY_LIMITATION",
    "Malignant compartments require final Part1-rebuild/external-CNV conservative support."
    if FINAL_MALIGNANCY
    else "Malignant compartments are Phase8 tumor-like candidates until external CNV consensus is finalized.",
)
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
MIN_CELLS = 20
MIN_PATIENTS = 5
SEED = 20260812
N_BOOT = 9999


RELATIONS = [
    ("H1_F1_to_SPP1", "f1_top100_weighted", "malignant_provisional", "spp1_state_no_spp1", "myeloid", 1),
    ("H1_TYMP_to_SPP1", "tymp_expression", "malignant_provisional", "spp1_state_no_spp1", "myeloid", 1),
    ("H2_F1_to_APOE_C1QC", "f1_top100_weighted", "malignant_provisional", "apoe_tam_expression_proxy", "myeloid", 1),
    ("H2_F1_to_C1QC", "f1_top100_weighted", "malignant_provisional", "c1qc_state", "myeloid", 1),
    ("H2_F1_to_HLAII", "f1_top100_weighted", "malignant_provisional", "hla_ii_supporting", "myeloid", 0),
    ("H3_F1_to_hypoxia", "f1_top100_weighted", "malignant_provisional", "hypoxia_angiogenesis", "malignant_provisional", 1),
    ("H3_TYMP_to_hypoxia", "tymp_expression", "malignant_provisional", "hypoxia_angiogenesis", "malignant_provisional", 1),
    ("H3_SPP1_to_myeloid_vascular", "spp1_state_no_spp1", "myeloid", "vascular_remodeling", "myeloid", 1),
    ("H4_SPP1_to_Tdysfunction", "spp1_state_no_spp1", "myeloid", "tcell_dysfunction", "t_nk", 1),
    ("H4_SPP1_to_cytotoxicity", "spp1_state_no_spp1", "myeloid", "cytotoxicity", "t_nk", -1),
    ("H5_F1_to_MTA", "f1_top100_weighted", "malignant_provisional", "polyamine_mta_context", "malignant_provisional", 1),
    ("H5_TYMP_to_MTA", "tymp_expression", "malignant_provisional", "polyamine_mta_context", "malignant_provisional", 1),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(f"Refusing existing output: {path}")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def atomic_text(text: str, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(f"Refusing existing output: {path}")
    part.write_text(text, encoding="utf-8")
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


def spearman_bootstrap(x: np.ndarray, y: np.ndarray, seed_offset: int) -> tuple[float, float, float, float]:
    observed = stats.spearmanr(x, y).statistic
    rng = np.random.default_rng(SEED + seed_offset)
    indices = rng.integers(0, len(x), size=(N_BOOT, len(x)))
    # Exact Spearman correlations for all bootstrap resamples, vectorized by
    # ranking each bootstrap row and applying row-wise Pearson correlation.
    x_rank = stats.rankdata(x[indices], axis=1)
    y_rank = stats.rankdata(y[indices], axis=1)
    x_centered = x_rank - x_rank.mean(axis=1, keepdims=True)
    y_centered = y_rank - y_rank.mean(axis=1, keepdims=True)
    numerator = np.sum(x_centered * y_centered, axis=1)
    denominator = np.sqrt(np.sum(x_centered ** 2, axis=1) * np.sum(y_centered ** 2, axis=1))
    values = np.divide(numerator, denominator, out=np.full(N_BOOT, np.nan), where=denominator > 0)
    valid = values[np.isfinite(values)]
    if not len(valid):
        return float(observed), np.nan, np.nan, 0.0
    return float(observed), float(np.quantile(valid, 0.025)), float(np.quantile(valid, 0.975)), float(len(valid) / N_BOOT)


def patient_wide(program: pd.DataFrame, value: str) -> pd.DataFrame:
    eligible = program.loc[program["n_cells"] >= MIN_CELLS].copy()
    eligible["key"] = eligible["signature_id"] + "__" + eligible["analysis_compartment"]
    return eligible.pivot(index=["dataset", "patient_key", "patient_id"], columns="key", values=value).reset_index()


def build_cross_sectional_patient_programs(samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build patient units from the prespecified biologically comparable samples.

    GSE178481 ecology analyses use explicit non-Normal (tumor) libraries. GSE264586
    cross-sectional analyses use Baseline only. Other datasets use all verified
    samples because they have one timepoint or a common surgical timepoint.
    """
    selected = samples.copy()
    include = pd.Series(True, index=selected.index)
    include &= ~((selected["dataset"] == "GSE178481") & selected["sample_id"].astype(str).str.contains("Normal", case=False, regex=False))
    include &= ~((selected["dataset"] == "GSE264586") & (selected["timepoint"] != "Baseline"))
    selected = selected.loc[include].copy()

    group_fields = ["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"]
    rows = []
    for keys, group in selected.groupby(group_fields, sort=True, dropna=False):
        weights = group["n_cells"].to_numpy(float)
        rows.append(
            dict(
                zip(group_fields, keys),
                unit="patient_cross_sectional_selected_samples",
                n_cells=int(weights.sum()),
                n_samples=int(group[["sample_id", "object_name"]].drop_duplicates().shape[0]),
                mean_primary=float(np.average(group["mean_primary"], weights=weights)),
                median_primary=np.nan,
                mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=weights)),
                median_sensitivity=np.nan,
            )
        )
    patient = pd.DataFrame(rows)
    audit = (
        selected[["dataset", "patient_key", "patient_id", "sample_id", "object_name", "timepoint", "treatment", "response"]]
        .drop_duplicates()
        .sort_values(["dataset", "patient_key", "sample_id"])
        .reset_index(drop=True)
    )
    audit["selection_rule"] = np.select(
        [audit["dataset"] == "GSE178481", audit["dataset"] == "GSE264586"],
        ["explicit tumor libraries only (sample_id does not contain Normal)", "Baseline only"],
        default="all verified samples",
    )
    return patient, audit


def analyze_associations(program: pd.DataFrame, method: str, seed_base: int) -> pd.DataFrame:
    value = "mean_primary" if method == "primary_gene_z_mean" else "mean_sensitivity"
    wide = patient_wide(program, value)
    rows = []
    for relation_index, (relation, x_sig, x_comp, y_sig, y_comp, expected) in enumerate(RELATIONS):
        x_col = f"{x_sig}__{x_comp}"
        y_col = f"{y_sig}__{y_comp}"
        subset = wide[["dataset", "patient_key", "patient_id", *[c for c in [x_col, y_col] if c in wide]]].dropna()
        n = len(subset) if x_col in subset and y_col in subset else 0
        if n >= MIN_PATIENTS and subset[x_col].nunique() > 1 and subset[y_col].nunique() > 1:
            x = subset[x_col].to_numpy(float)
            y = subset[y_col].to_numpy(float)
            test = stats.spearmanr(x, y)
            rho, ci_low, ci_high, boot_valid = spearman_bootstrap(x, y, seed_base + relation_index)
            p_value = float(test.pvalue)
            direction_match = bool(np.sign(rho) == expected) if expected else np.nan
            loo = []
            for i in range(n):
                kept = np.arange(n) != i
                if np.unique(x[kept]).size > 1 and np.unique(y[kept]).size > 1:
                    loo.append(float(stats.spearmanr(x[kept], y[kept]).statistic))
            loo_arr = np.array(loo, dtype=float)
            loo_min = float(np.min(loo_arr)) if len(loo_arr) else np.nan
            loo_max = float(np.max(loo_arr)) if len(loo_arr) else np.nan
            loo_same = float(np.mean(np.sign(loo_arr) == np.sign(rho))) if len(loo_arr) and rho != 0 else np.nan
            status = MALIGNANT_TEST_STATUS if "malignant" in x_comp or "malignant" in y_comp else "TESTED"
        else:
            rho = ci_low = ci_high = boot_valid = p_value = np.nan
            direction_match = np.nan
            loo_min = loo_max = loo_same = np.nan
            status = "NOT_TESTED_LT5_ELIGIBLE_PATIENTS"
        rows.append(
            {
                "dataset": program["dataset"].iloc[0],
                "method": method,
                "relation": relation,
                "x_signature": x_sig,
                "x_compartment": x_comp,
                "y_signature": y_sig,
                "y_compartment": y_comp,
                "expected_direction": expected,
                "n_patients": n,
                "spearman_rho": rho,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "bootstrap_valid_fraction": boot_valid,
                "p_value": p_value,
                "direction_matches_expected": direction_match,
                "loo_min_rho": loo_min,
                "loo_max_rho": loo_max,
                "loo_sign_stability_fraction": loo_same,
                "status": status,
            }
        )
    answer = pd.DataFrame(rows)
    answer["p_fdr_bh_within_dataset_method"] = bh(answer["p_value"])
    return answer


def subset_analysis(program: pd.DataFrame, dataset: str, subset_kind: str, members: set[str]) -> pd.DataFrame:
    frame = program.loc[program["patient_id"].astype(str).isin(members)].copy()
    if frame.empty:
        return pd.DataFrame()
    rows = []
    for method_index, method in enumerate(["primary_gene_z_mean", "sensitivity_mean_log_normalized"]):
        result = analyze_associations(frame, method, 5000 + method_index * 100)
        result.insert(1, "subset", subset_kind)
        result.insert(2, "subset_members", ";".join(sorted(members)))
        rows.append(result)
    return pd.concat(rows, ignore_index=True)


def paired_longitudinal(sample: pd.DataFrame, dataset: str) -> pd.DataFrame:
    frame = sample.loc[(sample["dataset"] == dataset) & (sample["n_cells"] >= MIN_CELLS)].copy()
    rows = []
    for (signature, compartment), group in frame.groupby(["signature_id", "analysis_compartment"], sort=True):
        baseline = group.loc[group["timepoint"] == "Baseline"].sort_values("sample_id").drop_duplicates("patient_key")
        for followup in ["Cycle2", "End_of_treatment"]:
            later = group.loc[group["timepoint"] == followup].sort_values("sample_id").drop_duplicates("patient_key")
            paired = baseline.merge(later, on=["dataset", "patient_key", "patient_id"], suffixes=("_baseline", "_followup"))
            for method, col in [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]:
                delta = paired[f"{col}_followup"].to_numpy(float) - paired[f"{col}_baseline"].to_numpy(float)
                n = len(delta)
                if n:
                    positive = int(np.sum(delta > 0))
                    negative = int(np.sum(delta < 0))
                    nonzero = positive + negative
                    sign_p = float(stats.binomtest(positive, nonzero, 0.5).pvalue) if nonzero else np.nan
                    wilcoxon_p = float(stats.wilcoxon(delta, alternative="two-sided", zero_method="wilcox").pvalue) if nonzero >= 5 else np.nan
                else:
                    positive = negative = nonzero = 0
                    sign_p = wilcoxon_p = np.nan
                rows.append(
                    {
                        "dataset": dataset,
                        "signature_id": signature,
                        "analysis_compartment": compartment,
                        "method": method,
                        "followup": followup,
                        "n_pairs": n,
                        "patient_keys": ";".join(paired["patient_key"].astype(str)),
                        "median_change": float(np.median(delta)) if n else np.nan,
                        "mean_change": float(np.mean(delta)) if n else np.nan,
                        "positive_pairs": positive,
                        "negative_pairs": negative,
                        "exact_sign_p_value": sign_p,
                        "wilcoxon_p_value": wilcoxon_p,
                        "status": "DESCRIPTIVE_LT5_PAIRS" if n < 5 else (MALIGNANT_TEST_STATUS if compartment == "malignant_provisional" else "TESTED"),
                    }
                )
    answer = pd.DataFrame(rows)
    answer["sign_p_fdr_bh"] = bh(answer["exact_sign_p_value"])
    answer["wilcoxon_p_fdr_bh"] = bh(answer["wilcoxon_p_value"])
    return answer


def response_analysis(sample: pd.DataFrame) -> pd.DataFrame:
    frame = sample.loc[(sample["dataset"] == "GSE294109") & (sample["n_cells"] >= MIN_CELLS)].copy()
    frame = frame.loc[frame["response"].isin(["CR", "PR", "No_effect"])]
    if frame.empty:
        return pd.DataFrame()
    patient_rows = []
    response_fields = ["dataset", "patient_key", "patient_id", "response", "signature_id", "analysis_compartment"]
    for keys, group in frame.groupby(response_fields, sort=True):
        weights = group["n_cells"].to_numpy(float)
        patient_rows.append(
            dict(
                zip(response_fields, keys),
                n_cells=int(weights.sum()),
                mean_primary=float(np.average(group["mean_primary"], weights=weights)),
                mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=weights)),
            )
        )
    patient = pd.DataFrame(patient_rows)
    patient["response_binary"] = np.where(patient["response"].isin(["CR", "PR"]), "benefit", "no_effect")
    rows = []
    for (signature, compartment), group in patient.groupby(["signature_id", "analysis_compartment"], sort=True):
        for method, col in [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]:
            benefit = group.loc[group["response_binary"] == "benefit", col].to_numpy(float)
            no_effect = group.loc[group["response_binary"] == "no_effect", col].to_numpy(float)
            tested = len(benefit) >= 3 and len(no_effect) >= 3
            mw = stats.mannwhitneyu(benefit, no_effect, alternative="two-sided") if tested else None
            rows.append(
                {
                    "dataset": "GSE294109",
                    "signature_id": signature,
                    "analysis_compartment": compartment,
                    "method": method,
                    "n_benefit_patients": len(benefit),
                    "n_no_effect_patients": len(no_effect),
                    "median_benefit": float(np.median(benefit)) if len(benefit) else np.nan,
                    "median_no_effect": float(np.median(no_effect)) if len(no_effect) else np.nan,
                    "median_difference_benefit_minus_no_effect": float(np.median(benefit) - np.median(no_effect)) if len(benefit) and len(no_effect) else np.nan,
                    "mann_whitney_u": float(mw.statistic) if mw else np.nan,
                    "p_value": float(mw.pvalue) if mw else np.nan,
                    "status": (MALIGNANT_TEST_STATUS if compartment == "malignant_provisional" else "TESTED") if tested else "DESCRIPTIVE_INSUFFICIENT_GROUP_SIZE",
                }
            )
    answer = pd.DataFrame(rows)
    answer["p_fdr_bh"] = bh(answer["p_value"])
    return answer


def malignant_reproducibility(programs_all_samples: pd.DataFrame) -> pd.DataFrame:
    rows = []
    f1 = programs_all_samples.loc[programs_all_samples["signature_id"] == "f1_top100_weighted"].copy()
    for dataset, dataset_frame in f1.groupby("dataset", sort=True):
        for method, value in [("primary_weighted_f1", "mean_primary"), ("sensitivity_unweighted_f1", "mean_sensitivity")]:
            malignant = dataset_frame.loc[
                (dataset_frame["analysis_compartment"] == "malignant_provisional") & (dataset_frame["n_cells"] >= MIN_CELLS),
                ["patient_key", "patient_id", "n_cells", value],
            ].rename(columns={"n_cells": "n_malignant_cells", value: "malignant_score"})
            normal = dataset_frame.loc[
                (dataset_frame["analysis_compartment"] == "normal_epithelial") & (dataset_frame["n_cells"] >= MIN_CELLS),
                ["patient_key", "patient_id", "n_cells", value],
            ].rename(columns={"n_cells": "n_normal_cells", value: "normal_score"})
            paired = malignant.merge(normal, on=["patient_key", "patient_id"], how="inner")
            delta = paired["malignant_score"].to_numpy(float) - paired["normal_score"].to_numpy(float)
            n = len(delta)
            positive = int((delta > 0).sum())
            negative = int((delta < 0).sum())
            nonzero = positive + negative
            sign_p = float(stats.binomtest(positive, nonzero, 0.5).pvalue) if nonzero else np.nan
            wilcoxon_p = float(stats.wilcoxon(delta, alternative="greater", zero_method="wilcox").pvalue) if nonzero >= MIN_PATIENTS else np.nan
            rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "malignant_patients_ge20_cells": int(malignant["patient_key"].nunique()),
                    "normal_patients_ge20_cells": int(normal["patient_key"].nunique()),
                    "paired_patients": n,
                    "paired_patient_keys": ";".join(paired["patient_key"].astype(str)),
                    "median_malignant_minus_normal": float(np.median(delta)) if n else np.nan,
                    "mean_malignant_minus_normal": float(np.mean(delta)) if n else np.nan,
                    "positive_pairs": positive,
                    "negative_pairs": negative,
                    "exact_sign_p_value": sign_p,
                    "one_sided_wilcoxon_p_value": wilcoxon_p,
                    "status": MALIGNANT_TEST_STATUS if n >= MIN_PATIENTS else ("DESCRIPTIVE_LT5_PAIRED_PATIENTS" if FINAL_MALIGNANCY else "DESCRIPTIVE_LT5_PAIRED_PATIENTS_PROVISIONAL_MALIGNANCY"),
                }
            )
    answer = pd.DataFrame(rows)
    answer["sign_p_fdr_bh"] = bh(answer["exact_sign_p_value"])
    answer["wilcoxon_p_fdr_bh"] = bh(answer["one_sided_wilcoxon_p_value"])
    return answer


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing existing output directory: {OUT}")
    OUT.mkdir(parents=False)

    all_programs_unselected = []
    all_samples = []
    input_rows = []
    for dataset_index, dataset in enumerate(DATASETS):
        patient_path = RUN / dataset / "tables" / "patient_program_summary_v1.tsv"
        sample_path = RUN / dataset / "tables" / "sample_program_summary_v1.tsv"
        validation_path = RUN / dataset / "status" / f"{dataset}_INDEPENDENT_VALIDATION_v1.tsv"
        for path in [patient_path, sample_path, validation_path]:
            if not path.is_file():
                raise FileNotFoundError(path)
            input_rows.append({"dataset": dataset, "input": str(path.relative_to(ROOT)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)})
        validation = pd.read_csv(validation_path, sep="\t")
        if not ((validation["check"] == "OVERALL") & (validation["result"] == "PASS")).any():
            raise RuntimeError(f"Dataset validation not PASS: {dataset}")
        program = pd.read_csv(patient_path, sep="\t")
        sample = pd.read_csv(sample_path, sep="\t")
        all_programs_unselected.append(program)
        all_samples.append(sample)

    programs_unselected = pd.concat(all_programs_unselected, ignore_index=True)
    samples = pd.concat(all_samples, ignore_index=True)
    programs, sample_selection_audit = build_cross_sectional_patient_programs(samples)
    all_associations = []
    for dataset_index, dataset in enumerate(DATASETS):
        program = programs.loc[programs["dataset"] == dataset]
        for method_index, method in enumerate(["primary_gene_z_mean", "sensitivity_mean_log_normalized"]):
            all_associations.append(analyze_associations(program, method, dataset_index * 1000 + method_index * 100))
    associations = pd.concat(all_associations, ignore_index=True)

    g178_old_patients = {"PR5", "PR6", "PR7", "PR8", "PR9"}
    g178_new_patients = {"BM1", "BM2", "PR1", "PR2", "PR3", "PR4"}
    subset_rows = [
        subset_analysis(programs.loc[programs["dataset"] == "GSE178481"], "GSE178481", "previously_included_5_patients_11_libraries", g178_old_patients),
        subset_analysis(programs.loc[programs["dataset"] == "GSE178481"], "GSE178481", "newly_added_6_patients_15_libraries", g178_new_patients),
    ]
    subset_results = pd.concat([x for x in subset_rows if not x.empty], ignore_index=True)

    longitudinal = paired_longitudinal(samples, "GSE264586")
    response = response_analysis(samples)
    reproducibility = malignant_reproducibility(programs_unselected)

    synthesis = (
        associations.loc[associations["method"] == "primary_gene_z_mean"]
        .groupby("relation", as_index=False)
        .agg(
            studies_tested=("spearman_rho", lambda x: int(x.notna().sum())),
            patients_total=("n_patients", "sum"),
            median_rho=("spearman_rho", "median"),
            positive_studies=("spearman_rho", lambda x: int((x > 0).sum())),
            negative_studies=("spearman_rho", lambda x: int((x < 0).sum())),
            fdr_lt_0_05_studies=("p_fdr_bh_within_dataset_method", lambda x: int((x < 0.05).sum())),
        )
    )
    expected_map = {row[0]: row[-1] for row in RELATIONS}
    synthesis["expected_direction"] = synthesis["relation"].map(expected_map)
    synthesis["studies_matching_expected"] = synthesis.apply(
        lambda row: row["positive_studies"] if row["expected_direction"] == 1 else (row["negative_studies"] if row["expected_direction"] == -1 else np.nan), axis=1
    )
    synthesis_status = "DESCRIPTIVE_CROSS_STUDY_SYNTHESIS" if FINAL_MALIGNANCY else "DESCRIPTIVE_CROSS_STUDY_SYNTHESIS_PROVISIONAL_MALIGNANCY"
    synthesis["status"] = np.where(synthesis["studies_tested"] >= 2, synthesis_status, "INSUFFICIENT_COMPARABLE_STUDIES")

    atomic_tsv(associations, OUT / f"within_dataset_patient_associations_{OUTPUT_VERSION}.tsv")
    atomic_tsv(subset_results, OUT / f"GSE178481_full_vs_added_subset_sensitivity_{OUTPUT_VERSION}.tsv")
    atomic_tsv(longitudinal, OUT / f"GSE264586_paired_longitudinal_{OUTPUT_VERSION}.tsv")
    atomic_tsv(response, OUT / f"GSE294109_response_stratified_{OUTPUT_VERSION}.tsv")
    atomic_tsv(reproducibility, OUT / f"malignant_F1_patient_reproducibility_{OUTPUT_VERSION}.tsv")
    atomic_tsv(sample_selection_audit, OUT / f"cross_sectional_sample_selection_audit_{OUTPUT_VERSION}.tsv")
    atomic_tsv(synthesis, OUT / f"cross_study_direction_synthesis_{OUTPUT_VERSION}.tsv")
    atomic_tsv(pd.DataFrame(input_rows), OUT / f"input_manifest_{OUTPUT_VERSION}.tsv")

    tested = associations["spearman_rho"].notna()
    status = {
        "status": PASS_STATUS,
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "datasets": DATASETS,
        "patients": int(programs["patient_key"].nunique()),
        "samples": int(samples[["dataset", "object_name"]].drop_duplicates().shape[0]),
        "association_rows": int(len(associations)),
        "association_rows_tested": int(tested.sum()),
        "minimum_cells": MIN_CELLS,
        "minimum_patients": MIN_PATIENTS,
        "bootstrap_replicates": N_BOOT,
        "seed": SEED,
        "limitations": [
            MALIGNANCY_LIMITATION,
            "GSE207493 old-versus-new 12/7 sensitivity is not reported until exact prior membership can be recovered; patient numbers were not guessed.",
            "GSE178481 old 5/new 6 patient membership is proven by the documented missing BM1/BM2/PR1-PR4 library rule.",
            "No cell-level pooled test was performed.",
        ],
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    status["cross_sectional_selection"] = {
        "GSE178481": "explicit tumor libraries only",
        "GSE264586": "Baseline only",
        "other_datasets": "all verified samples",
    }
    status["final_malignancy_gate"] = FINAL_MALIGNANCY
    status["malignancy_gate_label"] = MALIGNANCY_GATE_LABEL
    atomic_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", OUT / f"PATIENT_ANALYSIS_STATUS_{OUTPUT_VERSION}.json")

    manifest_name = f"FILE_MANIFEST_SHA256_{OUTPUT_VERSION}.tsv"
    files = sorted(path for path in OUT.iterdir() if path.name != manifest_name)
    manifest = pd.DataFrame([{"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)} for path in files])
    atomic_tsv(manifest, OUT / manifest_name)


if __name__ == "__main__":
    main()
