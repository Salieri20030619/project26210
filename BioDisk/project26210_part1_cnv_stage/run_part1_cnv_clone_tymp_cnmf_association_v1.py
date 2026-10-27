#!/usr/bin/env python3
"""Patient-level association of inferred CNV clones, TYMP, and rank-7 cNMF programs."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
CLONES = ROOT / "04_subclones_tymp/cnv_subclones_v1"
CLONE_VERIFICATION = ROOT / "04_subclones_tymp/cnv_subclones_verification_v1"
EXTENSION = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
CNMF = ROOT / "05_cnmf/full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz"
CNMF_SELECTION = ROOT / "05_cnmf/tymp_study_sensitivity_v1/selected_rank_v1.tsv"
OUT = ROOT / "04_subclones_tymp/cnv_clone_tymp_cnmf_association_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
THRESHOLDS = (0.225, 0.250, 0.275)
FACTORS = tuple(f"F{i}" for i in range(1, 8))
PRIMARY_MINIMUM = 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(p), dtype=float)
    result[order] = np.minimum(ranked, 1.0)
    return result


def signed_summary(values: pd.Series) -> dict[str, object]:
    x = values.dropna().to_numpy(float)
    nonzero = x[~np.isclose(x, 0)]
    if len(nonzero):
        wilcox_p = float(wilcoxon(nonzero, alternative="two-sided", zero_method="wilcox", method="auto").pvalue)
        positives = int(np.sum(nonzero > 0))
        sign_p = float(binomtest(positives, len(nonzero), 0.5, alternative="two-sided").pvalue)
    else:
        wilcox_p, positives, sign_p = 1.0, 0, 1.0
    return {
        "n_patients": int(len(x)),
        "n_nonzero_patient_deltas": int(len(nonzero)),
        "median_patient_tymp_delta": float(np.median(x)),
        "mean_patient_tymp_delta": float(np.mean(x)),
        "positive_patient_fraction": float(np.mean(x > 0)),
        "wilcoxon_two_sided_p": wilcox_p,
        "sign_test_two_sided_p": sign_p,
        "n_positive_nonzero": positives,
    }


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty association directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gate = subprocess.run([str(MEMORY_GATE), "cnv_clone_tymp_cnmf_association_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    if gate.returncode != 0 or "PASS" not in gate_text:
        raise RuntimeError(f"Memory gate blocked association: {gate_text}")
    if not (CLONE_VERIFICATION / "RUN_COMPLETE").is_file():
        raise RuntimeError("Formal CNV clone verification is required")

    selection = pd.read_csv(CNMF_SELECTION, sep="\t")
    if len(selection) != 1 or int(selection.iloc[0]["selected_rank"]) != 7 or selection.iloc[0]["selected_TYMP_factor"] != "F1":
        raise RuntimeError("Unexpected cNMF selection record")
    cells = pd.read_csv(CLONES / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    usage = pd.read_csv(CNMF, sep="\t", dtype={"exact_cell_id": str})
    grid = pd.read_csv(EXTENSION / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t")
    base = cells[["exact_cell_id", "study", "infercnv_unit", "infercnv_group", "TYMP_raw_count_v1", "TYMP_detected_v1", "TYMP_log1p_CP10K_v1", "threshold_sensitive_patient_v1"]].merge(
        usage, on="exact_cell_id", how="inner", validate="one_to_one"
    )
    if len(base) != 11950 or base["exact_cell_id"].duplicated().any():
        raise RuntimeError("Association cell input failed exact identity merge")
    if not np.allclose(base[list(FACTORS)].sum(axis=1), 1.0, atol=1e-6, rtol=0):
        raise RuntimeError("Rank-7 usage rows do not sum to one")

    clone_frames = []
    patient_contrasts = []
    patient_heterogeneity = []
    for threshold in THRESHOLDS:
        mapping = grid.loc[np.isclose(grid["distance_threshold"], threshold), ["infercnv_unit", "infercnv_group", "clone_id", "clone_cells"]]
        frame = base.merge(mapping, on=["infercnv_unit", "infercnv_group"], how="left", validate="many_to_one")
        if frame["clone_id"].isna().any() or len(frame) != 11950:
            raise RuntimeError(f"Threshold {threshold}: incomplete clone mapping")
        frame["primary_clone_ge20"] = frame["clone_cells"].ge(PRIMARY_MINIMUM)
        primary_counts = mapping.loc[mapping["clone_cells"].ge(PRIMARY_MINIMUM)].groupby("infercnv_unit")["clone_id"].nunique()
        eligible_units = set(primary_counts.loc[primary_counts >= 2].index)
        eligible = frame.loc[frame["primary_clone_ge20"] & frame["infercnv_unit"].isin(eligible_units)].copy()
        agg: dict[str, tuple[str, str]] = {
            "n_cells": ("exact_cell_id", "size"),
            "TYMP_median_log1p_CP10K": ("TYMP_log1p_CP10K_v1", "median"),
            "TYMP_mean_log1p_CP10K": ("TYMP_log1p_CP10K_v1", "mean"),
            "TYMP_detection_fraction": ("TYMP_detected_v1", "mean"),
            "threshold_sensitive_patient_v1": ("threshold_sensitive_patient_v1", "first"),
        }
        for factor in FACTORS:
            agg[f"{factor}_median_usage"] = (factor, "median")
            agg[f"{factor}_mean_usage"] = (factor, "mean")
        clone_summary = eligible.groupby(["study", "infercnv_unit", "clone_id"], observed=True).agg(**agg).reset_index()
        clone_summary.insert(0, "distance_threshold", threshold)
        clone_frames.append(clone_summary)

        for (study, unit), patient_clones in clone_summary.groupby(["study", "infercnv_unit"], sort=True, observed=True):
            patient_heterogeneity.append({
                "distance_threshold": threshold,
                "study": study,
                "infercnv_unit": unit,
                "n_primary_clones": int(len(patient_clones)),
                "n_primary_clone_cells": int(patient_clones["n_cells"].sum()),
                "TYMP_clone_median_range": float(patient_clones["TYMP_median_log1p_CP10K"].max() - patient_clones["TYMP_median_log1p_CP10K"].min()),
                **{f"{factor}_clone_median_range": float(patient_clones[f"{factor}_median_usage"].max() - patient_clones[f"{factor}_median_usage"].min()) for factor in FACTORS},
            })
            for factor in FACTORS:
                ordered = patient_clones.sort_values([f"{factor}_median_usage", "clone_id"], ascending=[False, True], kind="stable")
                high = ordered.iloc[0]
                others = ordered.iloc[1:]
                patient_contrasts.append({
                    "distance_threshold": threshold,
                    "study": study,
                    "infercnv_unit": unit,
                    "factor": factor,
                    "preselected_TYMP_factor": bool(factor == "F1"),
                    "n_primary_clones": int(len(ordered)),
                    "n_primary_clone_cells": int(ordered["n_cells"].sum()),
                    "factor_high_clone_id": high["clone_id"],
                    "factor_high_clone_median_usage": float(high[f"{factor}_median_usage"]),
                    "other_clones_unweighted_mean_median_usage": float(others[f"{factor}_median_usage"].mean()),
                    "factor_median_usage_separation": float(high[f"{factor}_median_usage"] - others[f"{factor}_median_usage"].mean()),
                    "factor_high_clone_TYMP_median": float(high["TYMP_median_log1p_CP10K"]),
                    "other_clones_unweighted_mean_TYMP_median": float(others["TYMP_median_log1p_CP10K"].mean()),
                    "patient_TYMP_delta_factor_high_vs_other_clones": float(high["TYMP_median_log1p_CP10K"] - others["TYMP_median_log1p_CP10K"].mean()),
                    "threshold_sensitive_patient_v1": bool(high["threshold_sensitive_patient_v1"]),
                })

    clone_summary_all = pd.concat(clone_frames, ignore_index=True)
    contrasts = pd.DataFrame(patient_contrasts)
    heterogeneity = pd.DataFrame(patient_heterogeneity)
    summary_rows = []
    for (threshold, factor), frame in contrasts.groupby(["distance_threshold", "factor"], sort=True, observed=True):
        summary_rows.append({"distance_threshold": threshold, "factor": factor, "preselected_TYMP_factor": bool(factor == "F1"), **signed_summary(frame["patient_TYMP_delta_factor_high_vs_other_clones"])})
    cohort_summary = pd.DataFrame(summary_rows)
    cohort_summary["wilcoxon_BH_FDR_within_threshold"] = cohort_summary.groupby("distance_threshold", group_keys=False)["wilcoxon_two_sided_p"].transform(bh)
    cohort_summary["sign_test_BH_FDR_within_threshold"] = cohort_summary.groupby("distance_threshold", group_keys=False)["sign_test_two_sided_p"].transform(bh)

    loo_rows = []
    f1 = contrasts.loc[contrasts["factor"].eq("F1")]
    for threshold, frame in f1.groupby("distance_threshold", sort=True, observed=True):
        studies = sorted(frame["study"].unique())
        for omitted in ["none", *studies]:
            retained = frame if omitted == "none" else frame.loc[~frame["study"].eq(omitted)]
            loo_rows.append({"distance_threshold": threshold, "omitted_study": omitted, **signed_summary(retained["patient_TYMP_delta_factor_high_vs_other_clones"])})
    loo = pd.DataFrame(loo_rows)
    study_rows = []
    selected_f1 = f1.loc[np.isclose(f1["distance_threshold"], 0.25)]
    for study, frame in selected_f1.groupby("study", sort=True, observed=True):
        study_rows.append({"study": study, **signed_summary(frame["patient_TYMP_delta_factor_high_vs_other_clones"])})
    study_summary = pd.DataFrame(study_rows)

    clone_summary_all.to_csv(OUT / "primary_clone_tymp_cnmf_summary_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    contrasts.to_csv(OUT / "patient_factor_high_clone_tymp_contrasts_v1.tsv", sep="\t", index=False)
    heterogeneity.to_csv(OUT / "patient_clone_program_heterogeneity_v1.tsv", sep="\t", index=False)
    cohort_summary.to_csv(OUT / "cohort_factor_contrast_summary_v1.tsv", sep="\t", index=False)
    loo.to_csv(OUT / "selected_F1_leave_one_study_out_v1.tsv", sep="\t", index=False)
    study_summary.to_csv(OUT / "selected_F1_study_stratified_v1.tsv", sep="\t", index=False)

    colors = {"0.225": "#4C78A8", "0.25": "#F58518", "0.275": "#54A24B"}
    fig, ax = plt.subplots(figsize=(9, 5.4))
    for threshold, frame in cohort_summary.groupby("distance_threshold"):
        x = np.arange(len(FACTORS)) + (threshold - 0.25) * 8
        ax.scatter(x, frame.set_index("factor").loc[list(FACTORS), "median_patient_tymp_delta"], s=55, label=f"threshold {threshold:.3f}", color=colors[str(threshold)])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(np.arange(len(FACTORS)), FACTORS)
    ax.set_ylabel("Patient TYMP median delta\n(factor-high clone minus other clones)")
    ax.set_title("CNV-clone / rank-7 program / TYMP patient-level concordance")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "factor_contrast_threshold_sensitivity_v1.png", dpi=180)
    plt.close(fig)

    selected_plot = selected_f1.sort_values("patient_TYMP_delta_factor_high_vs_other_clones")
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    y = np.arange(len(selected_plot))
    bar_colors = np.where(selected_plot["threshold_sensitive_patient_v1"], "#E45756", "#4C78A8")
    ax.barh(y, selected_plot["patient_TYMP_delta_factor_high_vs_other_clones"], color=bar_colors)
    ax.set_yticks(y, selected_plot["infercnv_unit"], fontsize=7)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("TYMP median delta: F1-high clone minus other primary clones")
    ax.set_title("Selected threshold 0.250; red = threshold-sensitive patient")
    fig.tight_layout()
    fig.savefig(OUT / "selected_F1_patient_tymp_deltas_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for idx, threshold in enumerate(THRESHOLDS):
        values = f1.loc[np.isclose(f1["distance_threshold"], threshold), "patient_TYMP_delta_factor_high_vs_other_clones"].to_numpy()
        jitter = np.linspace(-0.08, 0.08, len(values))
        ax.scatter(np.full(len(values), idx) + jitter, values, color=colors[str(threshold)], s=32, alpha=0.8)
        ax.plot([idx - 0.22, idx + 0.22], [np.median(values)] * 2, color="black", linewidth=2)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(range(3), [f"{value:.3f}" for value in THRESHOLDS])
    ax.set_xlabel("CNV profile distance threshold")
    ax.set_ylabel("Patient TYMP median delta for F1-high clone")
    ax.set_title("Mandatory CNV-threshold sensitivity")
    fig.tight_layout()
    fig.savefig(OUT / "selected_F1_threshold_patient_sensitivity_v1.png", dpi=180)
    plt.close(fig)

    formal_clones = clone_summary_all.loc[np.isclose(clone_summary_all["distance_threshold"], 0.25)].copy()
    formal_clones["clone_rank"] = formal_clones.groupby("infercnv_unit")["n_cells"].rank(method="first", ascending=False).astype(int)
    units = sorted(formal_clones["infercnv_unit"].unique())
    fig, ax = plt.subplots(figsize=(9.5, 7.0))
    for y_index, unit in enumerate(units):
        frame = formal_clones.loc[formal_clones["infercnv_unit"].eq(unit)]
        scatter = ax.scatter(frame["clone_rank"], np.full(len(frame), y_index), s=np.sqrt(frame["n_cells"]) * 14, c=frame["TYMP_median_log1p_CP10K"], cmap="viridis", vmin=0)
    ax.set_yticks(np.arange(len(units)), units, fontsize=7)
    ax.set_xlabel("Primary CNV clone rank by cell count within patient")
    ax.set_title("Formal primary inferred-CNV clones (size = cells; color = TYMP median)")
    fig.colorbar(scatter, ax=ax, label="TYMP median log1p(CP10K)")
    fig.tight_layout()
    fig.savefig(OUT / "formal_primary_clone_tymp_landscape_v1.png", dpi=180)
    plt.close(fig)

    metadata = {
        "schema_version": "part1_cnv_clone_tymp_cnmf_association_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_input_cells": 11950,
        "formal_threshold": 0.25,
        "mandatory_sensitivity_thresholds": [0.225, 0.275],
        "primary_clone_minimum_cells": 20,
        "selected_cnmf_rank": 7,
        "preselected_TYMP_factor": "F1",
        "statistical_unit": "patient",
        "primary_contrast": "per patient, TYMP median of factor-high primary clone minus unweighted mean TYMP medians of other primary clones",
        "interpretation_limits": [
            "cNMF F1 was preselected using TYMP evidence; association is concordance, not independent confirmation",
            "inferred CNV clones are not mutation-validated genetic clones",
            "no causal TYMP regulation claim",
            "F1 has a documented GSE207493 study-composition warning",
        ],
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({"file": [p.name for p in files], "size_bytes": [p.stat().st_size for p in files], "sha256": [sha256(p) for p in files]}).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(cohort_summary.to_string(index=False))
    print(loo.to_string(index=False))
    print("PART1_CNV_CLONE_TYMP_CNMF_ASSOCIATION_V1_PASS")


if __name__ == "__main__":
    main()
