#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT_ROOT / "results/part4_spatial"
OUTPUT_DIR = RESULT_ROOT / "advanced_patch_misty_cross_sample_v3"
INPUT_DIR = RESULT_ROOT / "advanced_patch_misty_inputs_v1"
MEMORY_GATE = PROJECT_ROOT / "check_memory_gate.sh"
SAMPLES = ("GLMF1", "GLMF2")
VIEWS = ("intra", "juxta.110", "para.200")


def memory_gate(stage: str) -> None:
    print(f"MEMORY_GATE stage={stage}", flush=True)
    subprocess.run([str(MEMORY_GATE)], check=True)


def sample_dir(sample: str) -> Path:
    return RESULT_ROOT / f"advanced_patch_misty_{sample.lower()}_v1"


def short_labels(feature_ids: list[str]) -> dict[str, str]:
    result = {}
    for feature in feature_ids:
        label = feature.replace("cell_", "").replace("program_", "P:").replace("_", " ")
        label = label.replace("APOE C1QC TAM", "APOE+C1QC TAM")
        label = label.replace("TYMP High Malignant", "TYMP-High Malig.")
        label = label.replace("TYMP Low Malignant", "TYMP-Low Malig.")
        label = label.replace("Other Macrophage", "Other Macro.")
        label = label.replace("Epithelial Normal", "Epi-Normal")
        label = label.replace("Tumor epithelial", "Tumor-epi")
        label = label.replace("TGFb response", "TGFb")
        label = label.replace("shared stromal", "Shared-stromal")
        result[feature] = label
    return result


def load_results() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    manifest = pd.read_csv(INPUT_DIR / "misty_feature_manifest_v1.tsv", sep="\t")
    feature_ids = manifest["feature_id"].tolist()
    performance_parts = []
    contribution_parts = []
    importance_parts = []
    for sample in SAMPLES:
        directory = sample_dir(sample)
        improvements = pd.read_csv(directory / f"{sample}_misty_improvement_stats_v1.tsv", sep="\t")
        improvement_wide = improvements.pivot(index="target", columns="measure", values="mean").reset_index()
        raw_performance = pd.read_csv(directory / "raw_misty/performance.txt", sep=r"\s+")
        improvement_wide = improvement_wide.merge(
            raw_performance[["target", "p.RMSE", "p.R2"]], on="target", validate="one_to_one"
        )
        improvement_wide.insert(0, "sample", sample)
        performance_parts.append(improvement_wide)

        contributions = pd.read_csv(directory / f"{sample}_misty_contribution_stats_v1.tsv", sep="\t")
        contributions.insert(0, "sample", sample)
        contribution_parts.append(contributions)

        importance = pd.read_csv(directory / f"{sample}_misty_importances_v1.tsv", sep="\t")
        importance = importance.loc[importance["Predictor"] != importance["Target"]].copy()
        importance["rank"] = importance.groupby("view")["Importance"].rank(method="min", ascending=False)
        importance["rank_fraction"] = importance["rank"] / importance.groupby("view")["Importance"].transform("count")
        importance_parts.append(importance)
    return (
        pd.concat(performance_parts, ignore_index=True),
        pd.concat(contribution_parts, ignore_index=True),
        pd.concat(importance_parts, ignore_index=True),
        feature_ids,
    )


def build_cross_tables(
    performance: pd.DataFrame,
    contributions: pd.DataFrame,
    importances: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    perf_wide = performance.pivot(index="target", columns="sample")
    perf_wide.columns = [f"{measure}_{sample}" for measure, sample in perf_wide.columns]
    perf_wide = perf_wide.reset_index()
    perf_wide["positive_gain_both_sections"] = (
        (perf_wide["gain.R2_GLMF1"] > 0)
        & (perf_wide["gain.R2_GLMF2"] > 0)
        & (perf_wide["gain.RMSE_GLMF1"] > 0)
        & (perf_wide["gain.RMSE_GLMF2"] > 0)
    )
    perf_wide["r2_gain_p_lt_0.05_both_sections"] = (
        (perf_wide["p.R2_GLMF1"] < 0.05) & (perf_wide["p.R2_GLMF2"] < 0.05)
    )
    perf_wide["rmse_gain_p_lt_0.05_both_sections"] = (
        (perf_wide["p.RMSE_GLMF1"] < 0.05) & (perf_wide["p.RMSE_GLMF2"] < 0.05)
    )

    contribution_wide = contributions.pivot(index=["target", "view"], columns="sample")
    contribution_wide.columns = [f"{measure}_{sample}" for measure, sample in contribution_wide.columns]
    contribution_wide = contribution_wide.reset_index()
    contribution_wide["mean_fraction"] = contribution_wide[["fraction_GLMF1", "fraction_GLMF2"]].mean(axis=1)
    contribution_wide["supported_p_lt_0.05_both_sections"] = (
        (contribution_wide["p.mean_GLMF1"] < 0.05) & (contribution_wide["p.mean_GLMF2"] < 0.05)
    )

    importance_wide = importances.pivot(
        index=["view", "Predictor", "Target"], columns="sample", values=["Importance", "rank", "rank_fraction"]
    )
    importance_wide.columns = [f"{measure}_{sample}" for measure, sample in importance_wide.columns]
    importance_wide = importance_wide.reset_index()
    importance_wide["mean_rank_fraction"] = importance_wide[
        ["rank_fraction_GLMF1", "rank_fraction_GLMF2"]
    ].mean(axis=1)
    importance_wide["replicated_top_decile"] = (
        (importance_wide["rank_fraction_GLMF1"] <= 0.10)
        & (importance_wide["rank_fraction_GLMF2"] <= 0.10)
        & (importance_wide["Importance_GLMF1"] > 0)
        & (importance_wide["Importance_GLMF2"] > 0)
    )
    importance_wide["interpretation_limit"] = (
        "relative multivariate predictive importance; not an effect direction or causal estimate"
    )
    return perf_wide, contribution_wide, importance_wide


def plot_performance(performance: pd.DataFrame, labels: dict[str, str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 8), constrained_layout=True)
    for axis, measure, title in zip(
        axes,
        ("gain.R2", "gain.RMSE"),
        ("R2 gain (percentage points)", "RMSE improvement (%)"),
    ):
        matrix = performance.pivot(index="target", columns="sample", values=measure)
        matrix = matrix.loc[list(labels)]
        values = matrix.to_numpy()
        cap = max(abs(np.nanpercentile(values, 2)), abs(np.nanpercentile(values, 98)), 1e-6)
        image = axis.imshow(values, cmap="coolwarm", vmin=-cap, vmax=cap, aspect="auto")
        axis.set_xticks(range(2), matrix.columns)
        axis.set_yticks(range(len(matrix)), [labels[item] for item in matrix.index])
        axis.set_title(title)
        for row in range(values.shape[0]):
            for column in range(values.shape[1]):
                axis.text(column, row, f"{values[row, column]:.1f}", ha="center", va="center", fontsize=7)
        fig.colorbar(image, ax=axis, fraction=0.04, pad=0.03)
    fig.suptitle("Cross-section MISTy multiview performance gains", fontsize=13)
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_performance_gains_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_performance_gains_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_contributions(contributions: pd.DataFrame, labels: dict[str, str]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), constrained_layout=True)
    colors = {"intra": "#4C78A8", "juxta.110": "#F58518", "para.200": "#54A24B"}
    titles = {"intra": "Intra", "juxta.110": "Juxta <=110 um", "para.200": "Para >110 um (l=200 um)"}
    program_targets = {item for item in labels if item.startswith("program_")}
    for axis, view in zip(axes, VIEWS):
        current = contributions.loc[contributions["view"] == view].pivot(
            index="target", columns="sample", values="fraction"
        )
        axis.scatter(current["GLMF1"], current["GLMF2"], s=28, color=colors[view], alpha=0.78)
        limit = max(float(current.to_numpy().max()) * 1.08, 0.1)
        axis.plot([0, limit], [0, limit], linestyle="--", linewidth=0.8, color="#999999")
        for target in current.index:
            if target in program_targets:
                axis.annotate(labels[target], (current.loc[target, "GLMF1"], current.loc[target, "GLMF2"]),
                    xytext=(3, 3), textcoords="offset points", fontsize=6.5)
        axis.set_xlim(0, limit)
        axis.set_ylim(0, limit)
        axis.set_title(titles[view])
        axis.set_xlabel("GLMF1 contribution fraction")
        axis.set_ylabel("GLMF2 contribution fraction")
    fig.suptitle("Cross-section comparison of relative MISTy view contributions", fontsize=13)
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_view_contributions_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_view_contributions_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_importance_scatter(importance: pd.DataFrame, labels: dict[str, str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    titles = {"juxta.110": "Juxta <=110 um", "para.200": "Para >110 um (l=200 um)"}
    for axis, view in zip(axes, ("juxta.110", "para.200")):
        current = importance.loc[importance["view"] == view].copy()
        axis.scatter(current["Importance_GLMF1"], current["Importance_GLMF2"],
            s=12, color="#B3B3B3", alpha=0.55, linewidth=0)
        replicated = current.loc[current["replicated_top_decile"]].sort_values("mean_rank_fraction")
        axis.scatter(replicated["Importance_GLMF1"], replicated["Importance_GLMF2"],
            s=34, color="#B40426", alpha=0.85, label=f"Top decile in both (n={len(replicated)})")
        for label_index, row in enumerate(replicated.head(5).itertuples(index=False)):
            label = f"{labels[row.Predictor]} -> {labels[row.Target]}"
            axis.annotate(label, (row.Importance_GLMF1, row.Importance_GLMF2),
                xytext=(0.46, 0.96 - 0.055 * label_index), textcoords="axes fraction",
                fontsize=6.5, ha="left", va="top",
                arrowprops={"arrowstyle": "-", "color": "#666666", "linewidth": 0.45})
        axis.axhline(0, color="#CCCCCC", linewidth=0.6)
        axis.axvline(0, color="#CCCCCC", linewidth=0.6)
        axis.set_title(titles[view])
        axis.set_xlabel("GLMF1 relative importance")
        axis.set_ylabel("GLMF2 relative importance")
        axis.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Cross-section replicated MISTy predictor-target importance", fontsize=13)
    fig.text(0.5, -0.02, "Importance is scaled within each target-view; arrows denote prediction, not causal direction.",
        ha="center", fontsize=9)
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_importance_replication_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "cross_sample_misty_importance_replication_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    memory_gate("misty_cross_sample_summary")
    performance, contributions, importances, feature_ids = load_results()
    if performance.shape[0] != 40 or contributions.shape[0] != 120 or importances.shape[0] != 2280:
        raise RuntimeError("Unexpected per-section MISTy result dimensions")
    perf_cross, contrib_cross, importance_cross = build_cross_tables(performance, contributions, importances)
    replicated = importance_cross.loc[importance_cross["replicated_top_decile"]].copy()
    replicated = replicated.sort_values(["view", "mean_rank_fraction", "Predictor", "Target"])

    perf_cross.to_csv(OUTPUT_DIR / "cross_sample_misty_performance_v3.tsv", sep="\t", index=False)
    contrib_cross.to_csv(OUTPUT_DIR / "cross_sample_misty_view_contributions_v3.tsv", sep="\t", index=False)
    importance_cross.to_csv(OUTPUT_DIR / "cross_sample_misty_importance_support_v3.tsv", sep="\t", index=False)
    replicated.to_csv(OUTPUT_DIR / "cross_sample_misty_replicated_top_decile_v3.tsv", sep="\t", index=False)

    labels = short_labels(feature_ids)
    plot_performance(performance, labels)
    plot_contributions(contributions, labels)
    plot_importance_scatter(importance_cross, labels)

    counts = replicated.groupby("view").size().reindex(VIEWS, fill_value=0).astype(int).to_dict()
    both_positive = int(perf_cross["positive_gain_both_sections"].sum())
    both_r2_sig = int(perf_cross["r2_gain_p_lt_0.05_both_sections"].sum())
    both_rmse_sig = int(perf_cross["rmse_gain_p_lt_0.05_both_sections"].sum())
    summary_lines = [
        "# Part4 MISTy cross-section summary v3",
        "",
        f"- Positive R2 and RMSE multiview gain in both sections: {both_positive}/20 targets.",
        f"- Nominal paired-fold p<0.05 for R2 gain in both sections: {both_r2_sig}/20 targets.",
        f"- Nominal paired-fold p<0.05 for RMSE gain in both sections: {both_rmse_sig}/20 targets.",
        f"- Replicated top-decile relative importance pairs: intra={counts['intra']}, "
        f"juxta={counts['juxta.110']}, para={counts['para.200']}.",
        "- Recurrent patterns include T/NK/CD8, stromal/pericyte/POSTN, macrophage/APOE-TAM, "
        "and malignant/epithelial programs across one or more spatial views.",
        "- MISTy importance is a multivariate predictive measure. It is not a signed effect, "
        "a linear coefficient, a ligand-receptor claim, or a causal estimate.",
    ]
    (OUTPUT_DIR / "misty_cross_sample_summary_v3.md").write_text("\n".join(summary_lines) + "\n")
    metadata = {
        "status": "complete",
        "n_targets": 20,
        "positive_gain_both_sections": both_positive,
        "nominal_r2_gain_p_lt_0.05_both_sections": both_r2_sig,
        "nominal_rmse_gain_p_lt_0.05_both_sections": both_rmse_sig,
        "replicated_top_decile_counts": counts,
        "replication_definition": (
            "predictor-target pair has positive scaled importance and ranks in the top 10% "
            "of all non-self pairs within the same view in each section"
        ),
        "interpretation_limit": (
            "MISTy importance is relative multivariate predictive importance, not effect direction or causality"
        ),
    }
    (OUTPUT_DIR / "cross_sample_run_metadata_v3.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n"
    )
    print(
        "MISTY_CROSS_SAMPLE_ALL_PASS "
        f"positive_both={both_positive} replicated_intra={counts['intra']} "
        f"replicated_juxta={counts['juxta.110']} replicated_para={counts['para.200']}"
    )


if __name__ == "__main__":
    main()
