#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
RESULT_ROOT = PROJECT_ROOT / "results/part4_spatial"
OUTPUT_DIR = RESULT_ROOT / "advanced_patch_commot_cross_sample_v3"
SAMPLES = ("GLMF1", "GLMF2")
EXPECTED_AXES = 12
EXPECTED_SUMMARY_ROWS = 23
EXPECTED_CELLTYPE_ROWS = 23 * 14 * 14
TOP_DECILE_RANK_CEILING = math.ceil(14 * 14 * 0.10)


def read_commot_v2(sample: str, stem: str) -> pd.DataFrame:
    path = RESULT_ROOT / f"advanced_patch_commot_{sample.lower()}_v2" / f"{sample}_{stem}_v2.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    if "sample" not in frame or set(frame["sample"]) != {sample}:
        raise RuntimeError(f"Invalid sample column in {path}")
    return frame


def read_offdiagonal_v3(sample: str, stem: str) -> pd.DataFrame:
    path = RESULT_ROOT / "advanced_patch_commot_offdiagonal_v3" / sample / f"{sample}_{stem}_v3.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    if "sample" not in frame or set(frame["sample"]) != {sample}:
        raise RuntimeError(f"Invalid sample column in {path}")
    return frame


def add_rank(frame: pd.DataFrame, value: str, name: str) -> pd.DataFrame:
    result = frame.copy()
    result[name] = result.groupby("scale")[value].rank(method="min", ascending=False).astype(int)
    result[f"{name}_fraction"] = result[name] / result.groupby("scale")[value].transform("size")
    return result


def merge_samples(frames: dict[str, pd.DataFrame], keys: list[str], values: list[str]) -> pd.DataFrame:
    merged = None
    for sample in SAMPLES:
        selected = frames[sample][keys + values].copy()
        selected = selected.rename(columns={column: f"{column}_{sample}" for column in values})
        merged = selected if merged is None else merged.merge(selected, on=keys, how="outer", validate="one_to_one")
    if merged is None or merged.isna().any(axis=None):
        raise RuntimeError(f"Incomplete cross-sample merge for keys={keys}")
    return merged


def primary_scale(row: pd.Series) -> bool:
    return (row["signaling_type"] == "Cell-Cell Contact" and row["scale"] == "contact_110um") or (
        row["signaling_type"] != "Cell-Cell Contact" and row["scale"] == "secreted_200um"
    )


def plot_transport_rank(transport: pd.DataFrame) -> None:
    secreted = transport[transport["scale"].isin(["secreted_200um", "secreted_500um"])].copy()
    axis_ids = secreted["axis_id"].drop_duplicates().tolist()
    colors = dict(zip(axis_ids, plt.colormaps["tab20"](np.linspace(0, 1, len(axis_ids)))))
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    handles = []
    for panel_index, (axis, scale) in enumerate(zip(axes, ["secreted_200um", "secreted_500um"])):
        data = secreted[secreted["scale"] == scale]
        for _, row in data.iterrows():
            handle = axis.scatter(row["offdiagonal_score_rank_GLMF1"], row["offdiagonal_score_rank_GLMF2"],
                                  s=65, color=colors[row["axis_id"]], edgecolor="black", linewidth=0.3,
                                  label=row["axis_id"])
            if panel_index == 0:
                handles.append(handle)
        axis.plot([1, 11], [1, 11], linestyle="--", color="grey", linewidth=1)
        axis.set(xlim=(0.5, 11.5), ylim=(11.5, 0.5), xticks=range(1, 12), yticks=range(1, 12),
                 xlabel="GLMF1 within-scale rank", ylabel="GLMF2 within-scale rank", title=scale)
    fig.suptitle("Targeted COMMOT off-diagonal transport-score ranks across two sections")
    fig.legend(handles=handles, labels=axis_ids, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_transport_rank_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_transport_rank_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_hypothesis(primary: pd.DataFrame) -> None:
    ordered = primary.sort_values(["replicated_proxy_support"], ascending=False).reset_index(drop=True)
    values = ordered[["offdiagonal_proxy_log2_enrichment_GLMF1", "offdiagonal_proxy_log2_enrichment_GLMF2"]].to_numpy()
    cap = max(float(np.nanpercentile(np.abs(values), 95)), 0.5)
    fig, axis = plt.subplots(figsize=(7, 8), constrained_layout=True)
    image = axis.imshow(values, cmap="coolwarm", vmin=-cap, vmax=cap, aspect="auto")
    axis.set_xticks([0, 1], ["GLMF1", "GLMF2"])
    labels = [f"{row.axis_id}\n{row.source_proxy} -> {row.target_proxy}" for row in ordered.itertuples()]
    axis.set_yticks(range(len(ordered)), labels, fontsize=8)
    axis.set_title("Prespecified primary-scale composition-weighted proxies")
    for row_index, row in ordered.iterrows():
        for column, sample in enumerate(SAMPLES):
            value = row[f"offdiagonal_proxy_log2_enrichment_{sample}"]
            rank = int(row[f"offdiagonal_proxy_mass_rank_among_196_{sample}"])
            marker = "*" if row["replicated_proxy_support"] else ""
            axis.text(column, row_index, f"{value:.2f}\nrank {rank}{marker}", ha="center", va="center", fontsize=7)
    fig.colorbar(image, ax=axis, fraction=0.04, label="off-diagonal proxy log2 enrichment")
    fig.text(0.5, -0.01, "* both sections: positive enrichment and rank <= 20 after self-spot removal; still not cell-resolved.",
             ha="center", fontsize=8)
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_hypothesis_proxy_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_hypothesis_proxy_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_replicated_edges(edges: pd.DataFrame) -> None:
    supported = edges[edges["replicated_top_decile_positive_enrichment"]].copy()
    supported["mean_mass_fraction"] = supported[["offdiagonal_weighted_mass_fraction_GLMF1", "offdiagonal_weighted_mass_fraction_GLMF2"]].mean(axis=1)
    supported = supported.nlargest(30, "mean_mass_fraction").sort_values("mean_mass_fraction")
    fig, axis = plt.subplots(figsize=(11, max(5, 0.28 * len(supported) + 2)), constrained_layout=True)
    if supported.empty:
        axis.text(0.5, 0.5, "No replicated top-decile positive-enrichment proxy edges", ha="center", va="center")
        axis.set_axis_off()
    else:
        labels = [f"{row.axis_id} [{row.scale}]\n{row.source_celltype} -> {row.target_celltype}"
                  for row in supported.itertuples()]
        y = np.arange(len(supported))
        axis.barh(y - 0.18, supported["offdiagonal_weighted_mass_fraction_GLMF1"], height=0.34, label="GLMF1", color="#2166ac")
        axis.barh(y + 0.18, supported["offdiagonal_weighted_mass_fraction_GLMF2"], height=0.34, label="GLMF2", color="#b2182b")
        axis.set_yticks(y, labels, fontsize=7)
        axis.set_xlabel("off-diagonal composition-weighted edge-mass fraction")
        axis.legend(frameon=False)
    axis.set_title("Replicated proxy edges: top 20 mass rank and positive enrichment in both sections")
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_replicated_proxy_edges_v3.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "commot_cross_sample_replicated_proxy_edges_v3.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    summaries = {sample: read_commot_v2(sample, "commot_transport_summary") for sample in SAMPLES}
    hypotheses = {sample: read_offdiagonal_v3(sample, "commot_offdiagonal_hypothesis_proxy_support") for sample in SAMPLES}
    edges = {sample: read_offdiagonal_v3(sample, "commot_offdiagonal_celltype_edge_mass_proxy") for sample in SAMPLES}
    for sample in SAMPLES:
        if len(summaries[sample]) != EXPECTED_SUMMARY_ROWS or len(hypotheses[sample]) != EXPECTED_SUMMARY_ROWS:
            raise RuntimeError(f"Unexpected summary row count for {sample}")
        if len(edges[sample]) != EXPECTED_CELLTYPE_ROWS:
            raise RuntimeError(f"Unexpected edge row count for {sample}: {len(edges[sample])}")
        summaries[sample] = add_rank(summaries[sample], "offdiagonal_score", "offdiagonal_score_rank")

    transport_values = [
        "total_score", "n_edges", "diagonal_share", "offdiagonal_score", "offdiagonal_edges",
        "weighted_mean_distance_um", "weighted_median_distance_um",
        "weighted_p90_distance_um", "score_share_within_110um", "score_share_within_200um",
        "offdiagonal_score_rank", "offdiagonal_score_rank_fraction",
    ]
    transport = merge_samples(summaries, ["scale", "axis_id", "ligand", "receptor", "signaling_type"], transport_values)
    transport.to_csv(OUTPUT_DIR / "commot_cross_sample_transport_summary_v3.tsv", sep="\t", index=False)

    concordance_rows = []
    for scale in ["secreted_200um", "secreted_500um"]:
        data = transport[transport["scale"] == scale]
        statistic, pvalue = spearmanr(data["offdiagonal_score_rank_GLMF1"], data["offdiagonal_score_rank_GLMF2"])
        concordance_rows.append({"scale": scale, "n_axes": len(data), "spearman_rho": statistic,
                                 "nominal_two_sided_p": pvalue,
                                 "interpretation_limit": "targeted-axis rank concordance across two sections; not patient-level inference"})
    concordance = pd.DataFrame(concordance_rows)
    concordance.to_csv(OUTPUT_DIR / "commot_cross_sample_transport_rank_concordance_v3.tsv", sep="\t", index=False)

    hypothesis_values = ["source_proxy", "target_proxy", "offdiagonal_proxy_weighted_mass_fraction",
                         "offdiagonal_proxy_expected_fraction", "offdiagonal_proxy_enrichment",
                         "offdiagonal_proxy_log2_enrichment", "offdiagonal_proxy_mass_rank_among_196"]
    hypothesis = merge_samples(hypotheses, ["scale", "axis_id", "ligand", "receptor"], hypothesis_values)
    for column in ["source_proxy", "target_proxy"]:
        if not (hypothesis[f"{column}_GLMF1"] == hypothesis[f"{column}_GLMF2"]).all():
            raise RuntimeError(f"Proxy definition mismatch: {column}")
        hypothesis[column] = hypothesis.pop(f"{column}_GLMF1")
        hypothesis.drop(columns=f"{column}_GLMF2", inplace=True)
    hypothesis["both_positive_proxy_enrichment"] = (
        (hypothesis["offdiagonal_proxy_log2_enrichment_GLMF1"] > 0)
        & (hypothesis["offdiagonal_proxy_log2_enrichment_GLMF2"] > 0)
    )
    hypothesis["both_top_decile_mass"] = (
        (hypothesis["offdiagonal_proxy_mass_rank_among_196_GLMF1"] <= TOP_DECILE_RANK_CEILING)
        & (hypothesis["offdiagonal_proxy_mass_rank_among_196_GLMF2"] <= TOP_DECILE_RANK_CEILING)
    )
    hypothesis["replicated_proxy_support"] = hypothesis["both_positive_proxy_enrichment"] & hypothesis["both_top_decile_mass"]
    hypothesis["replication_limit"] = "two-section recurrence only; not independent-patient replication or causal evidence"
    hypothesis.to_csv(OUTPUT_DIR / "commot_cross_sample_hypothesis_proxy_v3.tsv", sep="\t", index=False)

    signaling = summaries["GLMF1"][["axis_id", "signaling_type"]].drop_duplicates()
    primary = hypothesis.merge(signaling, on="axis_id", how="left", validate="many_to_one")
    primary = primary[primary.apply(primary_scale, axis=1)].copy()
    if len(primary) != EXPECTED_AXES:
        raise RuntimeError(f"Expected {EXPECTED_AXES} primary hypotheses, observed {len(primary)}")
    sensitivity = hypothesis[hypothesis["scale"] == "secreted_500um"][[
        "axis_id", "both_positive_proxy_enrichment", "both_top_decile_mass", "replicated_proxy_support"
    ]].rename(columns={
        "both_positive_proxy_enrichment": "sensitivity_500um_both_positive_enrichment",
        "both_top_decile_mass": "sensitivity_500um_both_top_decile_mass",
        "replicated_proxy_support": "sensitivity_500um_replicated_proxy_support",
    })
    primary = primary.merge(sensitivity, on="axis_id", how="left", validate="one_to_one")
    is_contact = primary["signaling_type"] == "Cell-Cell Contact"
    primary.loc[is_contact, ["sensitivity_500um_both_positive_enrichment", "sensitivity_500um_both_top_decile_mass",
                             "sensitivity_500um_replicated_proxy_support"]] = np.nan
    primary.to_csv(OUTPUT_DIR / "commot_cross_sample_primary_hypothesis_support_v3.tsv", sep="\t", index=False)

    edge_values = ["offdiagonal_weighted_mass_fraction", "offdiagonal_proxy_enrichment",
                   "offdiagonal_proxy_log2_enrichment", "offdiagonal_mass_rank_among_196"]
    edge_cross = merge_samples(edges, ["scale", "axis_id", "ligand", "receptor", "source_celltype", "target_celltype"], edge_values)
    edge_cross["replicated_top_decile_positive_enrichment"] = (
        (edge_cross["offdiagonal_mass_rank_among_196_GLMF1"] <= TOP_DECILE_RANK_CEILING)
        & (edge_cross["offdiagonal_mass_rank_among_196_GLMF2"] <= TOP_DECILE_RANK_CEILING)
        & (edge_cross["offdiagonal_proxy_log2_enrichment_GLMF1"] > 0)
        & (edge_cross["offdiagonal_proxy_log2_enrichment_GLMF2"] > 0)
    )
    edge_cross["replication_limit"] = "off-diagonal composition-weighted mixed-spot proxy in two sections; not cell-resolved attribution"
    edge_cross.to_csv(OUTPUT_DIR / "commot_cross_sample_celltype_edge_proxy_v3.tsv", sep="\t", index=False)
    edge_cross[edge_cross["replicated_top_decile_positive_enrichment"]].to_csv(
        OUTPUT_DIR / "commot_cross_sample_replicated_top_decile_proxy_edges_v3.tsv", sep="\t", index=False
    )

    plot_transport_rank(transport)
    plot_hypothesis(primary)
    plot_replicated_edges(edge_cross)
    metadata = {
        "status": "complete",
        "samples": list(SAMPLES),
        "n_targeted_axes": EXPECTED_AXES,
        "n_primary_hypotheses": len(primary),
        "top_decile_definition": f"mass rank <= ceil(196*0.10) = {TOP_DECILE_RANK_CEILING}",
        "primary_proxy_support_definition": "positive off-diagonal proxy log2 enrichment and mass rank <=20 in both sections",
        "n_primary_replicated_proxy_support": int(primary["replicated_proxy_support"].sum()),
        "n_replicated_top_decile_positive_enrichment_edges": int(edge_cross["replicated_top_decile_positive_enrichment"].sum()),
        "interpretation_limit": (
            "targeted COMMOT off-diagonal transport and composition-weighted mixed-spot proxy recurrence across two sections; "
            "not direct binding, cell-resolved communication, causality, or independent-patient replication"
        ),
    }
    (OUTPUT_DIR / "run_metadata_v3.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (OUTPUT_DIR / "status.json").write_text(json.dumps({"status": "complete"}, indent=2) + "\n")
    print(f"COMMOT_CROSS_SAMPLE_ALL_PASS primary_supported={metadata['n_primary_replicated_proxy_support']} "
          f"replicated_edges={metadata['n_replicated_top_decile_positive_enrichment_edges']}")


if __name__ == "__main__":
    main()
