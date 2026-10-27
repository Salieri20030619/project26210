#!/usr/bin/env python3
"""Patient-level synthesis of the TJ-RCC 12-sample advanced spatial analyses."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
GATE = PROJECT / "check_memory_gate.sh"
SAMPLES = ("H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T", "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T")
PROGRAMS = {
    "Tumor-epithelial": ("CA9", "EPCAM", "KRT8", "KRT18", "KRT19"),
    "POSTN-ECM": ("POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC"),
    "APOE-TAM": ("APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LPL"),
    "CD8-effector": ("CD8A", "CD8B", "NKG7", "CCL5", "PRF1", "GZMB"),
    "TGFb-response": ("TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF"),
    "TYMP-context": ("TYMP", "SPP1", "MIF", "CD74", "CXCR4", "APP"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def memory_gate(stage: str) -> str:
    result = subprocess.run([str(GATE)], check=True, capture_output=True, text=True)
    line = result.stdout.strip()
    print(f"MEMORY_GATE stage={stage} {line}", flush=True)
    return line


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    return series.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False}).astype(bool)


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    finite = np.isfinite(p)
    out = np.full(len(p), np.nan)
    if not finite.any():
        return out
    work = p[finite]
    order = np.argsort(work)
    ranked = work[order]
    adjusted = np.minimum.accumulate((ranked * len(work) / np.arange(1, len(work) + 1))[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[np.flatnonzero(finite)] = restored
    return out


def signed_rank(values: pd.Series) -> float:
    array = values.dropna().to_numpy(float)
    if len(array) < 5 or np.allclose(array, 0):
        return math.nan
    try:
        return float(wilcoxon(array, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return math.nan


def aggregate_numeric(frame: pd.DataFrame, groups: list[str], value: str, prefix: str) -> pd.DataFrame:
    result = frame.groupby(groups, dropna=False)[value].agg(
        n_patients="count", median="median", mean="mean",
        q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75),
        positive_n=lambda x: int((x > 0).sum()),
    ).reset_index()
    result = result.rename(columns={name: f"{prefix}_{name}" for name in result.columns if name not in groups})
    return result


def load_bayesspace(output: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    patient_rows, effect_rows, metadata = [], [], []
    for sample in SAMPLES:
        directory = ROOT / "bayesspace" / sample
        status = read_json(directory / "status.json")
        bayes = read_json(directory / "bayesspace_run_metadata.json")
        boundary_meta = read_json(directory / "boundary_run_metadata.json")
        if status.get("status") != "completed" or bayes.get("status") != "bayesspace_complete" or boundary_meta.get("status") != "completed":
            raise RuntimeError(f"{sample}: BayesSpace is not complete")
        spots = pd.read_csv(directory / f"{sample}_computational_boundary_spots.tsv", sep="\t")
        tumor = as_bool(spots["computational_tumor"]).to_numpy()
        boundary = as_bool(spots["computational_boundary_spot"]).to_numpy()
        patient_rows.append({
            "sample": sample, "n_spots": len(spots), "selected_q": int(bayes["selected_q"]),
            "n_tumor_like_spots": int(tumor.sum()), "tumor_like_fraction": float(tumor.mean()),
            "n_boundary_spots": int(boundary.sum()), "boundary_spot_fraction": float(boundary.mean()),
            "median_absolute_boundary_distance_um": float(np.median(np.abs(spots["signed_distance_um"]))),
            "microns_per_pixel": float(boundary_meta["microns_per_pixel"]),
        })
        for program in ("Tumor-epithelial", "POSTN-ECM", "APOE-TAM", "CD8-effector", "TGFb-response", "APOE_TAM_q05", "T_cell_q05"):
            values = spots[program].to_numpy(float)
            effect_rows.append({
                "sample": sample, "feature": program,
                "tumor_median": float(np.nanmedian(values[tumor])),
                "non_tumor_median": float(np.nanmedian(values[~tumor])),
                "tumor_minus_non_tumor_median": float(np.nanmedian(values[tumor]) - np.nanmedian(values[~tumor])),
            })
        metadata.append({"sample": sample, "bayesspace": bayes, "boundary": boundary_meta})
    patients = pd.DataFrame(patient_rows)
    effects = pd.DataFrame(effect_rows)
    patients.to_csv(output / "bayesspace_patient_summary_v1.tsv", sep="\t", index=False)
    effects.to_csv(output / "bayesspace_patient_program_effects_v1.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    sns.barplot(data=patients, x="sample", y="selected_q", color="#4C78A8", ax=axes[0])
    sns.barplot(data=patients, x="sample", y="tumor_like_fraction", color="#E45756", ax=axes[1])
    sns.barplot(data=patients, x="sample", y="boundary_spot_fraction", color="#72B7B2", ax=axes[2])
    for axis, title in zip(axes, ("Selected BayesSpace q", "Computational tumor-like fraction", "Boundary-spot fraction")):
        axis.set_title(title); axis.tick_params(axis="x", rotation=55); axis.set_xlabel("")
    fig.savefig(output / "bayesspace_patient_overview_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / "bayesspace_patient_overview_v1.pdf", bbox_inches="tight")
    plt.close(fig)

    matrix = effects.pivot(index="feature", columns="sample", values="tumor_minus_non_tumor_median").reindex(columns=SAMPLES)
    fig, ax = plt.subplots(figsize=(12, 5.6))
    cap = max(float(np.nanpercentile(np.abs(matrix), 95)), 0.1)
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-cap, vmax=cap, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Patient-level computational tumor-minus-non-tumor contrasts")
    fig.tight_layout()
    fig.savefig(output / "bayesspace_program_effects_v1.png", dpi=300)
    fig.savefig(output / "bayesspace_program_effects_v1.pdf")
    plt.close(fig)
    return patients, effects, metadata


def load_spatial_programs(output: Path) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    genes, modules, metadata = [], [], []
    for sample in SAMPLES:
        directory = ROOT / "spatial_programs" / sample
        status = read_json(directory / "status.json")
        spark = read_json(directory / "sparkx_run_metadata.json")
        hot = read_json(directory / "hotspot_run_metadata.json")
        if status.get("status") != "completed" or spark.get("status") != "completed" or hot.get("status") != "completed":
            raise RuntimeError(f"{sample}: spatial programs are not complete")
        frame = pd.read_csv(directory / f"{sample}_sparkx_hotspot_gene_statistics.tsv", sep="\t")
        frame.insert(0, "sample", sample)
        frame["consensus_rank_fraction"] = frame["consensus_rank"] / len(frame)
        frame["hotspot_positive_strict"] = (frame["hotspot_fdr"] < 0.05) & (frame["hotspot_C"] > 0)
        genes.append(frame)
        try:
            module = pd.read_csv(directory / f"{sample}_hotspot_module_summary.tsv", sep="\t")
        except pd.errors.EmptyDataError:
            module = pd.DataFrame()
        if not module.empty:
            module.insert(0, "sample", sample)
            modules.append(module)
        metadata.append({"sample": sample, "sparkx": spark, "hotspot": hot})
    long = pd.concat(genes, ignore_index=True)
    cohort = long.groupby("gene", as_index=False).agg(
        n_patients_tested=("sample", "nunique"),
        n_both_strict=("support_class", lambda x: int((x == "both_strict").sum())),
        n_sparkx_strict=("sparkx_by_fdr", lambda x: int((x < 0.05).sum())),
        n_hotspot_strict=("hotspot_positive_strict", "sum"),
        median_consensus_rank_fraction=("consensus_rank_fraction", "median"),
        q25_consensus_rank_fraction=("consensus_rank_fraction", lambda x: x.quantile(0.25)),
        q75_consensus_rank_fraction=("consensus_rank_fraction", lambda x: x.quantile(0.75)),
    )
    cohort = cohort.sort_values(["n_both_strict", "n_sparkx_strict", "n_hotspot_strict", "median_consensus_rank_fraction", "gene"], ascending=[False, False, False, True, True])
    cohort.to_csv(output / "spatial_program_gene_replication_v1.tsv", sep="\t", index=False)
    if modules:
        module_cohort = pd.concat(modules, ignore_index=True)
    else:
        module_cohort = pd.DataFrame(columns=[
            "sample", "module", "n_genes", "n_both_strict", "n_sparkx_strict",
            "n_hotspot_strict", "top_genes", "candidate_program_overlaps",
            "distance_spearman_rho_descriptive", "tumor_median_score",
            "nontumor_median_score", "tumor_minus_nontumor_median",
        ])
    module_cohort.to_csv(output / "hotspot_patient_module_summary_v1.tsv", sep="\t", index=False)

    overlay_rows = []
    lookup = cohort.set_index("gene")
    for program, members in PROGRAMS.items():
        for gene in members:
            row = {"program": program, "gene": gene, "present_in_cohort_table": gene in lookup.index}
            if gene in lookup.index:
                row.update(lookup.loc[gene].to_dict())
            overlay_rows.append(row)
    overlay = pd.DataFrame(overlay_rows)
    overlay.to_csv(output / "candidate_program_gene_replication_v1.tsv", sep="\t", index=False)

    plot = overlay.pivot(index="program", columns="gene", values="n_both_strict")
    order = [gene for members in PROGRAMS.values() for gene in members if gene in plot.columns]
    plot = plot.reindex(index=PROGRAMS, columns=list(dict.fromkeys(order)))
    fig, ax = plt.subplots(figsize=(16, 5.5))
    sns.heatmap(plot, cmap="Purples", vmin=0, vmax=len(SAMPLES), annot=True, fmt=".0f", linewidths=0.3, ax=ax)
    ax.set_title("SPARK-X/Hotspot dual-method replication across 12 patients")
    ax.set_xlabel("Candidate-overlay gene (discovery remained genome-wide)"); ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output / "spatial_program_candidate_replication_v1.png", dpi=300)
    fig.savefig(output / "spatial_program_candidate_replication_v1.pdf")
    plt.close(fig)

    top = cohort.loc[cohort["n_patients_tested"] >= 6].head(30).sort_values("n_both_strict")
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.barh(top["gene"], top["n_both_strict"], color="#6A3D9A")
    ax.set_xlim(0, len(SAMPLES)); ax.set_xlabel("Patients with dual-method strict support")
    ax.set_title("Most recurrent unbiased spatially variable genes")
    fig.tight_layout()
    fig.savefig(output / "spatial_program_top_recurrent_genes_v1.png", dpi=300)
    fig.savefig(output / "spatial_program_top_recurrent_genes_v1.pdf")
    plt.close(fig)
    return cohort, overlay, metadata


def load_misty(output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    improvements, contributions, importances, metadata = [], [], [], []
    for sample in SAMPLES:
        directory = ROOT / "misty" / "samples" / sample
        status = read_json(directory / "status.json")
        meta = read_json(directory / "run_metadata_v1.json")
        if status.get("status") != "complete" or meta.get("status") != "complete":
            raise RuntimeError(f"{sample}: MISTy is not complete")
        imp = pd.read_csv(directory / f"{sample}_misty_improvement_stats_v1.tsv", sep="\t")
        imp.insert(0, "sample", sample); improvements.append(imp)
        con = pd.read_csv(directory / f"{sample}_misty_contribution_stats_v1.tsv", sep="\t")
        con.insert(0, "sample", sample); contributions.append(con)
        importance = pd.read_csv(directory / f"{sample}_misty_importances_v1.tsv", sep="\t")
        importance = importance.loc[importance["Predictor"] != importance["Target"]].copy()
        importance["sample"] = sample
        importance["rank"] = importance.groupby("view")["Importance"].rank(method="min", ascending=False)
        importance["rank_fraction"] = importance["rank"] / importance.groupby("view")["Importance"].transform("count")
        importances.append(importance)
        metadata.append(meta)
    improvement_long = pd.concat(improvements, ignore_index=True)
    contribution_long = pd.concat(contributions, ignore_index=True)
    importance_long = pd.concat(importances, ignore_index=True)

    improvement = improvement_long.groupby(["target", "measure"])["mean"].agg(
        n_patients="count", median="median", mean="mean", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75),
        positive_n=lambda x: int((x > 0).sum()), wilcoxon_two_sided_p=signed_rank,
    ).reset_index()
    improvement["wilcoxon_bh_fdr_within_measure"] = improvement.groupby("measure")["wilcoxon_two_sided_p"].transform(bh)
    contribution = contribution_long.groupby(["target", "view"])["fraction"].agg(
        n_patients="count", median_fraction="median", mean_fraction="mean",
        q25_fraction=lambda x: x.quantile(.25), q75_fraction=lambda x: x.quantile(.75),
    ).reset_index()
    importance = importance_long.groupby(["view", "Predictor", "Target"]).agg(
        n_patients=("sample", "nunique"), positive_importance_n=("Importance", lambda x: int((x > 0).sum())),
        top_decile_n=("rank_fraction", lambda x: int((x <= .10).sum())),
        median_importance=("Importance", "median"), median_rank_fraction=("rank_fraction", "median"),
    ).reset_index().sort_values(["top_decile_n", "median_rank_fraction"], ascending=[False, True])
    improvement.to_csv(output / "misty_patient_performance_summary_v1.tsv", sep="\t", index=False)
    contribution.to_csv(output / "misty_patient_view_contribution_summary_v1.tsv", sep="\t", index=False)
    importance.to_csv(output / "misty_patient_importance_replication_v1.tsv", sep="\t", index=False)

    perf = improvement.loc[improvement["measure"].isin(["gain.R2", "gain.RMSE"])].pivot(index="target", columns="measure", values="median")
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), constrained_layout=True)
    for axis, measure in zip(axes, ("gain.R2", "gain.RMSE")):
        values = perf[[measure]]
        cap = max(float(np.nanpercentile(np.abs(values), 98)), .1)
        sns.heatmap(values, cmap="vlag", center=0, vmin=-cap, vmax=cap, annot=True, fmt=".2f", ax=axis, cbar=False)
        axis.set_title(f"Patient median {measure}"); axis.set_xlabel(""); axis.set_ylabel("")
    fig.savefig(output / "misty_patient_performance_v1.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / "misty_patient_performance_v1.pdf", bbox_inches="tight")
    plt.close(fig)

    view = contribution.pivot(index="target", columns="view", values="median_fraction").reindex(columns=["intra", "juxta.110", "para.200"])
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(view, cmap="viridis", vmin=0, vmax=1, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Median MISTy view contribution across patients"); ax.set_xlabel(""); ax.set_ylabel("")
    fig.tight_layout(); fig.savefig(output / "misty_patient_view_contributions_v1.png", dpi=300); fig.savefig(output / "misty_patient_view_contributions_v1.pdf"); plt.close(fig)
    return improvement, contribution, importance, metadata


def load_commot(output: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    manifests, summaries, hypotheses, edges, domains, metadata = [], [], [], [], [], []
    for sample in SAMPLES:
        directory = ROOT / "commot" / "samples" / sample
        status = read_json(directory / "status.json")
        meta = read_json(directory / "run_metadata_v2.json")
        if status.get("status") != "complete" or meta.get("status") != "complete":
            raise RuntimeError(f"{sample}: COMMOT is not complete")
        manifests.append(pd.read_csv(directory / f"{sample}_commot_candidate_manifest_v2.tsv", sep="\t"))
        summaries.append(pd.read_csv(directory / f"{sample}_commot_transport_summary_v2.tsv", sep="\t"))
        hypotheses.append(pd.read_csv(directory / f"{sample}_commot_hypothesis_proxy_support_v2.tsv", sep="\t"))
        edges.append(pd.read_csv(directory / f"{sample}_commot_celltype_edge_mass_proxy_v2.tsv", sep="\t"))
        domains.append(pd.read_csv(directory / f"{sample}_commot_computational_domain_edge_mass_v2.tsv", sep="\t"))
        metadata.append(meta)
    manifest = pd.concat(manifests, ignore_index=True)
    summary = pd.concat(summaries, ignore_index=True)
    hypothesis = pd.concat(hypotheses, ignore_index=True)
    edge = pd.concat(edges, ignore_index=True)
    domain = pd.concat(domains, ignore_index=True)
    summary["transport_rank_fraction_within_sample_scale"] = summary.groupby(["sample", "scale"])["total_score"].rank(method="average", ascending=False) / summary.groupby(["sample", "scale"])["axis_id"].transform("count")
    cohort_transport = summary.groupby(["scale", "axis_id", "ligand", "receptor", "signaling_type"], dropna=False).agg(
        n_patients=("sample", "nunique"), median_total_score=("total_score", "median"),
        q25_total_score=("total_score", lambda x: x.quantile(.25)), q75_total_score=("total_score", lambda x: x.quantile(.75)),
        median_transport_rank_fraction=("transport_rank_fraction_within_sample_scale", "median"),
        median_weighted_p90_distance_um=("weighted_p90_distance_um", "median"),
    ).reset_index()
    cohort_hypothesis = hypothesis.groupby(["scale", "axis_id", "ligand", "receptor", "source_proxy", "target_proxy"], dropna=False).agg(
        n_patients=("sample", "nunique"), median_proxy_log2_enrichment=("proxy_log2_enrichment", "median"),
        positive_enrichment_n=("proxy_log2_enrichment", lambda x: int((x > 0).sum())),
        top_decile_positive_n=("proxy_mass_rank_among_196", lambda x: int(((x <= 20) & (hypothesis.loc[x.index, "proxy_log2_enrichment"] > 0)).sum())),
        median_proxy_mass_rank=("proxy_mass_rank_among_196", "median"),
    ).reset_index()
    cohort_edge = edge.groupby(["scale", "axis_id", "source_celltype", "target_celltype"], dropna=False).agg(
        n_patients=("sample", "nunique"), median_proxy_log2_enrichment=("proxy_log2_enrichment", "median"),
        positive_enrichment_n=("proxy_log2_enrichment", lambda x: int((x > 0).sum())),
        top_decile_positive_n=("mass_rank_among_196", lambda x: int(((x <= 20) & (edge.loc[x.index, "proxy_log2_enrichment"] > 0)).sum())),
        median_mass_rank=("mass_rank_among_196", "median"),
    ).reset_index().sort_values(["top_decile_positive_n", "median_mass_rank"], ascending=[False, True])
    cohort_domain = domain.groupby(["scale", "axis_id", "source_region", "target_region"], dropna=False).agg(
        n_patients=("sample", "nunique"), median_edge_mass_fraction=("edge_mass_fraction", "median"),
        q25_edge_mass_fraction=("edge_mass_fraction", lambda x: x.quantile(.25)), q75_edge_mass_fraction=("edge_mass_fraction", lambda x: x.quantile(.75)),
    ).reset_index()
    availability = manifest.groupby("axis_id", as_index=False).agg(n_patients_available=("available", lambda x: int(as_bool(x).sum())), n_patients_total=("sample", "nunique"))
    availability["not_estimable_n"] = availability["n_patients_total"] - availability["n_patients_available"]
    availability.to_csv(output / "commot_axis_estimability_v1.tsv", sep="\t", index=False)
    cohort_transport.to_csv(output / "commot_patient_transport_summary_v1.tsv", sep="\t", index=False)
    cohort_hypothesis.to_csv(output / "commot_patient_hypothesis_proxy_summary_v1.tsv", sep="\t", index=False)
    cohort_edge.to_csv(output / "commot_patient_celltype_proxy_replication_v1.tsv", sep="\t", index=False)
    cohort_domain.to_csv(output / "commot_patient_domain_transport_v1.tsv", sep="\t", index=False)

    transport_plot = cohort_transport.pivot(index="axis_id", columns="scale", values="median_transport_rank_fraction")
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(transport_plot, cmap="mako_r", vmin=0, vmax=1, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Median within-patient COMMOT transport rank fraction"); ax.set_xlabel(""); ax.set_ylabel("")
    fig.tight_layout(); fig.savefig(output / "commot_patient_transport_rank_v1.png", dpi=300); fig.savefig(output / "commot_patient_transport_rank_v1.pdf"); plt.close(fig)
    hypothesis_plot = cohort_hypothesis.pivot(index="axis_id", columns="scale", values="median_proxy_log2_enrichment")
    fig, ax = plt.subplots(figsize=(9, 8))
    cap = max(float(np.nanpercentile(np.abs(hypothesis_plot), 95)), .25)
    sns.heatmap(hypothesis_plot, cmap="vlag", center=0, vmin=-cap, vmax=cap, annot=True, fmt=".2f", ax=ax)
    ax.set_title("Median prespecified source-target proxy enrichment"); ax.set_xlabel(""); ax.set_ylabel("")
    fig.tight_layout(); fig.savefig(output / "commot_patient_hypothesis_proxy_v1.png", dpi=300); fig.savefig(output / "commot_patient_hypothesis_proxy_v1.pdf"); plt.close(fig)
    return cohort_transport, cohort_hypothesis, availability, metadata


def main(output: Path) -> int:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    output.mkdir(parents=True)
    write_json(output / "status.json", {"status": "running", "stage": "patient_level_synthesis"})
    gates = [memory_gate("tjrcc_advanced_synthesis_start")]
    bayes, bayes_effects, bayes_meta = load_bayesspace(output)
    gates.append(memory_gate("tjrcc_advanced_synthesis_spatial_programs"))
    genes, overlay, spatial_meta = load_spatial_programs(output)
    gates.append(memory_gate("tjrcc_advanced_synthesis_misty"))
    misty_perf, misty_views, misty_importance, misty_meta = load_misty(output)
    gates.append(memory_gate("tjrcc_advanced_synthesis_commot"))
    transport, hypothesis, availability, commot_meta = load_commot(output)
    metadata = {
        "status": "completed", "n_patients": len(SAMPLES), "samples": list(SAMPLES),
        "n_spots": int(bayes["n_spots"].sum()),
        "bayesspace_selected_q_range": [int(bayes["selected_q"].min()), int(bayes["selected_q"].max())],
        "spatial_program_genes_tested_in_at_least_6_patients": int((genes["n_patients_tested"] >= 6).sum()),
        "misty_targets": int(misty_perf["target"].nunique()),
        "commot_axes_prespecified": int(availability.shape[0]),
        "commot_axes_estimable_all_patients": int((availability["n_patients_available"] == len(SAMPLES)).sum()),
        "source_metadata": {"bayesspace": bayes_meta, "spatial_programs": spatial_meta, "misty": misty_meta, "commot": commot_meta},
        "memory_gate_records": gates,
        "inference_boundary": "Patients are the independent biological units for cohort summaries. Section/spot-level tests and computational boundaries are descriptive and do not make spots independent replicates. COMMOT unavailable axes are not estimable, not biological negatives.",
        "method_scope": "BayesSpace/computational boundary; genome-wide SPARK-X and Hotspot; locked-feature MISTy; targeted 12-axis COMMOT.",
    }
    write_json(output / "run_metadata_v1.json", metadata)
    write_json(output / "status.json", {"status": "completed", "stage": "patient_level_synthesis"})
    print(f"TJRCC_ADVANCED_COHORT_SYNTHESIS_PASS patients={len(SAMPLES)} spots={metadata['n_spots']}", flush=True)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: summarize_tjrcc_advanced_cohort_v1.py OUTPUT_DIR")
    target = Path(sys.argv[1]).resolve()
    existed = target.exists()
    try:
        raise SystemExit(main(target))
    except Exception as exc:
        if not existed:
            target.mkdir(parents=True, exist_ok=True)
            write_json(target / "status.json", {"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
            (target / "failure_traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise
