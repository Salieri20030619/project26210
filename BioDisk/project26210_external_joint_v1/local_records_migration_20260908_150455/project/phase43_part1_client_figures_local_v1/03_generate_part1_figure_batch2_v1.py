from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from PIL import Image


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "source_remote_v1"
OUT = ROOT / "figures_v2"
QC = ROOT / "qc_v2"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
)


def panel_letter(ax, letter: str) -> None:
    ax.text(-0.05, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", ha="right")


def clean_axes(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#E1E1E1", linewidth=0.5, alpha=0.75)
    ax.set_axisbelow(True)


def save(fig, stem: str, width: float, height: float) -> list[Path]:
    fig.set_size_inches(width, height)
    outputs = []
    for extension, options in [("png", {"dpi": 300}), ("pdf", {})]:
        path = OUT / f"{stem}.{extension}"
        fig.savefig(path, bbox_inches="tight", facecolor="white", **options)
        outputs.append(path)
    plt.close(fig)
    return outputs


all_outputs: list[Path] = []

# Figure 4: cNMF rank choice, TYMP-associated program, and functional annotation.
rank_qc = pd.read_csv(SRC / "cnmf_rank_qc_v1.tsv", sep="\t")
factor = pd.read_csv(SRC / "cnmf_factor_evidence_v1.tsv", sep="\t")
rank7 = factor.loc[factor["rank"] == 7].sort_values("factor").copy()
ora_img = Image.open(SRC / "rank7_program_pathway_ORA_v1.png").convert("RGB")
de_img = Image.open(SRC / "F1_high_clone_top_DE_patient_heatmap_v1.png").convert("RGB")

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.28])
ax_a = fig.add_subplot(gs[0, 0])
ax_a.plot(rank_qc["rank"], rank_qc["mean_component_silhouette"], marker="o", color="#2B6F8E", label="Mean silhouette")
ax_a.plot(rank_qc["rank"], rank_qc["min_run_fraction"], marker="s", color="#D77A44", label="Minimum run coverage")
ax_a.axvspan(6.86, 7.14, color="#F1C75B", alpha=0.32, linewidth=0)
for row in rank_qc.itertuples(index=False):
    if bool(row.stability_gate):
        ax_a.scatter(row.rank, 0.05, marker="o", s=25, color="#278B6C", zorder=4)
    else:
        ax_a.scatter(row.rank, 0.05, marker="x", s=30, color="#A94442", zorder=4)
ax_a.set_xticks(rank_qc["rank"])
ax_a.set_ylim(0, 1.02)
ax_a.set_xlabel("cNMF rank")
ax_a.set_ylabel("Stability metric")
ax_a.legend(frameon=False, loc="lower left", fontsize=7.2)
clean_axes(ax_a)
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
colors = ["#D55E00" if value else "#7EA3B6" for value in rank7["selected_TYMP_candidate_within_rank"]]
ax_b.bar(rank7["factor"], rank7["median_patient_rho_TYMP"], color=colors, edgecolor="none")
ax_b.errorbar(
    np.arange(len(rank7)),
    rank7["median_patient_rho_TYMP"],
    yerr=np.vstack(
        [
            rank7["median_patient_rho_TYMP"] - rank7["q25_patient_rho_TYMP"],
            rank7["q75_patient_rho_TYMP"] - rank7["median_patient_rho_TYMP"],
        ]
    ),
    fmt="none",
    ecolor="#333333",
    elinewidth=0.9,
    capsize=2.5,
)
ax_b.axhline(0, color="#666666", linewidth=0.7)
ax_b.set_xlabel("Rank-7 program")
ax_b.set_ylabel("Patient-level median Spearman ρ with TYMP")
clean_axes(ax_b)
panel_letter(ax_b, "B")

ax_c = fig.add_subplot(gs[1, 0])
ax_c.imshow(ora_img)
ax_c.axis("off")
panel_letter(ax_c, "C")
ax_d = fig.add_subplot(gs[1, 1])
ax_d.imshow(de_img)
ax_d.axis("off")
panel_letter(ax_d, "D")
all_outputs += save(fig, "Fig04_cNMF_rank7_TYMP_program_v1", 12.0, 9.0)

# Figure 5: scDRS map, descriptive cell-level relation, and patient-level robustness.
scdrs_path = ROOT.parent / "phase25_part1_scdrs_rebuild_local_v1" / "server_results_sync_v1" / "analysis_v1"
scdrs = pd.read_csv(scdrs_path / "scdrs_factor_tymp_cell_source_v1.tsv.gz", sep="\t")
patient_assoc = pd.read_csv(scdrs_path / "patient_within_association_v1.tsv", sep="\t")
eligible = patient_assoc.loc[patient_assoc["eligible_primary"] == True].copy()
study_colors = {
    "GSE156632": "#0072B2",
    "GSE159115": "#E69F00",
    "GSE178481": "#009E73",
    "GSE207493": "#CC79A7",
    "GSE242299": "#6A51A3",
}

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 0.9])
ax_a = fig.add_subplot(gs[0, 0], facecolor="#080A0C")
order = np.argsort(scdrs["scdrs_norm_score"].to_numpy())
plot_a = ax_a.scatter(
    scdrs.loc[order, "UMAP1"],
    scdrs.loc[order, "UMAP2"],
    c=scdrs.loc[order, "scdrs_norm_score"],
    cmap="inferno",
    s=2.1,
    linewidths=0,
    rasterized=True,
)
ax_a.set_xlabel("UMAP1", color="white")
ax_a.set_ylabel("UMAP2", color="white")
ax_a.tick_params(colors="white")
for spine in ax_a.spines.values():
    spine.set_color("white")
cb = fig.colorbar(plot_a, ax=ax_a, fraction=0.04, pad=0.02)
cb.set_label("RCC scDRS normalized score")
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
hexbin = ax_b.hexbin(
    scdrs["F1"],
    scdrs["scdrs_norm_score"],
    gridsize=45,
    mincnt=1,
    cmap="Blues",
    bins="log",
)
fig.colorbar(hexbin, ax=ax_b, fraction=0.045, pad=0.02, label="log10(cell density)")
ax_b.text(
    0.04,
    0.96,
    "Pooled-cell Pearson r = 0.009\n(descriptive)",
    transform=ax_b.transAxes,
    ha="left",
    va="top",
    fontsize=8,
)
ax_b.set_xlabel("Factor-TYMP program usage (F1)")
ax_b.set_ylabel("RCC scDRS normalized score")
clean_axes(ax_b)
panel_letter(ax_b, "B")

ax_c = fig.add_subplot(gs[0, 2])
eligible = eligible.sort_values("spearman_rho_f1_vs_scdrs")
y = np.arange(len(eligible))
for study, subset in eligible.groupby("study"):
    idx = [eligible.index.get_loc(value) for value in subset.index]
    ax_c.scatter(
        subset["spearman_rho_f1_vs_scdrs"],
        idx,
        s=38,
        color=study_colors[study],
        edgecolor="white",
        linewidth=0.5,
        label=study,
    )
ax_c.axvline(0, color="#666666", linewidth=0.8)
ax_c.axvline(-0.0490078, color="#D55E00", linewidth=1.4)
ax_c.axvspan(-0.1291551, 0.0635501, color="#D55E00", alpha=0.13, linewidth=0)
ax_c.set_yticks(y, [value.replace("__", " · ") for value in eligible["curated_patient_id"]], fontsize=6.2)
ax_c.set_xlabel("Within-patient Spearman ρ")
ax_c.set_ylabel("")
ax_c.legend(frameon=False, fontsize=6.2, loc="lower right")
clean_axes(ax_c)
panel_letter(ax_c, "C")
all_outputs += save(fig, "Fig05_scDRS_Factor_TYMP_v1", 12.4, 5.2)

# Figure 6: pySCENIC dual-method/dual-endpoint driver evidence.
pyscenic_path = ROOT.parent / "phase26_part1_cnmf_pyscenic_local_v1" / "server_results_sync_v2" / "analysis_v3"
drivers = pd.read_csv(pyscenic_path / "f1_driver_consensus_v1.tsv", sep="\t")
top = drivers.head(15).copy()
endpoints = [
    ("GENIE3\ncontinuous", "median_patient_effect__GENIE3__continuous"),
    ("GRNBoost2\ncontinuous", "median_patient_effect__GRNBoost2__continuous"),
    ("GENIE3\nF1 group", "median_patient_effect_internal__GENIE3__f1_group"),
    ("GRNBoost2\nF1 group", "median_patient_effect_internal__GRNBoost2__f1_group"),
]
matrix = top[[column for _, column in endpoints]].apply(pd.to_numeric, errors="coerce").to_numpy()
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 0.72])
ax_a = fig.add_subplot(gs[0, 0])
limit = np.nanmax(np.abs(matrix))
image_a = ax_a.imshow(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), aspect="auto")
ax_a.set_xticks(np.arange(4), [label for label, _ in endpoints])
ax_a.set_yticks(np.arange(len(top)), top["tf"])
for row in range(matrix.shape[0]):
    for col in range(matrix.shape[1]):
        value = matrix[row, col]
        if np.isfinite(value):
            ax_a.text(col, row, f"{value:.3f}", ha="center", va="center", fontsize=5.7)
fig.colorbar(image_a, ax=ax_a, fraction=0.035, pad=0.02, label="Patient-level effect")
ax_a.set_xlabel("")
ax_a.set_ylabel("")
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
top_bar = top.sort_values("n_endpoints_passed", ascending=True)
bar_colors = ["#D55E00" if value >= 3 else "#7EA3B6" for value in top_bar["n_endpoints_passed"]]
ax_b.barh(top_bar["tf"], top_bar["n_endpoints_passed"], color=bar_colors, edgecolor="none")
ax_b.set_xlim(0, 4.2)
ax_b.set_xticks([0, 1, 2, 3, 4])
ax_b.set_xlabel("Evidence endpoints passed (of 4)")
ax_b.set_ylabel("")
clean_axes(ax_b)
ax_b.grid(axis="x", color="#E1E1E1", linewidth=0.5)
ax_b.grid(axis="y", visible=False)
panel_letter(ax_b, "B")
all_outputs += save(fig, "Fig06_pySCENIC_F1_driver_consensus_v1", 10.2, 6.4)

# Figure 7: MR-metaVIPER direction and F1-associated activity effects.
mv_path = ROOT.parent / "phase29_part1_metaviper_downstream_local_v1" / "server_results_sync_v1"
mr = pd.read_csv(mv_path / "m11_analysis_v1" / "m11_mr_effect_activity_points_v1.tsv", sep="\t")
effects = pd.read_csv(mv_path / "m12_analysis_v1" / "m12_patient_activity_effects_v1.tsv", sep="\t")
summary = pd.read_csv(mv_path / "m12_analysis_v1" / "m12_f1_activity_patient_summary_v1.tsv", sep="\t")
proteins = ["CASP9", "CHMP1A", "PRSS8", "PXN"]
endpoint_labels = {
    "continuous_F1_spearman": "Continuous F1",
    "F1_dominant_median_difference": "F1-dominant",
    "F1_Q4_minus_Q1_median_difference": "F1 Q4−Q1",
}

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 1.18, 1.15])
ax_a = fig.add_subplot(gs[0, 0])
for row in mr.itertuples(index=False):
    color = "#D55E00" if row.mr_tier == "Tier 1" else "#2B6F8E"
    ax_a.scatter(row.mr_beta_log_or, row.patient_equal_mean_activity, s=45 + 22 * -np.log10(row.mr_fdr), color=color, edgecolor="white", linewidth=0.7)
    ax_a.text(row.mr_beta_log_or, row.patient_equal_mean_activity, f"  {row.protein}", va="center", fontsize=7.2)
ax_a.axhline(0, color="#777777", linewidth=0.7)
ax_a.axvline(0, color="#777777", linewidth=0.7)
ax_a.set_xlabel("MR effect, log(odds ratio)")
ax_a.set_ylabel("Patient-equal mean metaVIPER activity")
clean_axes(ax_a)
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
matrix_b = summary.loc[summary["protein"].isin(proteins)].pivot(index="protein", columns="endpoint", values="median_patient_effect").reindex(index=proteins, columns=endpoint_labels.keys())
limit_b = np.nanmax(np.abs(matrix_b.to_numpy()))
image_b = ax_b.imshow(matrix_b, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit_b, vcenter=0, vmax=limit_b), aspect="auto")
ax_b.set_xticks(np.arange(3), [endpoint_labels[value] for value in matrix_b.columns], rotation=42, ha="right")
ax_b.set_yticks(np.arange(len(proteins)), proteins)
for i, protein in enumerate(proteins):
    for j, endpoint in enumerate(matrix_b.columns):
        value = matrix_b.loc[protein, endpoint]
        row = summary.loc[(summary["protein"] == protein) & (summary["endpoint"] == endpoint)]
        star = "*" if len(row) and float(row.iloc[0]["BH_FDR"]) < 0.05 else ""
        ax_b.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=7)
fig.colorbar(image_b, ax=ax_b, fraction=0.045, pad=0.02, label="Median patient effect")
panel_letter(ax_b, "B")

ax_c = fig.add_subplot(gs[0, 2])
endpoint = "F1_Q4_minus_Q1_median_difference"
subset = effects.loc[(effects["protein"].isin(proteins)) & (effects["endpoint"] == endpoint)].copy()
rng = np.random.default_rng(26210)
for i, protein in enumerate(proteins):
    values = subset.loc[subset["protein"] == protein, "effect"].dropna().to_numpy()
    violin = ax_c.violinplot([values], positions=[i], widths=0.72, showextrema=False)
    violin["bodies"][0].set_facecolor("#7EA3B6" if protein == "PXN" else "#D97A4A")
    violin["bodies"][0].set_edgecolor("none")
    violin["bodies"][0].set_alpha(0.42)
    ax_c.scatter(rng.normal(i, 0.045, len(values)), values, s=18, color="#3F3F3F", alpha=0.65, linewidths=0)
    ax_c.plot([i - 0.2, i + 0.2], [np.median(values), np.median(values)], color="#222222", linewidth=1.4)
ax_c.axhline(0, color="#666666", linewidth=0.7, linestyle="--")
ax_c.set_xticks(np.arange(len(proteins)), proteins, rotation=35, ha="right")
ax_c.set_ylabel("Patient activity effect, F1 Q4−Q1")
clean_axes(ax_c)
panel_letter(ax_c, "C")
all_outputs += save(fig, "Fig07_metaVIPER_MR_and_F1_activity_v1", 12.0, 5.0)

# Figure 8: requested enzyme/transporter synergy matrix and F1 endpoint matrix.
m17_path = ROOT.parent / "phase31_part1_m17_enzyme_matrix_local_v1" / "results_local_v1"
pairs = pd.read_csv(m17_path / "m17_gene_pair_correlation_summary_v1.tsv", sep="\t")
gene_effects = pd.read_csv(m17_path / "m17_gene_f1_summary_v1.tsv", sep="\t")
genes = ["TYMP", "TK1", "SLC29A1", "SLC29A2", "ODC1", "SRM", "SMS", "MTAP", "CHKA", "PCYT1A", "PCYT2", "CEPT1"]
corr = pd.DataFrame(np.nan, index=genes, columns=genes)
fdr = pd.DataFrame(np.nan, index=genes, columns=genes)
np.fill_diagonal(corr.values, 1.0)
np.fill_diagonal(fdr.values, 0.0)
for row in pairs.itertuples(index=False):
    if row.gene_a in genes and row.gene_b in genes:
        corr.loc[row.gene_a, row.gene_b] = row.median_patient_effect
        corr.loc[row.gene_b, row.gene_a] = row.median_patient_effect
        fdr.loc[row.gene_a, row.gene_b] = row.BH_FDR
        fdr.loc[row.gene_b, row.gene_a] = row.BH_FDR

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0])
ax_a = fig.add_subplot(gs[0, 0])
mask = np.triu(np.ones_like(corr, dtype=bool), 1)
matrix_a = np.ma.array(corr.to_numpy(), mask=mask)
image_a = ax_a.imshow(matrix_a, cmap="RdBu_r", vmin=-0.18, vmax=0.18)
ax_a.set_xticks(np.arange(len(genes)), genes, rotation=52, ha="right")
ax_a.set_yticks(np.arange(len(genes)), genes)
for i in range(len(genes)):
    for j in range(i + 1):
        if np.isfinite(corr.iloc[i, j]):
            label = f"{corr.iloc[i, j]:.2f}" if i != j else ""
            if i != j and fdr.iloc[i, j] < 0.05:
                label += "*"
            ax_a.text(j, i, label, ha="center", va="center", fontsize=5.7)
fig.colorbar(image_a, ax=ax_a, fraction=0.04, pad=0.02, label="Median within-patient Spearman ρ")
panel_letter(ax_a, "A")

ax_b = fig.add_subplot(gs[0, 1])
endpoint_order = ["F1_continuous", "F1_patient_quartile", "F1_dominant_group"]
effect_matrix = gene_effects.pivot(index="gene", columns="endpoint", values="median_patient_effect").reindex(index=genes, columns=endpoint_order)
limit_e = max(0.15, np.nanmax(np.abs(effect_matrix.to_numpy())))
image_b = ax_b.imshow(effect_matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit_e, vcenter=0, vmax=limit_e), aspect="auto")
ax_b.set_xticks(np.arange(3), ["Continuous F1", "F1 Q4−Q1", "F1-dominant"], rotation=38, ha="right")
ax_b.set_yticks(np.arange(len(genes)), genes)
for i, gene in enumerate(genes):
    for j, endpoint_name in enumerate(endpoint_order):
        value = effect_matrix.loc[gene, endpoint_name]
        row = gene_effects.loc[(gene_effects["gene"] == gene) & (gene_effects["endpoint"] == endpoint_name)]
        if np.isfinite(value):
            star = "*" if len(row) and bool(row.iloc[0]["statistical_support"]) else ""
            ax_b.text(j, i, f"{value:.2f}{star}", ha="center", va="center", fontsize=6)
fig.colorbar(image_b, ax=ax_b, fraction=0.05, pad=0.02, label="Median patient effect")
panel_letter(ax_b, "B")
all_outputs += save(fig, "Fig08_enzyme_transporter_synergy_matrix_v1", 11.0, 6.8)

# Figure 9: five-method malignant-cell trajectory panel in the requested small-multiple format.
trajectory_files = [
    ("Monocle3", "Monocle3_pseudotime_v1.png"),
    ("Slingshot", "Slingshot_pseudotime_v2.png"),
    ("CytoTRACE2", "CytoTRACE2_potency_v1.png"),
    ("scTour", "scTour_pseudotime_v1.png"),
    ("CellRank", "CellRank_fate_TYMP_state_v1.png"),
]
fig, axes = plt.subplots(2, 3, constrained_layout=True)
for ax, (method, filename), letter in zip(axes.flat, trajectory_files, ["A", "B", "C", "D", "E"]):
    ax.imshow(Image.open(SRC / filename).convert("RGB"))
    ax.axis("off")
    ax.set_title(method, pad=3)
    panel_letter(ax, letter)
axes.flat[-1].axis("off")
all_outputs += save(fig, "Fig09_malignant_trajectory_five_methods_v1", 12.0, 8.0)

# Hash inventory for this batch only.
rows = []
for path in all_outputs:
    rows.append(
        {
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    )
pd.DataFrame(rows).to_csv(QC / "batch2_figure_manifest_v1.tsv", sep="\t", index=False)
