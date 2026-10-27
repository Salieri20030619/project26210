from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


NAVY = "#23364D"
TEAL = "#0D8C82"
MAGENTA = "#B72C63"
GOLD = "#E99A2E"
MUTED = "#66717C"
LIGHT = "#F2F4F7"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(fig, path: Path):
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def heatmap(ax, matrix, row_labels, col_labels, cmap, vmin=None, vmax=None, fmt=".2f"):
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(col_labels)), col_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(row_labels)), row_labels, fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            if np.isfinite(val):
                span = max(abs(vmin or np.nanmin(matrix)), abs(vmax or np.nanmax(matrix)), 1e-9)
                color = "white" if abs(val) > 0.58 * span else NAVY
                ax.text(j, i, format(val, fmt), ha="center", va="center", fontsize=7.2, color=color)
    return im


def make_pseudobulk(source: Path, output: Path):
    df = pd.read_csv(source, sep="\t")
    pathways = [
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "HALLMARK_HYPOXIA",
        "HALLMARK_MTORC1_SIGNALING",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_GLYCOLYSIS",
        "HALLMARK_COMPLEMENT",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_ALLOGRAFT_REJECTION",
        "HALLMARK_ANGIOGENESIS",
        "HALLMARK_INFLAMMATORY_RESPONSE",
    ]
    row_names = {
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNFA / NFKB",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-gamma response",
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation",
        "HALLMARK_HYPOXIA": "Hypoxia",
        "HALLMARK_MTORC1_SIGNALING": "mTORC1 signaling",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
        "HALLMARK_GLYCOLYSIS": "Glycolysis",
        "HALLMARK_COMPLEMENT": "Complement",
        "HALLMARK_IL2_STAT5_SIGNALING": "IL2 / STAT5",
        "HALLMARK_ALLOGRAFT_REJECTION": "Allograft rejection",
        "HALLMARK_ANGIOGENESIS": "Angiogenesis",
        "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    }
    model_order = [
        ("GSE264586", "Cycle2_vs_Baseline", "CD4_naive_memory", "C2-BL CD4"),
        ("GSE264586", "Cycle2_vs_Baseline", "NK_cell", "C2-BL NK"),
        ("GSE264586", "End_of_treatment_vs_Baseline", "CD4_naive_memory", "EOT-BL CD4"),
        ("GSE264586", "End_of_treatment_vs_Baseline", "Macrophage_C1QC", "EOT-BL C1QC"),
        ("GSE264586", "End_of_treatment_vs_Baseline", "Monocyte_FCN1", "EOT-BL FCN1"),
        ("GSE264586", "End_of_treatment_vs_Baseline", "NK_cell", "EOT-BL NK"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "CD4_naive_memory", "Benefit CD4"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "CD8_cytotoxic", "Benefit CD8"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "Macrophage_C1QC", "Benefit C1QC"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "Monocyte_FCN1", "Benefit FCN1"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "NK_cell", "Benefit NK"),
        ("GSE294109", "ICB_benefit_vs_no_effect", "Treg", "Benefit Treg"),
    ]
    matrix = np.full((len(pathways), len(model_order)), np.nan)
    qmat = np.full_like(matrix, np.nan)
    for j, (dataset, contrast, annotation, _) in enumerate(model_order):
        sub = df[(df.dataset == dataset) & (df.contrast == contrast) & (df.fine_annotation == annotation)].set_index("pathway")
        for i, pathway in enumerate(pathways):
            if pathway in sub.index:
                row = sub.loc[pathway]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                matrix[i, j] = float(row.NES)
                qmat[i, j] = float(row.padj)
    vmax = np.nanmax(np.abs(matrix))
    cmap = LinearSegmentedColormap.from_list("nes", [TEAL, "#F7F8FA", MAGENTA])
    fig, ax = plt.subplots(figsize=(12.8, 6.7))
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(model_order)), [x[3] for x in model_order], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(pathways)), [row_names[x] for x in pathways], fontsize=8.5)
    ax.set_ylim(len(pathways) - 0.5, -1.4)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                if qmat[i, j] < 0.05:
                    ax.text(j, i, "•", ha="center", va="center", fontsize=13, color="white" if abs(matrix[i, j]) > 0.55 * vmax else NAVY, fontweight="bold")
    ax.axvline(5.5, color=NAVY, lw=1.1)
    ax.text(2.5, -0.82, "GSE264586 longitudinal", ha="center", va="bottom", fontsize=9, color=NAVY, fontweight="bold")
    ax.text(8.5, -0.82, "GSE294109 benefit vs no-effect", ha="center", va="bottom", fontsize=9, color=NAVY, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("NES (negative ←  → positive)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    ax.set_title("Patient-level immune pathway remodeling", loc="left", fontsize=13, color=NAVY, fontweight="bold", pad=18)
    ax.text(0, 1.01, "Dot = within-model FDR < 0.05", transform=ax.transAxes, fontsize=8, color=MUTED, va="bottom")
    fig.tight_layout()
    save(fig, output)


def make_forest(source: Path, output: Path):
    df = pd.read_csv(source, sep="\t")
    labels = [
        "E-MTAB PFS", "E-MTAB ORR", "CheckMate PFS interaction", "CheckMate nivolumab ORR",
        "CheckMate nivolumab PFS", "JAVELIN combo PFS", "JAVELIN PFS interaction",
    ]
    y = np.arange(len(df))[::-1]
    fig, ax = plt.subplots(figsize=(10.5, 5.3))
    colors = [TEAL if c == "E-MTAB-3218" else MAGENTA if "CheckMate" in c else GOLD for c in df.cohort]
    for yi, (_, row), color in zip(y, df.iterrows(), colors):
        ax.plot([row.ci_low, row.ci_high], [yi, yi], color=color, lw=2)
        ax.scatter(row.effect, yi, s=55, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        ax.text(max(df.ci_high) + 0.12, yi, f"P={row.p_value:.3g}  q={row.q_value:.3g}", va="center", fontsize=8, color=MUTED)
    ax.axvline(1, color=NAVY, lw=1.1)
    ax.set_yticks(y, labels, fontsize=8.5)
    ax.set_xscale("log")
    ax.set_xlim(0.28, 3.9)
    ax.set_xticks([0.3, 0.5, 1, 2, 3], ["0.3", "0.5", "1", "2", "3"])
    ax.set_xlabel("Hazard / odds / interaction ratio per SD", fontsize=9)
    ax.grid(axis="x", color="#DDE2E7", lw=0.7)
    ax.tick_params(axis="y", length=0)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.set_title("Independent clinical validation of the replicated Myeloid/DC core", loc="left", fontsize=13, color=NAVY, fontweight="bold", pad=14)
    ax.text(0.01, 1.01, "All seven preregistered models are FDR-negative", transform=ax.transAxes, fontsize=8.5, color=MUTED, va="bottom")
    fig.subplots_adjust(left=0.28, right=0.78, top=0.84, bottom=0.18)
    save(fig, output)


def clean_term(term: str) -> str:
    replacements = {
        "MHC protein complex": "MHC protein complex",
        "MHC class II protein complex": "MHC-II protein complex",
        "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES": "MHC-II antigen presentation",
        "antigen binding": "Antigen binding",
        "antigen processing and presentation of exogenous peptide antigen via MHC class II": "Exogenous antigen presentation",
        "complement activation": "Complement activation",
        "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES": "Complement / coagulation",
        "external encapsulating structure": "External matrix structure",
        "extracellular matrix": "Extracellular matrix",
        "humoral immune response": "Humoral immune response",
        "humoral immune response mediated by circulating immunoglobulin": "Circulating Ig humoral response",
    }
    return replacements.get(term, term[:54])


def make_ko(source: Path, output: Path):
    df = pd.read_csv(source, sep="\t")
    chosen = []
    for gene in ["TYMP", "MTAP"]:
        sub = df[df.knockout_gene == gene].head(5).copy()
        chosen.append(sub)
    plot = pd.concat(chosen, ignore_index=True)
    labels = [f"{g}: {clean_term(t)}" for g, t in zip(plot.knockout_gene, plot.term_name)]
    q1 = -np.log10(plot.GSE264586_q_value.astype(float).to_numpy())
    q2 = -np.log10(plot.GSE294109_q_value.astype(float).to_numpy())
    y = np.arange(len(plot))[::-1]
    fig, ax = plt.subplots(figsize=(11.4, 6.1))
    for yi, a, b in zip(y, q1, q2):
        ax.plot([a, b], [yi, yi], color="#CBD2D9", lw=2)
    ax.scatter(q1, y, color=TEAL, s=55, label="GSE264586", zorder=3)
    ax.scatter(q2, y, color=MAGENTA, marker="D", s=48, label="GSE294109", zorder=3)
    ax.axvline(-math.log10(0.05), color=NAVY, lw=1, ls="--")
    ax.set_yticks(y, labels, fontsize=8.2)
    ax.set_xlabel("-log10(FDR q)", fontsize=9)
    ax.grid(axis="x", color="#E1E5E9", lw=0.7)
    ax.tick_params(axis="y", length=0)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("Cross-cohort network-perturbation enrichment", loc="left", fontsize=13, color=NAVY, fontweight="bold", pad=14)
    ax.text(0.01, 1.01, "TYMP KO → antigen presentation; MTAP KO → complement / matrix / secretory programs", transform=ax.transAxes, fontsize=8.5, color=MUTED, va="bottom")
    fig.subplots_adjust(left=0.40, right=0.97, top=0.84, bottom=0.14)
    save(fig, output)


def make_identifiability(source: Path, output: Path):
    df = pd.read_csv(source, sep="\t")
    df = df[df.tier == "tissue_context"].copy()
    method_order = ["frozen13", "frozen13_anchor_excluded", "lodo13", "lodo3_anchor_excluded"]
    dataset_order = ["GSE207493", "GSE264586", "GSE294109"]
    df["method"] = pd.Categorical(df.method, method_order, ordered=True)
    df["dataset"] = pd.Categorical(df.dataset, dataset_order, ordered=True)
    df = df.sort_values(["method", "dataset"])
    rows = [f"{m.replace('_anchor_excluded','-no anchor')} | {d}" for m, d in zip(df.method.astype(str), df.dataset.astype(str))]
    recovery = df[["minimum_state_rho", "dominant_state_accuracy", "total_macrophage_rho"]].to_numpy(float)
    failure = df[["maximum_state_mae", "maximum_state_zero_fraction"]].to_numpy(float)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 6.7), gridspec_kw={"width_ratios": [1.45, 1]})
    cmap_rec = LinearSegmentedColormap.from_list("rec", ["#F7E8EF", "#F7F8FA", TEAL])
    cmap_fail = LinearSegmentedColormap.from_list("fail", [TEAL, "#F7F8FA", MAGENTA])
    heatmap(axes[0], recovery, rows, ["Min state rho", "Dominant accuracy", "Total macro rho"], cmap_rec, vmin=0, vmax=1)
    heatmap(axes[1], failure, rows, ["Max state MAE", "Max zero fraction"], cmap_fail, vmin=0, vmax=1)
    axes[1].set_yticklabels([])
    axes[0].set_title("Recovery metrics (higher is better)", loc="left", fontsize=11, color=NAVY, fontweight="bold")
    axes[1].set_title("Failure metrics (lower is better)", loc="left", fontsize=11, color=NAVY, fontweight="bold")
    fig.suptitle("Bulk fine-state identifiability benchmark", x=0.04, ha="left", fontsize=13, color=NAVY, fontweight="bold")
    fig.text(0.04, 0.925, "Fine-state recovery fails while coarse total-macrophage recovery remains usable", fontsize=8.5, color=MUTED)
    fig.subplots_adjust(left=0.28, right=0.98, top=0.84, bottom=0.12, wspace=0.18)
    save(fig, output)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 03_make_representative_figures_v1.py <project_root> <output_dir>")
    root = Path(sys.argv[1])
    output = Path(sys.argv[2])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")
    output.mkdir(parents=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Microsoft YaHei", "Arial", "DejaVu Sans"], "axes.unicode_minus": False})
    targets = {
        "immune_pathway_remodeling_heatmap.png": lambda p: make_pseudobulk(root / "phase62_part1_immune_pseudobulk_local_v1" / "results_v1" / "pseudobulk_hallmark_fgsea_v1.tsv", p),
        "clinical_validation_forest.png": lambda p: make_forest(root / "phase69_85_server_status_sync_local_v1" / "phase69" / "results_v1" / "R2_COMPOSITE_PRIMARY_MODELS_v1.tsv", p),
        "virtual_knockout_enrichment.png": lambda p: make_ko(root / "phase69_85_server_status_sync_local_v1" / "phase74" / "results_v1" / "phase74_top_two_external_terms_v1.tsv", p),
        "bulk_identifiability_diagnostic.png": lambda p: make_identifiability(root / "phase86_macrophage_state_identifiability_local_v1" / "results_v1" / "phase86_gate_decisions_v1.tsv", p),
    }
    records = []
    for name, fn in targets.items():
        target = output / name
        fn(target)
        records.append({"file": name, "bytes": target.stat().st_size, "sha256": sha256(target)})
    (output / "REPRESENTATIVE_FIGURE_MANIFEST_v1.json").write_text(json.dumps({"status": "PASS", "figures": records}, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "figures": len(records)}, indent=2))


if __name__ == "__main__":
    main()

