from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
PHASE = Path(__file__).resolve().parent
SRC = PHASE / "source_remote_v1"
OUT = PHASE / "figures_new_v1"
OUT.mkdir(parents=True, exist_ok=True)

PH35 = ROOT / "phase35_part1_c01_c04_cell_communication_local_v1" / "server_results_sync_v1"
PH36 = ROOT / "phase36_part1_c05_tcga_joint_survival_local_v1" / "server_results_sync_v1" / "analysis_v1"
PH25 = ROOT / "phase25_part1_scdrs_rebuild_local_v1" / "server_results_sync_v1" / "analysis_v1"
PH30 = ROOT / "phase30_part1_enzyme_function_axis_local_v1"

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.8,
})

BLUE = "#3E6FA8"
ORANGE = "#D9793D"
GREEN = "#4E8D73"
PURPLE = "#7D62A3"
RED = "#B84A4A"
GRAY = "#8B9098"
PALETTE = [BLUE, ORANGE, GREEN, PURPLE, RED, "#B69345"]


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel(ax, letter: str) -> None:
    ax.text(-0.09, 1.06, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def clean(ax, grid: bool = False) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#DADDE1", linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)


def black_umap(ax, x, y, value, title: str, cmap="inferno", s=1.2):
    ax.set_facecolor("#050505")
    finite = np.isfinite(value)
    if finite.any():
        lo, hi = np.nanpercentile(value[finite], [1, 99])
        if lo == hi:
            lo, hi = np.nanmin(value[finite]), np.nanmax(value[finite]) + 1e-9
    else:
        lo, hi = 0, 1
    ax.scatter(x[~finite], y[~finite], s=s, c="#292929", linewidths=0, rasterized=True)
    sc = ax.scatter(x[finite], y[finite], c=value[finite], s=s, cmap=cmap, vmin=lo, vmax=hi,
                    linewidths=0, rasterized=True)
    ax.set_title(title, pad=4)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    return sc


def fig10_regulon_networks():
    df = pd.read_csv(SRC / "pyscenic_selected_regulon_targets_v1.tsv", sep="\t")
    tfs = ["ILF2", "HOXD3", "HNF4A"]
    methods = ["GENIE3", "GRNBoost2"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5))
    for k, (ax, tf) in enumerate(zip(axes, tfs)):
        sub = df[(df.tf == tf) & (df.target != tf)].copy()
        ranks = []
        for method in methods:
            sm = sub[sub.method == method].nlargest(12, "weight").copy()
            sm["rank"] = np.arange(1, len(sm) + 1)
            ranks.append(sm)
        top = pd.concat(ranks, ignore_index=True)
        targets = list(dict.fromkeys(top.target.tolist()))[:20]
        n = len(targets)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
        pos = {g: (math.cos(a), math.sin(a)) for g, a in zip(targets, angles)}
        pos[tf] = (0.0, 0.0)
        maxw = max(top.weight.max(), 1e-12)
        for method, color, rad in [("GENIE3", BLUE, 0.12), ("GRNBoost2", ORANGE, -0.12)]:
            sm = top[(top.method == method) & (top.target.isin(targets))]
            for row in sm.itertuples():
                x, y = pos[row.target]
                patch = FancyArrowPatch((0, 0), (0.87*x, 0.87*y), arrowstyle="-|>", mutation_scale=7,
                                        linewidth=0.35 + 2.2 * row.weight / maxw, color=color,
                                        alpha=0.55, connectionstyle=f"arc3,rad={rad}")
                ax.add_patch(patch)
        shared = set(top[top.method == "GENIE3"].target) & set(top[top.method == "GRNBoost2"].target)
        for g in targets:
            x, y = pos[g]
            ax.scatter([x], [y], s=65 if g in shared else 40, c=GREEN if g in shared else "#D9DDE2",
                       edgecolors="#52575D", linewidths=0.5, zorder=4)
            ha = "left" if x > 0.2 else "right" if x < -0.2 else "center"
            va = "bottom" if y > 0.65 else "top" if y < -0.65 else "center"
            ax.text(1.10*x, 1.10*y, g, ha=ha, va=va, fontsize=7.2)
        ax.scatter([0], [0], s=850, c=PURPLE, edgecolors="white", linewidths=1.2, zorder=5)
        ax.text(0, 0, tf, ha="center", va="center", color="white", fontweight="bold", fontsize=8.2, zorder=6)
        ax.set_title(f"{tf} regulon")
        ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.axis("off")
        panel(ax, chr(65+k))
    handles = [mpl.lines.Line2D([0],[0], color=BLUE, lw=2, label="GENIE3"),
               mpl.lines.Line2D([0],[0], color=ORANGE, lw=2, label="GRNBoost2"),
               mpl.lines.Line2D([0],[0], marker="o", color="none", markerfacecolor=GREEN,
                                markeredgecolor="#52575D", label="shared top target")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Independent pySCENIC networks converge on selected F1-associated regulators", y=1.02, fontsize=13)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    save(fig, "Fig10_pySCENIC_selected_regulon_networks_v1")


def fig11_metaviper():
    act = pd.read_csv(SRC / "metaviper_selected_activity_v1.tsv", sep="\t")
    emb = pd.read_csv(PH25 / "scdrs_factor_tymp_cell_source_v1.tsv.gz", sep="\t",
                      usecols=["exact_cell_id", "UMAP1", "UMAP2", "F1"])
    d = emb.merge(act, on="exact_cell_id", how="inner", validate="one_to_one")
    genes = ["CASP9", "CHMP1A", "PRSS8", "PXN"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.7))
    for i, (ax, gene) in enumerate(zip(axes.flat, genes)):
        sc = black_umap(ax, d.UMAP1.values, d.UMAP2.values, d[f"{gene}_activity"].values,
                        f"{gene} protein activity", cmap="coolwarm", s=1.0)
        cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
        cb.set_label("metaVIPER activity", fontsize=7)
        cb.ax.tick_params(labelsize=7)
        panel(ax, chr(65+i))
    fig.suptitle("Candidate protein activity across malignant-cell states", y=0.995, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "Fig11_metaVIPER_candidate_activity_umaps_v1")


def fig12_metabolism():
    auc = pd.read_csv(SRC / "four_pathway_AUCell_cell_scores_v1.tsv.gz", sep="\t")
    emb = pd.read_csv(PH25 / "scdrs_factor_tymp_cell_source_v1.tsv.gz", sep="\t",
                      usecols=["exact_cell_id", "UMAP1", "UMAP2"])
    d = auc.merge(emb, on="exact_cell_id", how="inner", validate="one_to_one")
    corr = pd.read_csv(SRC / "four_pathway_AUCell_correlation_summary_v1.tsv", sep="\t")
    paired = pd.read_csv(SRC / "four_pathway_AUCell_paired_statistics_v1.tsv", sep="\t")
    cols = ["AUC_Pyrimidine", "AUC_CysMet", "AUC_ArgPro", "AUC_Glycerophospholipid"]
    labels = ["Pyrimidine", "Cysteine / methionine", "Arginine / proline", "Glycerophospholipid"]
    fig = plt.figure(figsize=(12.2, 8.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 0.92], hspace=0.36, wspace=0.28)
    for i, (col, label) in enumerate(zip(cols, labels)):
        ax = fig.add_subplot(gs[0, i])
        sc = black_umap(ax, d.UMAP1.values, d.UMAP2.values, d[col].values, label, s=0.8)
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=6)
        panel(ax, chr(65+i))
    ax = fig.add_subplot(gs[1, :2])
    piv = corr.pivot(index="pathway", columns="predictor", values="median_patient_rho").reindex(cols)
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax.set_xticks(range(2), ["F1", "TYMP"])
    ax.set_yticks(range(4), labels)
    for i in range(4):
        for j in range(2):
            row = corr[(corr.pathway == cols[i]) & (corr.predictor == piv.columns[j])].iloc[0]
            star = "***" if row.FDR < 0.001 else "**" if row.FDR < 0.01 else "*" if row.FDR < 0.05 else ""
            ax.text(j, i, f"{row.median_patient_rho:+.2f}{star}", ha="center", va="center",
                    color="white" if abs(row.median_patient_rho) > 0.14 else "#202225", fontsize=9)
    ax.set_title("Median within-patient correlation")
    panel(ax, "E")

    ax = fig.add_subplot(gs[1, 2:])
    x = np.arange(4)
    width = 0.34
    p1 = paired[paired.contrast == "TYMP-positive_vs_negative"].set_index("pathway").loc[cols]
    p2 = paired[paired.contrast == "F1-high_clone_vs_other"].set_index("pathway").loc[cols]
    ax.bar(x-width/2, p1.median_patient_delta, width, color=ORANGE, label="TYMP+ vs TYMP−")
    ax.bar(x+width/2, p2.median_patient_delta, width, color=BLUE, label="F1-high clone vs other")
    ax.axhline(0, color="#55585D", lw=0.8)
    for xx, row in zip(x-width/2, p1.itertuples()):
        if row.FDR < 0.05: ax.text(xx, row.median_patient_delta + 0.001, "*", ha="center", va="bottom", fontsize=12)
    for xx, row in zip(x+width/2, p2.itertuples()):
        if row.FDR < 0.05: ax.text(xx, row.median_patient_delta + 0.001, "*", ha="center", va="bottom", fontsize=12)
    ax.set_xticks(x, ["Pyr", "Cys/Met", "Arg/Pro", "Glycero-P"], rotation=0)
    ax.set_ylabel("Median patient activity difference")
    ax.set_title("Patient-paired metabolic contrasts")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    clean(ax, grid=True); panel(ax, "F")
    fig.suptitle("Four-pathway metabolic activity links TYMP and the F1 malignant program", y=0.995, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save(fig, "Fig12_metabolic_pathway_activity_and_contrasts_v1")


def _node_positions(names):
    ang = np.linspace(np.pi/2, np.pi/2 + 2*np.pi, len(names), endpoint=False)
    return {n: (np.cos(a), np.sin(a)) for n, a in zip(names, ang)}


def fig13_cellchat():
    cdir = PH35 / "c01_cellchat_integrated_v1"
    raw = pd.read_csv(cdir / "c01_cellchat_patient_interactions_v1.tsv.gz", sep="\t")
    lr = pd.read_csv(cdir / "c01_cellchat_lr_pair_summary_v1.tsv", sep="\t")
    endpoint = pd.read_csv(cdir / "c01_cellchat_endpoint_summary_v1.tsv", sep="\t")
    nodes = ["TYMP-high Malignant", "TYMP-low Malignant", "SPP1-state Macrophage",
             "C1QC-state Macrophage", "Other Macrophage", "Endothelial"]
    agg = raw.groupby(["patient", "source", "target"], as_index=False).prob.sum()
    med = agg.groupby(["source", "target"], as_index=False).prob.median()
    med = med[(med.source.isin(nodes)) & (med.target.isin(nodes)) & (med.source != med.target)]
    pos = _node_positions(nodes)
    fig = plt.figure(figsize=(13.2, 8.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], height_ratios=[1, 1], hspace=0.38, wspace=0.34)
    ax = fig.add_subplot(gs[:, 0])
    vmax = med.prob.max() if len(med) else 1
    node_colors = dict(zip(nodes, PALETTE))
    for row in med.sort_values("prob").itertuples():
        x1,y1=pos[row.source]; x2,y2=pos[row.target]
        patch = FancyArrowPatch((0.84*x1,0.84*y1),(0.84*x2,0.84*y2),arrowstyle="-|>",mutation_scale=7,
                                connectionstyle="arc3,rad=0.18",lw=0.25+4.2*row.prob/vmax,
                                color=node_colors[row.source],alpha=0.23+0.45*row.prob/vmax)
        ax.add_patch(patch)
    circle_labels = {
        "TYMP-high Malignant": "TYMP-high\nMalignant",
        "TYMP-low Malignant": "TYMP-low\nMalignant",
        "SPP1-state Macrophage": "SPP1-state\nMacrophage",
        "C1QC-state Macrophage": "C1QC-state\nMacrophage",
        "Other Macrophage": "Other\nMacrophage",
        "Endothelial": "Endothelial",
    }
    for n in nodes:
        x,y=pos[n]
        ax.scatter([x],[y],s=550,c=node_colors[n],edgecolor="white",lw=1.2,zorder=5)
        ax.text(1.30*x,1.30*y,circle_labels[n],ha="left" if x>0.2 else "right" if x<-0.2 else "center",
                va="bottom" if y>0.2 else "top",fontsize=8)
    ax.set_xlim(-1.62,1.62); ax.set_ylim(-1.48,1.48); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Median intercellular signaling network")
    panel(ax,"A")

    ax = fig.add_subplot(gs[0,1])
    outv = med.groupby("source").prob.sum().reindex(nodes, fill_value=0)
    inv = med.groupby("target").prob.sum().reindex(nodes, fill_value=0)
    short = {
        "TYMP-high Malignant": "TYMP-high Mal.", "TYMP-low Malignant": "TYMP-low Mal.",
        "SPP1-state Macrophage": "SPP1 Macro.", "C1QC-state Macrophage": "C1QC Macro.",
        "Other Macrophage": "Other Macro.", "Endothelial": "Endothelial",
    }
    offsets = {
        "TYMP-high Malignant": (8, 9), "TYMP-low Malignant": (8, -16),
        "SPP1-state Macrophage": (8, -2), "C1QC-state Macrophage": (-58, -2),
        "Other Macrophage": (8, -8), "Endothelial": (8, 3),
    }
    for n in nodes:
        ax.scatter(outv[n],inv[n],s=50+600*(outv[n]+inv[n])/(outv.add(inv).max()+1e-12),c=node_colors[n],
                   edgecolor="white",lw=0.7)
        ax.annotate(short[n], (outv[n],inv[n]), xytext=offsets[n], textcoords="offset points",fontsize=6.8)
    ax.set_xlabel("Outgoing strength"); ax.set_ylabel("Incoming strength"); ax.set_title("Sender–receiver roles")
    clean(ax,grid=True); panel(ax,"B")

    ax = fig.add_subplot(gs[0,2])
    ep = endpoint.copy()
    endpoint_labels = {
        "endothelial_to_malignant": "Endo. → Mal.",
        "macrophage_to_malignant": "Macro. → Mal.",
        "malignant_to_endothelial": "Mal. → Endo.",
        "malignant_to_macrophage": "Mal. → Macro.",
    }
    ep["label"] = ep.endpoint.map(endpoint_labels)
    colors=[GREEN if x else GRAY for x in ep.statistical_support]
    ax.barh(np.arange(len(ep)),ep.median_high_minus_low,color=colors)
    ax.set_yticks(np.arange(len(ep)),ep.label); ax.invert_yaxis(); ax.set_xlabel("Median high − low communication")
    ax.set_title("Patient-level endpoint shifts")
    for i,r in ep.iterrows(): ax.text(r.median_high_minus_low,i,f"  FDR={r.BH_FDR_4_endpoints:.3g}",va="center",fontsize=7)
    clean(ax,grid=True); panel(ax,"C")

    ax = fig.add_subplot(gs[1,1:])
    lr = lr.assign(score=-np.log10(lr.BH_FDR_within_direction_family.clip(lower=1e-300)))
    dirs = ["malignant_to_endothelial", "malignant_to_macrophage"]
    top = pd.concat([
        lr[lr.direction_family == direction]
          .sort_values(["BH_FDR_within_direction_family", "median_high_minus_low"], ascending=[True, False])
          .head(12)
        for direction in dirs
    ], ignore_index=True)
    ymap={d:i for i,d in enumerate(dirs)}
    top["dir_order"] = top.direction_family.map({d:i for i,d in enumerate(dirs)})
    top=top.sort_values(["dir_order","BH_FDR_within_direction_family"])
    x=np.arange(len(top))
    sc=ax.scatter(x,[ymap[v] for v in top.direction_family],s=25+70*top.score,
                  c=top.median_high_minus_low,cmap="RdBu_r",edgecolor="#444",lw=0.3)
    labels=[p if len(p)<=24 else p[:22]+"…" for p in top.pair_id]
    ax.set_xticks(x,labels,rotation=70,ha="right")
    direction_labels = {
        "malignant_to_endothelial": "Mal. → Endo.",
        "malignant_to_macrophage": "Mal. → Macro.",
    }
    ax.set_yticks(range(len(dirs)),[direction_labels[d] for d in dirs])
    ax.set_title("Top ligand–receptor pair shifts")
    fig.colorbar(sc,ax=ax,fraction=0.025,pad=0.02,label="Median high − low")
    clean(ax,grid=False); panel(ax,"D")
    fig.suptitle("TYMP-stratified malignant cells reshape multicellular communication",y=0.995,fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.97))
    save(fig,"Fig13_CellChat_communication_landscape_v1")


def fig14_nichenet_receptors():
    nndir = PH35 / "c03_nichenet_v1"
    c04dir = PH35 / "c02_c04_analysis_v1"
    pa = pd.read_csv(nndir / "c03_nichenet_patient_ligand_activity_v1.tsv",sep="\t")
    wt = pd.read_csv(nndir / "c03_nichenet_direct_target_weights_v1.tsv",sep="\t")
    um = pd.read_csv(c04dir / "c04_receptor_tymp_f1_umap_source_v1.tsv.gz",sep="\t")
    sm = pd.read_csv(c04dir / "c04_receptor_tymp_f1_summary_v1.tsv",sep="\t")
    fig=plt.figure(figsize=(12.2,9.0))
    gs=fig.add_gridspec(3,4,height_ratios=[1,1,0.95],hspace=0.42,wspace=0.32)
    ax=fig.add_subplot(gs[0,:2])
    rng=np.random.default_rng(26210)
    for i,(lig,color) in enumerate([("MST1",BLUE),("ICAM5",ORANGE)]):
        v=pa.loc[pa.ligand==lig,"complete_target_pearson"].dropna().values
        parts=ax.violinplot(v,positions=[i],showmedians=True,widths=0.65)
        for b in parts["bodies"]: b.set_facecolor(color); b.set_edgecolor(color); b.set_alpha(0.28)
        for key in ["cmedians","cbars","cmins","cmaxes"]:
            if key in parts: parts[key].set_color(color)
        ax.scatter(i+rng.normal(0,0.055,len(v)),v,s=22,c=color,edgecolor="white",lw=0.5,zorder=3)
    ax.axhline(0,color="#555",lw=0.8); ax.set_xticks([0,1],["MST1","ICAM5"])
    ax.set_ylabel("NicheNet ligand activity (Pearson)"); ax.set_title("Patient-level ligand activity")
    clean(ax,grid=True); panel(ax,"A")

    ax=fig.add_subplot(gs[0,2:])
    ligs=["MST1","ICAM5"]; targets=list(dict.fromkeys(wt.target.tolist()))
    for row in wt.itertuples():
        ax.scatter(targets.index(row.target),ligs.index(row.ligand),s=20+260*row.target_weight_percentile,
                   c=row.nichenet_weight,cmap="viridis",vmin=wt.nichenet_weight.min(),vmax=wt.nichenet_weight.max(),
                   edgecolor="#444",lw=0.3)
    ax.set_xticks(range(len(targets)),targets,rotation=45,ha="right"); ax.set_yticks(range(2),ligs)
    ax.set_title("Prioritized ligand–target links")
    norm=mpl.colors.Normalize(wt.nichenet_weight.min(),wt.nichenet_weight.max())
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap="viridis"),ax=ax,fraction=0.035,pad=0.03,label="NicheNet weight")
    clean(ax); panel(ax,"B")

    feats=["MST1R","ITGAL_ITGB2","PRTG","IGDCC4"]
    for i,feat in enumerate(feats):
        ax=fig.add_subplot(gs[1,i])
        sc=black_umap(ax,um.UMAP1.values,um.UMAP2.values,um[feat].values,feat,s=0.85)
        fig.colorbar(sc,ax=ax,fraction=0.045,pad=0.02).ax.tick_params(labelsize=6)
        panel(ax,chr(67+i))

    ax=fig.add_subplot(gs[2,:])
    outcomes=["TYMP_continuous","F1_continuous"]
    rfeats=["MST1R","ITGAL","ITGB2","ITGAL_ITGB2","CD209","PRTG","IGDCC4","ITGAM_ITGB2_sensitivity"]
    mat=sm.pivot(index="receptor_feature",columns="outcome",values="median_patient_rho").reindex(rfeats)[outcomes]
    im=ax.imshow(mat.values,cmap="RdBu_r",vmin=-0.12,vmax=0.12,aspect="auto")
    ax.set_xticks(range(2),["TYMP","F1"]); ax.set_yticks(range(len(rfeats)),[x.replace("_sensitivity"," (sens.)") for x in rfeats])
    for i,rf in enumerate(rfeats):
        for j,outcome in enumerate(outcomes):
            row=sm[(sm.receptor_feature==rf)&(sm.outcome==outcome)].iloc[0]
            txt=f"{row.median_patient_rho:+.2f}\nn={int(row.n_eligible_patients)}"
            if row.n_eligible_patients>=10 and row.BH_FDR_within_outcome8features<0.05: txt += "*"
            ax.text(j,i,txt,ha="center",va="center",fontsize=7,color="white" if abs(row.median_patient_rho)>0.07 else "#202225")
    ax.set_title("Within-patient receptor association (n < 10 shown directionally)")
    fig.colorbar(im,ax=ax,fraction=0.02,pad=0.02,label="Median Spearman rho")
    panel(ax,"G")
    fig.suptitle("NicheNet ligand activity and receptor context for the TYMP/F1 axis",y=0.995,fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.97))
    save(fig,"Fig14_NicheNet_ligand_receptor_evidence_v1")


def km_curve(time, event):
    order=np.argsort(time); time=np.asarray(time)[order]; event=np.asarray(event)[order].astype(int)
    event_times=np.unique(time[event==1]); surv=1.0; xs=[0.0]; ys=[1.0]
    for t in event_times:
        n=np.sum(time>=t); d=np.sum((time==t)&(event==1))
        xs.extend([t,t]); ys.extend([surv,surv*(1-d/n)]); surv*=1-d/n
    xs.append(float(np.nanmax(time))); ys.append(surv)
    return np.asarray(xs),np.asarray(ys)


def _survival_panel(ax,cohort,stats,axis_name,endpoint,letter):
    rcol=axis_name
    tcol=f"{endpoint}_time"; ecol=f"{endpoint}_event"
    d=cohort[["TYMP",rcol,tcol,ecol]].dropna().copy()
    med_t=d.TYMP.median(); med_r=d[rcol].median()
    d["group"]=np.where((d.TYMP>=med_t)&(d[rcol]>=med_r),"Both high",
                  np.where((d.TYMP<med_t)&(d[rcol]<med_r),"Both low","Single high"))
    colors={"Both low":BLUE,"Single high":GRAY,"Both high":RED}
    for group in ["Both low","Single high","Both high"]:
        z=d[d.group==group]; x,y=km_curve(z[tcol].values,z[ecol].values)
        ax.step(x,y,where="post",color=colors[group],lw=1.7,label=f"{group} (n={len(z)})")
    row=stats[(stats.receptor_axis==axis_name)&(stats.endpoint==endpoint)].iloc[0]
    ax.text(0.98,0.96,f"log-rank FDR = {row.BH_FDR_within_endpoint6axes:.2g}",transform=ax.transAxes,
            ha="right",va="top",fontsize=8)
    ax.set_ylim(0,1.03); ax.set_xlabel("Months"); ax.set_ylabel("Survival probability")
    ax.set_title(f"{axis_name} · {endpoint}"); clean(ax,grid=True); panel(ax,letter)


def fig15_survival():
    cohort=pd.read_csv(PH36/"c05_tcga_analysis_cohort_v1.tsv",sep="\t")
    stats=pd.read_csv(PH36/"c05_joint_group_logrank_summary_v1.tsv",sep="\t")
    fig,axes=plt.subplots(2,2,figsize=(10.2,7.8),sharey=True)
    items=[("PRTG","OS"),("PRTG","PFS"),("IGDCC4","OS"),("IGDCC4","PFS")]
    for i,(ax,(axis_name,endpoint)) in enumerate(zip(axes.flat,items)):
        _survival_panel(ax,cohort,stats,axis_name,endpoint,chr(65+i))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.5,-0.005))
    fig.suptitle("TCGA-KIRC joint TYMP–receptor survival stratification",y=0.995,fontsize=13)
    fig.tight_layout(rect=(0,0.06,1,0.97))
    save(fig,"Fig15_TCGA_TYMP_PRTG_IGDCC4_survival_v1")

    axes_names=["MST1R","ICAM5_ITGAL_ITGB2_primary","CD209","PRTG","IGDCC4","ICAM5_ITGAM_ITGB2_sensitivity"]
    fig,axes=plt.subplots(6,2,figsize=(10.2,18.0),sharey=True)
    k=0
    for i,axis_name in enumerate(axes_names):
        for j,endpoint in enumerate(["OS","PFS"]):
            _survival_panel(axes[i,j],cohort,stats,axis_name,endpoint,chr(65+k)); k+=1
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(0.5,0.006))
    fig.suptitle("TCGA-KIRC joint TYMP–receptor survival results across six receptor axes",y=0.997,fontsize=13)
    fig.tight_layout(rect=(0,0.025,1,0.987))
    save(fig,"FigS01_TCGA_TYMP_receptor_survival_all_axes_v1")


def main():
    fig10_regulon_networks()
    fig11_metaviper()
    fig12_metabolism()
    fig13_cellchat()
    fig14_nichenet_receptors()
    fig15_survival()
    print(f"generated 7 figure stems in {OUT}")


if __name__ == "__main__":
    main()
