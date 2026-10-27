#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import os
import random
import subprocess
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import cellrank as cr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.stats import spearmanr, wilcoxon
import torch
from cytotrace2_py.cytotrace2_py import cytotrace2
from sctour.train import Trainer


SEED = 26210
PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = PROJECT / "results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1"
INPUT, QC, FIG, TABLE, SOURCE, OBJECT = (OUT / x for x in ("input", "qc", "figures", "tables", "source_data", "objects"))
GATE = PROJECT / "check_memory_gate.sh"


def gate(stage: str) -> None:
    result = subprocess.run([str(GATE), stage], text=True, capture_output=True)
    with (QC / "memory_gate_log_v1.tsv").open("a", encoding="utf-8") as handle:
        handle.write(stage + "\t" + str(result.returncode) + "\t" + (result.stdout + result.stderr).strip().replace("\t", "|") + "\n")
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked {stage}")


def read_vector(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()


def scale01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float); out = np.full(x.shape, np.nan); ok = np.isfinite(x)
    if not ok.any(): return out
    lo, hi = x[ok].min(), x[ok].max(); out[ok] = 0.5 if hi == lo else (x[ok] - lo) / (hi - lo)
    return out


def rho(x: np.ndarray, y: np.ndarray, minimum: int = 20) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float); ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < minimum or np.std(x[ok]) == 0 or np.std(y[ok]) == 0: return np.nan
    return float(spearmanr(x[ok], y[ok]).statistic)


def bh(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float); out = np.full(p.shape, np.nan); ok = np.isfinite(p)
    if not ok.any(): return out
    v = p[ok]; order = np.argsort(v); ranked = v[order] * len(v) / np.arange(1, len(v) + 1); ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty_like(ranked); restored[order] = np.minimum(ranked, 1.0); out[ok] = restored; return out


def plot_embedding(umap: np.ndarray, values: np.ndarray, title: str, stem: str, cmap: str = "plasma") -> None:
    order = np.argsort(np.nan_to_num(values, nan=-np.inf)); fig, ax = plt.subplots(figsize=(4.4, 3.8))
    scatter = ax.scatter(umap[order, 0], umap[order, 1], c=values[order], s=1.0, lw=0, cmap=cmap, rasterized=True)
    ax.set_title(title, loc="left", fontsize=9); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False); fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.03)
    fig.tight_layout(); fig.savefig(FIG / f"{stem}.png", dpi=400, facecolor="white"); fig.savefig(FIG / f"{stem}.pdf"); plt.close(fig)


def main() -> None:
    if not (OUT / "R_METHODS_COMPLETE").exists() or (OUT / "PYTHON_METHODS_COMPLETE").exists():
        raise RuntimeError("R methods missing or Python methods already complete")
    random.seed(SEED); np.random.seed(SEED); torch.manual_seed(SEED); torch.set_num_threads(min(16, os.cpu_count() or 1))
    start = time.time()
    gate("part1_trajectory_python_input_v1")
    cells = read_vector(INPUT / "cells_v1.tsv"); genes = read_vector(INPUT / "genes_v1.tsv")
    matrix = mmread(INPUT / "raw_counts_genes_by_cells_v1.mtx.gz").T.tocsr().astype(np.float32)
    metadata = pd.read_csv(INPUT / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id").loc[cells]
    r_scores = pd.read_csv(SOURCE / "trajectory_r_cell_scores_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id").loc[cells]
    for column in [c for c in r_scores.columns if c.startswith("Monocle3_") or c.startswith("Slingshot_")]: metadata[column] = r_scores[column]
    latent = pd.read_csv(INPUT / "scvi_latent_v1.tsv.gz", sep="\t", index_col=0).loc[cells].to_numpy(np.float32)
    umap = pd.read_csv(INPUT / "scvi_umap_v1.tsv.gz", sep="\t", index_col=0).loc[cells].to_numpy(np.float32)
    if matrix.shape != (11950, 2500) or metadata.index.tolist() != cells: raise ValueError("Python trajectory input invariant failed")
    adata = ad.AnnData(X=matrix, obs=metadata.copy(), var=pd.DataFrame(index=genes)); adata.obsm["X_scVI"] = latent; adata.obsm["X_umap"] = umap
    adata.obs["n_genes_by_counts"] = np.asarray((matrix > 0).sum(axis=1)).ravel().astype(np.int32)
    adata.obs["total_counts"] = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)

    studies = metadata["study"].astype(str).unique().tolist()
    cyto_saved = SOURCE / "CytoTRACE2_predictions_by_study_v1.tsv.gz"
    if cyto_saved.exists() and cyto_saved.stat().st_size > 0:
        cyto = pd.read_csv(cyto_saved, sep="\t", index_col=0)
        if cyto.index.astype(str).tolist() != cells or cyto["study"].astype(str).nunique() != 5 or cyto["CytoTRACE2_Score"].notna().sum() != 11950:
            raise ValueError("Saved CytoTRACE2 checkpoint invariant failed")
        print("PART1_TRAJECTORY_REUSING_CYTOTRACE2_CHECKPOINT n_cells=11950 n_studies=5")
    else:
        gate("part1_trajectory_cytotrace2_load_v1")
        cyto_genes = read_vector(INPUT / "cytotrace_genes_v1.tsv")
        cyto_matrix = mmread(INPUT / "cytotrace_counts_genes_by_cells_v1.mtx.gz").T.tocsr().astype(np.float32)
        cyto_results = []
        cyto_input_dir = INPUT / "cytotrace_by_study"; cyto_output_dir = OUT / "python_output/cytotrace2_by_study"
        cyto_input_dir.mkdir(parents=True, exist_ok=True); cyto_output_dir.mkdir(parents=True, exist_ok=True)
        for study in studies:
            gate(f"part1_trajectory_cytotrace2_{study}_v1")
            idx = np.flatnonzero(metadata["study"].astype(str).to_numpy() == study); study_cells = [cells[i] for i in idx]
            study_file = cyto_input_dir / f"{study}_raw_counts_genes_by_cells_v1.tsv.gz"
            frame = pd.DataFrame(cyto_matrix[idx, :].T.toarray(), index=cyto_genes, columns=study_cells); frame.index.name = "gene"
            frame.to_csv(study_file, sep="\t", compression="gzip"); del frame
            result = cytotrace2(str(study_file), species="human", batch_size=len(study_cells), smooth_batch_size=min(1000, len(study_cells)), disable_parallelization=False, max_cores=8, seed=SEED, output_dir=str(cyto_output_dir / study), disable_plotting=True, disable_verbose=False)
            result = result.loc[study_cells].copy(); result["study"] = study; cyto_results.append(result)
        cyto = pd.concat(cyto_results).loc[cells]
        cyto.to_csv(cyto_saved, sep="\t", compression="gzip")
        del cyto_matrix
    for column in ("CytoTRACE2_Score", "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score"):
        adata.obs[column] = pd.to_numeric(cyto[column], errors="coerce").to_numpy()
    adata.obs["CytoTRACE2_Potency"] = pd.Categorical(cyto["CytoTRACE2_Potency"])

    gate("part1_trajectory_sctour_v1")
    sctour_scores_path = SOURCE / "sctour_cell_scores_v1.tsv.gz"
    sctour_latent_path = OBJECT / "sctour_latent_v1.npy"
    sctour_marker = OUT / "SCTOUR_COMPLETE"
    if sctour_marker.exists() and sctour_scores_path.exists() and sctour_latent_path.exists():
        saved = pd.read_csv(sctour_scores_path, sep="\t", dtype={"exact_cell_id": str}).set_index("exact_cell_id").loc[cells]
        mix_z = np.load(sctour_latent_path)
        if mix_z.shape != (11950, 5) or not bool(saved[["scTour_global_pseudotime_v1", "scTour_pseudotime_patient_v1"]].notna().all().all()):
            raise ValueError("Saved scTour checkpoint invariant failed")
        global_time = saved["scTour_global_pseudotime_v1"].to_numpy(float); patient_time = saved["scTour_pseudotime_patient_v1"].to_numpy(float)
        global_rho_before = float(saved["scTour_global_orientation_rho_before_v1"].iloc[0])
        print("PART1_TRAJECTORY_REUSING_SCTOUR_CHECKPOINT n_cells=11950")
    else:
        trainer = Trainer(adata, percent=0.85, n_latent=5, n_ode_hidden=25, n_vae_hidden=128, loss_mode="nb", nepoch=200, batch_size=1024, random_state=SEED, val_frac=0.10, use_gpu=False)
        trainer.train(); raw_time = np.asarray(trainer.get_time(), float); mix_z, _, _ = trainer.get_latentsp(batch_size=1024)
        anchor = pd.to_numeric(adata.obs["progression_anchor_v1"], errors="coerce").to_numpy()
        global_rho_before = rho(raw_time, anchor)
        if np.isfinite(global_rho_before) and global_rho_before < 0: raw_time = 1.0 - raw_time
        global_time = scale01(raw_time); patient_time = np.full(adata.n_obs, np.nan)
        for _, idx in adata.obs.groupby("curated_patient_id", observed=True, sort=False).indices.items(): patient_time[idx] = scale01(global_time[idx])
        trainer.save_model(str(OBJECT), "sctour_part1_seed26210_v1")
        losses = pd.DataFrame({"epoch": np.arange(1, len(trainer.log["train_loss"]) + 1), "train_loss": trainer.log["train_loss"], "validation_loss": trainer.log["validation_loss"]})
        losses.to_csv(QC / "sctour_training_loss_v1.tsv", sep="\t", index=False)
        np.save(sctour_latent_path, np.asarray(mix_z, np.float32))
        pd.DataFrame({"exact_cell_id": cells, "scTour_global_pseudotime_v1": global_time, "scTour_pseudotime_patient_v1": patient_time, "scTour_global_orientation_rho_before_v1": global_rho_before}).to_csv(sctour_scores_path, sep="\t", index=False, compression="gzip")
        sctour_marker.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
    adata.obs["scTour_global_pseudotime_v1"] = global_time; adata.obs["scTour_pseudotime_patient_v1"] = patient_time; adata.obsm["X_scTour"] = np.asarray(mix_z, np.float32)

    gate("part1_trajectory_cellrank_v1")
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_scVI", random_state=SEED)
    pk = cr.kernels.PseudotimeKernel(adata, time_key="scTour_pseudotime_patient_v1").compute_transition_matrix(threshold_scheme="soft", b=10.0, nu=0.5, n_jobs=8, show_progress_bar=True)
    ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix(); kernel = 0.8 * pk + 0.2 * ck; kernel.compute_transition_matrix()
    estimator = cr.estimators.GPCCA(kernel)
    terminal_cells = {}; used = set()
    for state in ("Hypoxia", "Cycle", "TYMP_state"):
        score = f"state_{state}_score_v1"; candidate = adata.obs.loc[adata.obs["malignant_state_v1"].astype(str) == state].sort_values(["progression_anchor_v1", score], ascending=False)
        chosen = [x for x in candidate.index if x not in used][:50]
        if len(chosen) < 20: chosen = [x for x in adata.obs.sort_values([score, "progression_anchor_v1"], ascending=False).index if x not in used][:50]
        terminal_cells[state] = chosen; used.update(chosen)
    estimator.set_terminal_states(terminal_cells)
    estimator.compute_fate_probabilities(solver="gmres", use_petsc=False, n_jobs=8, show_progress_bar=True, tol=1e-6)
    fate = np.asarray(estimator.fate_probabilities.X, float); fate_names = [str(x) for x in estimator.fate_probabilities.names]
    for j, name in enumerate(fate_names): adata.obs[f"CellRank_fate_{name}_v1"] = fate[:, j]
    adata.obs["CellRank_max_fate_probability_v1"] = fate.max(axis=1); adata.obs["CellRank_predicted_fate_v1"] = pd.Categorical(np.asarray(fate_names, object)[fate.argmax(axis=1)])
    pd.DataFrame(fate, index=cells, columns=fate_names).to_csv(SOURCE / "cellrank_fate_probabilities_v1.tsv.gz", sep="\t", compression="gzip")
    transition = kernel.transition_matrix.tocsr(); velocity = np.asarray(transition @ umap - umap, np.float32)
    pd.DataFrame({"exact_cell_id": cells, "UMAP1": umap[:,0], "UMAP2": umap[:,1], "dUMAP1": velocity[:,0], "dUMAP2": velocity[:,1]}).to_csv(SOURCE / "cellrank_pseudotime_transition_vectors_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    methods = ["Monocle3_pseudotime_patient_v1", "Slingshot_pseudotime_patient_v1", "CytoTRACE2_Score", "scTour_pseudotime_patient_v1"]
    targets = ["progression_anchor_v1", "TYMP_log1p_CP10K_v1", "F1", "infercnv_burden_over_threshold_v1", "AUC_Pyrimidine", "AUC_CysMet", "AUC_ArgPro", "AUC_Glycerophospholipid", "Pyrimidine_flux", "CysMet_flux", "ArgPro_flux", "Choline_proxy"]
    rows = []
    for method in methods:
        for target in targets:
            for patient, part in adata.obs.groupby("curated_patient_id", observed=True): rows.append({"method": method, "target": target, "curated_patient_id": patient, "n_cells": len(part), "rho": rho(part[method], part[target])})
    corr = pd.DataFrame(rows).dropna(subset=["rho"]); corr.to_csv(TABLE / "trajectory_all_methods_patient_correlations_v1.tsv", sep="\t", index=False)
    summary_rows = []
    for (method, target), part in corr.groupby(["method", "target"], observed=True):
        values = part["rho"].to_numpy(); p = 1.0 if np.allclose(values, 0) else float(wilcoxon(values).pvalue)
        summary_rows.append({"method": method, "target": target, "n_patients": len(values), "median_patient_rho": np.median(values), "fraction_positive": np.mean(values > 0), "signed_rank_PValue": p})
    summary = pd.DataFrame(summary_rows); summary["FDR"] = summary.groupby("method", group_keys=False)["signed_rank_PValue"].transform(lambda x: bh(x))
    summary.to_csv(TABLE / "trajectory_all_methods_correlation_summary_v1.tsv", sep="\t", index=False)
    output_cols = [c for c in adata.obs.columns if c.startswith(("Monocle3_", "Slingshot_", "CytoTRACE2_", "scTour_", "CellRank_"))] + ["curated_patient_id", "study", "sample_id", "progression_anchor_v1", "malignant_state_v1", "TYMP_log1p_CP10K_v1", "F1", "formal_inferred_cnv_clone_v1"]
    adata.obs[output_cols].to_csv(SOURCE / "trajectory_all_methods_cell_scores_v1.tsv.gz", sep="\t", compression="gzip")
    plot_embedding(umap, adata.obs["CytoTRACE2_Score"].to_numpy(float), "CytoTRACE2 potency", "CytoTRACE2_potency_v1", "viridis")
    plot_embedding(umap, patient_time, "scTour cross-sectional pseudotime", "scTour_pseudotime_v1")
    plot_embedding(umap, adata.obs["CellRank_max_fate_probability_v1"].to_numpy(float), "CellRank fate confidence", "CellRank_fate_confidence_v1", "viridis")
    for j, name in enumerate(fate_names): plot_embedding(umap, fate[:,j], f"CellRank fate: {name}", f"CellRank_fate_{name}_v1", "viridis")
    adata.write_h5ad(OBJECT / "part1_malignant_trajectory_all_methods_v1.h5ad", compression="gzip")
    metadata_out = {"schema_version": "part1_trajectory_python_methods_v1", "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"), "n_cells": adata.n_obs, "n_genes_sctour": adata.n_vars, "n_studies_cytotrace2": len(studies), "sctour_epochs": 200, "sctour_global_orientation_rho_before": global_rho_before, "cellrank_kernel": "0.8 pseudotime + 0.2 connectivity", "cellrank_terminal_states": fate_names, "has_velocity_layers": False, "interpretation": "pseudotime-informed cross-sectional fate probabilities;not RNA velocity or lineage tracing", "runtime_seconds": time.time()-start}
    (QC / "trajectory_python_metadata_v1.json").write_text(json.dumps(metadata_out, indent=2), encoding="utf-8")
    (OUT / "PYTHON_METHODS_COMPLETE").write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
    print(json.dumps(metadata_out, indent=2))


if __name__ == "__main__": main()
