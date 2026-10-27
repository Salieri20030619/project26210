#!/usr/bin/env python3
"""Build a fast patient/study-balanced reference-v2 for Part 4 sensitivity pilot."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse


INPUT = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/"
    "cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
)
USAGE = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/"
    "full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz"
)
OLD_SIGNATURE = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/ref/ref_signatures_matrix.csv"
)
MEMORY_GATE = Path("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/reference_v2")

SPP1_PROGRAM = ["SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
C1QC_PROGRAM = ["C1QA", "C1QB", "C1QC", "APOE"]
KEY_GENES = sorted(set(SPP1_PROGRAM + C1QC_PROGRAM + ["TYMP", "CD74", "TREM2", "VEGFA", "VEGFB"]))
MIN_MALIGNANT_PER_PATIENT = 20
TOP_MARKERS_PER_STATE = 100

BROAD_MAP = {
    "B cell": "B cell",
    "Dendritic cell": "Dendritic cell",
    "Endothelial": "Endothelial",
    "Mast cell": "Mast cell",
    "Monocyte": "Monocyte",
    "NK cell": "NK cell",
    "Pericyte/SMC": "SMC/Pericyte",
    "T cell": "T cell",
}


def gate(stage: str, rows: list[dict[str, object]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    rows.append({"stage": stage, "exit_code": result.returncode, "output": line})
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate failed at {stage}: {line}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    memory_rows: list[dict[str, object]] = []
    gate("before_load", memory_rows)

    old = pd.read_csv(OLD_SIGNATURE, index_col=0)
    old_genes = old.index.astype(str).tolist()
    obj = ad.read_h5ad(INPUT, backed="r")
    genes = [gene for gene in old_genes if gene in obj.var_names]
    for gene in KEY_GENES:
        if gene in obj.var_names and gene not in genes:
            genes.append(gene)
    if len(genes) < 12000:
        raise RuntimeError(f"Too few old-reference genes retained: {len(genes)}")

    required = [
        "cell_type_broad_v1", "eligible_for_subclone_selection_v1", "curated_patient_id",
        "study", "source_tissue", "total_counts_raw",
    ]
    obs = obj.obs[required].copy()
    obs.index = obs.index.astype(str)
    obs["patient"] = obs["curated_patient_id"].astype(str)
    obs["study_v1"] = obs["study"].astype(str)
    obs["broad"] = obs["cell_type_broad_v1"].astype(str)
    obs["state_v2"] = pd.Series(pd.NA, index=obs.index, dtype="string")

    usage = pd.read_csv(USAGE, sep="\t", usecols=["exact_cell_id", "F1"]).set_index("exact_cell_id")
    obs["F1"] = usage["F1"].reindex(obs.index).to_numpy()
    eligible = obs["eligible_for_subclone_selection_v1"].astype(bool)
    if int(eligible.sum()) != 11950 or int(obs.loc[eligible, "F1"].notna().sum()) != 11950:
        raise RuntimeError("Incomplete F1 mapping for frozen malignant cells")

    malignant_counts = obs.loc[eligible].groupby("patient", observed=True).size()
    malignant_patients = malignant_counts[malignant_counts >= MIN_MALIGNANT_PER_PATIENT].index
    malignant_mask = eligible & obs["patient"].isin(malignant_patients)
    f1_threshold = obs.loc[malignant_mask].groupby("patient", observed=True)["F1"].transform(lambda x: x.quantile(0.75))
    obs.loc[malignant_mask, "state_v2"] = np.where(
        obs.loc[malignant_mask, "F1"].to_numpy() >= f1_threshold.to_numpy(),
        "F1-high Malignant",
        "F1-low Malignant",
    )

    for broad, state in BROAD_MAP.items():
        mask = obs["broad"].eq(broad)
        obs.loc[mask, "state_v2"] = state
    epithelial_normal = obs["broad"].eq("Epithelial cell") & obs["source_tissue"].astype(str).eq("Normal")
    obs.loc[epithelial_normal, "state_v2"] = "Epithelial-Normal"

    gate("before_selected_counts", memory_rows)
    indices = [int(obj.var_names.get_loc(gene)) for gene in genes]
    counts = obj.layers["counts"][:, indices]
    if hasattr(counts, "to_memory"):
        counts = counts.to_memory()
    if not issparse(counts):
        counts = csr_matrix(np.asarray(counts))
    counts = counts.tocsr()
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise RuntimeError("Invalid raw counts")
    gene_index = {gene: i for i, gene in enumerate(genes)}

    macro = obs["broad"].eq("Macrophage").to_numpy()
    macro_pos = np.flatnonzero(macro)
    libs = pd.to_numeric(obs.iloc[macro_pos]["total_counts_raw"], errors="coerce").to_numpy(dtype=float)
    def program_score(program: list[str]) -> np.ndarray:
        matrix = counts[macro_pos][:, [gene_index[g] for g in program]].toarray().astype(float)
        normalized = np.log1p(matrix / np.maximum(libs[:, None], 1.0) * 1e4)
        return normalized.mean(axis=1)
    spp1_score = program_score(SPP1_PROGRAM)
    c1qc_score = program_score(C1QC_PROGRAM)
    macro_frame = pd.DataFrame({
        "patient": obs.iloc[macro_pos]["patient"].to_numpy(),
        "spp1": spp1_score,
        "c1qc": c1qc_score,
    }, index=obs.index[macro_pos])
    macro_frame["spp1_pct"] = macro_frame.groupby("patient", observed=True)["spp1"].rank(method="average", pct=True)
    macro_frame["c1qc_pct"] = macro_frame.groupby("patient", observed=True)["c1qc"].rank(method="average", pct=True)
    macro_state = np.full(len(macro_frame), "Other Macrophage", dtype=object)
    spp1_high = (macro_frame["spp1_pct"] >= 0.75) & (macro_frame["spp1_pct"] >= macro_frame["c1qc_pct"])
    c1qc_high = (macro_frame["c1qc_pct"] >= 0.75) & ~spp1_high
    macro_state[spp1_high.to_numpy()] = "SPP1-state Macrophage"
    macro_state[c1qc_high.to_numpy()] = "C1QC-state Macrophage"
    obs.loc[macro_frame.index, "state_v2"] = macro_state
    obs.loc[macro_frame.index, "SPP1_program_score_v2"] = macro_frame["spp1"].to_numpy()
    obs.loc[macro_frame.index, "C1QC_program_score_v2"] = macro_frame["c1qc"].to_numpy()

    included = obs["state_v2"].notna().to_numpy()
    included_pos = np.flatnonzero(included)
    included_obs = obs.iloc[included_pos].copy()
    included_obs["patient_state"] = included_obs["patient"] + "||" + included_obs["state_v2"].astype(str)
    patient_states = sorted(included_obs["patient_state"].unique())
    ps_code = pd.Categorical(included_obs["patient_state"], categories=patient_states).codes
    indicator = csr_matrix(
        (np.ones(len(included_pos)), (ps_code, np.arange(len(included_pos)))),
        shape=(len(patient_states), len(included_pos)),
    )
    pb_counts = indicator @ counts[included_pos]
    libsum_by_cell = pd.to_numeric(included_obs["total_counts_raw"], errors="coerce").to_numpy(dtype=float)
    pb_lib = np.asarray(indicator @ libsum_by_cell).reshape(-1)
    pb_cpm = pb_counts.multiply((1e6 / np.maximum(pb_lib, 1.0))[:, None]).toarray()
    ps_meta = pd.DataFrame({"patient_state": patient_states})
    ps_meta[["patient", "state_v2"]] = ps_meta["patient_state"].str.split(
        "||", n=1, expand=True, regex=False
    )
    patient_to_study = included_obs.drop_duplicates("patient").set_index("patient")["study_v1"]
    ps_meta["study"] = ps_meta["patient"].map(patient_to_study)

    signature_rows = []
    signatures = {}
    for state, state_rows in ps_meta.groupby("state_v2", observed=True):
        study_profiles = []
        for study, study_rows in state_rows.groupby("study", observed=True):
            study_profiles.append(pb_cpm[study_rows.index].mean(axis=0))
        signature = np.vstack(study_profiles).mean(axis=0)
        signatures[state] = signature
        signature_rows.append({
            "state_v2": state,
            "n_cells": int(included_obs["state_v2"].eq(state).sum()),
            "n_patients": int(state_rows["patient"].nunique()),
            "n_studies": int(state_rows["study"].nunique()),
            "median_cells_per_patient": float(
                included_obs.loc[included_obs["state_v2"].eq(state)].groupby("patient", observed=True).size().median()
            ),
        })
    signature = pd.DataFrame(signatures, index=genes)
    signature.to_csv(OUT / "reference_v2_signatures_study_balanced_cpm.tsv.gz", sep="\t")
    pd.DataFrame(signature_rows).sort_values("state_v2").to_csv(OUT / "reference_v2_state_qc.tsv", sep="\t", index=False)
    ps_counts = included_obs.groupby(["study_v1", "patient", "state_v2"], observed=True).size().rename("n_cells").reset_index()
    ps_counts.to_csv(OUT / "reference_v2_patient_state_counts.tsv", sep="\t", index=False)

    excluded_gene = pd.Series(False, index=signature.index)
    upper = signature.index.str.upper()
    excluded_gene |= upper.str.startswith(("MT-", "RPL", "RPS", "HBA", "HBB"))
    log_signature = np.log1p(signature)
    marker_rows = []
    for state in signature.columns:
        background = log_signature.drop(columns=state).mean(axis=1)
        specificity = log_signature[state] - background
        candidates = specificity.loc[~excluded_gene].nlargest(TOP_MARKERS_PER_STATE)
        for rank, (gene, score) in enumerate(candidates.items(), 1):
            marker_rows.append({"state_v2": state, "gene": gene, "rank": rank, "specificity_log1p_cpm": score})
    markers = pd.DataFrame(marker_rows)
    marker_genes = sorted(set(markers["gene"]) | set(KEY_GENES))
    marker_signature = signature.loc[signature.index.intersection(marker_genes)]
    marker_signature.to_csv(OUT / "reference_v2_marker_signatures_cpm.tsv", sep="\t")
    markers.to_csv(OUT / "reference_v2_marker_ranking.tsv", sep="\t", index=False)

    assignments = obs.loc[included, ["patient", "study_v1", "state_v2", "F1", "SPP1_program_score_v2", "C1QC_program_score_v2"]]
    assignments.index.name = "exact_cell_id"
    assignments.to_csv(OUT / "reference_v2_cell_assignments.tsv.gz", sep="\t")
    pd.DataFrame(memory_rows).to_csv(OUT / "memory_gates_v2.tsv", sep="\t", index=False)

    metadata = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input": str(INPUT),
        "input_verified_sha256": "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7",
        "rank7_usage": str(USAGE),
        "n_input_cells": int(obj.n_obs),
        "n_reference_cells": int(included.sum()),
        "n_reference_genes": int(len(genes)),
        "n_marker_genes": int(len(marker_signature)),
        "n_states": int(signature.shape[1]),
        "malignant_definition": "11,950 frozen malignant cells; patients >=20 cells; within-patient F1 upper quartile",
        "macrophage_definition": "within-patient continuous program percentiles; mutually exclusive upper-quartile SPP1/C1QC states",
        "signature_aggregation": "patient CPM -> study mean -> equal-weight study mean",
        "purpose": "fast NNLS sensitivity pilot, not replacement formal cell2location reference",
    }
    (OUT / "run_metadata_v2.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if getattr(obj, "file", None) is not None:
        obj.file.close()

    manifest_rows = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "file_manifest_v2.tsv":
            manifest_rows.append({
                "file": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    pd.DataFrame(manifest_rows).to_csv(OUT / "file_manifest_v2.tsv", sep="\t", index=False)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
