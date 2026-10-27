#!/usr/bin/env python3
"""Preflight patient-specific CNV clone consolidation from inferCNV HMM profiles."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
INFER = ROOT / "03_malignancy_cnv/infercnv_v2/patients"
META = ROOT / "05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz"
OUT = ROOT / "04_subclones_tymp/cnv_subclone_preflight_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
BIN_SIZE = 5_000_000
NEUTRAL_STATE = 3.0
DISTANCE_THRESHOLDS = (0.05, 0.10, 0.15, 0.20)
PRIMARY_MIN_CELLS = (10, 20, 30, 50)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE), stage], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    pd.DataFrame(records).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False)
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def chr_number(value: object) -> int | None:
    text = str(value).removeprefix("chr")
    return int(text) if text.isdigit() and 1 <= int(text) <= 22 else None


def profile_hash(profile: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(profile, dtype=np.float64).tobytes()).hexdigest()


def deterministic_clone_labels(raw_labels: np.ndarray, group_names: list[str], group_cells: np.ndarray) -> dict[int, int]:
    records = []
    for raw_label in sorted(np.unique(raw_labels)):
        idx = np.flatnonzero(raw_labels == raw_label)
        records.append((raw_label, int(group_cells[idx].sum()), "|".join(sorted(group_names[i] for i in idx))))
    records.sort(key=lambda x: (-x[1], x[2]))
    return {raw_label: order for order, (raw_label, _, _) in enumerate(records, start=1)}


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty preflight directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    gate("cnv_subclone_preflight_load", gates)

    meta = pd.read_csv(META, sep="\t", dtype={"exact_cell_id": str, "infercnv_unit": str, "infercnv_group": str})
    if len(meta) != 11950 or meta["exact_cell_id"].duplicated().any():
        raise RuntimeError("Unexpected strict cNMF/CNV metadata")
    group_counts = (
        meta.groupby(["infercnv_unit", "infercnv_group"], observed=True)
        .size().rename("n_selected_cells").reset_index()
    )
    if len(group_counts) != 264 or group_counts["n_selected_cells"].sum() != 11950:
        raise RuntimeError("Unexpected selected HMM group inventory")

    group_inventory: list[dict[str, object]] = []
    assignments: list[dict[str, object]] = []
    pairwise: list[dict[str, object]] = []
    threshold_patient: list[dict[str, object]] = []
    cell_assignments_by_threshold: dict[float, list[str]] = {value: [] for value in DISTANCE_THRESHOLDS}

    for unit, unit_groups in group_counts.groupby("infercnv_unit", sort=True, observed=True):
        gate(f"cnv_subclone_preflight_{unit}", gates)
        run_dir = INFER / unit / "infercnv_run"
        genes_path = run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat"
        pred_path = run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"
        if not genes_path.is_file() or not pred_path.is_file():
            raise RuntimeError(f"{unit}: missing inferCNV HMM input")
        genes = pd.read_csv(genes_path, sep="\t", index_col=0).reset_index(names="gene")
        if not {"chr", "start", "stop"}.issubset(genes.columns) or genes["gene"].duplicated().any():
            raise RuntimeError(f"{unit}: invalid genes-used table")
        all_gene_names = set(genes["gene"])
        genes["chr_number"] = genes["chr"].map(chr_number)
        genes = genes.loc[genes["chr_number"].notna()].copy()
        genes["chr_number"] = genes["chr_number"].astype(int)
        genes["midpoint"] = (pd.to_numeric(genes["start"]) + pd.to_numeric(genes["stop"])) / 2
        genes["bin_index"] = np.floor((genes["midpoint"] - 1) / BIN_SIZE).astype(int)
        genes["bin_key"] = genes["chr_number"].astype(str) + ":" + genes["bin_index"].astype(str)
        genes = genes.sort_values(["chr_number", "bin_index", "midpoint", "gene"], kind="stable").reset_index(drop=True)
        bin_codes, bin_names = pd.factorize(genes["bin_key"], sort=False)
        n_bins = len(bin_names)

        pred = pd.read_csv(pred_path, sep="\t", dtype={"cell_group_name": str, "gene": str, "state": float})
        required_pred = {"cell_group_name", "gene", "state"}
        if not required_pred.issubset(pred.columns):
            raise RuntimeError(f"{unit}: invalid predicted-CNV table columns")
        if pred[["cell_group_name", "gene"]].duplicated().any():
            raise RuntimeError(f"{unit}: duplicate group/gene predicted-CNV rows")
        if pred["state"].isna().any() or not pred["state"].isin(range(1, 7)).all():
            raise RuntimeError(f"{unit}: HMM states must be integers 1-6")
        unknown_genes = sorted(set(pred["gene"]) - all_gene_names)
        if unknown_genes:
            raise RuntimeError(f"{unit}: predicted genes absent from genes-used table: {unknown_genes[:10]}")
        pred = pred.loc[pred["cell_group_name"].isin(unit_groups["infercnv_group"])].copy()
        groups = unit_groups.sort_values("infercnv_group", kind="stable").reset_index(drop=True)
        group_names = groups["infercnv_group"].tolist()
        group_cells = groups["n_selected_cells"].to_numpy(int)
        if set(group_names) != set(pred["cell_group_name"].unique()):
            missing = sorted(set(group_names) - set(pred["cell_group_name"].unique()))
            raise RuntimeError(f"{unit}: selected groups missing HMM profiles: {missing}")
        gene_index = pd.Series(np.arange(len(genes)), index=genes["gene"])
        profiles = np.full((len(group_names), n_bins), NEUTRAL_STATE, dtype=np.float64)

        for group_idx, group_name in enumerate(group_names):
            altered = pred.loc[pred["cell_group_name"].eq(group_name), ["gene", "state"]]
            altered = altered.loc[altered["gene"].isin(gene_index.index)]
            state = np.full(len(genes), NEUTRAL_STATE, dtype=np.float64)
            state[gene_index.loc[altered["gene"]].to_numpy(int)] = altered["state"].to_numpy(float)
            bin_profile = pd.Series(state).groupby(bin_codes, sort=False).median().to_numpy(float)
            if len(bin_profile) != n_bins or np.any(~np.isfinite(bin_profile)):
                raise RuntimeError(f"{unit}/{group_name}: invalid binned profile")
            profiles[group_idx] = bin_profile
            shift = bin_profile - NEUTRAL_STATE
            group_inventory.append({
                "infercnv_unit": unit,
                "infercnv_group": group_name,
                "n_selected_cells": int(group_cells[group_idx]),
                "n_genes_used_autosomes": int(len(genes)),
                "n_5mb_bins": int(n_bins),
                "n_altered_5mb_bins": int(np.count_nonzero(shift)),
                "n_loss_5mb_bins": int(np.count_nonzero(shift < 0)),
                "n_gain_5mb_bins": int(np.count_nonzero(shift > 0)),
                "mean_abs_bin_state_shift": float(np.mean(np.abs(shift))),
                "profile_sha256": profile_hash(bin_profile),
            })

        if len(group_names) == 1:
            hamming = np.zeros((1, 1), dtype=float)
            l1 = np.zeros((1, 1), dtype=float)
            linkage_matrix = None
        else:
            hamming = squareform(pdist(profiles, metric="hamming"))
            l1 = squareform(pdist(profiles, metric="cityblock")) / n_bins
            linkage_matrix = linkage(squareform(hamming, checks=False), method="complete")
            for i in range(len(group_names)):
                for j in range(i + 1, len(group_names)):
                    pairwise.append({
                        "infercnv_unit": unit,
                        "group_a": group_names[i],
                        "group_b": group_names[j],
                        "group_a_cells": int(group_cells[i]),
                        "group_b_cells": int(group_cells[j]),
                        "bin_hamming_fraction": float(hamming[i, j]),
                        "mean_abs_bin_state_difference": float(l1[i, j]),
                    })

        meta_unit = meta.loc[meta["infercnv_unit"].eq(unit), ["exact_cell_id", "infercnv_group"]]
        for threshold in DISTANCE_THRESHOLDS:
            raw = np.ones(len(group_names), dtype=int) if linkage_matrix is None else fcluster(linkage_matrix, t=threshold, criterion="distance")
            relabel = deterministic_clone_labels(raw, group_names, group_cells)
            labels = np.asarray([relabel[value] for value in raw], dtype=int)
            clone_cells = {clone: int(group_cells[labels == clone].sum()) for clone in np.unique(labels)}
            threshold_code = f"T{int(round(threshold * 100)):02d}"
            mapping = dict(zip(group_names, labels, strict=True))
            cell_assignments_by_threshold[threshold].extend(
                f"{unit}__{threshold_code}_C{mapping[group]:02d}" for group in meta_unit["infercnv_group"]
            )
            for group_idx, group_name in enumerate(group_names):
                clone = int(labels[group_idx])
                assignments.append({
                    "distance_threshold": threshold,
                    "infercnv_unit": unit,
                    "infercnv_group": group_name,
                    "raw_group_cells": int(group_cells[group_idx]),
                    "clone_number": clone,
                    "clone_id": f"{unit}__{threshold_code}_C{clone:02d}",
                    "clone_cells": clone_cells[clone],
                })
            for minimum in PRIMARY_MIN_CELLS:
                primary = [value for value in clone_cells.values() if value >= minimum]
                threshold_patient.append({
                    "distance_threshold": threshold,
                    "primary_min_cells": minimum,
                    "infercnv_unit": unit,
                    "n_selected_cells": int(group_cells.sum()),
                    "n_raw_hmm_groups": int(len(group_names)),
                    "n_consolidated_clones": int(len(clone_cells)),
                    "n_primary_clones": int(len(primary)),
                    "n_cells_in_primary_clones": int(sum(primary)),
                    "fraction_cells_in_primary_clones": float(sum(primary) / group_cells.sum()),
                    "evaluable_two_primary_clones": bool(len(primary) >= 2),
                })

    inventory_df = pd.DataFrame(group_inventory)
    assignment_df = pd.DataFrame(assignments)
    pairwise_df = pd.DataFrame(pairwise)
    patient_df = pd.DataFrame(threshold_patient)
    if len(inventory_df) != 264 or inventory_df["n_selected_cells"].sum() != 11950:
        raise RuntimeError("Group inventory failed conservation")
    if len(assignment_df) != 264 * len(DISTANCE_THRESHOLDS):
        raise RuntimeError("Threshold assignment grid incomplete")

    cohort_df = (
        patient_df.groupby(["distance_threshold", "primary_min_cells"], observed=True)
        .agg(
            n_patients=("infercnv_unit", "size"),
            n_cells=("n_selected_cells", "sum"),
            n_raw_hmm_groups=("n_raw_hmm_groups", "sum"),
            n_consolidated_clones=("n_consolidated_clones", "sum"),
            n_primary_clones=("n_primary_clones", "sum"),
            n_cells_in_primary_clones=("n_cells_in_primary_clones", "sum"),
            n_evaluable_patients_two_primary_clones=("evaluable_two_primary_clones", "sum"),
        ).reset_index()
    )
    cohort_df["fraction_cells_in_primary_clones"] = cohort_df["n_cells_in_primary_clones"] / cohort_df["n_cells"]

    stability_rows = []
    for i, threshold_a in enumerate(DISTANCE_THRESHOLDS):
        labels_a = cell_assignments_by_threshold[threshold_a]
        if len(labels_a) != 11950:
            raise RuntimeError("Cell assignment order/count mismatch")
        for threshold_b in DISTANCE_THRESHOLDS[i + 1:]:
            labels_b = cell_assignments_by_threshold[threshold_b]
            stability_rows.append({
                "threshold_a": threshold_a,
                "threshold_b": threshold_b,
                "cell_weighted_adjusted_rand_index": float(adjusted_rand_score(labels_a, labels_b)),
            })
    stability_df = pd.DataFrame(stability_rows)

    inventory_df.to_csv(OUT / "raw_hmm_group_profile_inventory_v1.tsv", sep="\t", index=False)
    assignment_df.to_csv(OUT / "group_to_clone_threshold_grid_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pairwise_df.to_csv(OUT / "within_patient_group_profile_distances_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    patient_df.to_csv(OUT / "threshold_patient_summary_v1.tsv", sep="\t", index=False)
    cohort_df.to_csv(OUT / "threshold_cohort_summary_v1.tsv", sep="\t", index=False)
    stability_df.to_csv(OUT / "threshold_assignment_stability_v1.tsv", sep="\t", index=False)

    metadata = {
        "schema_version": "part1_cnv_subclone_preflight_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "purpose": "method and threshold preflight only; no formal clone selected",
        "n_cells": 11950,
        "n_patients": 26,
        "n_raw_hmm_groups": 264,
        "neutral_hmm_state": NEUTRAL_STATE,
        "autosomes": "chr1-chr22",
        "bin_size_bp": BIN_SIZE,
        "profile_aggregation": "median HMM state among genes in each 5 Mb bin; unreported genes set to neutral state 3",
        "distance": "fraction of 5 Mb bins with unequal median HMM state",
        "clustering": "patient-specific complete-linkage hierarchical clustering",
        "distance_thresholds": list(DISTANCE_THRESHOLDS),
        "primary_clone_minimum_cell_grid": list(PRIMARY_MIN_CELLS),
        "formal_clone_selected": False,
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    manifest = pd.DataFrame({
        "file": [path.name for path in files],
        "size_bytes": [path.stat().st_size for path in files],
        "sha256": [sha256(path) for path in files],
    })
    manifest.to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(cohort_df.to_string(index=False))
    print(stability_df.to_string(index=False))
    print("PART1_CNV_SUBCLONE_PREFLIGHT_V1_PASS")


if __name__ == "__main__":
    main()
