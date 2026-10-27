#!/usr/bin/env python3
"""Extend the verified patient-specific CNV consolidation threshold grid."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, silhouette_samples


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
SOURCE = ROOT / "04_subclones_tymp/cnv_subclone_preflight_v1"
VERIFICATION = ROOT / "04_subclones_tymp/cnv_subclone_preflight_verification_v1"
OUT = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
THRESHOLDS = tuple(round(value, 3) for value in np.arange(0.10, 0.5001, 0.025))
PRIMARY_MINIMUM = 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_labels(raw: np.ndarray, names: list[str], weights: np.ndarray) -> np.ndarray:
    clusters = []
    for label in np.unique(raw):
        idx = np.flatnonzero(raw == label)
        clusters.append((int(label), int(weights[idx].sum()), "|".join(sorted(names[i] for i in idx))))
    clusters.sort(key=lambda row: (-row[1], row[2]))
    mapping = {old: new for new, (old, _, _) in enumerate(clusters, start=1)}
    return np.asarray([mapping[int(value)] for value in raw], dtype=int)


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty threshold extension: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gate = subprocess.run([str(MEMORY_GATE), "cnv_subclone_threshold_extension_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    if gate.returncode != 0 or "PASS" not in gate_text:
        raise RuntimeError(f"Memory gate blocked threshold extension: {gate_text}")
    if not (VERIFICATION / "RUN_COMPLETE").is_file():
        raise RuntimeError("Verified preflight is required")

    inventory = pd.read_csv(SOURCE / "raw_hmm_group_profile_inventory_v1.tsv", sep="\t")
    distances = pd.read_csv(SOURCE / "within_patient_group_profile_distances_v1.tsv.gz", sep="\t")
    all_assignments: list[dict[str, object]] = []
    patient_rows: list[dict[str, object]] = []
    cell_labels: dict[float, list[str]] = {threshold: [] for threshold in THRESHOLDS}

    for unit, frame in inventory.groupby("infercnv_unit", sort=True, observed=True):
        groups = frame.sort_values("infercnv_group", kind="stable").reset_index(drop=True)
        names = groups["infercnv_group"].tolist()
        weights = groups["n_selected_cells"].to_numpy(int)
        matrix = np.zeros((len(names), len(names)), dtype=float)
        index = {name: idx for idx, name in enumerate(names)}
        unit_distances = distances.loc[distances["infercnv_unit"].eq(unit)]
        for row in unit_distances.itertuples(index=False):
            i, j = index[row.group_a], index[row.group_b]
            matrix[i, j] = matrix[j, i] = float(row.bin_hamming_fraction)
        tree = None if len(names) == 1 else linkage(squareform(matrix, checks=True), method="complete")

        for threshold in THRESHOLDS:
            raw = np.ones(len(names), dtype=int) if tree is None else fcluster(tree, t=threshold, criterion="distance")
            labels = deterministic_labels(raw, names, weights)
            code = f"T{int(round(threshold * 1000)):03d}"
            clone_cells = {label: int(weights[labels == label].sum()) for label in np.unique(labels)}
            for idx, group_name in enumerate(names):
                label = int(labels[idx])
                clone_id = f"{unit}__{code}_C{label:02d}"
                all_assignments.append({
                    "distance_threshold": threshold,
                    "infercnv_unit": unit,
                    "infercnv_group": group_name,
                    "raw_group_cells": int(weights[idx]),
                    "clone_number": label,
                    "clone_id": clone_id,
                    "clone_cells": clone_cells[label],
                })
                cell_labels[threshold].extend([clone_id] * int(weights[idx]))

            primary = [size for size in clone_cells.values() if size >= PRIMARY_MINIMUM]
            if 1 < len(np.unique(labels)) < len(labels):
                sil = silhouette_samples(matrix, labels, metric="precomputed")
                weighted_silhouette = float(np.average(sil, weights=weights))
            else:
                weighted_silhouette = np.nan
            diameters = []
            for label in np.unique(labels):
                idx = np.flatnonzero(labels == label)
                diameters.append(0.0 if len(idx) == 1 else float(matrix[np.ix_(idx, idx)].max()))
            patient_rows.append({
                "distance_threshold": threshold,
                "infercnv_unit": unit,
                "n_selected_cells": int(weights.sum()),
                "n_raw_hmm_groups": int(len(names)),
                "n_consolidated_clones": int(len(clone_cells)),
                "n_primary_clones_min20": int(len(primary)),
                "n_cells_in_primary_clones_min20": int(sum(primary)),
                "fraction_cells_in_primary_clones_min20": float(sum(primary) / weights.sum()),
                "evaluable_two_primary_clones_min20": bool(len(primary) >= 2),
                "cell_weighted_group_silhouette": weighted_silhouette,
                "maximum_clone_diameter": float(max(diameters)),
            })

    assignments = pd.DataFrame(all_assignments)
    patient = pd.DataFrame(patient_rows)
    if len(assignments) != len(inventory) * len(THRESHOLDS):
        raise RuntimeError("Extended assignment grid is incomplete")
    if any(len(cell_labels[threshold]) != 11950 for threshold in THRESHOLDS):
        raise RuntimeError("Cell-weighted threshold labels are incomplete")

    cohort = (
        patient.groupby("distance_threshold", observed=True)
        .agg(
            n_patients=("infercnv_unit", "size"),
            n_cells=("n_selected_cells", "sum"),
            n_consolidated_clones=("n_consolidated_clones", "sum"),
            n_primary_clones_min20=("n_primary_clones_min20", "sum"),
            n_cells_in_primary_clones_min20=("n_cells_in_primary_clones_min20", "sum"),
            n_evaluable_patients_two_primary_clones_min20=("evaluable_two_primary_clones_min20", "sum"),
            median_patient_cell_weighted_group_silhouette=("cell_weighted_group_silhouette", "median"),
            min_patient_cell_weighted_group_silhouette=("cell_weighted_group_silhouette", "min"),
            maximum_clone_diameter=("maximum_clone_diameter", "max"),
        ).reset_index()
    )
    cohort["fraction_cells_in_primary_clones_min20"] = cohort["n_cells_in_primary_clones_min20"] / cohort["n_cells"]

    stability_rows = []
    for i, threshold_a in enumerate(THRESHOLDS):
        for threshold_b in THRESHOLDS[i + 1:]:
            stability_rows.append({
                "threshold_a": threshold_a,
                "threshold_b": threshold_b,
                "threshold_difference": round(threshold_b - threshold_a, 3),
                "cell_weighted_adjusted_rand_index": float(adjusted_rand_score(cell_labels[threshold_a], cell_labels[threshold_b])),
                "adjacent_grid_thresholds": bool(abs(threshold_b - threshold_a - 0.025) < 1e-9),
            })
    stability = pd.DataFrame(stability_rows)

    assignments.to_csv(OUT / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    patient.to_csv(OUT / "extended_threshold_patient_summary_v1.tsv", sep="\t", index=False)
    cohort.to_csv(OUT / "extended_threshold_cohort_summary_v1.tsv", sep="\t", index=False)
    stability.to_csv(OUT / "extended_threshold_stability_v1.tsv", sep="\t", index=False)
    metadata = {
        "schema_version": "part1_cnv_subclone_threshold_extension_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "source_verified_preflight": str(VERIFICATION),
        "thresholds": list(THRESHOLDS),
        "threshold_step": 0.025,
        "primary_clone_minimum_cells": PRIMARY_MINIMUM,
        "purpose": "threshold sensitivity only; no formal clone selected",
        "formal_clone_selected": False,
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({
        "file": [path.name for path in files],
        "size_bytes": [path.stat().st_size for path in files],
        "sha256": [sha256(path) for path in files],
    }).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(cohort.to_string(index=False))
    print("\nAdjacent threshold stability:")
    print(stability.loc[stability["adjacent_grid_thresholds"]].to_string(index=False))
    print("PART1_CNV_SUBCLONE_THRESHOLD_EXTENSION_V1_PASS")


if __name__ == "__main__":
    main()
