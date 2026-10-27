#!/usr/bin/env python3
"""Independent verifier for the extended Part1 CNV consolidation threshold grid."""

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
from sklearn.metrics import adjusted_rand_score


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
PREFLIGHT = ROOT / "04_subclones_tymp/cnv_subclone_preflight_v1"
SOURCE = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_v1"
OUT = ROOT / "04_subclones_tymp/cnv_subclone_threshold_extension_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
THRESHOLDS = tuple(round(value, 3) for value in np.arange(0.10, 0.5001, 0.025))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty verification directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        rows.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
        if not passed:
            pd.DataFrame(rows).to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
            raise RuntimeError(f"Verification failed: {name}: {detail}")

    gate = subprocess.run([str(MEMORY_GATE), "cnv_subclone_threshold_extension_verification_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    check("memory_gate", gate.returncode == 0 and "PASS" in gate_text, gate_text)
    check("source_complete", (SOURCE / "RUN_COMPLETE").is_file(), str(SOURCE))

    manifest = pd.read_csv(SOURCE / "file_manifest_v1.tsv", sep="\t")
    for row in manifest.itertuples(index=False):
        path = SOURCE / row.file
        check(f"file_exists::{row.file}", path.is_file(), str(path))
        check(f"file_size::{row.file}", path.stat().st_size == row.size_bytes, f"observed={path.stat().st_size}; expected={row.size_bytes}")
        check(f"file_sha256::{row.file}", sha256(path) == row.sha256, row.sha256)

    metadata = json.loads((SOURCE / "run_metadata_v1.json").read_text(encoding="utf-8"))
    check("not_formal_selection", metadata["formal_clone_selected"] is False, str(metadata["formal_clone_selected"]))
    check("threshold_grid_metadata", np.allclose(metadata["thresholds"], THRESHOLDS), str(metadata["thresholds"]))

    inventory = pd.read_csv(PREFLIGHT / "raw_hmm_group_profile_inventory_v1.tsv", sep="\t")
    distances = pd.read_csv(PREFLIGHT / "within_patient_group_profile_distances_v1.tsv.gz", sep="\t")
    assignments = pd.read_csv(SOURCE / "group_to_clone_extended_grid_v1.tsv.gz", sep="\t")
    patient = pd.read_csv(SOURCE / "extended_threshold_patient_summary_v1.tsv", sep="\t")
    cohort = pd.read_csv(SOURCE / "extended_threshold_cohort_summary_v1.tsv", sep="\t")
    stability = pd.read_csv(SOURCE / "extended_threshold_stability_v1.tsv", sep="\t")
    check("assignment_rows", len(assignments) == 264 * len(THRESHOLDS), str(len(assignments)))
    check("assignment_unique", not assignments[["distance_threshold", "infercnv_unit", "infercnv_group"]].duplicated().any(), "threshold/unit/group")
    check("patient_rows", len(patient) == 26 * len(THRESHOLDS), str(len(patient)))
    check("cohort_rows", len(cohort) == len(THRESHOLDS), str(len(cohort)))
    check("stability_rows", len(stability) == len(THRESHOLDS) * (len(THRESHOLDS) - 1) // 2, str(len(stability)))

    cell_labels: dict[float, list[str]] = {threshold: [] for threshold in THRESHOLDS}
    max_partition_error = 0.0
    for unit, frame in inventory.groupby("infercnv_unit", sort=True, observed=True):
        groups = frame.sort_values("infercnv_group", kind="stable").reset_index(drop=True)
        names = groups["infercnv_group"].tolist()
        weights = groups["n_selected_cells"].to_numpy(int)
        matrix = np.zeros((len(names), len(names)), dtype=float)
        index = {name: idx for idx, name in enumerate(names)}
        for row in distances.loc[distances["infercnv_unit"].eq(unit)].itertuples(index=False):
            i, j = index[row.group_a], index[row.group_b]
            matrix[i, j] = matrix[j, i] = float(row.bin_hamming_fraction)
        tree = None if len(names) == 1 else linkage(squareform(matrix, checks=True), method="complete")
        for threshold in THRESHOLDS:
            expected = np.ones(len(names), dtype=int) if tree is None else fcluster(tree, t=threshold, criterion="distance")
            cut = assignments.loc[assignments["infercnv_unit"].eq(unit) & np.isclose(assignments["distance_threshold"], threshold)].set_index("infercnv_group").loc[names]
            observed = cut["clone_number"].to_numpy(int)
            partition_ari = adjusted_rand_score(expected, observed)
            max_partition_error = max(max_partition_error, abs(1.0 - partition_ari))
            check(f"partition::{unit}::{threshold}", abs(partition_ari - 1.0) <= 1e-12, f"ARI={partition_ari:.12g}")
            for clone_id, clone_frame in cut.groupby("clone_id"):
                idx = np.asarray([index[name] for name in clone_frame.index], dtype=int)
                expected_cells = int(weights[idx].sum())
                check(f"clone_cells::{unit}::{threshold}::{clone_id}", clone_frame["clone_cells"].nunique() == 1 and int(clone_frame["clone_cells"].iloc[0]) == expected_cells, str(expected_cells))
                diameter = 0.0 if len(idx) == 1 else float(matrix[np.ix_(idx, idx)].max())
                check(f"diameter::{unit}::{threshold}::{clone_id}", diameter <= threshold + 1e-12, f"diameter={diameter:.12g}")
            for name, weight, clone_id in zip(names, weights, cut["clone_id"], strict=True):
                cell_labels[threshold].extend([str(clone_id)] * int(weight))

            declared_patient = patient.loc[patient["infercnv_unit"].eq(unit) & np.isclose(patient["distance_threshold"], threshold)]
            check(f"patient_row::{unit}::{threshold}", len(declared_patient) == 1, str(len(declared_patient)))
            declared_patient = declared_patient.iloc[0]
            clone_sizes = cut.groupby("clone_id")["clone_cells"].first().astype(int)
            primary = clone_sizes.loc[clone_sizes >= 20]
            check(f"patient_clone_count::{unit}::{threshold}", int(declared_patient["n_consolidated_clones"]) == len(clone_sizes), str(len(clone_sizes)))
            check(f"patient_primary_count::{unit}::{threshold}", int(declared_patient["n_primary_clones_min20"]) == len(primary), str(len(primary)))
            check(f"patient_primary_cells::{unit}::{threshold}", int(declared_patient["n_cells_in_primary_clones_min20"]) == int(primary.sum()), str(int(primary.sum())))

    check("all_partitions_recomputed", max_partition_error <= 1e-12, f"max_1_minus_ari={max_partition_error:.3g}")
    for threshold in THRESHOLDS:
        check(f"cell_label_count::{threshold}", len(cell_labels[threshold]) == 11950, str(len(cell_labels[threshold])))
    for row in stability.itertuples(index=False):
        observed = adjusted_rand_score(cell_labels[float(row.threshold_a)], cell_labels[float(row.threshold_b)])
        check(f"stability::{row.threshold_a}::{row.threshold_b}", abs(observed - row.cell_weighted_adjusted_rand_index) <= 1e-12, f"ARI={observed:.12g}")
    for row in cohort.itertuples(index=False):
        subset = patient.loc[np.isclose(patient["distance_threshold"], row.distance_threshold)]
        check(f"cohort_cells::{row.distance_threshold}", int(subset["n_selected_cells"].sum()) == int(row.n_cells) == 11950, str(row.n_cells))
        check(f"cohort_clone_count::{row.distance_threshold}", int(subset["n_consolidated_clones"].sum()) == int(row.n_consolidated_clones), str(row.n_consolidated_clones))
        check(f"cohort_primary_count::{row.distance_threshold}", int(subset["n_primary_clones_min20"].sum()) == int(row.n_primary_clones_min20), str(row.n_primary_clones_min20))
        check(f"cohort_evaluable::{row.distance_threshold}", int(subset["evaluable_two_primary_clones_min20"].sum()) == int(row.n_evaluable_patients_two_primary_clones_min20), str(row.n_evaluable_patients_two_primary_clones_min20))

    results = pd.DataFrame(rows)
    results.to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "part1_cnv_subclone_threshold_extension_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "n_checks": int(len(results)),
        "n_pass": int(results["status"].eq("PASS").sum()),
        "n_fail": int(results["status"].eq("FAIL").sum()),
        "n_thresholds_recomputed": len(THRESHOLDS),
        "n_patient_threshold_partitions_recomputed": 26 * len(THRESHOLDS),
        "dependency_changes": False,
    }
    (OUT / "verification_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    files = sorted(path for path in OUT.iterdir() if path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"})
    pd.DataFrame({
        "file": [path.name for path in files],
        "size_bytes": [path.stat().st_size for path in files],
        "sha256": [sha256(path) for path in files],
    }).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("PART1_CNV_SUBCLONE_THRESHOLD_EXTENSION_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
