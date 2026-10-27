#!/usr/bin/env python3
"""Independent verification of Part1 CNV subclone profile preflight v1."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
INFER = ROOT / "03_malignancy_cnv/infercnv_v2/patients"
META = ROOT / "05_cnmf/input_v1/cnmf_cell_metadata_v1.tsv.gz"
SOURCE = ROOT / "04_subclones_tymp/cnv_subclone_preflight_v1"
OUT = ROOT / "04_subclones_tymp/cnv_subclone_preflight_verification_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
THRESHOLDS = (0.05, 0.10, 0.15, 0.20)
MINIMUMS = (10, 20, 30, 50)
BIN_SIZE = 5_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def autosome(value: object) -> int | None:
    text = str(value).removeprefix("chr")
    return int(text) if text.isdigit() and 1 <= int(text) <= 22 else None


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty verification directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})
        if not condition:
            pd.DataFrame(checks).to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
            raise RuntimeError(f"Verification failed: {name}: {detail}")

    gate = subprocess.run([str(MEMORY_GATE), "cnv_subclone_preflight_verification_v1"], text=True, capture_output=True)
    gate_text = " ".join((gate.stdout + " " + gate.stderr).split())
    (OUT / "memory_gate_v1.txt").write_text(gate_text + "\n", encoding="utf-8")
    check("memory_gate", gate.returncode == 0 and "PASS" in gate_text, gate_text)
    check("source_complete", (SOURCE / "RUN_COMPLETE").is_file(), str(SOURCE / "RUN_COMPLETE"))

    declared = pd.read_csv(SOURCE / "file_manifest_v1.tsv", sep="\t")
    check("source_manifest_schema", list(declared.columns) == ["file", "size_bytes", "sha256"], str(list(declared.columns)))
    for row in declared.itertuples(index=False):
        path = SOURCE / row.file
        check(f"source_file_exists::{row.file}", path.is_file(), str(path))
        check(f"source_size::{row.file}", path.stat().st_size == row.size_bytes, f"observed={path.stat().st_size}; expected={row.size_bytes}")
        check(f"source_sha256::{row.file}", sha256(path) == row.sha256, row.sha256)

    metadata = json.loads((SOURCE / "run_metadata_v1.json").read_text(encoding="utf-8"))
    check("preflight_not_formal", metadata["formal_clone_selected"] is False, str(metadata["formal_clone_selected"]))
    check("method_constants", metadata["bin_size_bp"] == BIN_SIZE and metadata["neutral_hmm_state"] == 3.0, str(metadata))

    meta = pd.read_csv(META, sep="\t", dtype={"exact_cell_id": str, "infercnv_unit": str, "infercnv_group": str})
    inventory = pd.read_csv(SOURCE / "raw_hmm_group_profile_inventory_v1.tsv", sep="\t")
    assignments = pd.read_csv(SOURCE / "group_to_clone_threshold_grid_v1.tsv.gz", sep="\t")
    distances = pd.read_csv(SOURCE / "within_patient_group_profile_distances_v1.tsv.gz", sep="\t")
    patient = pd.read_csv(SOURCE / "threshold_patient_summary_v1.tsv", sep="\t")
    cohort = pd.read_csv(SOURCE / "threshold_cohort_summary_v1.tsv", sep="\t")
    stability = pd.read_csv(SOURCE / "threshold_assignment_stability_v1.tsv", sep="\t")
    groups = meta.groupby(["infercnv_unit", "infercnv_group"], observed=True).size().rename("n_selected_cells").reset_index()

    check("cell_count", len(meta) == 11950 and not meta["exact_cell_id"].duplicated().any(), f"n={len(meta)}")
    check("patient_count", meta["infercnv_unit"].nunique() == 26, str(meta["infercnv_unit"].nunique()))
    check("raw_group_count", len(groups) == 264, str(len(groups)))
    check("inventory_unique", not inventory[["infercnv_unit", "infercnv_group"]].duplicated().any(), str(len(inventory)))
    merged_counts = groups.merge(inventory, on=["infercnv_unit", "infercnv_group"], how="outer", suffixes=("_meta", "_inventory"), indicator=True)
    check("inventory_group_identity", merged_counts["_merge"].eq("both").all(), merged_counts["_merge"].value_counts().to_dict().__str__())
    check("inventory_cell_conservation", merged_counts["n_selected_cells_meta"].equals(merged_counts["n_selected_cells_inventory"]), str(merged_counts[["n_selected_cells_meta", "n_selected_cells_inventory"]].sum().to_dict()))
    check("assignment_grid_rows", len(assignments) == 264 * len(THRESHOLDS), str(len(assignments)))
    check("assignment_grid_unique", not assignments[["distance_threshold", "infercnv_unit", "infercnv_group"]].duplicated().any(), "threshold/unit/group")
    check("distance_pair_count", len(distances) == sum(len(x) * (len(x) - 1) // 2 for _, x in groups.groupby("infercnv_unit")), str(len(distances)))
    check("distance_bounds", distances["bin_hamming_fraction"].between(0, 1).all(), str(distances["bin_hamming_fraction"].agg(["min", "max"]).to_dict()))

    max_distance_error = 0.0
    max_inventory_shift_error = 0.0
    for unit, unit_groups in groups.groupby("infercnv_unit", sort=True, observed=True):
        run_dir = INFER / unit / "infercnv_run"
        genes = pd.read_csv(run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat", sep="\t", index_col=0).reset_index(names="gene")
        genes["autosome"] = genes["chr"].map(autosome)
        genes = genes.loc[genes["autosome"].notna()].copy()
        genes["autosome"] = genes["autosome"].astype(int)
        genes["mid"] = (pd.to_numeric(genes["start"]) + pd.to_numeric(genes["stop"])) / 2
        genes["bin"] = genes["autosome"].astype(str) + ":" + np.floor((genes["mid"] - 1) / BIN_SIZE).astype(int).astype(str)
        genes = genes.sort_values(["autosome", "bin", "mid", "gene"], kind="stable").reset_index(drop=True)
        codes, names = pd.factorize(genes["bin"], sort=False)
        gene_to_row = pd.Series(np.arange(len(genes)), index=genes["gene"])
        pred = pd.read_csv(run_dir / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat", sep="\t")
        ordered = unit_groups.sort_values("infercnv_group", kind="stable").reset_index(drop=True)
        group_names = ordered["infercnv_group"].tolist()
        profiles = []
        for group_name in group_names:
            rows = pred.loc[pred["cell_group_name"].eq(group_name) & pred["gene"].isin(gene_to_row.index), ["gene", "state"]]
            gene_states = np.full(len(genes), 3.0)
            gene_states[gene_to_row.loc[rows["gene"]].to_numpy(int)] = rows["state"].to_numpy(float)
            profile = pd.Series(gene_states).groupby(codes, sort=False).median().to_numpy(float)
            profiles.append(profile)
            inv = inventory.loc[inventory["infercnv_unit"].eq(unit) & inventory["infercnv_group"].eq(group_name)].iloc[0]
            shift_error = abs(float(np.mean(np.abs(profile - 3.0))) - float(inv["mean_abs_bin_state_shift"]))
            max_inventory_shift_error = max(max_inventory_shift_error, shift_error)
            check(f"profile_bin_count::{unit}::{group_name}", len(profile) == len(names) == int(inv["n_5mb_bins"]), f"n={len(profile)}")

        profiles_array = np.asarray(profiles)
        matrix = np.zeros((len(group_names), len(group_names))) if len(group_names) == 1 else squareform(pdist(profiles_array, metric="hamming"))
        observed_pairs = distances.loc[distances["infercnv_unit"].eq(unit)]
        lookup = {(row.group_a, row.group_b): row.bin_hamming_fraction for row in observed_pairs.itertuples(index=False)}
        for i, group_a in enumerate(group_names):
            for j in range(i + 1, len(group_names)):
                group_b = group_names[j]
                error = abs(float(matrix[i, j]) - float(lookup[(group_a, group_b)]))
                max_distance_error = max(max_distance_error, error)

        unit_assignment = assignments.loc[assignments["infercnv_unit"].eq(unit)]
        cell_counts = ordered.set_index("infercnv_group")["n_selected_cells"]
        for threshold in THRESHOLDS:
            cut = unit_assignment.loc[np.isclose(unit_assignment["distance_threshold"], threshold)].set_index("infercnv_group")
            check(f"threshold_groups::{unit}::{threshold}", set(cut.index) == set(group_names), f"n={len(cut)}")
            expected_clone_cells = cut.groupby("clone_id").apply(lambda frame: int(cell_counts.loc[frame.index].sum()), include_groups=False)
            declared_clone_cells = cut.groupby("clone_id")["clone_cells"].first().astype(int)
            check(f"clone_cell_counts::{unit}::{threshold}", expected_clone_cells.equals(declared_clone_cells), f"sum={expected_clone_cells.sum()}")
            labels = cut.loc[group_names, "clone_id"].to_numpy()
            for clone in np.unique(labels):
                idx = np.flatnonzero(labels == clone)
                diameter = 0.0 if len(idx) == 1 else float(matrix[np.ix_(idx, idx)].max())
                check(f"complete_link_diameter::{unit}::{threshold}::{clone}", diameter <= threshold + 1e-12, f"diameter={diameter:.12g}")

    check("profile_mean_shift_recompute", max_inventory_shift_error <= 1e-12, f"max_abs_error={max_inventory_shift_error:.3g}")
    check("all_pairwise_distances_recomputed", max_distance_error <= 1e-12, f"max_abs_error={max_distance_error:.3g}")

    cell_labels: dict[float, list[str]] = {}
    for threshold in THRESHOLDS:
        cut = assignments.loc[np.isclose(assignments["distance_threshold"], threshold), ["infercnv_unit", "infercnv_group", "clone_id"]]
        joined = meta[["infercnv_unit", "infercnv_group"]].merge(cut, on=["infercnv_unit", "infercnv_group"], how="left", validate="many_to_one")
        check(f"cell_assignment_complete::{threshold}", joined["clone_id"].notna().all() and len(joined) == 11950, str(joined["clone_id"].isna().sum()))
        cell_labels[threshold] = joined["clone_id"].tolist()
    for row in stability.itertuples(index=False):
        ari = adjusted_rand_score(cell_labels[float(row.threshold_a)], cell_labels[float(row.threshold_b)])
        check(f"ari_recomputed::{row.threshold_a}::{row.threshold_b}", abs(ari - row.cell_weighted_adjusted_rand_index) <= 1e-12, f"observed={ari:.12g}")

    check("patient_grid_shape", len(patient) == 26 * len(THRESHOLDS) * len(MINIMUMS), str(len(patient)))
    check("cohort_grid_shape", len(cohort) == len(THRESHOLDS) * len(MINIMUMS), str(len(cohort)))
    for row in cohort.itertuples(index=False):
        subset = patient.loc[np.isclose(patient["distance_threshold"], row.distance_threshold) & patient["primary_min_cells"].eq(row.primary_min_cells)]
        check(f"cohort_cells::{row.distance_threshold}::{row.primary_min_cells}", int(subset["n_selected_cells"].sum()) == int(row.n_cells) == 11950, str(row.n_cells))
        check(f"cohort_clones::{row.distance_threshold}::{row.primary_min_cells}", int(subset["n_consolidated_clones"].sum()) == int(row.n_consolidated_clones), str(row.n_consolidated_clones))
        check(f"cohort_primary_cells::{row.distance_threshold}::{row.primary_min_cells}", int(subset["n_cells_in_primary_clones"].sum()) == int(row.n_cells_in_primary_clones), str(row.n_cells_in_primary_clones))

    results = pd.DataFrame(checks)
    results.to_csv(OUT / "validation_results_v1.tsv", sep="\t", index=False)
    summary = {
        "schema_version": "part1_cnv_subclone_preflight_verification_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "source": str(SOURCE),
        "n_checks": int(len(results)),
        "n_pass": int(results["status"].eq("PASS").sum()),
        "n_fail": int(results["status"].eq("FAIL").sum()),
        "recomputed_profiles": 264,
        "recomputed_pairwise_distances": int(len(distances)),
        "max_pairwise_distance_abs_error": max_distance_error,
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
    print("PART1_CNV_SUBCLONE_PREFLIGHT_VERIFICATION_V1_PASS")


if __name__ == "__main__":
    main()
