#!/usr/bin/env python3
"""Build patient-specific inferred-CNV distance trees with rank-7 cNMF annotations."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from collections import defaultdict, deque
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "input_snapshot_local_v1"
HMM = INPUT / "hmm_profiles"
OUT = ROOT / "results_local_v1"
BIN_SIZE = 5_000_000
NEUTRAL = 3.0
FACTORS = [f"F{i}" for i in range(1, 8)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, compression=compression, lineterminator="\n")
    part.replace(path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    part.replace(path)


def chr_number(value: object) -> int | None:
    text = str(value).removeprefix("chr")
    return int(text) if text.isdigit() and 1 <= int(text) <= 22 else None


def profile_hash(profile: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(profile, dtype=np.float64).tobytes()).hexdigest()


def weighted_median_states(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    result = np.empty(values.shape[1], dtype=np.float64)
    half = float(weights.sum()) / 2.0
    for col in range(values.shape[1]):
        order = np.argsort(values[:, col], kind="stable")
        cumulative = np.cumsum(weights[order])
        result[col] = values[order[np.searchsorted(cumulative, half, side="left")], col]
    return result


class UnionFind:
    def __init__(self, nodes: list[str]) -> None:
        self.parent = {node: node for node in nodes}

    def find(self, node: str) -> str:
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left: str, right: str) -> bool:
        root_left, root_right = self.find(left), self.find(right)
        if root_left == root_right:
            return False
        if root_left > root_right:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        return True


def deterministic_mst(profiles: dict[str, np.ndarray]) -> list[tuple[str, str, float]]:
    nodes = sorted(profiles)
    candidates: list[tuple[float, str, str]] = []
    for i, left in enumerate(nodes):
        for right in nodes[i + 1:]:
            distance = float(np.mean(profiles[left] != profiles[right]))
            candidates.append((distance, left, right))
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    uf = UnionFind(nodes)
    selected = []
    for distance, left, right in candidates:
        if uf.union(left, right):
            selected.append((left, right, distance))
            if len(selected) == len(nodes) - 1:
                break
    if len(selected) != len(nodes) - 1:
        raise RuntimeError("MST construction did not connect all nodes")
    return selected


def orient_tree(edges: list[tuple[str, str, float]], root: str) -> list[tuple[str, str, float, int]]:
    adjacency: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for left, right, distance in edges:
        adjacency[left].append((right, distance))
        adjacency[right].append((left, distance))
    oriented = []
    seen = {root}
    queue = deque([(root, 0)])
    while queue:
        parent, depth = queue.popleft()
        for child, distance in sorted(adjacency[parent], key=lambda item: (item[1], item[0])):
            if child in seen:
                continue
            seen.add(child)
            oriented.append((parent, child, distance, depth + 1))
            queue.append((child, depth + 1))
    if len(seen) != len(adjacency):
        raise RuntimeError("Rooted tree is disconnected")
    return oriented


def newick_from_oriented(edges: list[tuple[str, str, float, int]], root: str) -> str:
    children: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for parent, child, distance, _ in edges:
        children[parent].append((child, distance))

    def rec(node: str) -> str:
        label = node.replace(" ", "_").replace(":", "_").replace(";", "_")
        if not children[node]:
            return label
        parts = [f"{rec(child)}:{distance:.8f}" for child, distance in sorted(children[node])]
        return f"({','.join(parts)}){label}"

    return rec(root) + ";"


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing output directory: {OUT}")
    OUT.mkdir(parents=False)

    assignments = pd.read_csv(INPUT / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    group_map = pd.read_csv(INPUT / "formal_group_to_cnv_clone_v1.tsv.gz", sep="\t")
    patient = pd.read_csv(INPUT / "formal_patient_clone_summary_v1.tsv", sep="\t")
    usage = pd.read_csv(INPUT / "consensus_usage_rank7_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    source_meta = json.loads((INPUT / "run_metadata_v1.json").read_text(encoding="utf-8"))

    if len(assignments) != 11950 or assignments["exact_cell_id"].duplicated().any():
        raise RuntimeError("Formal clone assignments fail 11,950-cell identity gate")
    if len(usage) != 11950 or usage["exact_cell_id"].duplicated().any() or set(usage["exact_cell_id"]) != set(assignments["exact_cell_id"]):
        raise RuntimeError("Rank-7 usage and formal clone assignments do not share exact cell identities")
    if source_meta["selected_distance_threshold"] != 0.25 or source_meta["primary_clone_minimum_cells"] != 20:
        raise RuntimeError("Unexpected formal clone definition")
    eligible_units = patient.loc[patient["evaluable_within_patient_two_primary_clones_v1"].eq(True), "infercnv_unit"].astype(str).tolist()
    if len(eligible_units) != 14:
        raise RuntimeError(f"Expected 14 evaluable patients, observed {len(eligible_units)}")
    if set(path.name for path in HMM.iterdir() if path.is_dir()) != set(eligible_units):
        raise RuntimeError("HMM snapshot patient set differs from formal evaluable-patient set")

    input_rows = []
    for path in sorted(p for p in INPUT.rglob("*") if p.is_file()):
        input_rows.append({"relative_path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(input_rows), OUT / "input_manifest_v1.tsv")

    profile_records = []
    profile_store: dict[tuple[str, str], np.ndarray] = {}
    bins_store: dict[str, pd.DataFrame] = {}
    formal_subset = group_map.loc[group_map["infercnv_unit"].isin(eligible_units)].copy()
    for unit in sorted(eligible_units):
        genes_path = HMM / unit / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat"
        pred_path = HMM / unit / "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"
        genes = pd.read_csv(genes_path, sep="\t", index_col=0).reset_index(names="gene")
        all_gene_names = set(genes["gene"])
        genes["chr_number"] = genes["chr"].map(chr_number)
        genes = genes.loc[genes["chr_number"].notna()].copy()
        genes["chr_number"] = genes["chr_number"].astype(int)
        genes["midpoint"] = (pd.to_numeric(genes["start"]) + pd.to_numeric(genes["stop"])) / 2
        genes["bin_index"] = np.floor((genes["midpoint"] - 1) / BIN_SIZE).astype(int)
        genes["bin_key"] = genes["chr_number"].astype(str) + ":" + genes["bin_index"].astype(str)
        genes = genes.sort_values(["chr_number", "bin_index", "midpoint", "gene"], kind="stable").reset_index(drop=True)
        bin_codes, bin_names = pd.factorize(genes["bin_key"], sort=False)
        bin_meta = genes.groupby("bin_key", sort=False).agg(chr_number=("chr_number", "first"), bin_index=("bin_index", "first"), start_bp=("start", "min"), stop_bp=("stop", "max"), n_genes=("gene", "size")).reset_index()
        bins_store[unit] = bin_meta
        pred = pd.read_csv(pred_path, sep="\t", dtype={"cell_group_name": str, "gene": str, "state": float})
        if set(pred["gene"]) - all_gene_names or pred[["cell_group_name", "gene"]].duplicated().any() or not pred["state"].isin(range(1, 7)).all():
            raise RuntimeError(f"{unit}: invalid HMM prediction table")
        gene_index = pd.Series(np.arange(len(genes)), index=genes["gene"])
        expected_groups = formal_subset.loc[formal_subset["infercnv_unit"].eq(unit)].sort_values("infercnv_group")
        for row in expected_groups.itertuples(index=False):
            altered = pred.loc[pred["cell_group_name"].eq(row.infercnv_group), ["gene", "state"]]
            altered = altered.loc[altered["gene"].isin(gene_index.index)]
            state = np.full(len(genes), NEUTRAL, dtype=np.float64)
            state[gene_index.loc[altered["gene"]].to_numpy(int)] = altered["state"].to_numpy(float)
            profile = pd.Series(state).groupby(bin_codes, sort=False).median().to_numpy(float)
            observed_hash = profile_hash(profile)
            hash_match = observed_hash == row.profile_sha256
            if len(profile) != int(row.n_5mb_bins) or not hash_match:
                raise RuntimeError(f"{unit}/{row.infercnv_group}: profile reproduction failed")
            profile_store[(unit, row.infercnv_group)] = profile
            profile_records.append({
                "infercnv_unit": unit,
                "infercnv_group": row.infercnv_group,
                "n_5mb_bins": len(profile),
                "source_profile_sha256": row.profile_sha256,
                "reproduced_profile_sha256": observed_hash,
                "profile_hash_match": hash_match,
            })

    usage_values = usage[FACTORS].to_numpy(float)
    if np.any(~np.isfinite(usage_values)) or np.any(usage_values < 0):
        raise RuntimeError("Invalid cNMF usage values")
    usage["dominant_factor_v1"] = np.asarray(FACTORS)[np.argmax(usage_values, axis=1)]
    cells = assignments.merge(usage, on="exact_cell_id", how="left", validate="one_to_one")
    tree_cells = cells.loc[
        cells["infercnv_unit"].isin(eligible_units)
        & cells["clone_size_class_v1"].eq("primary_clone_ge20_v1")
    ].copy()

    clone_profiles: dict[tuple[str, str], np.ndarray] = {}
    node_rows = []
    edge_rows = []
    transition_rows = []
    pair_rows = []
    newick_rows = []
    patient_rows = []
    for unit in sorted(eligible_units):
        mapping = formal_subset.loc[
            formal_subset["infercnv_unit"].eq(unit)
            & formal_subset["clone_size_class_v1"].eq("primary_clone_ge20_v1")
        ].sort_values(["formal_inferred_cnv_clone_v1", "infercnv_group"])
        clone_ids = sorted(mapping["formal_inferred_cnv_clone_v1"].unique())
        if len(clone_ids) < 2:
            raise RuntimeError(f"{unit}: fewer than two primary clones")
        for clone_id in clone_ids:
            clone_groups = mapping.loc[mapping["formal_inferred_cnv_clone_v1"].eq(clone_id)]
            matrix = np.vstack([profile_store[(unit, group)] for group in clone_groups["infercnv_group"]])
            weights = clone_groups["raw_group_cells"].to_numpy(int)
            clone_profiles[(unit, clone_id)] = weighted_median_states(matrix, weights)

        profiles = {"NEUTRAL_ROOT": np.full(len(bins_store[unit]), NEUTRAL, dtype=np.float64)}
        profiles.update({clone_id: clone_profiles[(unit, clone_id)] for clone_id in clone_ids})
        undirected = deterministic_mst(profiles)
        oriented = orient_tree(undirected, "NEUTRAL_ROOT")
        newick_rows.append({"infercnv_unit": unit, "newick_v1": newick_from_oriented(oriented, "NEUTRAL_ROOT")})

        depth_map = {"NEUTRAL_ROOT": 0}
        path_map = {"NEUTRAL_ROOT": 0.0}
        for parent_id, child_id, distance, depth in oriented:
            depth_map[child_id] = depth
            path_map[child_id] = path_map[parent_id] + distance
            parent_profile, child_profile = profiles[parent_id], profiles[child_id]
            changed = parent_profile != child_profile
            pshift, cshift = parent_profile - NEUTRAL, child_profile - NEUTRAL
            edge_rows.append({
                "infercnv_unit": unit,
                "parent_node": parent_id,
                "child_node": child_id,
                "child_depth_from_neutral": depth,
                "bin_hamming_fraction": distance,
                "n_changed_5mb_bins": int(changed.sum()),
                "n_child_more_loss_like": int(np.sum(cshift < pshift)),
                "n_child_more_gain_like": int(np.sum(cshift > pshift)),
                "root_anchored_not_temporal": True,
            })
            bin_meta = bins_store[unit]
            for idx in np.flatnonzero(changed):
                transition_rows.append({
                    "infercnv_unit": unit,
                    "parent_node": parent_id,
                    "child_node": child_id,
                    "bin_key": bin_meta.iloc[idx]["bin_key"],
                    "chr_number": int(bin_meta.iloc[idx]["chr_number"]),
                    "bin_index": int(bin_meta.iloc[idx]["bin_index"]),
                    "parent_hmm_state": parent_profile[idx],
                    "child_hmm_state": child_profile[idx],
                    "state_delta_from_parent": child_profile[idx] - parent_profile[idx],
                })

        for i, left in enumerate(clone_ids):
            for right in clone_ids[i + 1:]:
                pair_rows.append({
                    "infercnv_unit": unit,
                    "clone_a": left,
                    "clone_b": right,
                    "bin_hamming_fraction": float(np.mean(profiles[left] != profiles[right])),
                    "mean_absolute_hmm_state_difference": float(np.mean(np.abs(profiles[left] - profiles[right]))),
                })

        unit_cells = tree_cells.loc[tree_cells["infercnv_unit"].eq(unit)]
        total_primary_cells = len(unit_cells)
        sensitive = bool(patient.loc[patient["infercnv_unit"].eq(unit), "threshold_sensitive_patient_v1"].iloc[0])
        node_rows.append({
            "infercnv_unit": unit,
            "node_id": "NEUTRAL_ROOT",
            "node_type": "neutral_reference_anchor",
            "clone_cells": 0,
            "clone_fraction_of_primary_cells": 0.0,
            "tree_depth_from_neutral": 0,
            "cumulative_hamming_from_neutral": 0.0,
            "dominant_cnmf_factor": "NA",
            "dominant_factor_cell_fraction": np.nan,
            "cnmf_usage_entropy": np.nan,
            "threshold_sensitive_patient_v1": sensitive,
            **{f"median_{factor}": np.nan for factor in FACTORS},
        })
        for clone_id in clone_ids:
            subset = unit_cells.loc[unit_cells["formal_inferred_cnv_clone_v1"].eq(clone_id)]
            if len(subset) < 20:
                raise RuntimeError(f"{unit}/{clone_id}: primary clone has fewer than 20 cells")
            medians = subset[FACTORS].median()
            dominant = str(medians.idxmax())
            cell_fraction = float(subset["dominant_factor_v1"].eq(dominant).mean())
            means = subset[FACTORS].mean().to_numpy(float)
            proportions = means / means.sum() if means.sum() > 0 else np.repeat(1 / len(means), len(means))
            entropy = float(-np.sum(np.where(proportions > 0, proportions * np.log(proportions), 0)) / math.log(len(proportions)))
            profile = profiles[clone_id]
            node_rows.append({
                "infercnv_unit": unit,
                "node_id": clone_id,
                "node_type": "formal_primary_inferred_cnv_clone",
                "clone_cells": len(subset),
                "clone_fraction_of_primary_cells": len(subset) / total_primary_cells,
                "tree_depth_from_neutral": depth_map[clone_id],
                "cumulative_hamming_from_neutral": path_map[clone_id],
                "dominant_cnmf_factor": dominant,
                "dominant_factor_cell_fraction": cell_fraction,
                "cnmf_usage_entropy": entropy,
                "threshold_sensitive_patient_v1": sensitive,
                "n_altered_5mb_bins_vs_neutral": int(np.sum(profile != NEUTRAL)),
                "n_loss_5mb_bins_vs_neutral": int(np.sum(profile < NEUTRAL)),
                "n_gain_5mb_bins_vs_neutral": int(np.sum(profile > NEUTRAL)),
                **{f"median_{factor}": float(medians[factor]) for factor in FACTORS},
            })
        factor_counts = pd.Series([row["dominant_cnmf_factor"] for row in node_rows if row["infercnv_unit"] == unit and row["node_type"].startswith("formal")]).value_counts()
        patient_rows.append({
            "infercnv_unit": unit,
            "n_primary_clone_nodes": len(clone_ids),
            "n_primary_clone_cells": total_primary_cells,
            "n_distinct_clone_dominant_factors": int(len(factor_counts)),
            "dominant_factor_clone_counts": ";".join(f"{factor}:{int(count)}" for factor, count in factor_counts.sort_index().items()),
            "total_tree_branch_length": float(sum(edge[2] for edge in oriented)),
            "maximum_tree_depth": max(edge[3] for edge in oriented),
            "minimum_root_to_clone_hamming": float(min(np.mean(profiles[clone] != profiles["NEUTRAL_ROOT"]) for clone in clone_ids)),
            "threshold_sensitive_patient_v1": sensitive,
            "interpretation": "root-anchored inferred-CNV distance tree; not a mutation-validated temporal phylogeny",
        })

    nodes = pd.DataFrame(node_rows)
    edges = pd.DataFrame(edge_rows)
    transitions = pd.DataFrame(transition_rows)
    pairs = pd.DataFrame(pair_rows)
    patient_summary = pd.DataFrame(patient_rows)
    newick = pd.DataFrame(newick_rows)
    reproduction = pd.DataFrame(profile_records)

    factor_summary = (
        nodes.loc[nodes["node_type"].eq("formal_primary_inferred_cnv_clone")]
        .groupby("dominant_cnmf_factor", observed=True)
        .agg(n_clone_nodes=("node_id", "size"), n_patients=("infercnv_unit", "nunique"), n_cells=("clone_cells", "sum"), median_node_purity=("dominant_factor_cell_fraction", "median"))
        .reset_index()
    )
    atomic_tsv(reproduction, OUT / "m03_hmm_profile_reproduction_v1.tsv")
    atomic_tsv(nodes, OUT / "m03_tree_nodes_cnmf_annotation_v1.tsv")
    atomic_tsv(edges, OUT / "m03_tree_edges_v1.tsv")
    atomic_tsv(transitions, OUT / "m03_tree_edge_bin_transitions_v1.tsv.gz", compression="gzip")
    atomic_tsv(pairs, OUT / "m03_within_patient_clone_pair_distances_v1.tsv")
    atomic_tsv(patient_summary, OUT / "m03_patient_tree_summary_v1.tsv")
    atomic_tsv(newick, OUT / "m03_patient_tree_newick_v1.tsv")
    atomic_tsv(factor_summary, OUT / "m03_clone_dominant_cnmf_factor_summary_v1.tsv")

    status = {
        "status": "PASS_M03_ANALYSIS_TABLES_COMPLETE_PRESENTATION_PENDING",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_formal_cells_source": len(assignments),
        "n_evaluable_patients": len(eligible_units),
        "n_primary_clone_nodes": int(nodes["node_type"].eq("formal_primary_inferred_cnv_clone").sum()),
        "n_neutral_root_nodes": int(nodes["node_type"].eq("neutral_reference_anchor").sum()),
        "n_tree_edges": len(edges),
        "n_profile_reproductions": len(reproduction),
        "n_profile_hash_matches": int(reproduction["profile_hash_match"].sum()),
        "n_threshold_sensitive_evaluable_patients": int(patient_summary["threshold_sensitive_patient_v1"].sum()),
        "tree_method": "deterministic patient-specific minimum spanning tree of clone-level weighted-median 5 Mb HMM states with neutral reference anchor",
        "interpretation_boundary": "inferred-CNV distance organization only; not mutation-validated temporal ancestry",
        "figures_created": False,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    atomic_json(status, OUT / "status.json")


if __name__ == "__main__":
    main()
