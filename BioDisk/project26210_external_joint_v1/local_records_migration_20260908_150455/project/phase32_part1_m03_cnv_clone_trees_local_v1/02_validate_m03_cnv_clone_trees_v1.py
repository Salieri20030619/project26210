#!/usr/bin/env python3
"""Independent structural and numerical validation for Phase 32 M03 tables."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_local_v1"
OUT = ROOT / "validation_local_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUT.exists():
        raise RuntimeError(f"Refusing to overwrite existing validation directory: {OUT}")
    OUT.mkdir(parents=False)
    checks: list[tuple[str, bool, str]] = []

    status = json.loads((RESULTS / "status.json").read_text(encoding="utf-8"))
    nodes = pd.read_csv(RESULTS / "m03_tree_nodes_cnmf_annotation_v1.tsv", sep="\t")
    edges = pd.read_csv(RESULTS / "m03_tree_edges_v1.tsv", sep="\t")
    reproduction = pd.read_csv(RESULTS / "m03_hmm_profile_reproduction_v1.tsv", sep="\t")
    patients = pd.read_csv(RESULTS / "m03_patient_tree_summary_v1.tsv", sep="\t")
    pairs = pd.read_csv(RESULTS / "m03_within_patient_clone_pair_distances_v1.tsv", sep="\t")
    newick = pd.read_csv(RESULTS / "m03_patient_tree_newick_v1.tsv", sep="\t")
    factors = pd.read_csv(RESULTS / "m03_clone_dominant_cnmf_factor_summary_v1.tsv", sep="\t")
    source_assign = pd.read_csv(ROOT / "input_snapshot_local_v1" / "formal_cell_cnv_clone_assignments_v1.tsv.gz", sep="\t")
    expected_tree_source = source_assign.loc[source_assign["eligible_for_primary_within_patient_clone_analysis_v1"].eq(True)]
    expected_clone_nodes = int(expected_tree_source["formal_inferred_cnv_clone_v1"].nunique())

    checks.append(("analysis_status", status["status"] == "PASS_M03_ANALYSIS_TABLES_COMPLETE_PRESENTATION_PENDING", status["status"]))
    checks.append(("source_cell_count", status["n_formal_cells_source"] == 11950 and len(source_assign) == 11950, str(len(source_assign))))
    checks.append(("patient_count", len(patients) == 14 and patients["infercnv_unit"].nunique() == 14, str(len(patients))))
    clone_nodes = nodes.loc[nodes["node_type"].eq("formal_primary_inferred_cnv_clone")]
    root_nodes = nodes.loc[nodes["node_type"].eq("neutral_reference_anchor")]
    checks.append(("primary_clone_node_count", len(clone_nodes) == expected_clone_nodes and clone_nodes["node_id"].nunique() == expected_clone_nodes, f"observed={len(clone_nodes)}; expected={expected_clone_nodes}"))
    checks.append(("one_root_per_patient", len(root_nodes) == 14 and root_nodes["infercnv_unit"].nunique() == 14, str(len(root_nodes))))
    checks.append(("edge_count", len(edges) == expected_clone_nodes, f"observed={len(edges)}; expected={expected_clone_nodes}"))
    checks.append(("profile_hash_reproduction", len(reproduction) > 0 and reproduction["profile_hash_match"].eq(True).all(), f"{int(reproduction['profile_hash_match'].sum())}/{len(reproduction)}"))
    checks.append(("profile_hash_unique_keys", not reproduction.duplicated(["infercnv_unit", "infercnv_group"]).any(), str(len(reproduction))))
    checks.append(("primary_minimum_20", clone_nodes["clone_cells"].ge(20).all(), str(int(clone_nodes["clone_cells"].min()))))
    checks.append(("node_fraction_conservation", all(np.isclose(frame["clone_fraction_of_primary_cells"].sum(), 1.0) for _, frame in clone_nodes.groupby("infercnv_unit")), "all patient sums=1"))
    checks.append(("cnmf_factor_labels", set(clone_nodes["dominant_cnmf_factor"]).issubset({f"F{i}" for i in range(1, 8)}), ",".join(sorted(clone_nodes["dominant_cnmf_factor"].unique()))))
    checks.append(("cnmf_median_finite", np.isfinite(clone_nodes[[f"median_F{i}" for i in range(1, 8)]].to_numpy(float)).all(), "all finite"))
    checks.append(("factor_summary_conservation", int(factors["n_clone_nodes"].sum()) == expected_clone_nodes and int(factors["n_cells"].sum()) == int(clone_nodes["clone_cells"].sum()), f"nodes={int(factors['n_clone_nodes'].sum())}; expected={expected_clone_nodes}"))
    checks.append(("pair_count", len(pairs) == int(sum(n * (n - 1) / 2 for n in patients["n_primary_clone_nodes"])), str(len(pairs))))
    checks.append(("distance_bounds", edges["bin_hamming_fraction"].between(0, 1).all() and pairs["bin_hamming_fraction"].between(0, 1).all(), "all within [0,1]"))
    checks.append(("newick_rows", len(newick) == 14 and newick["newick_v1"].str.endswith(";").all(), str(len(newick))))
    checks.append(("threshold_sensitive_count", int(patients["threshold_sensitive_patient_v1"].sum()) == 3, str(int(patients["threshold_sensitive_patient_v1"].sum()))))

    topology_ok = True
    topology_detail = []
    for unit, patient_nodes in nodes.groupby("infercnv_unit"):
        ids = set(patient_nodes["node_id"])
        patient_edges = edges.loc[edges["infercnv_unit"].eq(unit)]
        adjacency: dict[str, set[str]] = defaultdict(set)
        parent_count = defaultdict(int)
        for row in patient_edges.itertuples(index=False):
            adjacency[row.parent_node].add(row.child_node)
            adjacency[row.child_node].add(row.parent_node)
            parent_count[row.child_node] += 1
        seen = {"NEUTRAL_ROOT"}
        queue = deque(["NEUTRAL_ROOT"])
        while queue:
            current = queue.popleft()
            for nxt in adjacency[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        ok = len(patient_edges) == len(ids) - 1 and seen == ids and all(parent_count[node] == 1 for node in ids - {"NEUTRAL_ROOT"})
        topology_ok &= ok
        topology_detail.append(f"{unit}:{ok}")
    checks.append(("tree_topology_connected_acyclic_rooted", topology_ok, ";".join(topology_detail)))
    checks.append(("no_figures", status["figures_created"] is False and not any(ROOT.rglob("*.png")) and not any(ROOT.rglob("*.pdf")), "no PNG/PDF"))
    checks.append(("interpretation_guard", "not mutation-validated temporal ancestry" in status["interpretation_boundary"], status["interpretation_boundary"]))

    failed = [name for name, passed, _ in checks if not passed]
    check_path = OUT / "verification_checks_v1.tsv.part"
    with check_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "pass", "detail"])
        writer.writerows((name, str(passed), detail) for name, passed, detail in checks)
    check_path.replace(OUT / "verification_checks_v1.tsv")

    hashes = {path.name: sha256(path) for path in sorted(RESULTS.iterdir()) if path.is_file()}
    payload = {
        "status": "PASS_INDEPENDENT_VALIDATION" if not failed else "FAIL_INDEPENDENT_VALIDATION",
        "checks_passed": len(checks) - len(failed),
        "checks_total": len(checks),
        "failed_checks": failed,
        "analysis_file_sha256": hashes,
    }
    part = OUT / "status.json.part"
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    part.replace(OUT / "status.json")
    if failed:
        raise RuntimeError(f"Validation failed: {failed}")


if __name__ == "__main__":
    main()
