#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.stats import spearmanr, wilcoxon
from sklearn.neighbors import NearestNeighbors


DATASETS = ["GSE207493", "GSE264586", "GSE294109"]
METHODS = ["slingshot", "geodesic"]
TARGET_STATES = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
PRIMARY = {
    "P01_SPP1_minus_Monocyte_pseudotime": ("state", "SPP1-state Macrophage", "positive"),
    "P02_C1QC_minus_Monocyte_pseudotime": ("state", "C1QC-state Macrophage", "positive"),
    "P03_Other_minus_Monocyte_pseudotime": ("state", "Other Macrophage", "positive"),
    "P04_SPP1_anchor_excluded_marker_trend": ("program", "SPP1_MARKER20_ANCHOR_EXCLUDED", "positive"),
    "P05_C1QC_anchor_excluded_marker_trend": ("program", "C1QC_MARKER20_ANCHOR_EXCLUDED", "positive"),
    "P06_Other_anchor_excluded_marker_trend": ("program", "OTHER_MARKER20_ANCHOR_EXCLUDED", "positive"),
    "P07_Monocyte_program_trend": ("program", "MONOCYTE_PROGRAM", "negative"),
}
SECONDARY = {
    "S01_Immunosuppression_no_IL10_trend": ("program", "IMMUNOSUPPRESSION_NO_IL10", "two_sided"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    p = values.loc[valid].to_numpy(dtype=float)
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result.loc[valid] = restored
    return result


def safe_rho(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if len(x) < 3 or np.unique(x).size < 3 or np.unique(y).size < 3:
        return float("nan")
    return float(spearmanr(x, y).statistic)


def signed_rank(values: np.ndarray) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    nonzero = values[np.abs(values) > 1e-15]
    if not len(nonzero):
        return 0.0, 1.0
    test = wilcoxon(nonzero, zero_method="wilcox", correction=False, alternative="two-sided", method="auto")
    return float(test.statistic), float(test.pvalue)


def geodesic_pseudotime(frame: pd.DataFrame) -> tuple[np.ndarray, dict]:
    pc_columns = [f"PC_{index}" for index in range(1, 31)]
    values = frame[pc_columns].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError("Non-finite PCA")
    n_neighbors = min(31, len(frame))
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean", n_jobs=-1)
    model.fit(values)
    distances, indices = model.kneighbors(values)
    distances = distances[:, 1:]
    indices = indices[:, 1:]
    rows = np.repeat(np.arange(len(frame)), indices.shape[1])
    cols = indices.reshape(-1)
    data = distances.reshape(-1)
    graph = sparse.csr_matrix((data, (rows, cols)), shape=(len(frame), len(frame)))
    transpose = graph.T.tocsr()
    union = graph.maximum(transpose)
    mutual = graph.sign().multiply(transpose.sign())
    mutual_min = graph.minimum(transpose)
    graph = union - union.multiply(mutual) + mutual_min
    root = np.flatnonzero(frame["root_cell"].to_numpy(dtype=bool))
    if len(root) != 1:
        raise RuntimeError(f"Expected one root cell, found {len(root)}")
    distance = np.asarray(dijkstra(graph, directed=False, indices=int(root[0])), dtype=float)
    finite = np.isfinite(distance)
    if finite.sum() < 0.80 * len(distance):
        raise RuntimeError("Less than 80% finite geodesic pseudotime")
    minimum = float(distance[finite].min())
    maximum = float(distance[finite].max())
    if maximum <= minimum:
        raise RuntimeError("Degenerate geodesic pseudotime")
    result = np.full(len(distance), np.nan)
    result[finite] = (distance[finite] - minimum) / (maximum - minimum)
    return result, {
        "cells": len(frame),
        "neighbors": n_neighbors - 1,
        "undirected_edges": int(graph.nnz // 2),
        "finite_cells": int(finite.sum()),
        "disconnected_cells": int((~finite).sum()),
        "root_cell": str(frame.loc[root[0], "cell_id"]),
    }


def unit_effects(frame: pd.DataFrame, method: str) -> list[dict]:
    pseudotime = frame[f"{method}_pseudotime"].to_numpy(dtype=float)
    rows: list[dict] = []
    endpoint_registry = {**PRIMARY, **SECONDARY}
    for endpoint, (kind, target, direction) in endpoint_registry.items():
        tier = "PRIMARY" if endpoint in PRIMARY else "SECONDARY"
        effect = float("nan")
        n_target = 0
        n_monocyte = int((frame["phase8_state"] == "Monocyte").sum())
        if len(frame) >= 50:
            if kind == "state":
                target_mask = frame["phase8_state"].eq(target).to_numpy()
                monocyte_mask = frame["phase8_state"].eq("Monocyte").to_numpy()
                n_target = int(target_mask.sum())
                if n_target >= 10 and n_monocyte >= 10:
                    target_values = pseudotime[target_mask & np.isfinite(pseudotime)]
                    monocyte_values = pseudotime[monocyte_mask & np.isfinite(pseudotime)]
                    if len(target_values) >= 10 and len(monocyte_values) >= 10:
                        effect = float(np.median(target_values) - np.median(monocyte_values))
            else:
                program = frame[target].to_numpy(dtype=float)
                effect = safe_rho(pseudotime, program)
        rows.append({
            "dataset": frame["dataset"].iloc[0],
            "patient_id": frame["patient_id"].iloc[0],
            "sample_id": frame["sample_id"].iloc[0],
            "method": method,
            "endpoint_id": endpoint,
            "tier": tier,
            "expected_direction": direction,
            "n_cells": len(frame),
            "n_monocyte": n_monocyte,
            "n_target": n_target,
            "effect": effect,
        })
    return rows


def summarize(patient_effects: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    keys = ["dataset", "method", "endpoint_id", "tier", "expected_direction"]
    for key, frame in patient_effects.groupby(keys, sort=True):
        values = frame["patient_median_effect"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        statistic, p_value = signed_rank(values)
        positive = int((values > 0).sum())
        negative = int((values < 0).sum())
        rows.append({
            **dict(zip(keys, key)),
            "n_patients": len(values),
            "positive_patients": positive,
            "negative_patients": negative,
            "zero_patients": int((np.abs(values) <= 1e-15).sum()),
            "positive_fraction": positive / (positive + negative) if positive + negative else float("nan"),
            "median_patient_effect": float(np.median(values)) if len(values) else float("nan"),
            "wilcoxon_statistic": statistic,
            "p_value": p_value,
        })
    result = pd.DataFrame(rows)
    result["BH_q_within_dataset_method_primary7"] = np.nan
    for (dataset, method), index in result.loc[result["tier"] == "PRIMARY"].groupby(["dataset", "method"]).groups.items():
        result.loc[index, "BH_q_within_dataset_method_primary7"] = bh(result.loc[index, "p_value"])
    expected_positive = result["expected_direction"].eq("positive") & (result["median_patient_effect"] > 0)
    expected_negative = result["expected_direction"].eq("negative") & (result["median_patient_effect"] < 0)
    direction_ok = expected_positive | expected_negative | result["expected_direction"].eq("two_sided")
    result["formal_support"] = (
        result["tier"].eq("PRIMARY")
        & (result["n_patients"] >= 10)
        & (result["BH_q_within_dataset_method_primary7"] < 0.05)
        & direction_ok
    )
    result["small_n_nominal_support"] = (
        result["tier"].eq("PRIMARY")
        & result["n_patients"].between(5, 9)
        & (result["p_value"] < 0.05)
        & direction_ok
    )
    return result


def write_summary(path: Path, decision: dict, summaries: pd.DataFrame, concordance: pd.DataFrame) -> None:
    lines = [
        "# Phase 87 three-cohort myeloid trajectory",
        "",
        f"- Computational status: **{decision['status']}**.",
        f"- Scientific decision: **{decision['decision']}**.",
        f"- P01 all-three-cohort dual-method support: {decision['p01_all_three_dual_method']}.",
        f"- P04 two-external-cohort dual-method support: {decision['p04_two_external_dual_method']}.",
        f"- Method-concordance gate: {decision['method_concordance_gate']}.",
        "",
        "## Primary endpoint summary",
        "",
        "| dataset | method | endpoint | n | median | p | q | support |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    primary = summaries.loc[summaries["tier"] == "PRIMARY"].sort_values(["endpoint_id", "dataset", "method"])
    for row in primary.itertuples(index=False):
        lines.append(
            f"| {row.dataset} | {row.method} | {row.endpoint_id} | {row.n_patients} | "
            f"{row.median_patient_effect:.4f} | {row.p_value:.4g} | {row.BH_q_within_dataset_method_primary7:.4g} | {row.formal_support} |"
        )
    lines.extend([
        "",
        "## Method concordance",
        "",
        "| dataset | patients | median rho | positive fraction | gate |",
        "|---|---:|---:|---:|---|",
    ])
    for row in concordance.sort_values("dataset").itertuples(index=False):
        lines.append(f"| {row.dataset} | {row.n_patients} | {row.median_patient_rho:.3f} | {row.positive_fraction:.3f} | {row.concordance_gate_pass} |")
    lines.extend([
        "",
        "## Boundary",
        "",
        "These are cross-sectional trajectory orderings. They do not prove lineage conversion, treatment-induced transition, temporal ancestry or causality.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.phase_root.resolve()
    inputs = root / "trajectory_inputs_v1"
    output = root / "results_v1"
    part = root / "results_v1.part"
    if not inputs.is_dir() or output.exists() or part.exists():
        raise RuntimeError("Invalid Phase87 directory state")
    part.mkdir(parents=True)

    all_cells = []
    graph_qc = []
    effect_rows = []
    unit_concordance_rows = []
    for dataset in DATASETS:
        path = inputs / f"{dataset}_phase87_cells_v1.tsv.gz"
        frame = pd.read_csv(path, sep="\t")
        geodesic, qc = geodesic_pseudotime(frame)
        frame["geodesic_pseudotime"] = geodesic
        graph_qc.append({"dataset": dataset, **qc})
        for (_, _), unit in frame.groupby(["patient_id", "sample_id"], sort=True):
            for method in METHODS:
                effect_rows.extend(unit_effects(unit, method))
            correlation = safe_rho(unit["slingshot_pseudotime"].to_numpy(dtype=float), unit["geodesic_pseudotime"].to_numpy(dtype=float))
            unit_concordance_rows.append({
                "dataset": dataset,
                "patient_id": unit["patient_id"].iloc[0],
                "sample_id": unit["sample_id"].iloc[0],
                "n_cells": len(unit),
                "pseudotime_rho": correlation,
            })
        keep = [
            "cell_id", "dataset", "patient_id", "sample_id", "phase8_state", "phase7_cluster",
            "root_cluster", "root_cell", "slingshot_pseudotime", "geodesic_pseudotime",
            "SPP1_MARKER20_ANCHOR_EXCLUDED", "C1QC_MARKER20_ANCHOR_EXCLUDED",
            "OTHER_MARKER20_ANCHOR_EXCLUDED", "MONOCYTE_PROGRAM", "IMMUNOSUPPRESSION_NO_IL10",
        ]
        all_cells.append(frame[keep])

    cells = pd.concat(all_cells, ignore_index=True)
    unit_effects_frame = pd.DataFrame(effect_rows)
    patient_effects = (
        unit_effects_frame.dropna(subset=["effect"])
        .groupby(["dataset", "patient_id", "method", "endpoint_id", "tier", "expected_direction"], sort=True)
        .agg(patient_median_effect=("effect", "median"), n_samples=("sample_id", "nunique"), total_cells=("n_cells", "sum"))
        .reset_index()
    )
    summaries = summarize(patient_effects)

    unit_concordance = pd.DataFrame(unit_concordance_rows)
    patient_concordance = (
        unit_concordance.dropna(subset=["pseudotime_rho"])
        .groupby(["dataset", "patient_id"], sort=True)
        .agg(patient_median_rho=("pseudotime_rho", "median"), n_samples=("sample_id", "nunique"))
        .reset_index()
    )
    concordance_rows = []
    for dataset, frame in patient_concordance.groupby("dataset", sort=True):
        values = frame["patient_median_rho"].to_numpy(dtype=float)
        positive_fraction = float(np.mean(values > 0)) if len(values) else float("nan")
        median = float(np.median(values)) if len(values) else float("nan")
        concordance_rows.append({
            "dataset": dataset,
            "n_patients": len(values),
            "median_patient_rho": median,
            "positive_fraction": positive_fraction,
            "concordance_gate_pass": bool(median > 0 and positive_fraction >= 0.70),
        })
    concordance = pd.DataFrame(concordance_rows)

    def formally_supported(dataset: str, method: str, endpoint: str) -> bool:
        row = summaries.loc[(summaries["dataset"] == dataset) & (summaries["method"] == method) & (summaries["endpoint_id"] == endpoint)]
        return len(row) == 1 and bool(row["formal_support"].iloc[0])

    p01_all = all(formally_supported(dataset, method, "P01_SPP1_minus_Monocyte_pseudotime") for dataset in DATASETS for method in METHODS)
    p04_external = all(formally_supported(dataset, method, "P04_SPP1_anchor_excluded_marker_trend") for dataset in ["GSE264586", "GSE294109"] for method in METHODS)
    concordance_gate = bool(len(concordance) == 3 and concordance["concordance_gate_pass"].all())
    if p01_all and p04_external and concordance_gate:
        scientific_decision = "THREE_COHORT_MONOCYTE_TO_SPP1_ORDERING_REPLICATED"
    elif p01_all and concordance_gate:
        scientific_decision = "STATE_ORDERING_WITHOUT_ANCHOR_EXCLUDED_PROGRAM_REPLICATION"
    else:
        scientific_decision = "MONOCYTE_TO_SPP1_ORDERING_NOT_THREE_COHORT_REPLICATED"

    decision = {
        "status": "PASS_COMPUTATION_COMPLETE",
        "decision": scientific_decision,
        "p01_all_three_dual_method": p01_all,
        "p04_two_external_dual_method": p04_external,
        "method_concordance_gate": concordance_gate,
        "datasets": DATASETS,
        "trajectory_cells": int(len(cells)),
        "patients": int(cells[["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "patient_sample_units": int(cells[["dataset", "patient_id", "sample_id"]].drop_duplicates().shape[0]),
        "preregistration_sha256": sha256(root / "PHASE87_PREREGISTRATION_v1.md"),
        "marker_table_sha256": sha256(root / "inputs_v1" / "phase83_marker_members_v1.tsv"),
    }

    cells.to_csv(part / "phase87_cell_pseudotime_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(graph_qc).to_csv(part / "phase87_geodesic_graph_qc_v1.tsv", sep="\t", index=False)
    unit_effects_frame.to_csv(part / "phase87_patient_sample_effects_v1.tsv", sep="\t", index=False)
    patient_effects.to_csv(part / "phase87_patient_effects_v1.tsv", sep="\t", index=False)
    summaries.to_csv(part / "phase87_dataset_method_endpoint_summary_v1.tsv", sep="\t", index=False)
    unit_concordance.to_csv(part / "phase87_unit_method_concordance_v1.tsv", sep="\t", index=False)
    patient_concordance.to_csv(part / "phase87_patient_method_concordance_v1.tsv", sep="\t", index=False)
    concordance.to_csv(part / "phase87_dataset_method_concordance_summary_v1.tsv", sep="\t", index=False)
    with (part / "PHASE87_DECISION_v1.json").open("w", encoding="utf-8") as handle:
        json.dump(decision, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    write_summary(part / "PHASE87_TRAJECTORY_SUMMARY_v1.md", decision, summaries, concordance)
    (part / "PHASE87_COMPUTATION_PASS_v1").write_text("PASS_COMPUTATION_COMPLETE\n", encoding="utf-8")
    manifest_rows = []
    for path in sorted(part.iterdir()):
        if path.name == "FILE_MANIFEST_SHA256_v1.tsv":
            continue
        manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(part / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    os.replace(part, output)
    print(json.dumps(decision, ensure_ascii=False))


if __name__ == "__main__":
    main()
