#!/usr/bin/env python3
"""Rebuild Phase 17 abundance/composition/niche tables, separately per section."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.spatial import cKDTree
import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


EXPECTED_STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]
EXPECTED_EPOCH = 10_000
K_RANGE = range(3, 9)
K_NEIGHBORS = 6
SEED = 20260813


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    part = Path(str(path) + ".part")
    with part.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(part, path)


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = Path(str(path) + ".part")
    if path.suffix == ".gz":
        with gzip.open(part, "wt", encoding="utf-8", newline="") as handle:
            frame.to_csv(handle, sep="\t", index=False, lineterminator="\n")
    else:
        frame.to_csv(part, sep="\t", index=False, lineterminator="\n")
    os.replace(part, path)


def clean_state(value: str) -> str:
    prefix = "q05cell_abundance_w_sf_"
    value = str(value)
    return value[len(prefix):] if value.startswith(prefix) else value


def discover_samples(root: Path) -> list[dict]:
    records = []
    cohorts = [
        ("GLMF", root / "glmf12_cell2location_tymp_continuous_v1"),
        ("TJRCC", root / "tjrcc12_cell2location_tymp_continuous_v1"),
    ]
    for cohort, cohort_root in cohorts:
        for status_path in sorted(cohort_root.glob("*/status.json")):
            status = json.loads(status_path.read_text(encoding="utf-8"))
            sample = str(status.get("sample", status_path.parent.name))
            if status.get("status") != "completed" or int(status.get("epoch", -1)) != EXPECTED_EPOCH:
                raise RuntimeError(f"Completion gate failed: {status_path}")
            h5ad = Path(status["output_h5ad"])
            abundance = status_path.parent / "cell_abundance_q05_tymp_continuous_v1.tsv.gz"
            if not h5ad.is_file() or not abundance.is_file():
                raise FileNotFoundError(f"Missing Phase 17 output for {sample}")
            if h5ad.stat().st_size != int(status["output_size_bytes"]):
                raise RuntimeError(f"H5AD size mismatch for {sample}")
            actual_hash = sha256(h5ad)
            if actual_hash != status["output_sha256"]:
                raise RuntimeError(f"H5AD SHA256 mismatch for {sample}")
            records.append({
                "cohort": cohort, "sample": sample, "status_path": status_path,
                "h5ad": h5ad, "abundance": abundance, "status": status,
                "h5ad_sha256": actual_hash,
            })
    if len(records) != 14:
        raise RuntimeError(f"Expected 14 completed samples, found {len(records)}")
    if len({r["sample"] for r in records}) != 14:
        raise RuntimeError("Sample names are not globally unique")
    return records


def load_sample(record: dict) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    abundance = pd.read_csv(record["abundance"], sep="\t", index_col=0)
    abundance.index = abundance.index.astype(str)
    abundance.columns = [clean_state(x) for x in abundance.columns]
    if abundance.columns.tolist() != EXPECTED_STATES:
        raise RuntimeError(f"Unexpected state names/order for {record['sample']}: {abundance.columns.tolist()}")
    values = abundance.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError(f"Invalid abundance values for {record['sample']}")
    totals = values.sum(axis=1)
    if (totals <= 0).any():
        raise RuntimeError(f"Non-positive total abundance for {record['sample']}")

    obj = ad.read_h5ad(record["h5ad"], backed="r")
    try:
        if obj.n_obs != len(abundance) or obj.n_vars != 8830:
            raise RuntimeError(f"Unexpected H5AD shape for {record['sample']}: {obj.shape}")
        if obj.obs_names.astype(str).tolist() != abundance.index.tolist():
            raise RuntimeError(f"H5AD/TSV spot order mismatch for {record['sample']}")
        if "spatial" not in obj.obsm or "q05_cell_abundance_w_sf" not in obj.obsm:
            raise RuntimeError(f"Required obsm missing for {record['sample']}")
        coordinates = np.asarray(obj.obsm["spatial"], dtype=float)[:, :2]
        h5_values = np.asarray(obj.obsm["q05_cell_abundance_w_sf"], dtype=float)
        if h5_values.shape != values.shape or not np.allclose(h5_values, values, rtol=2e-6, atol=2e-7):
            raise RuntimeError(f"H5AD/TSV abundance mismatch for {record['sample']}")
    finally:
        obj.file.close()
    if coordinates.shape != (len(abundance), 2) or not np.isfinite(coordinates).all():
        raise RuntimeError(f"Invalid spatial coordinates for {record['sample']}")
    return abundance, values, coordinates


def neighborhood_graph(coordinates: np.ndarray) -> sparse.csr_matrix:
    k = min(K_NEIGHBORS, len(coordinates) - 1)
    neighbors = cKDTree(coordinates).query(coordinates, k=k + 1)[1][:, 1:]
    rows = np.repeat(np.arange(len(coordinates)), k)
    cols = neighbors.reshape(-1)
    graph = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(coordinates),) * 2)
    graph = graph.maximum(graph.T)
    graph.data[:] = 1.0
    return graph


def select_niches(features: np.ndarray, sample_index: int) -> tuple[np.ndarray, pd.DataFrame]:
    seed = SEED + sample_index * 1000
    minimum = max(25, int(np.ceil(0.01 * len(features))))
    rows = []
    labels = {}
    for k in K_RANGE:
        model = KMeans(n_clusters=k, n_init=50, max_iter=1000, random_state=seed + k)
        current = model.fit_predict(features)
        labels[k] = current
        if len(features) > 2500:
            rng = np.random.default_rng(seed + k)
            idx = np.sort(rng.choice(len(features), 2500, replace=False))
            silhouette = float(silhouette_score(features[idx], current[idx]))
        else:
            silhouette = float(silhouette_score(features, current))
        smallest = int(np.bincount(current, minlength=k).min())
        rows.append({"k": k, "silhouette": silhouette, "minimum_niche_spots": smallest,
                     "minimum_required": minimum, "valid": smallest >= minimum})
    table = pd.DataFrame(rows)
    valid = table[table["valid"]].sort_values(["silhouette", "k"], ascending=[False, True])
    if valid.empty:
        raise RuntimeError("No K in 3..8 passed the minimum niche-size gate")
    selected = int(valid.iloc[0]["k"])
    table["selected"] = table["k"].eq(selected)
    raw = labels[selected]
    old, counts = np.unique(raw, return_counts=True)
    order = old[np.argsort(-counts, kind="stable")]
    mapping = {int(old_id): new_id + 1 for new_id, old_id in enumerate(order)}
    relabeled = np.array([mapping[int(x)] for x in raw], dtype=int)
    return relabeled, table


def process_sample(record: dict, sample_index: int) -> dict[str, pd.DataFrame]:
    abundance, values, coordinates = load_sample(record)
    totals = values.sum(axis=1)
    fractions = values / totals[:, None]
    graph = neighborhood_graph(coordinates)
    graph_self = graph + sparse.eye(len(abundance), format="csr")
    neighborhood = np.asarray(graph_self @ values)
    neighborhood_composition = neighborhood / neighborhood.sum(axis=1, keepdims=True)
    features = np.sqrt(neighborhood_composition)
    niche_id, k_table = select_niches(features, sample_index)
    sample = record["sample"]
    cohort = record["cohort"]
    niche = np.array([f"{sample}-N{x}" for x in niche_id], dtype=object)
    dominant_idx = values.argmax(axis=1)

    wide = pd.DataFrame({
        "cohort": cohort, "sample": sample, "spot_id": abundance.index,
        "x": coordinates[:, 0], "y": coordinates[:, 1], "total_q05_abundance": totals,
        "dominant_state": np.asarray(EXPECTED_STATES)[dominant_idx], "niche": niche,
    })
    for j, state in enumerate(EXPECTED_STATES):
        wide[f"fraction__{state}"] = fractions[:, j]

    state_rows = []
    for j, state in enumerate(EXPECTED_STATES):
        state_rows.append({
            "cohort": cohort, "sample": sample, "state": state, "n_spots": len(abundance),
            "total_q05_abundance": float(values[:, j].sum()),
            "mean_q05_abundance": float(values[:, j].mean()),
            "median_q05_abundance": float(np.median(values[:, j])),
            "mean_fraction": float(fractions[:, j].mean()),
            "median_fraction": float(np.median(fractions[:, j])),
            "dominant_spots": int((dominant_idx == j).sum()),
        })

    global_comp = neighborhood_composition.mean(axis=0)
    niche_summary = []
    niche_comp = []
    for current in range(1, int(niche_id.max()) + 1):
        mask = niche_id == current
        label = f"{sample}-N{current}"
        mean_comp = neighborhood_composition[mask].mean(axis=0)
        enrichment = np.log2((mean_comp + 1e-8) / (global_comp + 1e-8))
        top = np.argsort(-enrichment)[:3]
        niche_summary.append({
            "cohort": cohort, "sample": sample, "niche": label,
            "n_spots": int(mask.sum()), "spot_fraction": float(mask.mean()),
            "top1_state": EXPECTED_STATES[top[0]], "top1_log2_enrichment": float(enrichment[top[0]]),
            "top2_state": EXPECTED_STATES[top[1]], "top2_log2_enrichment": float(enrichment[top[1]]),
            "top3_state": EXPECTED_STATES[top[2]], "top3_log2_enrichment": float(enrichment[top[2]]),
        })
        for j, state in enumerate(EXPECTED_STATES):
            niche_comp.append({
                "cohort": cohort, "sample": sample, "niche": label, "state": state,
                "mean_neighborhood_fraction": float(mean_comp[j]),
                "log2_enrichment": float(enrichment[j]),
                "mean_spot_q05_abundance": float(values[mask, j].mean()),
            })

    edges = sparse.triu(graph, k=1).tocoo()
    adjacency_rows = []
    counts = np.bincount(niche_id, minlength=int(niche_id.max()) + 1)[1:]
    frequency = counts / counts.sum()
    for a in range(1, len(counts) + 1):
        for b in range(a, len(counts) + 1):
            observed = int(np.sum(((niche_id[edges.row] == a) & (niche_id[edges.col] == b)) |
                                  ((niche_id[edges.row] == b) & (niche_id[edges.col] == a))))
            expected_fraction = frequency[a - 1] ** 2 if a == b else 2 * frequency[a - 1] * frequency[b - 1]
            expected = float(len(edges.row) * expected_fraction)
            adjacency_rows.append({
                "cohort": cohort, "sample": sample, "niche_1": f"{sample}-N{a}",
                "niche_2": f"{sample}-N{b}", "observed_edges": observed,
                "expected_edges_independent_labels": expected,
                "observed_over_expected": float(observed / expected) if expected > 0 else np.nan,
                "interpretation": "descriptive_not_inferential",
            })

    k_table.insert(0, "sample", sample)
    k_table.insert(0, "cohort", cohort)
    provenance = pd.DataFrame([{
        "cohort": cohort, "sample": sample, "status_path": str(record["status_path"]),
        "input_h5ad": str(record["h5ad"]), "input_h5ad_size_bytes": record["h5ad"].stat().st_size,
        "input_h5ad_sha256": record["h5ad_sha256"], "input_abundance_tsv": str(record["abundance"]),
        "input_abundance_sha256": sha256(record["abundance"]), "n_spots": len(abundance),
        "n_genes": 8830, "n_states": len(EXPECTED_STATES), "epoch": EXPECTED_EPOCH,
        "median_symmetric_knn_degree": float(np.median(np.asarray(graph.sum(axis=1)).ravel())),
        "selected_k": int(k_table.loc[k_table["selected"], "k"].iloc[0]),
    }])
    return {
        "wide": wide, "state_summary": pd.DataFrame(state_rows),
        "niche_summary": pd.DataFrame(niche_summary), "niche_composition": pd.DataFrame(niche_comp),
        "adjacency": pd.DataFrame(adjacency_rows), "k_selection": k_table, "provenance": provenance,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase17-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    phase17 = args.phase17_root.resolve(strict=True)
    output = args.output_root.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output}")
    output.mkdir(parents=True)
    tables = output / "tables"
    tables.mkdir()
    atomic_json(output / "RUNNING.json", {"status": "RUNNING", "started_at": now_iso()})

    records = discover_samples(phase17)
    collected: dict[str, list[pd.DataFrame]] = {key: [] for key in (
        "wide", "state_summary", "niche_summary", "niche_composition", "adjacency", "k_selection", "provenance"
    )}
    for index, record in enumerate(records):
        print(f"[{now_iso()}] {record['cohort']} {record['sample']}", flush=True)
        result = process_sample(record, index)
        for key, value in result.items():
            collected[key].append(value)

    outputs = {
        "spot_composition_wide.tsv.gz": pd.concat(collected["wide"], ignore_index=True),
        "sample_state_summary.tsv": pd.concat(collected["state_summary"], ignore_index=True),
        "niche_summary.tsv": pd.concat(collected["niche_summary"], ignore_index=True),
        "niche_state_composition.tsv": pd.concat(collected["niche_composition"], ignore_index=True),
        "niche_adjacency_descriptive.tsv": pd.concat(collected["adjacency"], ignore_index=True),
        "niche_k_selection.tsv": pd.concat(collected["k_selection"], ignore_index=True),
        "input_provenance.tsv": pd.concat(collected["provenance"], ignore_index=True),
    }
    for name, frame in outputs.items():
        atomic_tsv(frame, tables / name)

    software = {
        "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__,
        "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__,
        "script_sha256": sha256(Path(__file__)),
    }
    atomic_json(output / "software_versions.json", software)
    manifest_rows = []
    for path in sorted(p for p in output.rglob("*") if p.is_file() and p.name != "RUNNING.json"):
        manifest_rows.append({"file": path.relative_to(output).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(pd.DataFrame(manifest_rows), output / "OUTPUT_MANIFEST.tsv")
    status = {
        "schema_version": "phase18_composition_niche_v1", "status": "PASS_ANALYSIS_COMPLETE",
        "completed_at": now_iso(), "n_sections": 14, "n_glmf": 2, "n_tjrcc": 12,
        "n_spots": int(len(outputs["spot_composition_wide.tsv.gz"])), "n_states": 13,
        "n_niches": int(len(outputs["niche_summary.tsv"])),
        "method": "section-specific symmetric 6-NN plus self neighborhood composition; Hellinger transform; KMeans K=3..8",
        "inference_policy": "descriptive only; sections are biological units; spots are not replicates",
    }
    atomic_json(output / "PHASE18_COMPOSITION_NICHE_STATUS_v1.json", status)
    (output / "RUNNING.json").unlink()
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
