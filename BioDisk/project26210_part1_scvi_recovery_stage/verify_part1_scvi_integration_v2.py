#!/usr/bin/env python3
"""Independent verification v2 for the Part 1 scVI integration v1 output.

The analytical run used ``logger=False`` with scvi-tools 1.4.2.  In that
configuration scvi-tools intentionally leaves ``model.history`` empty, so the
v1 verifier's requirement for a populated model-history table was impossible
even though Lightning completed all 31 epochs.  This verifier preserves that
failed v1 audit and independently closes the training-completion check from
the immutable runtime log.  It does not modify or retrain the model.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT = PROJECT / "results/part1_rebuild_v1/01_raw_counts_qc/qc_v1/part1_raw_counts_qc_filtered_v1.h5ad"
OLD = PROJECT / "object/news_scVI_adata.h5ad"
ROOT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/scvi_integration_v1"
OUTPUT = ROOT / "part1_scvi_integrated_v1.h5ad"
VERIFY = ROOT / "verification_v2"
RUNTIME_LOG = PROJECT / "results/part1_rebuild_v1/logs/scvi_integration_v1.log"
RUN_SCRIPT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/run_part1_scvi_integration_v1.py"
GATE = PROJECT / "check_memory_gate.sh"


def sha256(path: Path, block: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(block):
            h.update(chunk)
    return h.hexdigest()


def dataset_hash(ds: h5py.Dataset, block: int = 10_000_000) -> str:
    h = hashlib.sha256()
    for start in range(0, ds.shape[0], block):
        h.update(np.asarray(ds[start : start + block]).tobytes())
    return h.hexdigest()


def main() -> None:
    VERIFY.mkdir(parents=True, exist_ok=True)
    gate = subprocess.run([str(GATE)], text=True, capture_output=True)
    gate_line = (gate.stdout or gate.stderr).strip()
    print(gate_line, flush=True)
    if gate.returncode != 0:
        raise RuntimeError("memory gate blocked verification")

    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})
        print(f"CHECK {name} passed={bool(passed)} detail={detail}", flush=True)

    check("run_complete", (ROOT / "RUN_COMPLETE").exists(), ROOT / "RUN_COMPLETE")
    check("output_h5ad_exists", OUTPUT.exists() and OUTPUT.stat().st_size > 0, OUTPUT.stat().st_size if OUTPUT.exists() else 0)
    metadata = json.loads((ROOT / "run_metadata_v1.json").read_text())
    status = json.loads((ROOT / "status.json").read_text())
    check("metadata_passed", metadata.get("status") == "passed", metadata.get("status"))
    check("status_completed", status.get("status") == "completed", status.get("status"))
    check("frozen_scope", metadata.get("scope") == "frozen existing 75 libraries; no external data added", metadata.get("scope"))
    check("legacy_embeddings_not_reused", metadata.get("legacy_embeddings_reused") is False, metadata.get("legacy_embeddings_reused"))
    check("legacy_annotations_not_included", metadata.get("legacy_annotations_included") is False, metadata.get("legacy_annotations_included"))

    inp = ad.read_h5ad(INPUT, backed="r")
    out = ad.read_h5ad(OUTPUT, backed="r")
    check("shape", out.shape == (258120, 37795), out.shape)
    check("cell_ids_exact", np.array_equal(inp.obs_names.to_numpy(), out.obs_names.to_numpy()), f"n={out.n_obs}")
    check("gene_ids_exact", np.array_equal(inp.var_names.to_numpy(), out.var_names.to_numpy()), f"n={out.n_vars}")
    check("cell_ids_unique", out.obs_names.is_unique, out.obs_names.nunique())
    check("gene_ids_unique", out.var_names.is_unique, out.var_names.nunique())
    check("library_count", out.obs["sample_id"].nunique() == 75, out.obs["sample_id"].nunique())
    check("patient_count", out.obs["curated_patient_id"].nunique() == 40, out.obs["curated_patient_id"].nunique())
    check("study_count", out.obs["study"].nunique() == 5, out.obs["study"].nunique())
    check("disease_count", out.obs["Disease"].nunique() == 2, out.obs["Disease"].value_counts().to_dict())
    check("all_cells_kept", bool(out.obs["analysis_keep"].all()), int(out.obs["analysis_keep"].sum()))
    check("no_final_annotation_leakage", "cell_type" not in out.obs.columns, list(out.obs.columns))

    hvg = out.var["highly_variable"].to_numpy()
    check("hvg_count", int(hvg.sum()) == 8000, int(hvg.sum()))
    hvg_table = pd.read_csv(ROOT / "tables/highly_variable_genes_v1.tsv.gz", sep="\t")
    check("hvg_table_rows", len(hvg_table) == 37795, len(hvg_table))
    check("hvg_table_exact", np.array_equal(hvg_table["gene"].astype(str), out.var_names.astype(str)), "gene order")
    check("hvg_batch_key", metadata["hvg"] == {"flavor": "seurat_v3", "n_top_genes": 8000, "batch_key": "curated_patient_id"}, metadata["hvg"])

    batch = np.asarray(out.obs["_scvi_batch"]).astype(int)
    mapping = pd.DataFrame({"patient": out.obs["curated_patient_id"].astype(str), "batch": batch}).drop_duplicates()
    check("scvi_batch_codes", len(np.unique(batch)) == 40 and mapping.shape[0] == 40, f"codes={len(np.unique(batch))};mapping={mapping.shape[0]}")
    check("scvi_latent_shape", out.obsm["X_scVI"].shape == (258120, 30), out.obsm["X_scVI"].shape)
    z = np.asarray(out.obsm["X_scVI"])
    check("scvi_latent_finite", np.isfinite(z).all(), f"min={z.min():.6g};max={z.max():.6g}")
    check("pca_shape", out.obsm["X_pca"].shape == (258120, 50), out.obsm["X_pca"].shape)
    check("pca_finite", np.isfinite(np.asarray(out.obsm["X_pca"])).all(), "X_pca")
    check("umap_shape", out.obsm["X_umap"].shape == (258120, 2), out.obsm["X_umap"].shape)
    check("umap_finite", np.isfinite(np.asarray(out.obsm["X_umap"])).all(), "X_umap")
    params = out.uns["neighbors"]["params"]
    check("neighbor_parameters", int(params["n_neighbors"]) == 15 and params["use_rep"] == "X_scVI" and params["metric"] == "euclidean", params)

    expected_res = {"leiden01": 0.1, "leiden02": 0.2, "leiden03": 0.3, "leiden04": 0.4}
    for key, resolution in expected_res.items():
        no_na = not out.obs[key].isna().any()
        n_clusters = out.obs[key].nunique()
        stored = float(out.uns[key]["params"]["resolution"])
        check(f"{key}_complete", no_na and n_clusters >= 2 and stored == resolution, f"clusters={n_clusters};resolution={stored}")

    with h5py.File(INPUT, "r") as fi, h5py.File(OUTPUT, "r") as fo:
        check("counts_data_exact", dataset_hash(fi["X/data"]) == dataset_hash(fo["layers/counts/data"]), f"input_nnz={fi['X/data'].shape[0]};output_nnz={fo['layers/counts/data'].shape[0]}")
        check("counts_indices_exact", dataset_hash(fi["X/indices"]) == dataset_hash(fo["layers/counts/indices"]), "indices sha256")
        check("counts_indptr_exact", dataset_hash(fi["X/indptr"]) == dataset_hash(fo["layers/counts/indptr"]), "indptr sha256")
        check("normalized_structure_indices", dataset_hash(fo["X/indices"]) == dataset_hash(fo["layers/counts/indices"]), "X vs counts indices")
        check("normalized_structure_indptr", dataset_hash(fo["X/indptr"]) == dataset_hash(fo["layers/counts/indptr"]), "X vs counts indptr")
        count_data = fo["layers/counts/data"]
        x_data = fo["X/data"]
        finite_nonnegative_integer = True
        normalized_finite_nonnegative = True
        for start in range(0, count_data.shape[0], 10_000_000):
            c = np.asarray(count_data[start : start + 10_000_000])
            x = np.asarray(x_data[start : start + 10_000_000])
            finite_nonnegative_integer &= bool(np.isfinite(c).all() and np.all(c >= 0) and np.equal(c, np.floor(c)).all())
            normalized_finite_nonnegative &= bool(np.isfinite(x).all() and np.all(x >= 0))
        check("counts_numeric_integrity", finite_nonnegative_integer, "finite/nonnegative/integer")
        check("normalized_numeric_integrity", normalized_finite_nonnegative, "finite/nonnegative")

    rng = np.random.default_rng(20260802)
    rows = np.sort(rng.choice(out.n_obs, size=1000, replace=False))
    sample = out[rows, :].to_memory()
    counts = sample.layers["counts"].tocsr().astype(np.float64)
    observed = sample.X.tocsr().astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    expected = counts.copy()
    expected.data *= np.repeat(1e4 / totals, np.diff(expected.indptr))
    np.log1p(expected.data, out=expected.data)
    diff = observed - expected
    max_abs = float(np.max(np.abs(diff.data))) if diff.nnz else 0.0
    check("normalization_recomputed", max_abs < 2e-5, f"sampled_rows=1000;max_abs={max_abs:.3g}")

    old = ad.read_h5ad(OLD, backed="r")
    old_pos = old.obs_names.get_indexer(out.obs_names[rows])
    legacy_z = np.asarray(old.obsm["X_scVI"])[old_pos]
    new_z = z[rows]
    check("legacy_latent_not_copied", not np.array_equal(legacy_z, new_z), f"sampled_rows=1000;max_abs={np.max(np.abs(legacy_z-new_z)):.6g}")
    old.file.close()

    history = pd.read_csv(ROOT / "tables/scvi_training_history_v1.tsv", sep="\t")
    run_script_text = RUN_SCRIPT.read_text(encoding="utf-8")
    check(
        "empty_model_history_explained_by_logger_false",
        len(history) == 0 and "logger=False" in run_script_text,
        f"history_rows={len(history)};logger_false={'logger=False' in run_script_text}",
    )

    runtime_text = RUNTIME_LOG.read_text(encoding="utf-8", errors="replace").replace("\r", "\n")
    completed_epochs: dict[int, float] = {}
    for line in runtime_text.splitlines():
        match = re.search(
            r"Epoch\s+(\d+)/31:.*?\|\s*(\d+)/31\s*\[.*?train_loss=([0-9.eE+\-]+)",
            line,
        )
        if match and int(match.group(1)) == int(match.group(2)):
            completed_epochs[int(match.group(1))] = float(match.group(3))
    progress = pd.DataFrame(
        {
            "epoch": sorted(completed_epochs),
            "displayed_train_loss": [completed_epochs[i] for i in sorted(completed_epochs)],
            "evidence_source": str(RUNTIME_LOG),
            "value_precision_note": "terminal progress-bar display; rounded; not model.history",
        }
    )
    progress.to_csv(VERIFY / "scvi_training_progress_from_runtime_v2.tsv", sep="\t", index=False)
    check(
        "runtime_log_all_31_epochs_complete",
        set(completed_epochs) == set(range(1, 32)),
        f"epochs={sorted(completed_epochs)}",
    )
    check(
        "runtime_log_trainer_reached_max_epochs",
        "Trainer.fit` stopped: `max_epochs=31` reached." in runtime_text,
        "max_epochs=31 terminal marker",
    )
    check(
        "runtime_log_terminal_analysis_pass",
        "PART1_SCVI_INTEGRATION_V1_PASS cells=258120 genes=37795 hvg=8000" in runtime_text,
        "terminal analysis marker",
    )
    check("model_saved", (ROOT / "model/scvi_model_v1/model.pt").exists(), ROOT / "model/scvi_model_v1/model.pt")
    gates = pd.read_csv(ROOT / "tables/memory_gate_v1.tsv", sep="\t")
    check("all_memory_gates_passed", len(gates) >= 10 and (gates["returncode"] == 0).all(), f"records={len(gates)}")

    cluster_table = pd.read_csv(ROOT / "tables/leiden_cluster_counts_v1.tsv", sep="\t")
    sums = cluster_table.groupby("cluster_key")["n_cells"].sum().to_dict()
    check("cluster_table_sums", set(sums) == set(expected_res) and all(v == 258120 for v in sums.values()), sums)
    embedding = pd.read_csv(ROOT / "tables/cell_embedding_clusters_v1.tsv.gz", sep="\t", nrows=5)
    required = {"cell_id", "UMAP_1", "UMAP_2", "scVI_1", "scVI_30", *expected_res.keys()}
    check("embedding_table_schema", required.issubset(embedding.columns), sorted(required - set(embedding.columns)))

    inp.file.close()
    out.file.close()

    qc = pd.DataFrame(checks)
    qc.to_csv(VERIFY / "verification_qc_v1.tsv", sep="\t", index=False)
    failures = qc.loc[~qc["passed"]]
    summary = {
        "status": "passed" if failures.empty else "failed",
        "verified_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "checks_passed": int(qc["passed"].sum()),
        "checks_total": int(len(qc)),
        "memory_gate": gate_line,
        "output_h5ad_bytes": OUTPUT.stat().st_size,
    }

    manifest_rows = []
    for path in sorted(
        p
        for p in ROOT.rglob("*")
        if p.is_file()
        and VERIFY not in p.parents
        and p.name != ".part1_scvi_integrated_v1.h5ad.tmp"
    ):
        manifest_rows.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(VERIFY / "file_manifest_v1.tsv", sep="\t", index=False)
    summary["files_manifested"] = len(manifest_rows)
    summary["verification_version"] = "v2"
    summary["analytical_output_retrained_or_modified"] = False
    summary["v1_failure_retained"] = True
    summary["training_completion_evidence"] = {
        "runtime_log": str(RUNTIME_LOG),
        "epochs_complete": sorted(completed_epochs),
        "trainer_terminal_marker": "max_epochs=31 reached",
    }
    (VERIFY / "verification_summary_v2.json").write_text(json.dumps(summary, indent=2) + "\n")
    (VERIFY / "training_history_logging_audit_v2.md").write_text(
        "# scVI training-history logging audit v2\n\n"
        "- The analytical script used `logger=False` with scvi-tools 1.4.2.\n"
        "- In this configuration, scvi-tools leaves `model.history` empty; therefore the v1 history-table check failed.\n"
        "- The immutable runtime log records completion of epochs 1 through 31, the Lightning `max_epochs=31 reached` marker, and the final analytical PASS marker.\n"
        "- Displayed loss values in the v2 progress table are rounded terminal progress-bar values, not reconstructed full-precision model history.\n"
        "- No model retraining and no analytical-output modification were performed for v2 verification.\n"
        "- The failed v1 verification directory is retained as audit evidence.\n",
        encoding="utf-8",
    )
    if not failures.empty:
        raise SystemExit(f"verification failed: {failures[['check','detail']].to_dict('records')}")
    print(f"PART1_SCVI_INTEGRATION_VERIFICATION_V2_ALL_PASS checks={len(qc)} files={len(manifest_rows)}", flush=True)


if __name__ == "__main__":
    main()
