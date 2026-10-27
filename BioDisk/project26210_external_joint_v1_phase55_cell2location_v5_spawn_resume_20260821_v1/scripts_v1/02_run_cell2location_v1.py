#!/usr/bin/env python3
"""CPU cell2location reference training and independent mapping of 24 GSE175540 sections."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing
import os
import platform
import traceback
from datetime import datetime
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy import sparse
import torch
import scvi
from cell2location.models import RegressionModel, Cell2location
from cell2location.utils.filtering import filter_genes


SEED = 26210


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path, gzip_output: bool = False) -> None:
    temporary = path.with_name(path.name + ".part")
    frame.to_csv(
        temporary, sep="\t", index=False, na_rep="NA", encoding="utf-8",
        lineterminator="\n", compression="gzip" if gzip_output else None,
    )
    os.replace(temporary, path)


def set_upper_unique_var_names(adata: ad.AnnData) -> None:
    adata.var_names = pd.Index(adata.var_names.astype(str).str.upper())
    adata.var_names_make_unique()


def as_bool(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(bool)
    return series.astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy(bool)


def balanced_reference(adata: ad.AnnData, label_key: str, batch_key: str) -> tuple[ad.AnnData, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    keep_indices = []
    rows = []
    labels = sorted(adata.obs[label_key].astype(str).unique())
    for label in labels:
        label_indices = np.flatnonzero(adata.obs[label_key].astype(str).to_numpy() == label)
        selected_label = []
        batches = sorted(adata.obs.iloc[label_indices][batch_key].astype(str).unique())
        for batch in batches:
            batch_mask = adata.obs.iloc[label_indices][batch_key].astype(str).to_numpy() == batch
            candidates = label_indices[batch_mask]
            take = min(3000, len(candidates))
            chosen = rng.choice(candidates, size=take, replace=False) if take < len(candidates) else candidates
            selected_label.extend(chosen.tolist())
            rows.append({"cell_type": label, "batch": batch, "available": len(candidates), "selected_before_label_cap": take})
        if len(selected_label) > 15000:
            selected_label = rng.choice(np.asarray(selected_label), size=15000, replace=False).tolist()
        keep_indices.extend(selected_label)
    selected = adata[np.sort(np.asarray(keep_indices, dtype=int))].copy()
    selection = pd.DataFrame(rows)
    final_counts = selected.obs.groupby([label_key, batch_key], observed=True).size().rename("final_selected").reset_index()
    selection = selection.merge(final_counts, left_on=["cell_type", "batch"], right_on=[label_key, batch_key], how="left")
    selection["final_selected"] = selection["final_selected"].fillna(0).astype(int)
    return selected, selection


def extract_reference_signatures(adata_ref: ad.AnnData) -> pd.DataFrame:
    key = "means_per_cluster_mu_fg"
    if key not in adata_ref.varm:
        raise KeyError(f"Missing {key} in reference posterior")
    value = adata_ref.varm[key]
    factors = [str(x) for x in adata_ref.uns["mod"]["factor_names"]]
    if isinstance(value, pd.DataFrame):
        frame = value.copy()
        expected = [f"means_per_cluster_mu_fg_{name}" for name in factors]
        if all(column in frame.columns for column in expected):
            frame = frame[expected]
            frame.columns = factors
        elif frame.shape[1] == len(factors):
            frame.columns = factors
        else:
            raise ValueError("Reference signature column count does not match factor names")
        frame.index = adata_ref.var_names
        return frame
    array = np.asarray(value)
    if array.shape != (adata_ref.n_vars, len(factors)):
        raise ValueError(f"Unexpected reference signature shape {array.shape}")
    return pd.DataFrame(array, index=adata_ref.var_names, columns=factors)


def prepare_and_train_reference(reference_path: Path, output_root: Path) -> Path:
    final = output_root / "reference_model_v1"
    part = output_root / "reference_model_v1.part"
    signature_path = final / "reference_cell_state_signatures_v1.tsv.gz"
    if (final / "REFERENCE_MODEL_PASS_v1").exists() and signature_path.exists():
        return signature_path
    if final.exists() or part.exists():
        raise FileExistsError("Existing incomplete reference-model output found")
    part.mkdir(parents=True)
    print(f"[{now()}] Loading reference {reference_path}", flush=True)
    source = sc.read_h5ad(reference_path)
    label_key = "cell_type_broad_v1"
    batch_key = "study"
    required = {label_key, batch_key, "annotation_ambiguous_v1"}
    if not required.issubset(source.obs.columns) or "counts" not in source.layers:
        raise RuntimeError("Reference lacks required labels, batch, ambiguity, or counts layer")
    keep = ~as_bool(source.obs["annotation_ambiguous_v1"])
    keep &= ~source.obs[label_key].astype(str).isin(["Epithelial/Myeloid mixed"]).to_numpy()
    if "analysis_keep" in source.obs:
        keep &= as_bool(source.obs["analysis_keep"])
    source = source[keep].copy()
    source.obs[label_key] = source.obs[label_key].astype(str)
    source.obs[batch_key] = source.obs[batch_key].astype(str)
    source.X = sparse.csr_matrix(source.layers["counts"])
    set_upper_unique_var_names(source)
    if source.X.data.size and (
        not np.isfinite(source.X.data).all() or source.X.data.min() < 0
        or not np.equal(source.X.data, np.floor(source.X.data)).all()
    ):
        raise ValueError("Reference counts fail nonnegative finite integer checks")
    source.layers["counts"] = source.X.copy()
    reference, selection = balanced_reference(source, label_key, batch_key)
    del source
    selected_genes = filter_genes(
        reference, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
    )
    reference = reference[:, selected_genes].copy()
    reference.X = reference.layers["counts"].copy()
    atomic_tsv(selection, part / "reference_balancing_v1.tsv")
    label_counts = reference.obs[label_key].value_counts().rename_axis("cell_type").reset_index(name="cells")
    atomic_tsv(label_counts, part / "reference_cell_type_counts_v1.tsv")
    reference.write_h5ad(part / "reference_balanced_counts_v1.h5ad", compression="gzip")

    scvi.settings.seed = SEED
    torch.set_num_threads(8)
    RegressionModel.setup_anndata(reference, layer="counts", batch_key=batch_key, labels_key=label_key)
    model = RegressionModel(reference)
    print(f"[{now()}] Training reference model cells={reference.n_obs} genes={reference.n_vars}", flush=True)
    model.train(
        max_epochs=250, batch_size=2500, train_size=1, lr=0.002,
        accelerator="cpu", enable_progress_bar=False,
    )
    model.save(part / "model", overwrite=True)
    reference = model.export_posterior(
        reference,
        sample_kwargs={"num_samples": 1000, "batch_size": 2500},
    )
    signatures = extract_reference_signatures(reference)
    sig_out = signatures.copy()
    sig_out.insert(0, "gene", sig_out.index)
    atomic_tsv(sig_out, part / "reference_cell_state_signatures_v1.tsv.gz", gzip_output=True)
    reference.write_h5ad(part / "reference_posterior_v1.h5ad", compression="gzip")
    metadata = {
        "status": "PASS", "completed_at": now(), "source_reference": str(reference_path),
        "source_sha256": sha256(reference_path), "cells": reference.n_obs, "genes": reference.n_vars,
        "cell_types": list(signatures.columns), "epochs": 250, "device": "cpu",
    }
    (part / "reference_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (part / "REFERENCE_MODEL_PASS_v1").write_text(
        f"status=PASS\ncells={reference.n_obs}\ngenes={reference.n_vars}\ncell_types={len(signatures.columns)}\n",
        encoding="utf-8",
    )
    os.replace(part, final)
    print(f"[{now()}] Reference model PASS", flush=True)
    return final / "reference_cell_state_signatures_v1.tsv.gz"


def abundance_frame(adata_vis: ad.AnnData, factors: list[str]) -> pd.DataFrame:
    result = pd.DataFrame(index=adata_vis.obs_names)
    for slot, prefix in [
        ("means_cell_abundance_w_sf", "mean"),
        ("q05_cell_abundance_w_sf", "q05"),
        ("q95_cell_abundance_w_sf", "q95"),
    ]:
        if slot not in adata_vis.obsm:
            continue
        value = adata_vis.obsm[slot]
        if isinstance(value, pd.DataFrame):
            array = value.to_numpy()
        else:
            array = np.asarray(value)
        if array.shape[1] != len(factors):
            raise ValueError(f"Abundance slot {slot} has shape {array.shape}")
        for index, factor in enumerate(factors):
            result[f"{prefix}__{factor}"] = array[:, index]
    return result


def map_one_section(section_dir_text: str, signature_path_text: str, output_sections_text: str) -> dict:
    torch.set_num_threads(4)
    scvi.settings.seed = SEED
    section_dir = Path(section_dir_text)
    sample = section_dir.name
    final = Path(output_sections_text) / sample
    part = Path(output_sections_text) / f"{sample}.part"
    if (final / "DECONVOLUTION_PASS_v1").exists():
        return json.loads((final / "deconvolution_summary_v1.json").read_text(encoding="utf-8"))
    if final.exists() or part.exists():
        raise FileExistsError(f"Existing incomplete deconvolution output for {sample}")
    part.mkdir(parents=True)
    print(f"[{now()}] Worker started {sample}", flush=True)
    started = now()
    try:
        signatures = pd.read_csv(signature_path_text, sep="\t", compression="gzip").set_index("gene")
        signatures.index = signatures.index.astype(str).str.upper()
        spatial_path = section_dir / "processed_section_v1.h5ad"
        spatial = sc.read_h5ad(spatial_path)
        if "counts" not in spatial.layers:
            raise RuntimeError("Spatial object lacks counts layer")
        spatial.X = sparse.csr_matrix(spatial.layers["counts"])
        set_upper_unique_var_names(spatial)
        common = spatial.var_names.intersection(signatures.index)
        if len(common) < 5000:
            raise RuntimeError(f"Only {len(common)} shared genes")
        spatial = spatial[:, common].copy()
        spatial.layers["counts"] = spatial.X.copy()
        cell_state = signatures.loc[common].copy()
        factors = list(cell_state.columns)
        spatial.obs["section_batch"] = sample
        Cell2location.setup_anndata(spatial, layer="counts", batch_key="section_batch")
        model = Cell2location(
            spatial, cell_state_df=cell_state,
            N_cells_per_location=10, detection_alpha=20,
        )
        print(f"[{now()}] Training {sample} spots={spatial.n_obs} genes={spatial.n_vars}", flush=True)
        model.train(
            max_epochs=20000, batch_size=None, train_size=1, lr=0.002,
            accelerator="cpu", enable_progress_bar=False,
        )
        model.save(part / "model", overwrite=True)
        spatial = model.export_posterior(
            spatial,
            sample_kwargs={"num_samples": 1000, "batch_size": spatial.n_obs},
        )
        abundance = abundance_frame(spatial, factors)
        metadata_columns = [column for column in ["array_row", "array_col", "pxl_row", "pxl_col", "TLS_annotation", "leiden_0_5"] if column in spatial.obs]
        table = spatial.obs[metadata_columns].copy()
        table.insert(0, "barcode", table.index)
        table = table.join(abundance)
        atomic_tsv(table, part / "spot_cell_abundance_v1.tsv.gz", gzip_output=True)
        spatial.write_h5ad(part / "cell2location_posterior_v1.h5ad", compression="gzip")
        summary = {
            "status": "PASS", "sample_id": sample, "started_at": started, "completed_at": now(),
            "spots": spatial.n_obs, "shared_genes": spatial.n_vars, "cell_types": len(factors),
            "epochs": 20000, "source_spatial_sha256": sha256(spatial_path),
        }
        (part / "deconvolution_summary_v1.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        (part / "DECONVOLUTION_PASS_v1").write_text(
            f"status=PASS\nsample_id={sample}\nspots={spatial.n_obs}\nshared_genes={spatial.n_vars}\n",
            encoding="utf-8",
        )
        os.replace(part, final)
        return summary
    except Exception:
        (part / "DECONVOLUTION_FAILURE_v1.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--spatial-sections", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    signature_path = prepare_and_train_reference(args.reference, args.output_root)
    sections = sorted(path for path in args.spatial_sections.iterdir() if path.is_dir() and path.name.startswith("GSM"))
    if len(sections) != 24:
        raise RuntimeError(f"Expected 24 spatial section directories, found {len(sections)}")
    output_sections = args.output_root / "sections_v1"
    output_sections.mkdir(exist_ok=True)
    summaries = []
    failures = []
    # PyTorch/scvi initialize native thread pools while training the reference
    # model. Forking that initialized parent can leave every child blocked on
    # an inherited futex. Fresh spawn interpreters avoid those inherited locks
    # while preserving the same parameters and independent section jobs.
    spawn_context = multiprocessing.get_context("spawn")
    print(
        f"[{now()}] Starting spatial pool workers={args.workers} start_method=spawn",
        flush=True,
    )
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers, mp_context=spawn_context
    ) as pool:
        futures = {
            pool.submit(map_one_section, str(section), str(signature_path), str(output_sections)): section.name
            for section in sections
        }
        for future in concurrent.futures.as_completed(futures):
            sample = futures[future]
            try:
                result = future.result()
                summaries.append(result)
                print(f"[{now()}] PASS {sample}", flush=True)
            except Exception as exc:
                failures.append({"sample_id": sample, "error": repr(exc)})
                print(f"[{now()}] FAIL {sample}: {exc!r}", flush=True)
    status = pd.DataFrame(summaries).sort_values("sample_id") if summaries else pd.DataFrame()
    atomic_tsv(status, args.output_root / "section_deconvolution_status_v1.tsv")
    atomic_tsv(pd.DataFrame(failures, columns=["sample_id", "error"]), args.output_root / "section_deconvolution_failures_v1.tsv")
    checks = pd.DataFrame([
        {"check": "reference_model_pass", "status": "PASS" if (args.output_root / "reference_model_v1" / "REFERENCE_MODEL_PASS_v1").exists() else "FAIL", "observed": (args.output_root / "reference_model_v1" / "REFERENCE_MODEL_PASS_v1").exists(), "expected": True},
        {"check": "spatial_sections_expected", "status": "PASS" if len(sections) == 24 else "FAIL", "observed": len(sections), "expected": 24},
        {"check": "spatial_sections_passed", "status": "PASS" if len(summaries) == 24 else "FAIL", "observed": len(summaries), "expected": 24},
        {"check": "spatial_sections_failed", "status": "PASS" if len(failures) == 0 else "FAIL", "observed": len(failures), "expected": 0},
    ])
    atomic_tsv(checks, args.output_root / "validation_checks_v1.tsv")
    final_status = "PASS" if (checks["status"] == "PASS").all() else "PARTIAL"
    metadata = {
        "status": final_status, "generated_at": now(), "workers": args.workers,
        "multiprocessing_start_method": "spawn",
        "reference": str(args.reference), "reference_sha256": sha256(args.reference),
        "sections_passed": len(summaries), "sections_failed": len(failures),
        "python": platform.python_version(), "torch": torch.__version__, "cell2location": "0.1.5",
        "scvi_tools": scvi.__version__, "anndata": ad.__version__, "scanpy": sc.__version__,
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "device": "cpu",
    }
    (args.output_root / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    marker = args.output_root / ("PHASE55_CELL2LOCATION_PASS_v1" if final_status == "PASS" else "PHASE55_CELL2LOCATION_PARTIAL_v1")
    marker.write_text(f"status={final_status}\nsections_passed={len(summaries)}\nsections_failed={len(failures)}\n", encoding="utf-8")
    return 0 if final_status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
