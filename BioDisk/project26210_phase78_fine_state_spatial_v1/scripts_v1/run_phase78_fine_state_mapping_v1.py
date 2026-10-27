#!/usr/bin/env python3
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
import scipy
from scipy import sparse
import torch
import scvi
from cell2location.models import Cell2location

SEED = 26210
EXPECTED_STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]


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


def load_signatures(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    frame.index = frame.index.astype(str).str.upper()
    frame = frame.loc[~frame.index.duplicated()].copy()
    if list(frame.columns) != EXPECTED_STATES:
        raise RuntimeError(f"Unexpected reference states: {list(frame.columns)}")
    values = frame.to_numpy(dtype=float)
    if frame.shape[0] < 12000 or not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError(f"Invalid reference signature matrix: {frame.shape}")
    return frame


def prepare_spatial(path: Path, signatures: pd.DataFrame) -> ad.AnnData:
    spatial = ad.read_h5ad(path)
    if "counts" not in spatial.layers:
        raise RuntimeError(f"Missing counts layer: {path}")
    spatial.X = sparse.csr_matrix(spatial.layers["counts"])
    spatial.var_names = pd.Index(spatial.var_names.astype(str).str.upper())
    spatial.var_names_make_unique()
    common = spatial.var_names.intersection(signatures.index)
    if len(common) < 8000:
        raise RuntimeError(f"Only {len(common)} shared genes: {path}")
    spatial = spatial[:, common].copy()
    spatial.layers["counts"] = sparse.csr_matrix(spatial.X.copy())
    return spatial


def abundance_frame(spatial: ad.AnnData, factors: list[str]) -> pd.DataFrame:
    result = pd.DataFrame({"spot_id": spatial.obs_names.astype(str)})
    for slot, prefix in [
        ("means_cell_abundance_w_sf", "mean"),
        ("q05_cell_abundance_w_sf", "q05"),
        ("q95_cell_abundance_w_sf", "q95"),
    ]:
        if slot not in spatial.obsm:
            raise RuntimeError(f"Missing posterior slot {slot}")
        value = spatial.obsm[slot]
        array = value.to_numpy() if isinstance(value, pd.DataFrame) else np.asarray(value)
        if array.shape != (spatial.n_obs, len(factors)):
            raise RuntimeError(f"Unexpected {slot} shape {array.shape}")
        for index, factor in enumerate(factors):
            result[f"{prefix}__{factor}"] = array[:, index]
    return result


def map_one(
    sample: str,
    section_path_text: str,
    signature_path_text: str,
    output_sections_text: str,
    sample_index: int,
) -> dict:
    torch.set_num_threads(4)
    scvi.settings.seed = SEED + sample_index
    np.random.seed(SEED + sample_index)
    section_path = Path(section_path_text)
    signatures = load_signatures(Path(signature_path_text))
    final = Path(output_sections_text) / sample
    part = final.with_name(final.name + ".part")
    if (final / "MAPPING_PASS_v1").exists():
        return json.loads((final / "mapping_summary_v1.json").read_text(encoding="utf-8"))
    if final.exists() or part.exists():
        raise FileExistsError(f"Existing incomplete output for {sample}")
    part.mkdir(parents=True)
    started = now()
    try:
        spatial = prepare_spatial(section_path, signatures)
        common = spatial.var_names
        cell_state = signatures.loc[common].copy()
        factors = list(cell_state.columns)
        spatial.obs["section_batch"] = sample
        Cell2location.setup_anndata(spatial, layer="counts", batch_key="section_batch")
        model = Cell2location(
            spatial,
            cell_state_df=cell_state,
            N_cells_per_location=10,
            detection_alpha=20,
        )
        print(f"[{now()}] training {sample} spots={spatial.n_obs} genes={spatial.n_vars}", flush=True)
        model.train(
            max_epochs=10000,
            batch_size=None,
            train_size=1,
            lr=0.002,
            accelerator="cpu",
            enable_progress_bar=False,
        )
        model.save(part / "model", overwrite=True)
        spatial = model.export_posterior(
            spatial,
            sample_kwargs={"num_samples": 1000, "batch_size": spatial.n_obs},
        )
        abund = abundance_frame(spatial, factors)
        obs = spatial.obs.copy()
        obs.insert(0, "spot_id", spatial.obs_names.astype(str))
        atomic_tsv(abund, part / "fine_state_abundance_v1.tsv.gz", gzip_output=True)
        atomic_tsv(obs.reset_index(drop=True), part / "spot_metadata_v1.tsv.gz", gzip_output=True)
        spatial.write_h5ad(part / "fine_state_posterior_v1.h5ad", compression="gzip")
        summary = {
            "status": "PASS",
            "sample": sample,
            "started_at": started,
            "completed_at": now(),
            "spots": int(spatial.n_obs),
            "genes": int(spatial.n_vars),
            "states": factors,
            "seed": SEED + sample_index,
            "input": str(section_path),
            "input_sha256": sha256(section_path),
        }
        (part / "mapping_summary_v1.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        (part / "MAPPING_PASS_v1").write_text("status=PASS\n", encoding="utf-8")
        os.replace(part, final)
        print(f"[{now()}] PASS {sample}", flush=True)
        return summary
    except Exception:
        (part / "TRACEBACK_v1.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", required=True, type=Path)
    parser.add_argument("--section-root", required=True, type=Path)
    parser.add_argument("--section-status", required=True, type=Path)
    parser.add_argument("--preregistration", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    signatures = load_signatures(args.signature)
    status = pd.read_csv(args.section_status, sep="\t")
    status = status.loc[status["status"].eq("PASS")].copy()
    if len(status) != 24 or status["sample_id"].nunique() != 24:
        raise RuntimeError(f"Expected 24 passed sections, found {len(status)}")
    if set(status["preparation"]) != {"FFPE", "frozen"}:
        raise RuntimeError("Expected FFPE and frozen preparations")
    records = []
    for row in status.itertuples(index=False):
        path = args.section_root / str(row.sample_id) / "processed_section_v1.h5ad"
        if not path.is_file():
            raise FileNotFoundError(path)
        probe = ad.read_h5ad(path, backed="r")
        records.append({
            "sample_id": str(row.sample_id),
            "preparation": str(row.preparation),
            "spots": int(probe.n_obs),
            "input_genes": int(probe.n_vars),
            "path": str(path),
        })
        probe.file.close()

    args.output_root.mkdir(parents=True, exist_ok=True)
    output_sections = args.output_root / "sections_v1"
    output_sections.mkdir(exist_ok=True)
    atomic_tsv(pd.DataFrame(records), args.output_root / "preflight_sections_v1.tsv")
    metadata = {
        "status": "RUNNING",
        "started_at": now(),
        "workers": args.workers,
        "signature": str(args.signature),
        "signature_sha256": sha256(args.signature),
        "signature_shape": list(signatures.shape),
        "preregistration": str(args.preregistration),
        "preregistration_sha256": sha256(args.preregistration),
        "parameters": {
            "N_cells_per_location": 10,
            "detection_alpha": 20,
            "max_epochs": 10000,
            "posterior_samples": 1000,
            "threads_per_worker": 4,
            "seed_base": SEED,
        },
        "environment": {
            "python": platform.python_version(),
            "torch": torch.__version__,
            "scvi": scvi.__version__,
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    metadata_path = args.output_root / "run_metadata_v1.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    failures = []
    completed = []
    context = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers, mp_context=context
    ) as pool:
        futures = {}
        for index, record in enumerate(records):
            future = pool.submit(
                map_one,
                record["sample_id"],
                record["path"],
                str(args.signature),
                str(output_sections),
                index,
            )
            futures[future] = record["sample_id"]
        for future in concurrent.futures.as_completed(futures):
            sample = futures[future]
            try:
                completed.append(future.result())
            except Exception as error:
                failures.append({"sample_id": sample, "error": repr(error)})
                print(f"[{now()}] FAIL {sample}: {error!r}", flush=True)

    atomic_tsv(pd.DataFrame(completed), args.output_root / "mapping_completed_v1.tsv")
    atomic_tsv(pd.DataFrame(failures), args.output_root / "mapping_failures_v1.tsv")
    metadata["completed_at"] = now()
    metadata["sections_passed"] = len(completed)
    metadata["sections_failed"] = len(failures)
    metadata["status"] = "PASS" if len(completed) == 24 and not failures else "FAIL"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    marker = "PHASE78_MAPPING_PASS_v1" if metadata["status"] == "PASS" else "PHASE78_MAPPING_FAIL_v1"
    (args.output_root / marker).write_text(
        f"status={metadata['status']}\nsections_passed={len(completed)}\nsections_failed={len(failures)}\n",
        encoding="utf-8",
    )
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
