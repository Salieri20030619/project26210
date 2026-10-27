#!/usr/bin/env python3
"""TJ-RCC three-case cell2location compatibility pilot.

Uses the existing 14-state reference signature matrix. It never edits source
Visium files or the GLMF1/GLMF2 outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
SOURCE_ROOT = (
    ROOT
    / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/extracted/ST_visium"
)
OUTPUT_ROOT = (
    ROOT
    / "external/ccRCC_spatial/TJ_RCC_figshare_24599295/"
    "cell2location_pilot_v2"
)
SIGNATURES = ROOT / "results/part4_spatial/ref/ref_signatures_matrix.csv"
SAMPLES = {
    "R51_T": "high",
    "X49_T": "medium",
    "X98_T": "low",
}
SEED = 20260724
N_CELLS_PER_LOCATION = 30
DETECTION_ALPHA = 20
MAX_EPOCHS = 10000
POSTERIOR_SAMPLES = 1000


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def locate_outs(sample: str) -> Path:
    sample_root = SOURCE_ROOT / sample
    candidates = [sample_root / "outs", sample_root]
    for candidate in candidates:
        if (
            (candidate / "filtered_feature_bc_matrix.h5").is_file()
            and (candidate / "spatial").is_dir()
        ):
            return candidate
    raise FileNotFoundError(f"No Space Ranger output found for {sample}")


def load_sample(sample: str):
    outs = locate_outs(sample)
    adata = sc.read_visium(
        outs,
        count_file="filtered_feature_bc_matrix.h5",
        load_images=True,
    )
    adata.var["gene_symbol_original"] = adata.var_names.astype(str)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample
    adata.var["mt"] = adata.var["gene_symbol_original"].str.upper().str.startswith(
        "MT-"
    )
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    return adata, outs


def load_signatures() -> pd.DataFrame:
    signatures = pd.read_csv(SIGNATURES, index_col=0)
    prefix = "means_per_cluster_mu_fg_"
    signatures.columns = [
        x[len(prefix) :] if x.startswith(prefix) else x for x in signatures.columns
    ]
    if signatures.index.has_duplicates or signatures.columns.has_duplicates:
        raise ValueError("Reference signature matrix has duplicate identifiers")
    if signatures.shape[1] != 14:
        raise ValueError(f"Expected 14 reference states, found {signatures.shape[1]}")
    values = signatures.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Reference signature matrix contains invalid values")
    return signatures


def preflight_one(sample: str, signatures: pd.DataFrame) -> dict:
    adata, outs = load_sample(sample)
    input_h5 = outs / "filtered_feature_bc_matrix.h5"
    original_symbols = adata.var["gene_symbol_original"].astype(str)
    shared = adata.var_names.intersection(signatures.index)
    target_states = [
        "TYMP-High Malignant",
        "TYMP-Low Malignant",
        "SPP1+ TAM",
        "APOE+C1QC+ TAM",
    ]
    return {
        "sample": sample,
        "complexity_group": SAMPLES[sample],
        "space_ranger_outs": str(outs),
        "input_h5": str(input_h5),
        "input_h5_size_bytes": input_h5.stat().st_size,
        "input_h5_sha256": sha256(input_h5),
        "spots_before_qc": int(adata.n_obs),
        "genes_before_intersection": int(adata.n_vars),
        "duplicate_gene_symbols_before_make_unique": int(
            original_symbols.duplicated(keep=False).sum()
        ),
        "shared_reference_genes": int(len(shared)),
        "reference_genes": int(signatures.shape[0]),
        "shared_reference_fraction": float(len(shared) / signatures.shape[0]),
        "reference_states": int(signatures.shape[1]),
        "key_reference_states_present": all(x in signatures.columns for x in target_states),
        "median_umi": float(adata.obs["total_counts"].median()),
        "median_detected_genes": float(adata.obs["n_genes_by_counts"].median()),
        "median_mito_fraction": float(adata.obs["pct_counts_mt"].median() / 100),
    }


def write_preflight(records: list[dict], signatures: pd.DataFrame) -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(OUTPUT_ROOT / "preflight_samples.tsv", sep="\t", index=False)
    metadata = {
        "created_at": now(),
        "status": "preflight_passed",
        "selected_samples": SAMPLES,
        "selection_basis": (
            "R51_T highest median UMI; X49_T intermediate complexity; "
            "X98_T lowest median UMI among 12 tumor cases"
        ),
        "reference_signatures": str(SIGNATURES),
        "reference_signatures_size_bytes": SIGNATURES.stat().st_size,
        "reference_signatures_sha256": sha256(SIGNATURES),
        "reference_shape": list(signatures.shape),
        "parameters": {
            "N_cells_per_location": N_CELLS_PER_LOCATION,
            "detection_alpha": DETECTION_ALPHA,
            "max_epochs": MAX_EPOCHS,
            "posterior_samples": POSTERIOR_SAMPLES,
            "seed": SEED,
            "accelerator": "cpu",
            "execution": "sequential",
        },
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
        },
    }
    (OUTPUT_ROOT / "run_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def run_sample(sample: str, signatures: pd.DataFrame) -> None:
    from cell2location.models import Cell2location
    import cell2location

    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.set_num_threads(32)

    sample_dir = OUTPUT_ROOT / f"{SAMPLES[sample]}_{sample}"
    if sample_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {sample_dir}")
    sample_dir.mkdir(parents=True)
    status_path = sample_dir / "status.json"

    def write_status(status: str, **extra) -> None:
        payload = {
            "sample": sample,
            "complexity_group": SAMPLES[sample],
            "status": status,
            "updated_at": now(),
            **extra,
        }
        status_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

    write_status("loading_input")
    adata, outs = load_sample(sample)
    input_h5 = outs / "filtered_feature_bc_matrix.h5"
    input_hash = sha256(input_h5)

    keep = (
        (adata.obs["total_counts"] >= 1000)
        & (adata.obs["n_genes_by_counts"] >= 300)
        & (adata.obs["pct_counts_mt"] < 20)
    )
    before = int(adata.n_obs)
    adata = adata[keep].copy()
    after = int(adata.n_obs)
    genes = adata.var_names.intersection(signatures.index)
    if len(genes) < 10000:
        raise ValueError(f"Too few shared reference genes: {len(genes)}")
    adata = adata[:, genes].copy()
    cell_state_df = signatures.loc[genes].copy()
    adata.uns["tjrcc_cell2location_pilot"] = {
        "sample": sample,
        "complexity_group": SAMPLES[sample],
        "input_h5": str(input_h5),
        "input_h5_sha256": input_hash,
        "spots_before_qc": before,
        "spots_after_qc": after,
        "shared_reference_genes": int(len(genes)),
        "seed": SEED,
    }
    write_status(
        "training",
        spots_before_qc=before,
        spots_after_qc=after,
        shared_reference_genes=int(len(genes)),
        input_h5_sha256=input_hash,
        cell2location_version=cell2location.__version__,
    )

    Cell2location.setup_anndata(adata)
    model = Cell2location(
        adata,
        cell_state_df=cell_state_df,
        N_cells_per_location=N_CELLS_PER_LOCATION,
        detection_alpha=DETECTION_ALPHA,
    )
    model.train(
        max_epochs=MAX_EPOCHS,
        batch_size=None,
        train_size=1,
        accelerator="cpu",
    )
    write_status("exporting_posterior")
    adata = model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": POSTERIOR_SAMPLES,
            "batch_size": min(1024, adata.n_obs),
            "accelerator": "cpu",
        },
    )
    model.save(sample_dir / "model", overwrite=False)
    output_h5ad = sample_dir / f"{sample}_cell2location_pilot.h5ad"
    adata.write(output_h5ad)

    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key not in adata.obsm:
        raise KeyError(f"Missing expected posterior abundance: {abundance_key}")
    abundance = adata.obsm[abundance_key]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(
            abundance, index=adata.obs_names, columns=cell_state_df.columns
        )
    abundance.to_csv(sample_dir / "cell_abundance_q05.tsv", sep="\t")
    summary = pd.DataFrame(
        {
            "mean_q05": abundance.mean(axis=0),
            "median_q05": abundance.median(axis=0),
            "nonzero_fraction": (abundance > 0).mean(axis=0),
        }
    )
    summary.to_csv(sample_dir / "cell_abundance_summary.tsv", sep="\t")

    if sha256(input_h5) != input_hash:
        raise RuntimeError("Input H5 hash changed during mapping")
    write_status(
        "completed",
        spots_before_qc=before,
        spots_after_qc=after,
        shared_reference_genes=int(len(genes)),
        output_h5ad=str(output_h5ad),
        output_h5ad_size_bytes=output_h5ad.stat().st_size,
        output_h5ad_sha256=sha256(output_h5ad),
        input_h5_sha256=input_hash,
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--sample", choices=list(SAMPLES))
    return parser.parse_args()


def main():
    args = parse_args()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    signatures = load_signatures()
    records = [preflight_one(sample, signatures) for sample in SAMPLES]
    write_preflight(records, signatures)
    if args.preflight_only:
        print(pd.DataFrame(records).to_string(index=False))
        return
    if not args.sample:
        raise SystemExit("--sample is required unless --preflight-only is used")
    run_sample(args.sample, signatures)


if __name__ == "__main__":
    main()
