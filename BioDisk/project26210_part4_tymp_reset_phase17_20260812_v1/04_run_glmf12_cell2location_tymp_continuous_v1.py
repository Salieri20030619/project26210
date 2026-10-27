#!/usr/bin/env python3
"""Phase 17 total-malignant Cell2location remapping for GLMF1 and GLMF2."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import random
import subprocess
import sys
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.sparse import issparse


ROOT = Path("/mnt/BioDisk/wangmingyuan/project26210_part4_tymp_reset_phase17_20260812_v1")
PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
SIGNATURES = (
    ROOT / "reference_tymp_continuous_v1/reference_tymp_continuous_v1_signatures_study_balanced_mean_counts.tsv.gz"
)
GENE_TEMPLATE = PROJECT / "results/part4_spatial/cell2loc/GLMF1_cell2loc.h5ad"
OUTPUT_ROOT = ROOT / "glmf12_cell2location_tymp_continuous_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
SAMPLES = {
    "GLMF1": PROJECT / "results/part4_spatial/mapping/GLMF1_QC.h5ad",
    "GLMF2": PROJECT / "results/part4_spatial/mapping/GLMF2_QC.h5ad",
}
SEED = 20260805
N_CELLS_PER_LOCATION = 30
DETECTION_ALPHA = 20
MAX_EPOCHS = 10000
POSTERIOR_SAMPLES = 1000
CHECKPOINT_EVERY_N_EPOCHS = 1000
THREADS_PER_SAMPLE = 32


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def memory_gate(stage: str, sample: str) -> None:
    resource_dir = OUTPUT_ROOT / "resource_logs"
    resource_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    message = (result.stdout or result.stderr).strip().replace("\n", " | ")
    with (resource_dir / f"{sample}.tsv").open("a", encoding="utf-8") as handle:
        handle.write(f"{now()}\t{stage}\t{message}\texit_code={result.returncode}\n")
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked {sample} at {stage}: {message}")


def load_signatures() -> pd.DataFrame:
    signatures = pd.read_csv(SIGNATURES, sep="\t", index_col=0)
    if signatures.shape != (12758, 13):
        raise RuntimeError(f"Unexpected reference shape: {signatures.shape}")
    if signatures.index.has_duplicates or signatures.columns.has_duplicates:
        raise RuntimeError("Duplicate reference identifiers")
    values = signatures.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError("Invalid reference values")
    return signatures


def load_gene_template() -> pd.Index:
    template = ad.read_h5ad(GENE_TEMPLATE, backed="r")
    genes = template.var_names.astype(str).copy()
    template.file.close()
    if len(genes) != 8830 or genes.has_duplicates:
        raise RuntimeError(f"Unexpected fixed comparison gene set: {len(genes)}")
    return genes


def preflight_one(sample: str, signatures: pd.DataFrame, genes: pd.Index) -> dict[str, object]:
    memory_gate("preflight", sample)
    input_path = SAMPLES[sample]
    obj = ad.read_h5ad(input_path, backed="r")
    missing_input = genes.difference(obj.var_names)
    missing_reference = genes.difference(signatures.index)
    shared = genes.intersection(obj.var_names).intersection(signatures.index)
    x = obj.X[: min(100, obj.n_obs), [int(obj.var_names.get_loc(g)) for g in shared[:1000]]]
    if issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=float)
    record = {
        "sample": sample,
        "input_h5ad": str(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "input_sha256": sha256(input_path),
        "n_spots": int(obj.n_obs),
        "n_input_genes": int(obj.n_vars),
        "fixed_model_genes": int(len(genes)),
        "shared_model_genes": int(len(shared)),
        "missing_input_genes": int(len(missing_input)),
        "missing_reference_genes": int(len(missing_reference)),
        "sample_counts_nonnegative": bool(np.min(x) >= 0),
        "sample_counts_integer_like": bool(np.allclose(x, np.rint(x))),
        "spatial_shape": "x".join(map(str, np.asarray(obj.obsm["spatial"]).shape)),
    }
    obj.file.close()
    if len(shared) != len(genes) or not record["sample_counts_nonnegative"] or not record["sample_counts_integer_like"]:
        raise RuntimeError(f"Preflight failed: {record}")
    return record


def write_preflight(records: list[dict[str, object]], signatures: pd.DataFrame) -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(OUTPUT_ROOT / "preflight_glmf12_tymp_continuous_v1.tsv", sep="\t", index=False)
    metadata = {
        "status": "preflight_passed",
        "created_at": now(),
        "samples": list(SAMPLES),
        "reference_signatures": str(SIGNATURES),
        "reference_sha256": sha256(SIGNATURES),
        "reference_shape": list(signatures.shape),
        "gene_set": "fixed 8,830-gene set from the prior GLMF1/GLMF2 cell2location mapping",
        "parameters": {
            "N_cells_per_location": N_CELLS_PER_LOCATION,
            "detection_alpha": DETECTION_ALPHA,
            "max_epochs": MAX_EPOCHS,
            "posterior_samples": POSTERIOR_SAMPLES,
            "checkpoint_every_n_epochs": CHECKPOINT_EVERY_N_EPOCHS,
            "threads_per_sample": THREADS_PER_SAMPLE,
            "execution": "GLMF1 and GLMF2 as independent parallel CPU jobs",
            "seed": SEED,
        },
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
        },
    }
    (OUTPUT_ROOT / "run_metadata_tymp_continuous_v1.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def run_sample(sample: str, signatures: pd.DataFrame, genes: pd.Index) -> None:
    import cell2location
    from cell2location.models import Cell2location
    from lightning.pytorch.callbacks import Callback

    sample_seed = SEED + (1 if sample == "GLMF2" else 0)
    random.seed(sample_seed)
    np.random.seed(sample_seed)
    torch.manual_seed(sample_seed)
    torch.set_num_threads(THREADS_PER_SAMPLE)

    sample_dir = OUTPUT_ROOT / sample
    if sample_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {sample_dir}")
    sample_dir.mkdir(parents=True)
    status_path = sample_dir / "status.json"
    milestone_path = sample_dir / "training_milestones.tsv"
    checkpoint_dir = sample_dir / "checkpoints"
    checkpoint_dir.mkdir()

    def write_status(status: str, **extra: object) -> None:
        payload = {"sample": sample, "status": status, "updated_at": now(), **extra}
        temporary = status_path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, status_path)

    try:
        memory_gate("run_start", sample)
        input_path = SAMPLES[sample]
        input_hash = sha256(input_path)
        write_status("loading_input", input_sha256=input_hash)
        adata = ad.read_h5ad(input_path)
        adata = adata[:, genes].copy()
        x = adata.X
        values = x.data if issparse(x) else np.asarray(x).reshape(-1)
        if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values)):
            raise RuntimeError("Spatial X is not valid raw integer-like counts")
        cell_state_df = signatures.loc[genes].copy()
        adata.uns["reference_tymp_continuous_v1_cell2location"] = {
            "sample": sample,
            "input_h5ad": str(input_path),
            "input_sha256": input_hash,
            "reference": str(SIGNATURES),
            "reference_sha256": sha256(SIGNATURES),
            "aggregation": "patient mean raw counts -> study mean -> equal-weight study mean",
            "seed": sample_seed,
            "n_model_genes": len(genes),
        }

        Cell2location.setup_anndata(adata)
        model = Cell2location(
            adata,
            cell_state_df=cell_state_df,
            N_cells_per_location=N_CELLS_PER_LOCATION,
            detection_alpha=DETECTION_ALPHA,
        )

        class MilestoneRecorder(Callback):
            def on_train_epoch_end(self, trainer, pl_module) -> None:
                epoch = int(trainer.current_epoch) + 1
                if epoch % CHECKPOINT_EVERY_N_EPOCHS:
                    return
                metrics = {
                    key: float(value.detach().cpu())
                    for key, value in trainer.callback_metrics.items()
                    if hasattr(value, "detach") and value.numel() == 1
                }
                checkpoint = checkpoint_dir / f"epoch={epoch:05d}.ckpt"
                trainer.save_checkpoint(checkpoint)
                with milestone_path.open("a", encoding="utf-8") as handle:
                    handle.write(f"{now()}\t{epoch}\t{checkpoint}\t{json.dumps(metrics, sort_keys=True)}\n")
                write_status(
                    "training",
                    epoch=epoch,
                    metrics=metrics,
                    checkpoint=str(checkpoint),
                    n_spots=adata.n_obs,
                    n_genes=adata.n_vars,
                )

        milestone_path.write_text("recorded_at\tepoch\tcheckpoint\tmetrics_json\n", encoding="utf-8")
        memory_gate("training_start", sample)
        write_status(
            "training", epoch=0, n_spots=adata.n_obs, n_genes=adata.n_vars,
            max_epochs=MAX_EPOCHS, cell2location_version=cell2location.__version__,
        )
        model.train(
            max_epochs=MAX_EPOCHS,
            batch_size=None,
            train_size=1,
            accelerator="cpu",
            callbacks=[MilestoneRecorder()],
            enable_checkpointing=False,
            enable_progress_bar=False,
        )

        memory_gate("posterior_export_start", sample)
        write_status("exporting_posterior", epoch=MAX_EPOCHS)
        adata = model.export_posterior(
            adata,
            sample_kwargs={
                "num_samples": POSTERIOR_SAMPLES,
                "batch_size": min(1024, adata.n_obs),
                "accelerator": "cpu",
            },
        )
        model.save(sample_dir / "model", overwrite=False)
        output_h5ad = sample_dir / f"{sample}_tymp_continuous_v1_cell2location.h5ad"
        adata.write(output_h5ad)

        abundance_key = "q05_cell_abundance_w_sf"
        if abundance_key not in adata.obsm:
            raise KeyError(f"Missing posterior key: {abundance_key}")
        abundance = adata.obsm[abundance_key]
        if not isinstance(abundance, pd.DataFrame):
            abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=cell_state_df.columns)
        abundance.to_csv(sample_dir / "cell_abundance_q05_tymp_continuous_v1.tsv.gz", sep="\t")
        pd.DataFrame({
            "mean_q05": abundance.mean(axis=0),
            "median_q05": abundance.median(axis=0),
            "nonzero_fraction": (abundance > 0).mean(axis=0),
        }).to_csv(sample_dir / "cell_abundance_summary_tymp_continuous_v1.tsv", sep="\t")

        if sha256(input_path) != input_hash:
            raise RuntimeError("Input H5AD changed during mapping")
        write_status(
            "completed",
            epoch=MAX_EPOCHS,
            output_h5ad=str(output_h5ad),
            output_size_bytes=output_h5ad.stat().st_size,
            output_sha256=sha256(output_h5ad),
            input_sha256=input_hash,
        )
    except Exception as error:
        write_status(
            "failed",
            error_type=type(error).__name__,
            error=str(error),
            traceback=traceback.format_exc(),
        )
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--sample", choices=list(SAMPLES))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    signatures = load_signatures()
    genes = load_gene_template()
    if args.preflight_only:
        records = [preflight_one(sample, signatures, genes) for sample in SAMPLES]
        write_preflight(records, signatures)
        print(pd.DataFrame(records).to_string(index=False))
        return
    if not args.sample:
        raise SystemExit("--sample is required unless --preflight-only is used")
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    run_sample(args.sample, signatures, genes)


if __name__ == "__main__":
    main()
