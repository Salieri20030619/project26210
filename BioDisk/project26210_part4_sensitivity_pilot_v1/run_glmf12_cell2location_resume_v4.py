#!/usr/bin/env python3
"""Resume the interrupted GLMF1/GLMF2 reference-v3 mappings without overwriting v3."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import re
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


SEED = 20260805
N_CELLS_PER_LOCATION = 30
DETECTION_ALPHA = 20
MAX_EPOCHS = 10000
POSTERIOR_SAMPLES = 1000
CHECKPOINT_EVERY_N_EPOCHS = 1000
THREADS_PER_SAMPLE = 32
EXPECTED_REFERENCE_SHAPE = (12758, 14)
EXPECTED_MODEL_GENES = 8830
SAMPLE_RELATIVE_INPUTS = {
    "GLMF1": "results/part4_spatial/mapping/GLMF1_QC.h5ad",
    "GLMF2": "results/part4_spatial/mapping/GLMF2_QC.h5ad",
}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_name(f"{path.name}.{os.getpid()}.part")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    atomic_text(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def checkpoint_epoch(path: Path) -> int:
    match = re.fullmatch(r"epoch=(\d{5})\.ckpt", path.name)
    if not match:
        raise ValueError(f"Unexpected checkpoint name: {path}")
    return int(match.group(1))


def latest_checkpoint(directory: Path) -> Path | None:
    checkpoints = list(directory.glob("epoch=*.ckpt")) if directory.is_dir() else []
    return max(checkpoints, key=checkpoint_epoch) if checkpoints else None


def validate_checkpoint(path: Path) -> dict[str, object]:
    payload = torch.load(path, map_location="cpu", weights_only=False)
    recorded_epoch = int(payload.get("epoch", -1)) + 1
    named_epoch = checkpoint_epoch(path)
    if recorded_epoch != named_epoch:
        raise RuntimeError(
            f"Checkpoint epoch mismatch: filename={named_epoch}, payload={recorded_epoch}"
        )
    if not payload.get("state_dict"):
        raise RuntimeError(f"Checkpoint has no state_dict: {path}")
    if not payload.get("optimizer_states"):
        raise RuntimeError(f"Checkpoint has no optimizer state: {path}")
    result = {
        "path": str(path),
        "epoch": named_epoch,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "global_step": int(payload.get("global_step", -1)),
        "optimizer_state_count": len(payload["optimizer_states"]),
    }
    del payload
    return result


def memory_gate(project_root: Path, output_root: Path, stage: str, sample: str) -> None:
    resource_dir = output_root / "resource_logs"
    resource_dir.mkdir(parents=True, exist_ok=True)
    gate = project_root / "check_memory_gate.sh"
    result = subprocess.run([str(gate)], capture_output=True, text=True, check=False)
    message = (result.stdout or result.stderr).strip().replace("\n", " | ")
    with (resource_dir / f"{sample}.tsv").open("a", encoding="utf-8") as handle:
        handle.write(f"{now()}\t{stage}\t{message}\texit_code={result.returncode}\n")
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked {sample} at {stage}: {message}")


def load_signatures(work_root: Path) -> tuple[Path, pd.DataFrame]:
    path = work_root / (
        "cell2location_reference_v3/"
        "reference_v3_signatures_study_balanced_mean_counts.tsv.gz"
    )
    signatures = pd.read_csv(path, sep="\t", index_col=0)
    if signatures.shape != EXPECTED_REFERENCE_SHAPE:
        raise RuntimeError(f"Unexpected reference shape: {signatures.shape}")
    if signatures.index.has_duplicates or signatures.columns.has_duplicates:
        raise RuntimeError("Duplicate reference identifiers")
    values = signatures.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise RuntimeError("Invalid reference values")
    return path, signatures


def load_gene_template(project_root: Path) -> pd.Index:
    path = project_root / "results/part4_spatial/cell2loc/GLMF1_cell2loc.h5ad"
    template = ad.read_h5ad(path, backed="r")
    genes = template.var_names.astype(str).copy()
    template.file.close()
    if len(genes) != EXPECTED_MODEL_GENES or genes.has_duplicates:
        raise RuntimeError(f"Unexpected fixed comparison gene set: {len(genes)}")
    return genes


def output_root(work_root: Path) -> Path:
    return work_root / "glmf12_cell2location_reference_v4_resume"


def source_root(work_root: Path) -> Path:
    return work_root / "glmf12_cell2location_reference_v3"


def preflight(project_root: Path, work_root: Path) -> None:
    out = output_root(work_root)
    out.mkdir(parents=True, exist_ok=True)
    metadata_path = out / "run_metadata_v4.json"
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if existing.get("status") != "preflight_passed":
            raise RuntimeError(f"Existing preflight metadata is not PASS: {metadata_path}")
        print(json.dumps(existing, ensure_ascii=False, indent=2))
        return

    signatures_path, signatures = load_signatures(work_root)
    genes = load_gene_template(project_root)
    records: list[dict[str, object]] = []
    for sample, relative in SAMPLE_RELATIVE_INPUTS.items():
        input_path = project_root / relative
        obj = ad.read_h5ad(input_path, backed="r")
        missing_input = genes.difference(obj.var_names)
        missing_reference = genes.difference(signatures.index)
        obj.file.close()
        source = latest_checkpoint(source_root(work_root) / sample / "checkpoints")
        if source is None:
            raise FileNotFoundError(f"No source checkpoint for {sample}")
        source_info = validate_checkpoint(source)
        if missing_input.size or missing_reference.size:
            raise RuntimeError(f"Gene mismatch for {sample}")
        records.append(
            {
                "sample": sample,
                "input_h5ad": str(input_path),
                "input_size_bytes": input_path.stat().st_size,
                "input_sha256": sha256(input_path),
                "n_model_genes": len(genes),
                "source_checkpoint": source_info,
            }
        )

    metadata = {
        "status": "preflight_passed",
        "created_at": now(),
        "purpose": "resume interrupted reference-v3 mappings in a non-overwriting v4 directory",
        "project_root": str(project_root),
        "work_root": str(work_root),
        "reference_signatures": str(signatures_path),
        "reference_sha256": sha256(signatures_path),
        "reference_shape": list(signatures.shape),
        "samples": records,
        "parameters": {
            "N_cells_per_location": N_CELLS_PER_LOCATION,
            "detection_alpha": DETECTION_ALPHA,
            "max_epochs_total": MAX_EPOCHS,
            "posterior_samples": POSTERIOR_SAMPLES,
            "checkpoint_every_n_epochs": CHECKPOINT_EVERY_N_EPOCHS,
            "threads_per_sample": THREADS_PER_SAMPLE,
            "execution": "serial GLMF1 then GLMF2",
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
    atomic_json(metadata_path, metadata)
    pd.json_normalize(records).to_csv(out / "preflight_glmf12_v4.tsv", sep="\t", index=False)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


def choose_resume_checkpoint(work_root: Path, sample: str, sample_dir: Path) -> Path:
    current = latest_checkpoint(sample_dir / "checkpoints")
    source = latest_checkpoint(source_root(work_root) / sample / "checkpoints")
    candidates = [path for path in (current, source) if path is not None]
    if not candidates:
        raise FileNotFoundError(f"No checkpoint available for {sample}")
    return max(candidates, key=checkpoint_epoch)


def run_sample(project_root: Path, work_root: Path, sample: str) -> None:
    import cell2location
    from cell2location.models import Cell2location
    from lightning.pytorch.callbacks import Callback

    out = output_root(work_root)
    metadata_path = out / "run_metadata_v4.json"
    if not metadata_path.exists():
        raise RuntimeError("Run preflight before sample training")

    sample_seed = SEED + (1 if sample == "GLMF2" else 0)
    random.seed(sample_seed)
    np.random.seed(sample_seed)
    torch.manual_seed(sample_seed)
    torch.set_num_threads(THREADS_PER_SAMPLE)

    sample_dir = out / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = sample_dir / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    status_path = sample_dir / "status.json"
    milestone_path = sample_dir / "training_milestones.tsv"
    if not milestone_path.exists():
        atomic_text(milestone_path, "recorded_at\tepoch\tcheckpoint\tmetrics_json\n")

    if status_path.exists():
        existing_status = json.loads(status_path.read_text(encoding="utf-8"))
        if existing_status.get("status") == "completed":
            print(f"{sample} already completed; no action taken")
            return

    forbidden = [
        sample_dir / f"{sample}_reference_v4_cell2location.h5ad",
        sample_dir / f"{sample}_reference_v4_cell2location.h5ad.part",
        sample_dir / "model",
        sample_dir / "model.part",
    ]
    if any(path.exists() for path in forbidden):
        raise FileExistsError(f"Refusing ambiguous export state: {[str(p) for p in forbidden if p.exists()]}")

    def write_status(status: str, **extra: object) -> None:
        atomic_json(status_path, {"sample": sample, "status": status, "updated_at": now(), **extra})

    try:
        signatures_path, signatures = load_signatures(work_root)
        genes = load_gene_template(project_root)
        resume_path = choose_resume_checkpoint(work_root, sample, sample_dir)
        resume_info = validate_checkpoint(resume_path)
        if int(resume_info["epoch"]) >= MAX_EPOCHS:
            resume_label = "training_complete_checkpoint"
        else:
            resume_label = "training_resume_checkpoint"

        memory_gate(project_root, out, "run_start", sample)
        input_path = project_root / SAMPLE_RELATIVE_INPUTS[sample]
        input_hash = sha256(input_path)
        write_status("loading_input", resume_checkpoint=resume_info, input_sha256=input_hash)
        adata = ad.read_h5ad(input_path)
        adata = adata[:, genes].copy()
        values = adata.X.data if issparse(adata.X) else np.asarray(adata.X).reshape(-1)
        if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values)):
            raise RuntimeError("Spatial X is not valid raw integer-like counts")

        cell_state_df = signatures.loc[genes].copy()
        adata.uns["reference_v4_cell2location_resume"] = {
            "sample": sample,
            "input_h5ad": str(input_path),
            "input_sha256": input_hash,
            "reference": str(signatures_path),
            "reference_sha256": sha256(signatures_path),
            "source_checkpoint": resume_info,
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
                checkpoint = checkpoint_dir / f"epoch={epoch:05d}.ckpt"
                if checkpoint.exists():
                    return
                metrics = {
                    key: float(value.detach().cpu())
                    for key, value in trainer.callback_metrics.items()
                    if hasattr(value, "detach") and value.numel() == 1
                }
                temporary = checkpoint.with_name(f"{checkpoint.name}.{os.getpid()}.part")
                trainer.save_checkpoint(temporary)
                os.replace(temporary, checkpoint)
                with milestone_path.open("a", encoding="utf-8") as handle:
                    handle.write(
                        f"{now()}\t{epoch}\t{checkpoint}\t"
                        f"{json.dumps(metrics, sort_keys=True)}\n"
                    )
                write_status(
                    "training",
                    epoch=epoch,
                    metrics=metrics,
                    checkpoint=str(checkpoint),
                    resumed_from=resume_info,
                    n_spots=adata.n_obs,
                    n_genes=adata.n_vars,
                )

        memory_gate(project_root, out, resume_label, sample)
        write_status(
            "training",
            epoch=int(resume_info["epoch"]),
            resumed_from=resume_info,
            max_epochs_total=MAX_EPOCHS,
            n_spots=adata.n_obs,
            n_genes=adata.n_vars,
            cell2location_version=cell2location.__version__,
        )
        model.train(
            max_epochs=MAX_EPOCHS,
            batch_size=None,
            train_size=1,
            accelerator="cpu",
            callbacks=[MilestoneRecorder()],
            enable_checkpointing=False,
            enable_progress_bar=False,
            ckpt_path=str(resume_path),
        )

        memory_gate(project_root, out, "posterior_export_start", sample)
        write_status("exporting_posterior", epoch=MAX_EPOCHS, resumed_from=resume_info)
        adata = model.export_posterior(
            adata,
            sample_kwargs={
                "num_samples": POSTERIOR_SAMPLES,
                "batch_size": min(1024, adata.n_obs),
                "accelerator": "cpu",
            },
        )

        model_part = sample_dir / "model.part"
        model_final = sample_dir / "model"
        model.save(model_part, overwrite=False)
        os.replace(model_part, model_final)

        output_h5ad = sample_dir / f"{sample}_reference_v4_cell2location.h5ad"
        output_part = sample_dir / f"{sample}_reference_v4_cell2location.h5ad.part"
        adata.write_h5ad(output_part)
        os.replace(output_part, output_h5ad)

        abundance_key = "q05_cell_abundance_w_sf"
        if abundance_key not in adata.obsm:
            raise KeyError(f"Missing posterior key: {abundance_key}")
        abundance = adata.obsm[abundance_key]
        if not isinstance(abundance, pd.DataFrame):
            abundance = pd.DataFrame(abundance, index=adata.obs_names, columns=cell_state_df.columns)
        abundance_values = abundance.to_numpy(dtype=float)
        if abundance.shape != (adata.n_obs, cell_state_df.shape[1]):
            raise RuntimeError(f"Unexpected abundance shape: {abundance.shape}")
        if not np.isfinite(abundance_values).all() or (abundance_values < 0).any():
            raise RuntimeError("Invalid q05 abundance values")

        abundance_path = sample_dir / "cell_abundance_q05_v4.tsv.gz"
        abundance_part = sample_dir / "cell_abundance_q05_v4.tsv.gz.part"
        abundance.to_csv(abundance_part, sep="\t", compression="gzip")
        os.replace(abundance_part, abundance_path)

        summary_path = sample_dir / "cell_abundance_summary_v4.tsv"
        summary_part = sample_dir / "cell_abundance_summary_v4.tsv.part"
        pd.DataFrame(
            {
                "mean_q05": abundance.mean(axis=0),
                "median_q05": abundance.median(axis=0),
                "nonzero_fraction": (abundance > 0).mean(axis=0),
            }
        ).to_csv(summary_part, sep="\t")
        os.replace(summary_part, summary_path)

        check = ad.read_h5ad(output_h5ad, backed="r")
        validation = {
            "n_spots": int(check.n_obs),
            "n_genes": int(check.n_vars),
            "q05_shape": list(check.obsm[abundance_key].shape),
            "factor_names": list(map(str, check.uns["mod"]["factor_names"])),
        }
        check.file.close()
        if validation["n_spots"] != adata.n_obs or validation["n_genes"] != EXPECTED_MODEL_GENES:
            raise RuntimeError(f"Written H5AD validation failed: {validation}")
        if len(validation["factor_names"]) != EXPECTED_REFERENCE_SHAPE[1]:
            raise RuntimeError(f"Factor validation failed: {validation}")
        if sha256(input_path) != input_hash:
            raise RuntimeError("Input H5AD changed during mapping")

        deliverables = [output_h5ad, abundance_path, summary_path]
        manifest_rows = [
            {"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in deliverables
        ]
        manifest_path = sample_dir / "file_manifest_v4.tsv"
        manifest_part = sample_dir / "file_manifest_v4.tsv.part"
        pd.DataFrame(manifest_rows).to_csv(manifest_part, sep="\t", index=False)
        os.replace(manifest_part, manifest_path)

        write_status(
            "completed",
            epoch=MAX_EPOCHS,
            resumed_from=resume_info,
            output_h5ad=str(output_h5ad),
            output_size_bytes=output_h5ad.stat().st_size,
            output_sha256=sha256(output_h5ad),
            input_sha256=input_hash,
            validation=validation,
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
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--sample", choices=list(SAMPLE_RELATIVE_INPUTS))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve(strict=True)
    work_root = args.work_root.resolve(strict=True)
    if args.preflight_only:
        preflight(project_root, work_root)
        return
    if not args.sample:
        raise SystemExit("--sample is required unless --preflight-only is used")
    run_sample(project_root, work_root, args.sample)


if __name__ == "__main__":
    main()
