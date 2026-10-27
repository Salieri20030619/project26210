#!/usr/bin/env python3
"""Resume 12 interrupted TJRCC mappings as non-overwriting auxiliary validation."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import os
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch


EXPECTED_SAMPLES = (
    "H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T",
    "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T",
)
OUTPUT_NAME = "tjrcc12_cell2location_reference_v4r1_aux_resume"
SOURCE_NAME = "tjrcc12_cell2location_reference_v3"
BASE_ENGINE = "run_glmf12_cell2location_resume_v4r1.py"
MAX_CONCURRENT_SAMPLES = 4
THREADS_PER_SAMPLE = 12
SEED_BASE = 20260815


def load_base(work_root: Path):
    path = work_root / BASE_ENGINE
    spec = importlib.util.spec_from_file_location("project26210_cell2loc_resume_base", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load recovery engine: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def input_paths(project_root: Path) -> dict[str, Path]:
    manifest_path = project_root / (
        "external/ccRCC_spatial/TJ_RCC_figshare_24599295/"
        "analysis_cell2location_full_v1/input_manifest.tsv"
    )
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    paths = {
        str(row["sample"]): Path(str(row["path"])).resolve()
        for _, row in manifest.iterrows()
    }
    if set(paths) != set(EXPECTED_SAMPLES) or len(paths) != 12:
        raise RuntimeError(f"Unexpected TJRCC manifest samples: {sorted(paths)}")
    return {sample: paths[sample] for sample in EXPECTED_SAMPLES}


def configure_base(base, project_root: Path, work_root: Path, paths: dict[str, Path]) -> None:
    base.SAMPLE_RELATIVE_INPUTS = paths
    base.THREADS_PER_SAMPLE = THREADS_PER_SAMPLE
    base.output_root = lambda _: work_root / OUTPUT_NAME
    base.source_root = lambda _: work_root / SOURCE_NAME

    # TJRCC inputs are older mapped H5ADs.  For training, retain only raw X,
    # spot/feature metadata and spatial geometry so no old posterior or scvi
    # registry can enter the rebuilt mapping.  Backed preflight/final reloads
    # still use the original reader unchanged.
    original_read = ad.read_h5ad
    input_set = {path.resolve() for path in paths.values()}
    genes = base.load_gene_template(project_root)

    def clean_read(filename, *args, **kwargs):
        path = Path(filename).resolve()
        backed = kwargs.get("backed")
        obj = original_read(filename, *args, **kwargs)
        if path not in input_set or backed is not None:
            return obj
        selected = obj[:, genes]
        clean = ad.AnnData(
            X=selected.X.copy(),
            obs=obj.obs.copy(),
            var=selected.var.copy(),
        )
        if "spatial" not in obj.obsm:
            raise KeyError(f"Missing spatial coordinates: {path}")
        clean.obsm["spatial"] = np.asarray(obj.obsm["spatial"]).copy()
        if "spatial" in obj.uns:
            clean.uns["spatial"] = copy.deepcopy(obj.uns["spatial"])
        return clean

    base.ad.read_h5ad = clean_read


def preflight(base, project_root: Path, work_root: Path, paths: dict[str, Path]) -> None:
    out = work_root / OUTPUT_NAME
    out.mkdir(parents=True, exist_ok=True)
    # The tested base recovery engine uses this filename internally even in
    # its non-overwriting v4r1 output directory.
    metadata_path = out / "run_metadata_v4.json"
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if existing.get("status") != "preflight_passed":
            raise RuntimeError("Existing preflight metadata is not PASS")
        print(json.dumps(existing, ensure_ascii=False, indent=2))
        return

    signatures_path, signatures = base.load_signatures(work_root)
    genes = base.load_gene_template(project_root)
    source_root = work_root / SOURCE_NAME
    prior = pd.read_csv(source_root / "preflight_tjrcc12_v3.tsv", sep="\t", dtype=str)
    prior = prior.set_index("sample", drop=False)
    records: list[dict[str, object]] = []
    for sample, path in paths.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(path)
        obj = ad.read_h5ad(path, backed="r")
        missing = genes.difference(obj.var_names)
        n_spots, n_input_genes = int(obj.n_obs), int(obj.n_vars)
        spatial_shape = list(np.asarray(obj.obsm["spatial"]).shape)
        obj.file.close()
        if missing.size or spatial_shape != [n_spots, 2]:
            raise RuntimeError(f"Input structure failed for {sample}")
        current_sha = base.sha256(path)
        if sample not in prior.index or current_sha != str(prior.loc[sample, "input_sha256"]):
            raise RuntimeError(f"Input hash changed for {sample}")
        checkpoint = base.latest_checkpoint(source_root / sample / "checkpoints")
        if checkpoint is None:
            raise FileNotFoundError(f"No source checkpoint for {sample}")
        checkpoint_info = base.validate_checkpoint(checkpoint)
        source_status = json.loads((source_root / sample / "status.json").read_text(encoding="utf-8"))
        if int(source_status.get("epoch", -1)) != int(checkpoint_info["epoch"]):
            raise RuntimeError(f"Source status/checkpoint epoch mismatch for {sample}")
        records.append({
            "sample": sample,
            "input_h5ad": str(path),
            "input_size_bytes": path.stat().st_size,
            "input_sha256": current_sha,
            "n_spots": n_spots,
            "n_input_genes": n_input_genes,
            "n_model_genes": len(genes),
            "spatial_shape": "x".join(map(str, spatial_shape)),
            "source_checkpoint": str(checkpoint),
            "source_checkpoint_epoch": int(checkpoint_info["epoch"]),
            "source_checkpoint_size_bytes": int(checkpoint_info["size_bytes"]),
            "source_checkpoint_sha256": str(checkpoint_info["sha256"]),
            "source_checkpoint_global_step": int(checkpoint_info["global_step"]),
            "source_optimizer_state_count": int(checkpoint_info["optimizer_state_count"]),
        })

    if len(records) != 12 or sum(int(row["n_spots"]) for row in records) != 34609:
        raise RuntimeError("TJRCC12 sample/spot conservation gate failed")
    metadata = {
        "status": "preflight_passed",
        "created_at": base.now(),
        "purpose": "auxiliary 12-sample spatial validation; does not replace GLMF1/GLMF2",
        "project_root": str(project_root),
        "work_root": str(work_root),
        "source_root": str(source_root),
        "output_root": str(out),
        "reference_signatures": str(signatures_path),
        "reference_sha256": base.sha256(signatures_path),
        "reference_shape": list(signatures.shape),
        "samples": records,
        "parameters": {
            "N_cells_per_location": base.N_CELLS_PER_LOCATION,
            "detection_alpha": base.DETECTION_ALPHA,
            "max_epochs_total": base.MAX_EPOCHS,
            "posterior_samples": base.POSTERIOR_SAMPLES,
            "checkpoint_every_n_epochs": base.CHECKPOINT_EVERY_N_EPOCHS,
            "threads_per_sample": THREADS_PER_SAMPLE,
            "maximum_concurrent_samples": MAX_CONCURRENT_SAMPLES,
            "seed_base": SEED_BASE,
        },
        "input_policy": "fresh AnnData from raw X, obs, var, spatial obsm and spatial uns only",
        "environment": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "torch": torch.__version__,
            "cuda_available": torch.cuda.is_available(),
        },
        "scripts": {
            "wrapper": str(Path(__file__).resolve()),
            "wrapper_sha256": base.sha256(Path(__file__).resolve()),
            "base_engine": str((work_root / BASE_ENGINE).resolve()),
            "base_engine_sha256": base.sha256((work_root / BASE_ENGINE).resolve()),
        },
    }
    base.atomic_json(metadata_path, metadata)
    table_path = out / "preflight_tjrcc12_v4r1.tsv"
    table_part = out / "preflight_tjrcc12_v4r1.tsv.part"
    pd.DataFrame(records).to_csv(table_part, sep="\t", index=False)
    os.replace(table_part, table_path)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--sample", choices=EXPECTED_SAMPLES)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve(strict=True)
    work_root = args.work_root.resolve(strict=True)
    paths = input_paths(project_root)
    base = load_base(work_root)
    configure_base(base, project_root, work_root, paths)
    if args.preflight_only:
        preflight(base, project_root, work_root, paths)
        return
    if not args.sample:
        raise SystemExit("--sample is required unless --preflight-only is used")
    metadata_path = work_root / OUTPUT_NAME / "run_metadata_v4.json"
    if not metadata_path.exists():
        raise RuntimeError("Run preflight before training")
    base.SEED = SEED_BASE + EXPECTED_SAMPLES.index(args.sample)
    base.run_sample(project_root, work_root, args.sample)


if __name__ == "__main__":
    main()
