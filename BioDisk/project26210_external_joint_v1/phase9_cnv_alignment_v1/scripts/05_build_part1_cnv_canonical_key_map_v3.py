#!/usr/bin/env python3
"""Validated Part 1 CNV transfer for identical cells using unique GSM+10x keys."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd

BARCODE_RE = re.compile(r"([ACGTN]+-[0-9]+)$", re.IGNORECASE)


def core_barcode(value: object) -> str | None:
    match = BARCODE_RE.search(str(value))
    return match.group(1).upper() if match else None


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, final: Path, compression: str | None = None) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    if part.exists():
        raise FileExistsError(f"Incomplete file exists: {part}")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", compression=compression)
    os.replace(part, final)


def atomic_text(text: str, final: Path) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    part.write_text(text, encoding="utf-8")
    os.replace(part, final)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--crosswalk", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad).resolve()
    crosswalk_path = Path(args.crosswalk).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "cell_map": out_dir / "part1_cnv_identical_cell_map_v3.tsv.gz",
        "coverage": out_dir / "part1_cnv_coverage_by_patient_v3.tsv",
        "state_cross": out_dir / "phase8_state_by_part1_cnv_status_v3.tsv",
        "status": out_dir / "PHASE9_PART1_CNV_CANONICAL_KEY_STATUS",
        "metadata": out_dir / "run_metadata_v3.json",
    }
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("No-overwrite gate: one or more v3 outputs already exist")

    ext = pd.read_csv(crosswalk_path, sep="\t", dtype=str, compression="infer")
    required_ext = {
        "cell_id", "dataset", "patient_id", "sample_id", "src__gsm", "barcode_original",
        "phase8_spatial_reference_state", "phase8_spatial_mapping_status",
        "phase8_f1_projection_score", "phase8_f1_percentile_patient",
        "phase8_f1_eligible_patient",
    }
    missing = required_ext.difference(ext.columns)
    if missing:
        raise RuntimeError(f"Missing external columns: {sorted(missing)}")
    if len(ext) != 585119 or ext["cell_id"].duplicated().any():
        raise RuntimeError("External count/uniqueness gate failed")

    ext["core_10x_barcode"] = ext["barcode_original"].map(core_barcode)
    ext["canonical_key"] = ext["src__gsm"].astype("string") + "|" + ext["core_10x_barcode"].astype("string")
    ext_valid = ext["src__gsm"].notna() & ext["core_10x_barcode"].notna()
    if ext.loc[ext_valid, "canonical_key"].duplicated(keep=False).any():
        raise RuntimeError("External canonical keys are not unique")

    a = ad.read_h5ad(h5ad_path, backed="r")
    obs = a.obs.copy()
    obs["part1_obs_name"] = obs.index.astype(str)
    obs["core_10x_barcode"] = obs["part1_obs_name"].map(core_barcode)
    obs["canonical_key"] = obs["geo_accession"].astype("string") + "|" + obs["core_10x_barcode"].astype("string")
    obs_valid = obs["geo_accession"].notna() & obs["core_10x_barcode"].notna()
    if obs.loc[obs_valid, "canonical_key"].duplicated(keep=False).any():
        raise RuntimeError("Part 1 canonical keys are not unique")
    obs_by_key = obs.loc[obs_valid].set_index("canonical_key", drop=False)

    ext["part1_identical_cell_match"] = ext_valid & ext["canonical_key"].isin(obs_by_key.index)
    matched = ext.loc[ext["part1_identical_cell_match"]].copy()
    cnv_fields = [
        "part1_obs_name", "study", "curated_patient_id", "source_patient_id",
        "geo_accession", "cell_type_broad_v1", "cnv_assessed_v1",
        "cnv_analysis_scope_v1", "infercnv_call_v1", "copykat_call_v1",
        "malignancy_status_v1", "malignant_conservative_v1",
        "malignancy_evidence_v1", "infercnv_unit_v1", "copykat_unit_v1",
        "included_in_primary_malignancy_summary_v1",
        "eligible_for_subclone_selection_v1",
    ]
    aligned = obs_by_key.loc[matched["canonical_key"], cnv_fields].copy()
    aligned.index = matched.index
    for field in cnv_fields:
        matched[f"part1__{field}"] = aligned[field].astype(str)

    dataset_ok = matched["dataset"].eq(matched["part1__study"])
    patient_ok = (matched["dataset"] + "__" + matched["patient_id"]).eq(matched["part1__curated_patient_id"])
    gsm_ok = matched["src__gsm"].eq(matched["part1__geo_accession"])
    if not (dataset_ok.all() and patient_ok.all() and gsm_ok.all()):
        raise RuntimeError(
            f"Identity gate failed: dataset={dataset_ok.sum()}/{len(matched)}, "
            f"patient={patient_ok.sum()}/{len(matched)}, gsm={gsm_ok.sum()}/{len(matched)}"
        )
    if matched["part1__part1_obs_name"].duplicated().any():
        raise RuntimeError("A Part 1 cell was mapped more than once")

    ext["phase8_tumor_like"] = ext["phase8_spatial_reference_state"].astype(str).str.startswith("F1-")
    ext["matched_cnv_assessed"] = False
    ext["matched_conservative_malignant"] = False
    assessed = matched["part1__cnv_assessed_v1"].str.lower().eq("true")
    malignant = matched["part1__malignant_conservative_v1"].str.lower().eq("true")
    ext.loc[matched.index, "matched_cnv_assessed"] = assessed.to_numpy()
    ext.loc[matched.index, "matched_conservative_malignant"] = malignant.to_numpy()

    coverage = (
        ext.groupby(["dataset", "patient_id"], observed=True, dropna=False)
        .agg(
            n_external_cells=("cell_id", "size"),
            n_identical_part1_cells=("part1_identical_cell_match", "sum"),
            n_phase8_tumor_like=("phase8_tumor_like", "sum"),
            n_phase8_tumor_like_identical_part1=(
                "part1_identical_cell_match",
                lambda s: int((s & ext.loc[s.index, "phase8_tumor_like"]).sum()),
            ),
            n_part1_cnv_assessed=("matched_cnv_assessed", "sum"),
            n_part1_conservative_malignant=("matched_conservative_malignant", "sum"),
        )
        .reset_index()
    )
    coverage["fraction_external_identical_part1"] = coverage["n_identical_part1_cells"] / coverage["n_external_cells"]
    coverage["fraction_tumor_like_identical_part1"] = (
        coverage["n_phase8_tumor_like_identical_part1"] / coverage["n_phase8_tumor_like"].replace(0, pd.NA)
    )

    state_cross = (
        matched.groupby(
            ["dataset", "phase8_spatial_reference_state", "part1__malignancy_status_v1", "part1__malignant_conservative_v1"],
            observed=True,
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
    )

    map_columns = [
        "cell_id", "dataset", "patient_id", "sample_id", "src__gsm",
        "core_10x_barcode", "canonical_key", "phase8_spatial_reference_state",
        "phase8_spatial_mapping_status", "phase8_f1_projection_score",
        "phase8_f1_percentile_patient", "phase8_f1_eligible_patient",
    ] + [f"part1__{field}" for field in cnv_fields]
    atomic_tsv(matched[map_columns], paths["cell_map"], compression="gzip")
    atomic_tsv(coverage, paths["coverage"])
    atomic_tsv(state_cross, paths["state_cross"])

    counts = {
        "external_cells": int(len(ext)),
        "external_valid_unique_keys": int(ext_valid.sum()),
        "part1_valid_unique_keys": int(obs_valid.sum()),
        "identical_part1_cells": int(len(matched)),
        "identity_dataset_pass": int(dataset_ok.sum()),
        "identity_patient_pass": int(patient_ok.sum()),
        "identity_gsm_pass": int(gsm_ok.sum()),
        "matched_datasets": sorted(matched["dataset"].unique().tolist()),
        "matched_patients": int(matched[["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "matched_phase8_tumor_like": int(matched["phase8_spatial_reference_state"].str.startswith("F1-").sum()),
        "matched_cnv_assessed": int(assessed.sum()),
        "matched_conservative_malignant": int(malignant.sum()),
    }
    metadata = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": "unique exact canonical key: GEO accession/GSM plus terminal 10x barcode; dataset, patient and GSM reconciliation required",
        "python": platform.python_version(),
        "anndata": importlib.metadata.version("anndata"),
        "pandas": pd.__version__,
        "inputs": {
            "h5ad": {"path": str(h5ad_path), "bytes": h5ad_path.stat().st_size, "sha256": file_sha256(h5ad_path)},
            "crosswalk": {"path": str(crosswalk_path), "bytes": crosswalk_path.stat().st_size, "sha256": file_sha256(crosswalk_path)},
        },
        "counts": counts,
        "outputs": {
            key: {"path": str(paths[key]), "bytes": paths[key].stat().st_size, "sha256": file_sha256(paths[key])}
            for key in ("cell_map", "coverage", "state_cross")
        },
    }
    atomic_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", paths["metadata"])
    status_lines = ["PASS"] + [f"{key}={value}" for key, value in counts.items()]
    status_lines += [f"{key}_sha256={value['sha256']}" for key, value in metadata["outputs"].items()]
    atomic_text("\n".join(status_lines) + "\n", paths["status"])
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
