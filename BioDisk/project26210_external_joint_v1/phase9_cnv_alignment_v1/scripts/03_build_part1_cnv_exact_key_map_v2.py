#!/usr/bin/env python3
"""Build a validated cell-level Part 1 CNV map using an exact GSM + 10x-barcode key."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


BARCODE_RE = re.compile(r"([ACGTN]+-[0-9]+)$", re.IGNORECASE)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_tsv(frame: pd.DataFrame, final: Path, compression: str | None = None) -> None:
    if final.exists():
        raise FileExistsError(f"No-overwrite gate: {final}")
    part = final.with_name(final.name + ".part")
    if part.exists():
        raise FileExistsError(f"Incomplete file already exists: {part}")
    frame.to_csv(part, sep="\t", index=False, compression=compression, na_rep="NA")
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
    map_path = out_dir / "part1_cnv_exact_cell_map_v2.tsv.gz"
    summary_path = out_dir / "part1_cnv_exact_key_coverage_v2.tsv"
    status_path = out_dir / "PHASE9_PART1_CNV_EXACT_KEY_STATUS"
    metadata_path = out_dir / "run_metadata_v2.json"
    for path in (map_path, summary_path, status_path, metadata_path):
        if path.exists():
            raise FileExistsError(f"No-overwrite gate: {path}")

    ext = pd.read_csv(crosswalk_path, sep="\t", dtype=str, compression="infer")
    required_ext = {
        "cell_id", "dataset", "patient_id", "src__gsm", "barcode_original",
        "phase8_spatial_reference_state", "phase8_spatial_mapping_status",
        "phase8_f1_projection_score", "phase8_f1_percentile_patient",
        "phase8_f1_eligible_patient",
    }
    missing_ext = required_ext.difference(ext.columns)
    if missing_ext:
        raise RuntimeError(f"Crosswalk missing columns: {sorted(missing_ext)}")
    if len(ext) != 585119 or ext["cell_id"].duplicated().any():
        raise RuntimeError("External cell count or cell_id uniqueness gate failed")

    def core_barcode(value: str) -> str | None:
        match = BARCODE_RE.search(str(value))
        return match.group(1).upper() if match else None

    ext["core_10x_barcode"] = ext["barcode_original"].map(core_barcode)
    ext["part1_exact_key"] = ext["src__gsm"].astype("string") + "_" + ext["core_10x_barcode"].astype("string")
    valid_key = ext["src__gsm"].notna() & ext["core_10x_barcode"].notna()
    duplicate_valid_keys = ext.loc[valid_key, "part1_exact_key"].duplicated(keep=False)
    if duplicate_valid_keys.any():
        examples = ext.loc[valid_key].loc[duplicate_valid_keys, "part1_exact_key"].head().tolist()
        raise RuntimeError(f"External GSM+barcode key is not unique: {examples}")

    a = ad.read_h5ad(h5ad_path, backed="r")
    obs = a.obs.copy()
    obs.index = obs.index.astype(str)
    if obs.index.duplicated().any():
        raise RuntimeError("Part 1 obs_names are not unique")
    required_obs = {
        "study", "curated_patient_id", "geo_accession", "cnv_assessed_v1",
        "infercnv_call_v1", "copykat_call_v1", "malignancy_status_v1",
        "malignant_conservative_v1", "malignancy_evidence_v1",
        "eligible_for_subclone_selection_v1",
    }
    missing_obs = required_obs.difference(obs.columns)
    if missing_obs:
        raise RuntimeError(f"Part 1 obs missing columns: {sorted(missing_obs)}")

    ext["part1_match"] = valid_key & ext["part1_exact_key"].isin(obs.index)
    matched_ext = ext.loc[ext["part1_match"]].copy()
    obs_cols = [
        "study", "curated_patient_id", "source_patient_id", "geo_accession",
        "cell_type_broad_v1", "cnv_assessed_v1", "cnv_analysis_scope_v1",
        "infercnv_call_v1", "copykat_call_v1", "malignancy_status_v1",
        "malignant_conservative_v1", "malignancy_evidence_v1",
        "infercnv_unit_v1", "copykat_unit_v1",
        "included_in_primary_malignancy_summary_v1",
        "eligible_for_subclone_selection_v1",
    ]
    aligned_obs = obs.loc[matched_ext["part1_exact_key"], obs_cols].reset_index(names="part1_obs_name")
    aligned_obs.index = matched_ext.index
    for col in obs_cols + ["part1_obs_name"]:
        matched_ext[f"part1__{col}"] = aligned_obs[col].astype(str)

    dataset_ok = matched_ext["dataset"].astype(str).eq(matched_ext["part1__study"].astype(str))
    patient_expected = matched_ext["dataset"].astype(str) + "__" + matched_ext["patient_id"].astype(str)
    patient_ok = patient_expected.eq(matched_ext["part1__curated_patient_id"].astype(str))
    gsm_ok = matched_ext["src__gsm"].astype(str).eq(matched_ext["part1__geo_accession"].astype(str))
    if not dataset_ok.all() or not patient_ok.all() or not gsm_ok.all():
        raise RuntimeError(
            f"Identity reconciliation failed: dataset={dataset_ok.mean():.6f}, "
            f"patient={patient_ok.mean():.6f}, gsm={gsm_ok.mean():.6f}"
        )

    tumor_like = ext["phase8_spatial_reference_state"].astype(str).str.contains("F1-", regex=False)
    ext["phase8_tumor_like"] = tumor_like
    ext["matched_cnv_assessed"] = False
    ext["matched_conservative_malignant"] = False
    if len(matched_ext):
        assessed = matched_ext["part1__cnv_assessed_v1"].str.lower().eq("true")
        malignant = matched_ext["part1__malignant_conservative_v1"].str.lower().eq("true")
        ext.loc[matched_ext.index, "matched_cnv_assessed"] = assessed.to_numpy()
        ext.loc[matched_ext.index, "matched_conservative_malignant"] = malignant.to_numpy()

    summary = (
        ext.groupby(["dataset", "patient_id"], observed=True, dropna=False)
        .agg(
            n_external_cells=("cell_id", "size"),
            n_valid_exact_keys=("part1_exact_key", lambda s: int(s.notna().sum())),
            n_part1_exact_matches=("part1_match", "sum"),
            n_phase8_tumor_like=("phase8_tumor_like", "sum"),
            n_phase8_tumor_like_part1_matches=("part1_match", lambda s: int((s & ext.loc[s.index, "phase8_tumor_like"]).sum())),
            n_part1_cnv_assessed=("matched_cnv_assessed", "sum"),
            n_part1_conservative_malignant=("matched_conservative_malignant", "sum"),
        )
        .reset_index()
    )
    summary["fraction_external_exact_match"] = summary["n_part1_exact_matches"] / summary["n_external_cells"]
    summary["fraction_tumor_like_exact_match"] = summary["n_phase8_tumor_like_part1_matches"] / summary["n_phase8_tumor_like"].replace(0, pd.NA)

    map_cols = [
        "cell_id", "dataset", "patient_id", "src__gsm", "core_10x_barcode",
        "part1_exact_key", "phase8_spatial_reference_state", "phase8_spatial_mapping_status",
        "phase8_f1_projection_score", "phase8_f1_percentile_patient",
        "phase8_f1_eligible_patient",
    ] + [f"part1__{c}" for c in ["part1_obs_name"] + obs_cols]
    atomic_tsv(matched_ext[map_cols], map_path, compression="gzip")
    atomic_tsv(summary, summary_path)

    counts = {
        "external_cells": int(len(ext)),
        "valid_exact_keys": int(valid_key.sum()),
        "part1_exact_matches": int(ext["part1_match"].sum()),
        "matched_datasets": sorted(matched_ext["dataset"].unique().tolist()),
        "matched_patients": int(matched_ext[["dataset", "patient_id"]].drop_duplicates().shape[0]),
        "matched_cnv_assessed": int(ext["matched_cnv_assessed"].sum()),
        "matched_conservative_malignant": int(ext["matched_conservative_malignant"].sum()),
        "matched_phase8_tumor_like": int((ext["part1_match"] & ext["phase8_tumor_like"]).sum()),
    }
    metadata = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "h5ad": str(h5ad_path),
        "h5ad_bytes": h5ad_path.stat().st_size,
        "h5ad_sha256": sha256(h5ad_path),
        "crosswalk": str(crosswalk_path),
        "crosswalk_bytes": crosswalk_path.stat().st_size,
        "crosswalk_sha256": sha256(crosswalk_path),
        "counts": counts,
        "outputs": {
            "cell_map": {"path": str(map_path), "bytes": map_path.stat().st_size, "sha256": sha256(map_path)},
            "coverage": {"path": str(summary_path), "bytes": summary_path.stat().st_size, "sha256": sha256(summary_path)},
        },
    }
    metadata_part = metadata_path.with_name(metadata_path.name + ".part")
    metadata_part.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(metadata_part, metadata_path)

    status_lines = ["PASS"] + [f"{key}={value}" for key, value in counts.items()]
    status_lines += [
        f"cell_map_sha256={metadata['outputs']['cell_map']['sha256']}",
        f"coverage_sha256={metadata['outputs']['coverage']['sha256']}",
    ]
    status_part = status_path.with_name(status_path.name + ".part")
    status_part.write_text("\n".join(status_lines) + "\n", encoding="utf-8")
    os.replace(status_part, status_path)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
