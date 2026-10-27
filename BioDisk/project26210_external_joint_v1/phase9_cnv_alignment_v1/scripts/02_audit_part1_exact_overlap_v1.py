#!/usr/bin/env python3
"""Read-only exact overlap audit between the Phase 8 crosswalk and Part 1 CNV H5AD."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--crosswalk")
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    a = ad.read_h5ad(h5ad_path, backed="r")
    obs = a.obs
    print(f"H5AD={h5ad_path}")
    print(f"SHAPE={a.n_obs}x{a.n_vars}")
    print("OBS_COLUMNS=" + "|".join(map(str, obs.columns)))
    print("OBS_NAMES_HEAD=" + "|".join(map(str, obs.index[:10])))
    for col in obs.columns:
        low = str(col).lower()
        if any(token in low for token in ("dataset", "patient", "sample", "barcode", "cell_id", "malignan", "cnv", "copykat", "infercnv")):
            values = obs[col].astype(str)
            print(f"COLUMN={col}\tN_UNIQUE={values.nunique(dropna=False)}\tHEAD=" + "|".join(values.head(5)))

    if not args.crosswalk:
        return

    ext = pd.read_csv(args.crosswalk, sep="\t", dtype=str, compression="infer")
    required = {"cell_id", "dataset", "patient_id", "barcode_original", "object_name", "phase8_spatial_reference_state"}
    missing = required.difference(ext.columns)
    if missing:
        raise RuntimeError(f"Crosswalk missing columns: {sorted(missing)}")
    if ext["cell_id"].duplicated().any():
        raise RuntimeError("External crosswalk cell_id is not unique")

    part1_ids = pd.Index(obs.index.astype(str))
    exact = ext["cell_id"].isin(part1_ids)
    ext = ext.assign(exact_obs_name_match=exact)
    print(f"EXTERNAL_CELLS={len(ext)}")
    print(f"EXACT_OBS_NAME_MATCH={int(exact.sum())}")
    print("EXACT_BY_DATASET")
    print(ext.groupby("dataset", observed=True)["exact_obs_name_match"].agg(["size", "sum"]).to_csv(sep="\t"), end="")
    print("EXACT_BY_DATASET_PATIENT")
    print(ext.groupby(["dataset", "patient_id"], observed=True)["exact_obs_name_match"].agg(["size", "sum"]).to_csv(sep="\t"), end="")

    cnv_columns = [
        c for c in (
            "cnv_assessed_v1", "infercnv_call_v1", "copykat_call_v1",
            "malignancy_status_v1", "malignant_conservative_v1",
            "malignancy_evidence_v1", "eligible_for_subclone_selection_v1",
        ) if c in obs.columns
    ]
    if exact.any():
        aligned = obs.loc[ext.loc[exact, "cell_id"], cnv_columns].copy()
        print("EXACT_MATCH_CNV_COUNTS")
        for col in cnv_columns:
            print(f"FIELD={col}")
            print(aligned[col].astype(str).value_counts(dropna=False).to_csv(sep="\t"), end="")


if __name__ == "__main__":
    main()
