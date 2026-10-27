#!/usr/bin/env python3
import argparse
import re

import anndata as ad
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", required=True)
parser.add_argument("--crosswalk", required=True)
args = parser.parse_args()

barcode_re = re.compile(r"([ACGTN]+-[0-9]+)$", re.IGNORECASE)
core = lambda value: (barcode_re.search(str(value)).group(1).upper() if barcode_re.search(str(value)) else None)

a = ad.read_h5ad(args.h5ad, backed="r")
obs = a.obs.copy()
obs["obs_name"] = obs.index.astype(str)
obs["core_barcode"] = obs["obs_name"].map(core)
ext = pd.read_csv(
    args.crosswalk,
    sep="\t",
    compression="infer",
    dtype=str,
    usecols=["cell_id", "dataset", "patient_id", "src__gsm", "barcode_original"],
)
ext["core_barcode"] = ext["barcode_original"].map(core)

for dataset in ("GSE178481", "GSE207493"):
    left = ext.loc[ext["dataset"].eq(dataset)].copy()
    right = obs.loc[obs["study"].astype(str).eq(dataset)].copy()
    left_key = left["src__gsm"].astype(str) + "|" + left["core_barcode"].astype(str)
    right_key = right["geo_accession"].astype(str) + "|" + right["core_barcode"].astype(str)
    print(f"DATASET={dataset}")
    print(f"EXT={len(left)} PART1={len(right)}")
    print(f"EXT_KEY_DUP={int(left_key.duplicated(keep=False).sum())} PART1_KEY_DUP={int(right_key.duplicated(keep=False).sum())}")
    print(f"KEY_INTERSECTION={int(left_key.isin(set(right_key)).sum())}")
    print("EXT_HEAD")
    print(left[["cell_id", "patient_id", "src__gsm", "barcode_original", "core_barcode"]].head(5).to_csv(sep="\t", index=False), end="")
    print("PART1_HEAD")
    print(right[["obs_name", "curated_patient_id", "geo_accession", "core_barcode"]].head(5).to_csv(sep="\t", index=False), end="")
