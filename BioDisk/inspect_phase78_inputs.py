#!/usr/bin/env python3
import gzip
from collections import Counter
from pathlib import Path

import h5py

REF = Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad")
META = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2")
RAW = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_phase60_v1/01_qc/cohorts")

with h5py.File(REF, "r") as handle:
    print("REFERENCE_OBS_KEYS")
    print("\n".join(sorted(handle["obs"].keys())))
    print("REFERENCE_VAR_KEYS")
    print("\n".join(sorted(handle["var"].keys())))

for path in sorted(META.glob("GSE*_metadata_v1.tsv.gz")):
    counts = Counter()
    patients = Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").replace('"', '').split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in handle:
            row = line.rstrip("\n").replace('"', '').split("\t")
            if row[ix["phase7_broad_annotation"]] == "Myeloid":
                fine = row[ix["phase7_fine_annotation"]]
                counts[fine] += 1
                patients[(fine, row[ix["patient_key"]])] += 1
    print(f"METADATA {path.name}")
    for label, n in counts.most_common():
        n_patients = sum(1 for (lab, _), m in patients.items() if lab == label and m > 0)
        print(f"{label}\t{n}\tpatients={n_patients}")

for name in ["GSE264586.qc_filtered_raw_v1.h5ad", "GSE294109.qc_filtered_raw_v1.h5ad"]:
    path = RAW / name
    with h5py.File(path, "r") as handle:
        print(f"RAW {name} X_type={type(handle['X']).__name__}")
        print("obs=" + ",".join(sorted(handle["obs"].keys())))
        print("var=" + ",".join(sorted(handle["var"].keys())))
        print("layers=" + ",".join(sorted(handle.get("layers", {}).keys())))
