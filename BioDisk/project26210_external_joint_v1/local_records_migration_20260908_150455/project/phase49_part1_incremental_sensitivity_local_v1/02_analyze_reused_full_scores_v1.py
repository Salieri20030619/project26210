#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE_ROOT = ROOT / "phase14_patient_program_validation_local_v1"
RUN = SOURCE_ROOT / "runs" / "20260812_082000_provisional_v3"
SOURCE_ANALYSIS = SOURCE_ROOT / "03_run_patient_level_analysis_v1.py"
SOURCE_GSE178_SENSITIVITY = SOURCE_ROOT / "patient_analysis_v3" / "GSE178481_full_vs_added_subset_sensitivity_v3.tsv"
FINAL = HERE / "reused_full_scores_v2"
PART = HERE / "reused_full_scores_v2.part"

OLD = {
    "GSE178481": {"PR5", "PR6", "PR7", "PR8", "PR9"},
    "GSE207493": {"RCC86", "RCC87", "RCC94", "RCC96", "RCC99", "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113"},
}
ADDED = {
    "GSE178481": {"BM1", "BM2", "PR1", "PR2", "PR3", "PR4"},
    "GSE207493": {"RCC81", "RCC84", "RCC114", "RCC115", "RCC116", "RCC119", "RCC120"},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


if FINAL.exists() or PART.exists():
    raise FileExistsError("Refusing existing v1 output directory")
PART.mkdir()

spec = importlib.util.spec_from_file_location("phase14_patient", SOURCE_ANALYSIS)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

association_frames = []
composition_frames = []
reproduction_rows = []
for dataset in ("GSE178481", "GSE207493"):
    table_root = RUN / dataset / "tables"
    samples = pd.read_csv(table_root / "sample_program_summary_v1.tsv", sep="\t")
    patient_program, selection_audit = module.build_cross_sectional_patient_programs(samples)
    expected_all = OLD[dataset] | ADDED[dataset]
    observed_all = set(patient_program["patient_id"].astype(str))
    expected_eligible = expected_all - ({"PR8"} if dataset == "GSE178481" else set())
    if observed_all != expected_eligible:
        raise RuntimeError(f"{dataset} eligible patient set mismatch: observed={sorted(observed_all)} expected={sorted(expected_eligible)}")

    subsets = [
        ("old_subset", OLD[dataset]),
        ("added_subset", ADDED[dataset]),
        ("expanded_full", expected_all),
    ]
    for label, members in subsets:
        result = module.subset_analysis(patient_program, dataset, label, members)
        result.insert(3, "scope_members", len(members))
        association_frames.append(result)

    counts = pd.read_csv(table_root / "patient_compartment_counts_v1.tsv", sep="\t")
    counts["patient_id"] = counts["patient_id"].astype(str)
    for label, members in subsets:
        sub = counts.loc[counts["patient_id"].isin(members)].copy()
        grouped = sub.groupby("analysis_compartment", as_index=False)["n_cells"].sum()
        grouped.insert(0, "dataset", dataset)
        grouped.insert(1, "subset", label)
        grouped.insert(2, "subset_members", len(members))
        grouped["fraction_subset_cells"] = grouped["n_cells"] / grouped["n_cells"].sum()
        composition_frames.append(grouped)

    if dataset == "GSE178481":
        previous = pd.read_csv(SOURCE_GSE178_SENSITIVITY, sep="\t")
        mapping = {
            "previously_included_5_patients_11_libraries": "old_subset",
            "newly_added_6_patients_15_libraries": "added_subset",
        }
        current = pd.concat(association_frames, ignore_index=True)
        current = current.loc[(current["dataset"] == dataset) & current["subset"].isin(mapping.values())]
        for old_label, new_label in mapping.items():
            a = previous.loc[previous["subset"] == old_label].sort_values(["method", "relation"]).reset_index(drop=True)
            b = current.loc[current["subset"] == new_label].sort_values(["method", "relation"]).reset_index(drop=True)
            for field in ("n_patients", "spearman_rho", "p_value", "p_fdr_bh_within_dataset_method", "status"):
                if field in ("status",):
                    equal = a[field].fillna("NA").astype(str).equals(b[field].fillna("NA").astype(str))
                    maximum = 0.0 if equal else np.nan
                else:
                    av = pd.to_numeric(a[field], errors="coerce").to_numpy(float)
                    bv = pd.to_numeric(b[field], errors="coerce").to_numpy(float)
                    equal = bool(np.allclose(av, bv, equal_nan=True, atol=1e-12, rtol=0))
                    maximum = float(np.nanmax(np.abs(av - bv))) if np.any(np.isfinite(av - bv)) else 0.0
                reproduction_rows.append({"subset": new_label, "field": field, "status": "PASS" if equal else "FAIL", "maximum_absolute_difference": maximum})

associations = pd.concat(association_frames, ignore_index=True)
composition = pd.concat(composition_frames, ignore_index=True)
reproduction = pd.DataFrame(reproduction_rows)
if (reproduction["status"] != "PASS").any():
    raise RuntimeError("GSE178481 frozen sensitivity reproduction failed")

key = associations.loc[associations["relation"].isin(["H1_TYMP_to_SPP1", "H3_TYMP_to_hypoxia", "H5_TYMP_to_MTA", "H1_F1_to_SPP1", "H3_F1_to_hypoxia"])].copy()
key["direction"] = np.select([key["spearman_rho"] > 0, key["spearman_rho"] < 0], ["positive", "negative"], default="not_estimable_or_zero")

atomic_tsv(associations, PART / "incremental_subset_associations_v1.tsv")
atomic_tsv(key, PART / "key_relation_subset_comparison_v1.tsv")
atomic_tsv(composition, PART / "incremental_subset_composition_v1.tsv")
atomic_tsv(reproduction, PART / "gse178481_reproduction_validation_v1.tsv")

metadata = {
    "status": "PASS",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "datasets": ["GSE178481", "GSE207493"],
    "subsets": ["old_subset", "added_subset", "expanded_full"],
    "relations_per_dataset_method_subset": 12,
    "bootstrap_replicates": module.N_BOOT,
    "minimum_cells_per_patient_compartment": module.MIN_CELLS,
    "minimum_patients": module.MIN_PATIENTS,
    "source_analysis_sha256": sha256(SOURCE_ANALYSIS),
    "source_gse178_sensitivity_sha256": sha256(SOURCE_GSE178_SENSITIVITY),
    "gse178_reproduction_checks": int(len(reproduction)),
    "gse178_reproduction_pass": int((reproduction["status"] == "PASS").sum()),
}
(PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
(PART / "INCREMENTAL_REUSED_SCORES_PASS_v1").write_text(
    "status=PASS\n"
    "datasets=GSE178481|GSE207493\n"
    "subsets=old_subset|added_subset|expanded_full\n"
    f"association_rows={len(associations)}\n"
    f"composition_rows={len(composition)}\n"
    f"gse178_reproduction_checks={len(reproduction)}\n"
    "normalization_performed=0\n"
    "integration_performed=0\n"
    "clustering_performed=0\n",
    encoding="utf-8",
)
os.replace(PART, FINAL)
print(f"PASS associations={len(associations)} composition={len(composition)}")
