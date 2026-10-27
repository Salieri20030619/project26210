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
RUN = SOURCE_ROOT / "runs" / "20260812_final_cnv_gated_v6"
SOURCE_ANALYSIS = SOURCE_ROOT / "03_run_patient_level_analysis_v1.py"
SOURCE_GSE178 = SOURCE_ROOT / "patient_analysis_final_cnv_v6" / "GSE178481_full_vs_added_subset_sensitivity_v6.tsv"
SOURCE_FULL = SOURCE_ROOT / "patient_analysis_final_cnv_v6" / "within_dataset_patient_associations_v6.tsv"
FINAL = HERE / "final_cnv_incremental_v3"
PART = HERE / "final_cnv_incremental_v3.part"

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
    raise FileExistsError("Refusing existing final-CNV incremental output directory")
PART.mkdir()

spec = importlib.util.spec_from_file_location("phase14_patient", SOURCE_ANALYSIS)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

frames: list[pd.DataFrame] = []
composition_frames: list[pd.DataFrame] = []
validation_rows: list[dict[str, object]] = []

for dataset in ("GSE178481", "GSE207493"):
    table_root = RUN / dataset / "tables"
    samples = pd.read_csv(table_root / "sample_program_summary_v1.tsv", sep="\t")
    patient_program, _ = module.build_cross_sectional_patient_programs(samples)
    expected_scope = OLD[dataset] | ADDED[dataset]
    # PR8 contributes only normal tissue and is intentionally absent from the
    # cross-sectional tumour-selected patient program table.
    expected = expected_scope - ({"PR8"} if dataset == "GSE178481" else set())
    observed = set(patient_program["patient_id"].astype(str))
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    validation_rows.append(
        {
            "dataset": dataset,
            "check": "patient_scope_available",
            "status": "PASS" if not missing and not extra else "FAIL",
            "detail": f"observed={len(observed)};expected={len(expected)};missing={';'.join(missing) or 'none'};extra={';'.join(extra) or 'none'}",
        }
    )
    if missing or extra:
        raise RuntimeError(f"{dataset} final-CNV patient scope mismatch")

    subsets = [
        ("old_subset", OLD[dataset]),
        ("added_subset", ADDED[dataset]),
        ("expanded_full", expected_scope),
    ]
    for label, members in subsets:
        result = module.subset_analysis(patient_program, dataset, label, members)
        result.insert(3, "scope_members", len(members))
        frames.append(result)

    counts = pd.read_csv(table_root / "final_cnv_malignancy_counts_v4.tsv", sep="\t")
    id_column = "patient_id" if "patient_id" in counts.columns else None
    if id_column:
        counts[id_column] = counts[id_column].astype(str)
        for label, members in subsets:
            sub = counts.loc[counts[id_column].isin(members)].copy()
            numeric = [c for c in sub.columns if c != id_column and pd.api.types.is_numeric_dtype(sub[c])]
            row = {"dataset": dataset, "subset": label, "scope_members": len(members), "rows": len(sub)}
            row.update({c: sub[c].sum() for c in numeric})
            composition_frames.append(pd.DataFrame([row]))

associations = pd.concat(frames, ignore_index=True)
validation = pd.DataFrame(validation_rows)

# Reproduce the frozen final-CNV GSE178481 old/added subset table exactly.
frozen = pd.read_csv(SOURCE_GSE178, sep="\t")
mapping = {
    "previously_included_5_patients_11_libraries": "old_subset",
    "newly_added_6_patients_15_libraries": "added_subset",
}
for frozen_label, current_label in mapping.items():
    a = frozen.loc[frozen["subset"] == frozen_label].sort_values(["method", "relation"]).reset_index(drop=True)
    b = associations.loc[(associations["dataset"] == "GSE178481") & (associations["subset"] == current_label)].sort_values(["method", "relation"]).reset_index(drop=True)
    for field in ("n_patients", "spearman_rho", "p_value", "p_fdr_bh_within_dataset_method", "status"):
        if field == "status":
            av = a[field].fillna("NA").astype(str).replace({"TESTED": "TESTED_PROVISIONAL_MALIGNANCY"})
            bv = b[field].fillna("NA").astype(str).replace({"TESTED": "TESTED_PROVISIONAL_MALIGNANCY"})
            passed = av.equals(bv)
            maximum = 0.0 if passed else np.nan
        else:
            av = pd.to_numeric(a[field], errors="coerce").to_numpy(float)
            bv = pd.to_numeric(b[field], errors="coerce").to_numpy(float)
            passed = bool(np.allclose(av, bv, equal_nan=True, atol=1e-12, rtol=0))
            maximum = float(np.nanmax(np.abs(av - bv))) if np.any(np.isfinite(av - bv)) else 0.0
        validation_rows.append(
            {
                "dataset": "GSE178481",
                "check": f"reproduce_{current_label}_{field}",
                "status": "PASS" if passed else "FAIL",
                "detail": f"maximum_absolute_difference={maximum}",
            }
        )

# Reproduce the frozen full-cohort final-CNV rows for both datasets.
frozen_full = pd.read_csv(SOURCE_FULL, sep="\t")
for dataset in ("GSE178481", "GSE207493"):
    a = frozen_full.loc[frozen_full["dataset"] == dataset].sort_values(["method", "relation"]).reset_index(drop=True)
    b = associations.loc[(associations["dataset"] == dataset) & (associations["subset"] == "expanded_full")].sort_values(["method", "relation"]).reset_index(drop=True)
    for field in ("n_patients", "spearman_rho", "p_value", "p_fdr_bh_within_dataset_method", "status"):
        if field == "status":
            av = a[field].fillna("NA").astype(str).replace({"TESTED": "TESTED_PROVISIONAL_MALIGNANCY"})
            bv = b[field].fillna("NA").astype(str).replace({"TESTED": "TESTED_PROVISIONAL_MALIGNANCY"})
            passed = av.equals(bv)
            maximum = 0.0 if passed else np.nan
        else:
            av = pd.to_numeric(a[field], errors="coerce").to_numpy(float)
            bv = pd.to_numeric(b[field], errors="coerce").to_numpy(float)
            passed = bool(np.allclose(av, bv, equal_nan=True, atol=1e-12, rtol=0))
            maximum = float(np.nanmax(np.abs(av - bv))) if np.any(np.isfinite(av - bv)) else 0.0
        validation_rows.append(
            {
                "dataset": dataset,
                "check": f"reproduce_expanded_full_{field}",
                "status": "PASS" if passed else "FAIL",
                "detail": f"maximum_absolute_difference={maximum}",
            }
        )

validation = pd.DataFrame(validation_rows)
if (validation["status"] != "PASS").any():
    raise RuntimeError("Final-CNV incremental validation failed")

key_relations = ["H1_F1_to_SPP1", "H1_TYMP_to_SPP1", "H3_F1_to_hypoxia", "H3_TYMP_to_hypoxia", "H5_TYMP_to_MTA"]
key = associations.loc[associations["relation"].isin(key_relations)].copy()
key["direction"] = np.select(
    [key["spearman_rho"] > 0, key["spearman_rho"] < 0],
    ["positive", "negative"],
    default="not_estimable_or_zero",
)

atomic_tsv(associations, PART / "final_cnv_subset_associations_v1.tsv")
atomic_tsv(key, PART / "final_cnv_key_relation_comparison_v1.tsv")
atomic_tsv(validation, PART / "final_cnv_validation_v1.tsv")
if composition_frames:
    atomic_tsv(pd.concat(composition_frames, ignore_index=True, sort=False), PART / "final_cnv_count_summary_v1.tsv")

metadata = {
    "status": "PASS",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "datasets": ["GSE178481", "GSE207493"],
    "subsets": ["old_subset", "added_subset", "expanded_full"],
    "source_analysis_sha256": sha256(SOURCE_ANALYSIS),
    "source_gse178_final_sensitivity_sha256": sha256(SOURCE_GSE178),
    "source_full_final_associations_sha256": sha256(SOURCE_FULL),
    "validation_checks": len(validation),
    "validation_pass": int((validation["status"] == "PASS").sum()),
}
(PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
(PART / "FINAL_CNV_INCREMENTAL_PASS_v1").write_text(
    "status=PASS\n"
    "datasets=GSE178481|GSE207493\n"
    "subsets=old_subset|added_subset|expanded_full\n"
    f"association_rows={len(associations)}\n"
    f"validation_checks={len(validation)}\n"
    "normalization_performed=0\n"
    "integration_performed=0\n"
    "clustering_performed=0\n",
    encoding="utf-8",
)
os.replace(PART, FINAL)
print(f"PASS associations={len(associations)} validation={len(validation)}")
