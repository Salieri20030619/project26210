#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCOPE = HERE / "scope_v1" / "incremental_library_scope_v1.tsv"
PROGRAM_DIR = HERE / "reused_full_scores_v2"
CNV_DIR = HERE / "final_cnv_incremental_v3"
GSE242_DIR = HERE / "server_results_v1" / "gse242299_added_normals_v1"
FINAL = HERE / "incremental_evidence_summary_v2"
PART = HERE / "incremental_evidence_summary_v2.part"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    temp = path.with_name(path.name + ".part")
    frame.to_csv(temp, sep="\t", index=False, na_rep="NA", encoding="utf-8", lineterminator="\n")
    os.replace(temp, path)


if FINAL.exists() or PART.exists():
    raise FileExistsError("Refusing existing integrated Phase 49 output directory")
PART.mkdir()

scope = pd.read_csv(SCOPE, sep="\t")
program = pd.read_csv(PROGRAM_DIR / "incremental_subset_associations_v1.tsv", sep="\t")
cnv = pd.read_csv(CNV_DIR / "final_cnv_subset_associations_v1.tsv", sep="\t")
composition = pd.read_csv(GSE242_DIR / "added_library_celltype_composition_v1.tsv", sep="\t")
library = pd.read_csv(GSE242_DIR / "library_tymp_f1_summary_v1.tsv", sep="\t")

checks: list[dict[str, object]] = []


def check(name: str, passed: bool, observed: object, expected: object) -> None:
    checks.append(
        {
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
        }
    )


check("scope_rows", len(scope) == 99, len(scope), 99)
check("scope_added_rows", int((scope["comparison_group"] == "added24").sum()) == 24, int((scope["comparison_group"] == "added24").sum()), 24)
check("program_association_rows", len(program) == 144, len(program), 144)
check("final_cnv_association_rows", len(cnv) == 144, len(cnv), 144)
check("program_pass_marker", (PROGRAM_DIR / "INCREMENTAL_REUSED_SCORES_PASS_v1").exists(), "present", "present")
check("final_cnv_pass_marker", (CNV_DIR / "FINAL_CNV_INCREMENTAL_PASS_v1").exists(), "present", "present")
check("gse242_pass_marker", (GSE242_DIR / "GSE242299_INCREMENTAL_PASS_v1").exists(), "present", "present")
check("gse242_library_count", library["library_id"].nunique() == 31, library["library_id"].nunique(), 31)
check("gse242_added_library_ids", set(library.loc[library["comparison_group"] == "added2", "library_id"]) == {"N3_1", "N4_1"}, ";".join(sorted(library.loc[library["comparison_group"] == "added2", "library_id"].unique())), "N3_1;N4_1")
check("gse242_f1_gene_coverage", int(library.shape[0]) == 31, int(library.shape[0]), 31)

# Tissue-matched GSE242299 comparison: old healthy libraries versus the two added healthy libraries.
healthy_library = library.loc[library["tissue"] == "Healthy"].copy()
score_columns = ["tymp_detection_fraction", "tymp_mean_log1p_cp10k", "f1_primary_mean", "f1_sensitivity_mean"]
healthy_rows: list[dict[str, object]] = []
for group, frame in healthy_library.groupby("comparison_group", sort=False):
    cells = int(frame["cells"].sum())
    row: dict[str, object] = {
        "comparison_group": group,
        "libraries": int(frame["library_id"].nunique()),
        "cells": cells,
        "library_ids": ";".join(sorted(frame["library_id"].astype(str))),
    }
    for column in score_columns:
        row[column + "_cell_weighted"] = float(np.average(frame[column], weights=frame["cells"]))
        row[column + "_library_median"] = float(frame[column].median())
    healthy_rows.append(row)
healthy_summary = pd.DataFrame(healthy_rows).sort_values("comparison_group").reset_index(drop=True)
check("gse242_old_healthy_libraries", int(healthy_summary.loc[healthy_summary["comparison_group"] == "old29", "libraries"].iloc[0]) == 10, int(healthy_summary.loc[healthy_summary["comparison_group"] == "old29", "libraries"].iloc[0]), 10)
check("gse242_added_healthy_libraries", int(healthy_summary.loc[healthy_summary["comparison_group"] == "added2", "libraries"].iloc[0]) == 2, int(healthy_summary.loc[healthy_summary["comparison_group"] == "added2", "libraries"].iloc[0]), 2)

healthy_comp = composition.loc[composition["tissue"] == "Healthy"].groupby(["comparison_group", "broad_cell_type"], as_index=False)["cells"].sum()
healthy_comp["fraction_within_group"] = healthy_comp["cells"] / healthy_comp.groupby("comparison_group")["cells"].transform("sum")
healthy_comp = healthy_comp.sort_values(["comparison_group", "broad_cell_type"]).reset_index(drop=True)

# Long-form evidence table keeps program-direction and strict final-CNV tiers separate.
key_relations = ["H1_F1_to_SPP1", "H1_TYMP_to_SPP1", "H3_F1_to_hypoxia", "H3_TYMP_to_hypoxia", "H5_TYMP_to_MTA"]
columns = [
    "dataset", "subset", "method", "relation", "n_patients", "spearman_rho", "p_value",
    "p_fdr_bh_within_dataset_method", "status",
]
program_key = program.loc[program["relation"].isin(key_relations), columns].copy()
program_key.insert(0, "evidence_tier", "program_direction")
cnv_key = cnv.loc[cnv["relation"].isin(key_relations), columns].copy()
cnv_key.insert(0, "evidence_tier", "final_cnv_gate")
key_table = pd.concat([program_key, cnv_key], ignore_index=True)
key_table["direction"] = np.select(
    [key_table["spearman_rho"] > 0, key_table["spearman_rho"] < 0],
    ["positive", "negative"],
    default="not_estimable_or_zero",
)

# Explicit decision rows used in the status report.
decision = pd.DataFrame(
    [
        {
            "dataset": "GSE178481",
            "increment": "15 libraries / 6 scoped patients",
            "local_result": "Expanded program-direction TYMP-SPP1: n=8, rho=0.928571, FDR=0.010356 in both score methods",
            "strict_result": "Final-CNV malignant eligibility: old n=1, added n=1, expanded n=2; core malignant relations not estimable",
            "decision": "directional_support_strengthened_but_not_strict_malignant_validation",
        },
        {
            "dataset": "GSE207493",
            "increment": "7 libraries / 7 patients",
            "local_result": "Expanded program-direction F1-SPP1 and TYMP-SPP1 remain weakly negative; no FDR-significant core relation",
            "strict_result": "Expanded final-CNV: F1-SPP1 rho=-0.554545 and TYMP-SPP1 rho=-0.409091; neither significant",
            "decision": "core_conclusion_unchanged; added_subset_has_three_method_specific_myeloid_relations",
        },
        {
            "dataset": "GSE242299",
            "increment": "2 healthy libraries / 1663 cells",
            "local_result": "Added healthy TYMP detection=0.032471 versus 0.048513 in old healthy cells; F1 primary mean=0.106740 versus 0.108030",
            "strict_result": "Added libraries are healthy; only 5/1663 cells carry the source Tumor label and no malignant test was run",
            "decision": "normal_reference_extended_without_changing_tumour_axis_conclusion",
        },
        {
            "dataset": "ALL_AFFECTED_DATASETS",
            "increment": "24 libraries / 13 patients",
            "local_result": "Existing full-cohort scores plus targeted H5AD comparison cover every omitted library",
            "strict_result": "No new strict malignant-axis positive result appears after the final-CNV gate",
            "decision": "full_reanalysis_not_required_for_current_scope_sensitivity_question",
        },
    ]
)

validation = pd.DataFrame(checks)
if (validation["status"] != "PASS").any():
    atomic_tsv(validation, PART / "incremental_validation_checks_v1.tsv")
    raise RuntimeError("Integrated Phase 49 validation failed")

atomic_tsv(healthy_summary, PART / "gse242299_healthy_matched_summary_v1.tsv")
atomic_tsv(healthy_comp, PART / "gse242299_healthy_composition_v1.tsv")
atomic_tsv(key_table, PART / "key_relation_evidence_tiers_v1.tsv")
atomic_tsv(decision, PART / "incremental_decision_summary_v1.tsv")
atomic_tsv(validation, PART / "incremental_validation_checks_v1.tsv")

report = """# Phase 49 omitted-library incremental sensitivity status

## Status

**PASS — LOCAL_COMPARISON_COMPLETE; FULL_REANALYSIS_NOT_REQUIRED_FOR_CURRENT_QUESTION**

The previous 75-library object and the complete 99-library scope differ by 24 libraries and 13 patients. The omitted material is confined to GSE178481 (15 libraries), GSE207493 (7 libraries) and GSE242299 (2 healthy libraries). All omitted libraries are covered here using frozen full-cohort program results, frozen final-CNV-gated results, or a targeted raw-count calculation for the two GSE242299 healthy libraries. No five-dataset normalization, integration, clustering or annotation was performed.

## Findings

### GSE178481

At the program-direction level, the six newly scoped patients strengthen the TYMP–SPP1 relationship. In the expanded eligible cohort, TYMP–SPP1 is positive in both scoring methods (n=8, Spearman rho=0.929, p=0.000863, within-method FDR=0.0104). The added subset also shows a positive F1–SPP1 relationship (n=5, rho=1.000, FDR<0.001).

The final CNV gate leaves one eligible old patient and one eligible added patient, for two patients in the expanded cohort. Malignant-cell F1/TYMP relations therefore remain non-estimable under the strict gate. The GSE178481 increment strengthens directional support but does not create strict malignant validation.

### GSE207493

The added seven patients do not reverse the expanded-cohort core result. At the program-direction level, expanded F1–SPP1 and TYMP–SPP1 remain weakly negative. Under the final CNV gate, expanded F1–SPP1 is rho=-0.555 and TYMP–SPP1 is rho=-0.409; neither is significant.

The added subset has three FDR-significant relations in the sensitivity scoring method (n=7, |rho|=0.821, FDR=0.0234): SPP1–myeloid vascular is positive, SPP1–T-cell dysfunction is negative and SPP1–cytotoxicity is positive. These are method-specific secondary immune findings and do not change the malignant F1/TYMP conclusion.

### GSE242299

The two omitted libraries are N3_1 and N4_1, both healthy tissue, with 1,663 cells in total. Compared with the ten previously included healthy libraries, cell-weighted TYMP detection is lower (3.25% versus 4.85%), TYMP mean log1p(CP10K) is lower (0.0828 versus 0.1047), and the primary F1 mean is nearly unchanged (0.1067 versus 0.1080). Broad cell-type composition remains comparable. Five of 1,663 added cells carry the source dataset's Tumor label; this 0.30% fraction is comparable to 0.27% in the previously included healthy cells and is not used as malignant evidence.

## Decision

The local comparison answers the cohort-omission question without a complete five-dataset rebuild. The omitted libraries strengthen one GSE178481 directional result, add three method-specific GSE207493 immune associations, and extend the GSE242299 normal reference. They do not produce a new strict malignant-axis positive result and do not change the overall interpretation of the current Part 1 analysis.

A complete rebuild is only needed if the revised client plan requires a new unified cell atlas, new clustering/annotation, or re-estimation of every CNV and cNMF result from raw counts. It is not required for the present omitted-library sensitivity check.

## Reproducibility

- Scope validation: 99 total libraries, 24 added libraries, 13 added patients.
- Reused program-result reproduction: 10/10 checks passed.
- Final-CNV subset reproduction and scope validation: 22/22 checks passed.
- GSE242299 targeted H5AD analysis: 31/31 libraries recognized; N3_1 and N4_1 analyzed from validated integer raw counts.
- Integrated validation: {validation_pass}/{validation_total} checks passed.
""".format(validation_pass=int((validation["status"] == "PASS").sum()), validation_total=len(validation))
(PART / "PHASE49_INCREMENTAL_SENSITIVITY_STATUS_v1.md").write_text(report, encoding="utf-8")
(PART / "PHASE49_INCREMENTAL_SENSITIVITY_PASS_v1").write_text(
    "status=PASS\n"
    "decision=FULL_REANALYSIS_NOT_REQUIRED_FOR_CURRENT_SCOPE_SENSITIVITY_QUESTION\n"
    "old_libraries=75\nexpanded_libraries=99\nadded_libraries=24\nadded_patients=13\n"
    f"validation_checks={len(validation)}\n"
    "normalization_scope=GSE242299_TARGETED_SCORE_ONLY\n"
    "integration_performed=0\nclustering_performed=0\nannotation_performed=0\n",
    encoding="utf-8",
)

deliverables = sorted(p for p in PART.iterdir() if p.is_file())
manifest_rows = [
    {"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
    for p in deliverables
]
atomic_tsv(pd.DataFrame(manifest_rows), PART / "FILE_MANIFEST_SHA256_v1.tsv")
metadata = {
    "status": "PASS",
    "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    "decision": "FULL_REANALYSIS_NOT_REQUIRED_FOR_CURRENT_SCOPE_SENSITIVITY_QUESTION",
    "old_libraries": 75,
    "expanded_libraries": 99,
    "added_libraries": 24,
    "added_patients": 13,
    "input_sha256": {
        "scope": sha256(SCOPE),
        "program_associations": sha256(PROGRAM_DIR / "incremental_subset_associations_v1.tsv"),
        "final_cnv_associations": sha256(CNV_DIR / "final_cnv_subset_associations_v1.tsv"),
        "gse242_library_summary": sha256(GSE242_DIR / "library_tymp_f1_summary_v1.tsv"),
        "gse242_composition": sha256(GSE242_DIR / "added_library_celltype_composition_v1.tsv"),
    },
    "validation_checks": len(validation),
    "validation_pass": int((validation["status"] == "PASS").sum()),
}
(PART / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
os.replace(PART, FINAL)
print(f"PASS validation={len(validation)} decision={metadata['decision']}")
