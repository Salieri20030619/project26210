#!/usr/bin/env python3
"""Complete Phase 49 with the frozen Phase 14 expanded-malignancy layer."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PHASE14 = ROOT / "phase14_patient_program_validation_local_v1"
PHASE49 = ROOT / "phase49_part1_incremental_sensitivity_local_v1"
EXPANDED_RUN = PHASE14 / "runs" / "20260812_expanded_malignancy_v2"
SOURCE_ANALYSIS = PHASE14 / "03_run_patient_level_analysis_v1.py"
SOURCE_GSE178_EXPANDED = (
    PHASE14
    / "patient_analysis_expanded_malignancy_v1"
    / "GSE178481_full_vs_added_subset_sensitivity_expanded_v1.tsv"
)
SOURCE_FULL_EXPANDED = (
    PHASE14
    / "patient_analysis_expanded_malignancy_v1"
    / "within_dataset_patient_associations_expanded_v1.tsv"
)
SOURCE_PROVISIONAL = PHASE49 / "reused_full_scores_v2" / "incremental_subset_associations_v1.tsv"
SOURCE_STRICT = PHASE49 / "final_cnv_incremental_v3" / "final_cnv_subset_associations_v1.tsv"
FINAL = HERE / "analysis_local_v1"
PART = HERE / "analysis_local_v1.part"

OLD = {
    "GSE178481": {"PR5", "PR6", "PR7", "PR8", "PR9"},
    "GSE207493": {
        "RCC86", "RCC87", "RCC94", "RCC96", "RCC99", "RCC100",
        "RCC101", "RCC103", "RCC104", "RCC106", "RCC112", "RCC113",
    },
}
ADDED = {
    "GSE178481": {"BM1", "BM2", "PR1", "PR2", "PR3", "PR4"},
    "GSE207493": {"RCC81", "RCC84", "RCC114", "RCC115", "RCC116", "RCC119", "RCC120"},
}
LAYERS = {
    "provisional_tumor_like": "tumor-like epithelial annotation; no positive CNV call required",
    "expanded_any_positive_cnv": "tumor-like epithelial annotation plus inferCNV-high or CopyKAT-aneuploid",
    "strict_conservative_cnv": "tumor-like epithelial annotation plus frozen conservative CNV support",
}
CORE_RELATIONS = ["H1_F1_to_SPP1", "H1_TYMP_to_SPP1"]
KEY_RELATIONS = [
    "H1_F1_to_SPP1", "H1_TYMP_to_SPP1", "H3_F1_to_hypoxia",
    "H3_TYMP_to_hypoxia", "H5_TYMP_to_MTA",
]


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


def normalize_status(values: pd.Series) -> pd.Series:
    return (
        values.fillna("NA")
        .astype(str)
        .replace(
            {
                "TESTED_PROVISIONAL_MALIGNANCY": "TESTED_POSTHOC_EXPANDED_MALIGNANCY",
                "TESTED": "TESTED_POSTHOC_EXPANDED_MALIGNANCY",
            }
        )
    )


def compare_frames(
    frozen: pd.DataFrame,
    current: pd.DataFrame,
    dataset: str,
    label: str,
    rows: list[dict[str, object]],
) -> None:
    frozen = frozen.sort_values(["method", "relation"]).reset_index(drop=True)
    current = current.sort_values(["method", "relation"]).reset_index(drop=True)
    same_keys = frozen[["method", "relation"]].equals(current[["method", "relation"]])
    rows.append(
        {
            "check": f"{dataset}_{label}_row_keys",
            "status": "PASS" if same_keys else "FAIL",
            "observed": len(current),
            "expected": len(frozen),
        }
    )
    if not same_keys:
        return
    for field in ("n_patients", "spearman_rho", "p_value", "p_fdr_bh_within_dataset_method", "status"):
        if field == "status":
            passed = normalize_status(frozen[field]).equals(normalize_status(current[field]))
            detail = "exact_after_expanded_status_normalization" if passed else "status_mismatch"
        else:
            old = pd.to_numeric(frozen[field], errors="coerce").to_numpy(float)
            new = pd.to_numeric(current[field], errors="coerce").to_numpy(float)
            passed = bool(np.allclose(old, new, equal_nan=True, atol=1e-12, rtol=0))
            delta = np.abs(old - new)
            maximum = float(np.nanmax(delta)) if np.any(np.isfinite(delta)) else 0.0
            detail = f"maximum_absolute_difference={maximum}"
        rows.append(
            {
                "check": f"{dataset}_{label}_reproduce_{field}",
                "status": "PASS" if passed else "FAIL",
                "observed": detail,
                "expected": "frozen_phase14_expanded_result",
            }
        )


def reporting_tier(row: pd.Series) -> str:
    n = int(row["n_patients"])
    p = pd.to_numeric(pd.Series([row["p_value"]]), errors="coerce").iloc[0]
    q = pd.to_numeric(pd.Series([row["p_fdr_bh_within_dataset_method"]]), errors="coerce").iloc[0]
    if n < 5:
        return "direction_only_lt5"
    if n < 10:
        return "small_n_significant_reportable" if np.isfinite(p) and p < 0.05 else "direction_only_5_to_9"
    return "fdr_significant_ge10" if np.isfinite(q) and q < 0.05 else "tested_not_fdr_significant_ge10"


def get_row(frame: pd.DataFrame, dataset: str, layer: str, method: str, relation: str) -> pd.Series:
    selected = frame.loc[
        (frame["dataset"] == dataset)
        & (frame["subset"] == "expanded_full")
        & (frame["evidence_layer"] == layer)
        & (frame["method"] == method)
        & (frame["relation"] == relation)
    ]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one summary row, found {len(selected)}")
    return selected.iloc[0]


def fmt(value: object, digits: int = 3) -> str:
    number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return "NA" if not np.isfinite(number) else f"{number:.{digits}f}"


def main() -> int:
    if FINAL.exists() or PART.exists():
        raise FileExistsError(f"Refusing to overwrite {FINAL} or {PART}")
    required = [
        SOURCE_ANALYSIS, SOURCE_GSE178_EXPANDED, SOURCE_FULL_EXPANDED,
        SOURCE_PROVISIONAL, SOURCE_STRICT,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen input(s): " + "; ".join(missing))
    PART.mkdir(parents=True)

    spec = importlib.util.spec_from_file_location("phase14_patient", SOURCE_ANALYSIS)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    expanded_frames: list[pd.DataFrame] = []
    validation: list[dict[str, object]] = []
    for dataset in ("GSE178481", "GSE207493"):
        source = EXPANDED_RUN / dataset / "tables" / "sample_program_summary_v1.tsv"
        samples = pd.read_csv(source, sep="\t")
        patient_program, _ = module.build_cross_sectional_patient_programs(samples)
        expected_scope = OLD[dataset] | ADDED[dataset]
        expected_observed = expected_scope - ({"PR8"} if dataset == "GSE178481" else set())
        observed = set(patient_program["patient_id"].astype(str))
        validation.append(
            {
                "check": f"{dataset}_expanded_patient_scope",
                "status": "PASS" if observed == expected_observed else "FAIL",
                "observed": ";".join(sorted(observed)),
                "expected": ";".join(sorted(expected_observed)),
            }
        )
        for subset, members in (
            ("old_subset", OLD[dataset]),
            ("added_subset", ADDED[dataset]),
            ("expanded_full", expected_scope),
        ):
            result = module.subset_analysis(patient_program, dataset, subset, members)
            result.insert(3, "scope_members", len(members))
            result["status"] = result["status"].replace(
                {"TESTED_PROVISIONAL_MALIGNANCY": "TESTED_POSTHOC_EXPANDED_MALIGNANCY"}
            )
            expanded_frames.append(result)

    expanded = pd.concat(expanded_frames, ignore_index=True)
    frozen_gse178 = pd.read_csv(SOURCE_GSE178_EXPANDED, sep="\t")
    for frozen_label, current_label in {
        "previously_included_5_patients_11_libraries": "old_subset",
        "newly_added_6_patients_15_libraries": "added_subset",
    }.items():
        compare_frames(
            frozen_gse178.loc[frozen_gse178["subset"] == frozen_label],
            expanded.loc[(expanded["dataset"] == "GSE178481") & (expanded["subset"] == current_label)],
            "GSE178481",
            current_label,
            validation,
        )

    frozen_full = pd.read_csv(SOURCE_FULL_EXPANDED, sep="\t")
    for dataset in ("GSE178481", "GSE207493"):
        compare_frames(
            frozen_full.loc[frozen_full["dataset"] == dataset],
            expanded.loc[(expanded["dataset"] == dataset) & (expanded["subset"] == "expanded_full")],
            dataset,
            "expanded_full",
            validation,
        )

    expanded.insert(4, "evidence_layer", "expanded_any_positive_cnv")
    expanded.insert(5, "malignancy_definition", LAYERS["expanded_any_positive_cnv"])
    provisional = pd.read_csv(SOURCE_PROVISIONAL, sep="\t")
    provisional.insert(4, "evidence_layer", "provisional_tumor_like")
    provisional.insert(5, "malignancy_definition", LAYERS["provisional_tumor_like"])
    strict = pd.read_csv(SOURCE_STRICT, sep="\t")
    strict.insert(4, "evidence_layer", "strict_conservative_cnv")
    strict.insert(5, "malignancy_definition", LAYERS["strict_conservative_cnv"])

    expected_rows_per_layer = 2 * 3 * 2 * len(module.RELATIONS)
    for label, frame in (
        ("provisional_tumor_like", provisional),
        ("expanded_any_positive_cnv", expanded),
        ("strict_conservative_cnv", strict),
    ):
        keys = ["dataset", "subset", "method", "relation"]
        passed = len(frame) == expected_rows_per_layer and not frame.duplicated(keys).any()
        validation.append(
            {
                "check": f"{label}_complete_unique_grid",
                "status": "PASS" if passed else "FAIL",
                "observed": f"rows={len(frame)};duplicates={int(frame.duplicated(keys).sum())}",
                "expected": f"rows={expected_rows_per_layer};duplicates=0",
            }
        )

    combined = pd.concat([strict, expanded, provisional], ignore_index=True, sort=False)
    order = {"strict_conservative_cnv": 1, "expanded_any_positive_cnv": 2, "provisional_tumor_like": 3}
    combined["layer_order"] = combined["evidence_layer"].map(order)
    combined = combined.sort_values(
        ["dataset", "subset", "method", "relation", "layer_order"]
    ).reset_index(drop=True)
    combined["direction"] = np.select(
        [combined["spearman_rho"] > 0, combined["spearman_rho"] < 0],
        ["positive", "negative"],
        default="not_estimable_or_zero",
    )
    combined["reporting_tier"] = combined.apply(reporting_tier, axis=1)

    patient_pivot = combined.pivot(
        index=["dataset", "subset", "method", "relation"],
        columns="evidence_layer",
        values="n_patients",
    ).reset_index()
    nesting = (
        patient_pivot["strict_conservative_cnv"]
        <= patient_pivot["expanded_any_positive_cnv"]
    ) & (
        patient_pivot["expanded_any_positive_cnv"]
        <= patient_pivot["provisional_tumor_like"]
    )
    validation.append(
        {
            "check": "patient_eligibility_layer_nesting",
            "status": "PASS" if bool(nesting.all()) else "FAIL",
            "observed": f"passing_rows={int(nesting.sum())};rows={len(nesting)}",
            "expected": f"passing_rows={len(nesting)};strict<=expanded<=provisional",
        }
    )

    validation_frame = pd.DataFrame(validation)
    if (validation_frame["status"] != "PASS").any():
        atomic_tsv(validation_frame, PART / "validation_checks_v1.tsv")
        raise RuntimeError("Phase52 validation failed")

    key = combined.loc[combined["relation"].isin(KEY_RELATIONS)].copy()
    core = combined.loc[combined["relation"].isin(CORE_RELATIONS)].copy()
    atomic_tsv(expanded, PART / "expanded_malignancy_subset_associations_v1.tsv")
    atomic_tsv(combined, PART / "three_tier_incremental_associations_v1.tsv")
    atomic_tsv(key, PART / "three_tier_key_relations_v1.tsv")
    atomic_tsv(core, PART / "core_spp1_three_tier_v1.tsv")
    atomic_tsv(patient_pivot, PART / "patient_eligibility_by_layer_v1.tsv")
    atomic_tsv(validation_frame, PART / "validation_checks_v1.tsv")

    g178_strict = get_row(combined, "GSE178481", "strict_conservative_cnv", "primary_gene_z_mean", "H1_TYMP_to_SPP1")
    g178_expanded = get_row(combined, "GSE178481", "expanded_any_positive_cnv", "primary_gene_z_mean", "H1_TYMP_to_SPP1")
    g178_provisional = get_row(combined, "GSE178481", "provisional_tumor_like", "primary_gene_z_mean", "H1_TYMP_to_SPP1")
    g207_strict = get_row(combined, "GSE207493", "strict_conservative_cnv", "primary_gene_z_mean", "H1_TYMP_to_SPP1")
    g207_expanded = get_row(combined, "GSE207493", "expanded_any_positive_cnv", "primary_gene_z_mean", "H1_TYMP_to_SPP1")
    g207_provisional = get_row(combined, "GSE207493", "provisional_tumor_like", "primary_gene_z_mean", "H1_TYMP_to_SPP1")

    status_text = f"""# Phase 52 Part 1 incremental three-tier completion

## Status

**PASS — the previously frozen expanded-malignancy layer is now included in the exact old/added/full incremental comparison.**

No malignancy threshold was changed. The analysis reuses the Phase 14 definition: tumor-like epithelial annotation plus at least one positive CNV method, with Part 1 rebuild calls taking priority. No normalization, integration, clustering or annotation was run.

## GSE178481 TYMP–SPP1

- Strict conservative CNV: n={int(g178_strict['n_patients'])}; not estimable.
- Expanded any-positive-CNV: n={int(g178_expanded['n_patients'])}, rho={fmt(g178_expanded['spearman_rho'])}, p={fmt(g178_expanded['p_value'])}, within-method FDR={fmt(g178_expanded['p_fdr_bh_within_dataset_method'])}.
- Provisional tumor-like: n={int(g178_provisional['n_patients'])}, rho={fmt(g178_provisional['spearman_rho'])}, p={fmt(g178_provisional['p_value'], 6)}, within-method FDR={fmt(g178_provisional['p_fdr_bh_within_dataset_method'])}.

The positive direction is retained in the expanded CNV-supported layer but reaches significance only in the broad tumor-like layer. Under the frozen small-sample rule, the expanded n=5 result is directional only, whereas the provisional n=8 significant result is reportable with its layer stated explicitly.

## GSE207493 TYMP–SPP1

- Strict conservative CNV: n={int(g207_strict['n_patients'])}, rho={fmt(g207_strict['spearman_rho'])}, p={fmt(g207_strict['p_value'])}.
- Expanded any-positive-CNV: n={int(g207_expanded['n_patients'])}, rho={fmt(g207_expanded['spearman_rho'])}, p={fmt(g207_expanded['p_value'])}.
- Provisional tumor-like: n={int(g207_provisional['n_patients'])}, rho={fmt(g207_provisional['spearman_rho'])}, p={fmt(g207_provisional['p_value'])}.

All three layers remain negative or near-null. Malignant expansion improves patient coverage but does not create a positive TYMP–SPP1 result in this cohort.

## Interpretation

The completed three-tier comparison resolves the threshold question without creating a new gate. GSE178481 provides tier-dependent positive support; GSE207493 provides consistent non-positive cohort evidence. The contrast should be reported as biological/cohort heterogeneity, not as a single universally replicated malignant-axis effect.

## Validation

- Frozen Phase 14 expanded results were reproduced for GSE178481 old and added subsets and for both full cohorts.
- Each layer contains {expected_rows_per_layer} unique association rows; the combined table contains {len(combined)} rows.
- Patient eligibility is nested as strict <= expanded <= provisional for all {len(patient_pivot)} dataset/subset/method/relation combinations.
- Validation checks: {len(validation_frame)}/{len(validation_frame)} PASS.
"""
    (PART / "PHASE52_INCREMENTAL_THREE_TIER_STATUS_v1.md").write_text(status_text, encoding="utf-8")

    metadata = {
        "schema_version": "phase52_part1_incremental_three_tier_v1",
        "status": "PASS",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "datasets": ["GSE178481", "GSE207493"],
        "subsets": ["old_subset", "added_subset", "expanded_full"],
        "evidence_layers": LAYERS,
        "rows_per_layer": expected_rows_per_layer,
        "combined_rows": len(combined),
        "validation_checks": len(validation_frame),
        "validation_pass": int((validation_frame["status"] == "PASS").sum()),
        "minimum_cells_per_patient_compartment": module.MIN_CELLS,
        "minimum_patients_for_correlation_test": module.MIN_PATIENTS,
        "bootstrap_replicates": module.N_BOOT,
        "inputs_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in required},
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "normalization_performed": False,
        "integration_performed": False,
        "clustering_performed": False,
        "annotation_performed": False,
    }
    (PART / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (PART / "PHASE52_INCREMENTAL_THREE_TIER_PASS_v1").write_text(
        "status=PASS\n"
        f"expanded_association_rows={len(expanded)}\n"
        f"combined_association_rows={len(combined)}\n"
        f"validation_checks={len(validation_frame)}\n"
        "threshold_changed=0\n"
        "normalization_performed=0\n"
        "integration_performed=0\n"
        "clustering_performed=0\n",
        encoding="utf-8",
    )
    os.replace(PART, FINAL)
    print(f"PASS expanded={len(expanded)} combined={len(combined)} validation={len(validation_frame)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
