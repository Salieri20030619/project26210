#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


EXPECTED_LIBRARIES = {"GSE264586": 19, "GSE294109": 25}
EXPECTED_PATIENTS = {"GSE264586": 12, "GSE294109": 15}
REQUIRED_OBS = {
    "patient_id", "library_id", "response_binary", "timepoint", "leiden06_preliminary_v1",
    "cell_type_broad_preliminary_v1", "F1_projection_score_v1", "TYMP_log1p_CP10K",
    "hif1a_downstream_5gene_score_v1", "ODC1_log1p_CP10K", "SRM_log1p_CP10K",
    "PCYT2_log1p_CP10K",
}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--integration-root", type=Path, required=True)
    parser.add_argument("--qc-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []

    def add(name: str, passed: bool, observed, expected) -> None:
        rows.append({"check": name, "pass": bool(passed), "observed": str(observed), "expected": str(expected)})

    for cohort in EXPECTED_LIBRARIES:
        status_path = args.integration_root / f"{cohort}.integration_status_v1.json"
        integrated_path = args.integration_root / f"{cohort}.integrated_annotated_preliminary_v1.h5ad"
        source_path = args.qc_root / "cohorts" / f"{cohort}.qc_filtered_raw_v1.h5ad"
        add(f"{cohort}_status_exists", status_path.is_file(), status_path.is_file(), True)
        add(f"{cohort}_h5ad_exists", integrated_path.is_file(), integrated_path.is_file(), True)
        if not (status_path.is_file() and integrated_path.is_file() and source_path.is_file()):
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        add(f"{cohort}_status_pass", status.get("status") == "PASS", status.get("status"), "PASS")
        obj = ad.read_h5ad(integrated_path, backed="r")
        source = ad.read_h5ad(source_path, backed="r")
        add(f"{cohort}_cell_conservation", obj.n_obs == source.n_obs, obj.n_obs, source.n_obs)
        add(f"{cohort}_gene_conservation", obj.n_vars == source.n_vars, obj.n_vars, source.n_vars)
        add(f"{cohort}_library_count", obj.obs["library_id"].nunique() == EXPECTED_LIBRARIES[cohort], obj.obs["library_id"].nunique(), EXPECTED_LIBRARIES[cohort])
        add(f"{cohort}_patient_count", obj.obs["patient_id"].nunique() == EXPECTED_PATIENTS[cohort], obj.obs["patient_id"].nunique(), EXPECTED_PATIENTS[cohort])
        missing = sorted(REQUIRED_OBS.difference(obj.obs.columns))
        add(f"{cohort}_required_obs", not missing, missing, [])
        add(f"{cohort}_counts_layer", "counts" in obj.layers, list(obj.layers.keys()), "contains counts")
        add(f"{cohort}_umap_present", "X_umap" in obj.obsm, list(obj.obsm.keys()), "contains X_umap")
        if "X_umap" in obj.obsm:
            coords = np.asarray(obj.obsm["X_umap"])
            add(f"{cohort}_umap_finite", np.isfinite(coords).all() and coords.shape == (obj.n_obs, 2), coords.shape, (obj.n_obs, 2))
        add(f"{cohort}_clusters_positive", obj.obs["leiden06_preliminary_v1"].nunique() >= 5, obj.obs["leiden06_preliminary_v1"].nunique(), ">=5")
        add(f"{cohort}_epithelial_positive", obj.obs["cell_type_broad_preliminary_v1"].eq("Epithelial").sum() > 0, obj.obs["cell_type_broad_preliminary_v1"].eq("Epithelial").sum(), ">0")
        f1_meta = obj.uns.get("F1_projection_definition_v1", {})
        add(f"{cohort}_f1_coverage", float(f1_meta.get("coverage", 0)) >= 0.60, f1_meta.get("coverage"), ">=0.60")
        boundary = obj.uns.get("phase60_preliminary_annotation_boundary_v1", {})
        add(f"{cohort}_preliminary_boundary", boundary.get("malignant_gate_complete") is False and boundary.get("response_screen_formal") is False, boundary, "malignant_gate_complete=false; response_screen_formal=false")
        for suffix in ["preliminary_patient_endpoint_summary_v1.tsv", "preliminary_patient_cell_composition_v1.tsv", "preliminary_response_screen_v1.tsv", "preliminary_cluster_marker_evidence_v1.tsv", "marker_panel_coverage_v1.tsv", "preliminary_broad_annotation_umap_v1.png"]:
            path = args.integration_root / f"{cohort}.{suffix}"
            add(f"{cohort}_{suffix}_exists", path.is_file() and path.stat().st_size > 0, path.stat().st_size if path.is_file() else 0, ">0")
        obj.file.close(); source.file.close()

    results = pd.DataFrame(rows)
    results.to_csv(args.output / "integration_verification_checks_v1.tsv", sep="\t", index=False)
    failed = results.loc[~results["pass"]]
    payload = {"status": "PASS" if failed.empty else "FAIL", "verified_at": now(),
               "n_checks": len(results), "n_failed": len(failed), "failed_checks": failed["check"].tolist()}
    (args.output / "integration_verification_summary_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, ensure_ascii=False), flush=True)
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

