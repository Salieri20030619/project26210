#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


EXPECTED_INPUT = {"GSE264586": 95004, "GSE294109": 145316}
EXPECTED_LIBRARIES = {"GSE264586": 19, "GSE294109": 25}
EXPECTED_PATIENTS = {"GSE264586": 12, "GSE294109": 15}


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--qc-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    checks: list[dict] = []

    def check(name: str, passed: bool, observed, expected) -> None:
        checks.append({"check": name, "pass": bool(passed), "observed": str(observed), "expected": str(expected)})

    check("manifest_rows", len(manifest) == 44, len(manifest), 44)
    summaries = []
    for _, row in manifest.iterrows():
        library = row["library_id"]
        checkpoint = args.qc_root / "checkpoints" / f"{library}.json"
        h5ad = args.qc_root / "libraries" / f"{library}.filtered_v1.h5ad"
        cell_qc = args.qc_root / "cell_qc" / f"{library}.cell_qc_v1.tsv.gz"
        check(f"{library}_checkpoint_exists", checkpoint.is_file(), checkpoint.is_file(), True)
        check(f"{library}_h5ad_exists", h5ad.is_file(), h5ad.is_file(), True)
        check(f"{library}_cell_qc_exists", cell_qc.is_file(), cell_qc.is_file(), True)
        if not (checkpoint.is_file() and h5ad.is_file() and cell_qc.is_file()):
            continue
        payload = json.loads(checkpoint.read_text(encoding="utf-8"))
        table = pd.read_csv(cell_qc, sep="\t", compression="gzip", index_col=0)
        obj = ad.read_h5ad(h5ad, backed="r")
        check(f"{library}_status", payload.get("status") == "PASS", payload.get("status"), "PASS")
        check(f"{library}_input_conserved", len(table) == payload["n_input"], len(table), payload["n_input"])
        check(f"{library}_retained_conserved", obj.n_obs == payload["n_retained"], obj.n_obs, payload["n_retained"])
        check(f"{library}_retained_positive", obj.n_obs > 0, obj.n_obs, ">0")
        check(f"{library}_counts_layer", "counts" in obj.layers, list(obj.layers.keys()), "contains counts")
        check(f"{library}_unique_cells", not obj.obs_names.duplicated().any(), obj.obs_names.duplicated().sum(), 0)
        obj.file.close()
        summaries.append(payload)

    summary = pd.DataFrame(summaries)
    for cohort in EXPECTED_INPUT:
        subset = summary.loc[summary["cohort"].eq(cohort)]
        check(f"{cohort}_library_count", len(subset) == EXPECTED_LIBRARIES[cohort], len(subset), EXPECTED_LIBRARIES[cohort])
        check(f"{cohort}_input_barcodes", int(subset["n_input"].sum()) == EXPECTED_INPUT[cohort], int(subset["n_input"].sum()), EXPECTED_INPUT[cohort])
        cohort_path = args.qc_root / "cohorts" / f"{cohort}.qc_filtered_raw_v1.h5ad"
        check(f"{cohort}_combined_exists", cohort_path.is_file(), cohort_path.is_file(), True)
        if cohort_path.is_file():
            obj = ad.read_h5ad(cohort_path, backed="r")
            check(f"{cohort}_combined_cells", obj.n_obs == int(subset["n_retained"].sum()), obj.n_obs, int(subset["n_retained"].sum()))
            check(f"{cohort}_combined_libraries", obj.obs["library_id"].nunique() == EXPECTED_LIBRARIES[cohort], obj.obs["library_id"].nunique(), EXPECTED_LIBRARIES[cohort])
            check(f"{cohort}_combined_patients", obj.obs["patient_id"].nunique() == EXPECTED_PATIENTS[cohort], obj.obs["patient_id"].nunique(), EXPECTED_PATIENTS[cohort])
            check(f"{cohort}_combined_unique_cells", not obj.obs_names.duplicated().any(), obj.obs_names.duplicated().sum(), 0)
            obj.file.close()

    results = pd.DataFrame(checks)
    results.to_csv(args.output / "qc_verification_checks_v1.tsv", sep="\t", index=False)
    failed = results.loc[~results["pass"]]
    payload = {
        "status": "PASS" if failed.empty else "FAIL",
        "verified_at": now(), "n_checks": int(len(results)), "n_failed": int(len(failed)),
        "failed_checks": failed["check"].tolist(),
    }
    (args.output / "qc_verification_summary_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, ensure_ascii=False), flush=True)
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
