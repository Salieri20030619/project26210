#!/usr/bin/env python3
"""Independent post-run validation for the Phase 55 cell2location v2 outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_h5ad_shape(path: Path) -> tuple[int, int, set[str]]:
    obj = ad.read_h5ad(path, backed="r")
    try:
        return int(obj.n_obs), int(obj.n_vars), set(obj.obsm.keys())
    finally:
        if obj.file is not None:
            obj.file.close()


def add_check(checks: list[dict], name: str, observed, expected, passed: bool) -> None:
    checks.append(
        {
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": str(observed),
            "expected": str(expected),
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--server-root", required=True, type=Path)
    args = parser.parse_args()

    server_root = args.server_root.resolve()
    run_root = server_root / "run_v1"
    if not (server_root / "PHASE55_BACKGROUND_PASS_v1").exists():
        print(json.dumps({"status": "INCOMPLETE", "server_root": str(server_root)}, sort_keys=True))
        return 75

    final = server_root / "postrun_validation_v2"
    part = server_root / "postrun_validation_v2.part"
    if final.exists() or part.exists():
        raise FileExistsError("Existing post-run validation output found")
    part.mkdir(parents=True)

    checks: list[dict] = []
    section_rows: list[dict] = []
    manifest_rows: list[dict] = []

    reference_dir = run_root / "reference_model_v1"
    signature_path = reference_dir / "reference_cell_state_signatures_v1.tsv.gz"
    reference_metadata_path = reference_dir / "reference_metadata_v1.json"
    reference_posterior_path = reference_dir / "reference_posterior_v1.h5ad"
    reference_pass = reference_dir / "REFERENCE_MODEL_PASS_v1"
    add_check(checks, "reference_pass_marker", reference_pass.exists(), True, reference_pass.exists())
    for path in [signature_path, reference_metadata_path, reference_posterior_path]:
        add_check(checks, f"reference_file_{path.name}", path.exists(), True, path.exists())

    signatures = pd.read_csv(signature_path, sep="\t", compression="gzip")
    if "gene" not in signatures.columns:
        raise KeyError("Reference signature table lacks gene column")
    factors = [str(column) for column in signatures.columns if column != "gene"]
    add_check(checks, "reference_signature_genes", len(signatures), ">=1000", len(signatures) >= 1000)
    add_check(checks, "reference_signature_cell_types", len(factors), ">=2", len(factors) >= 2)
    add_check(checks, "reference_signature_unique_genes", signatures["gene"].nunique(), len(signatures), signatures["gene"].nunique() == len(signatures))
    signature_values = signatures[factors].to_numpy(dtype=float)
    add_check(checks, "reference_signature_finite", np.isfinite(signature_values).all(), True, bool(np.isfinite(signature_values).all()))
    add_check(checks, "reference_signature_nonnegative", float(np.nanmin(signature_values)), ">=0", bool(np.nanmin(signature_values) >= 0))

    reference_meta = json.loads(reference_metadata_path.read_text(encoding="utf-8"))
    reference_obs, reference_vars, _ = read_h5ad_shape(reference_posterior_path)
    add_check(checks, "reference_posterior_cells", reference_obs, int(reference_meta["cells"]), reference_obs == int(reference_meta["cells"]))
    add_check(checks, "reference_posterior_genes", reference_vars, int(reference_meta["genes"]), reference_vars == int(reference_meta["genes"]))

    status_path = run_root / "section_deconvolution_status_v1.tsv"
    failures_path = run_root / "section_deconvolution_failures_v1.tsv"
    internal_checks_path = run_root / "validation_checks_v1.tsv"
    status = pd.read_csv(status_path, sep="\t")
    failures = pd.read_csv(failures_path, sep="\t")
    internal_checks = pd.read_csv(internal_checks_path, sep="\t")
    add_check(checks, "section_status_rows", len(status), 24, len(status) == 24)
    add_check(checks, "section_status_unique_samples", status["sample_id"].nunique(), 24, status["sample_id"].nunique() == 24)
    add_check(checks, "section_failure_rows", len(failures), 0, len(failures) == 0)
    add_check(checks, "internal_checks_all_pass", bool((internal_checks["status"] == "PASS").all()), True, bool((internal_checks["status"] == "PASS").all()))

    sections_root = run_root / "sections_v1"
    section_dirs = sorted(path for path in sections_root.iterdir() if path.is_dir() and not path.name.endswith(".part"))
    add_check(checks, "section_output_directories", len(section_dirs), 24, len(section_dirs) == 24)
    expected_obsm = {"means_cell_abundance_w_sf", "q05_cell_abundance_w_sf", "q95_cell_abundance_w_sf"}

    key_paths = [signature_path, reference_metadata_path, reference_posterior_path, status_path, failures_path, internal_checks_path]
    for section_dir in section_dirs:
        sample = section_dir.name
        pass_path = section_dir / "DECONVOLUTION_PASS_v1"
        summary_path = section_dir / "deconvolution_summary_v1.json"
        abundance_path = section_dir / "spot_cell_abundance_v1.tsv.gz"
        posterior_path = section_dir / "cell2location_posterior_v1.h5ad"
        required_exist = all(path.exists() for path in [pass_path, summary_path, abundance_path, posterior_path])
        if not required_exist:
            section_rows.append({"sample_id": sample, "status": "FAIL", "failure": "missing_required_output"})
            continue

        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        abundance = pd.read_csv(abundance_path, sep="\t", compression="gzip")
        posterior_obs, posterior_vars, obsm_keys = read_h5ad_shape(posterior_path)
        mean_columns = [f"mean__{factor}" for factor in factors]
        q05_columns = [f"q05__{factor}" for factor in factors]
        q95_columns = [f"q95__{factor}" for factor in factors]
        abundance_columns_ok = all(column in abundance.columns for column in mean_columns + q05_columns + q95_columns)
        numeric = abundance[mean_columns + q05_columns + q95_columns].to_numpy(dtype=float) if abundance_columns_ok else np.empty((0, 0))
        finite = bool(np.isfinite(numeric).all()) if numeric.size else False
        nonnegative = bool(np.nanmin(numeric) >= 0) if numeric.size else False
        ordered = True
        if abundance_columns_ok:
            for factor in factors:
                ordered &= bool((abundance[f"q05__{factor}"] <= abundance[f"mean__{factor}"] + 1e-8).all())
                ordered &= bool((abundance[f"mean__{factor}"] <= abundance[f"q95__{factor}"] + 1e-8).all())
        else:
            ordered = False

        tests = {
            "pass_marker": pass_path.exists(),
            "summary_status": summary.get("status") == "PASS",
            "summary_sample": summary.get("sample_id") == sample,
            "rows_match_summary": len(abundance) == int(summary.get("spots", -1)),
            "posterior_obs_match": posterior_obs == int(summary.get("spots", -1)),
            "posterior_vars_match": posterior_vars == int(summary.get("shared_genes", -1)),
            "barcodes_unique": "barcode" in abundance.columns and abundance["barcode"].nunique() == len(abundance),
            "abundance_columns": abundance_columns_ok,
            "abundance_finite": finite,
            "abundance_nonnegative": nonnegative,
            "posterior_quantiles_ordered": ordered,
            "posterior_obsm_slots": expected_obsm.issubset(obsm_keys),
            "cell_type_count": int(summary.get("cell_types", -1)) == len(factors),
        }
        passed = all(tests.values())
        section_rows.append(
            {
                "sample_id": sample,
                "status": "PASS" if passed else "FAIL",
                "spots": len(abundance),
                "shared_genes": posterior_vars,
                "cell_types": len(factors),
                "failed_checks": ";".join(name for name, value in tests.items() if not value),
            }
        )
        key_paths.extend([summary_path, abundance_path, posterior_path])

    section_frame = pd.DataFrame(section_rows).sort_values("sample_id")
    section_frame.to_csv(part / "section_validation_v2.tsv", sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    section_passed = int((section_frame["status"] == "PASS").sum()) if len(section_frame) else 0
    add_check(checks, "independent_section_validations", section_passed, 24, section_passed == 24)

    for path in sorted(set(key_paths)):
        manifest_rows.append(
            {
                "relative_path": path.relative_to(server_root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(part / "KEY_FILE_MANIFEST_SHA256_v2.tsv", sep="\t", index=False, encoding="utf-8", lineterminator="\n")

    check_frame = pd.DataFrame(checks)
    check_frame.to_csv(part / "validation_checks_v2.tsv", sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    overall_pass = bool((check_frame["status"] == "PASS").all())
    summary = {
        "status": "PASS" if overall_pass else "FAIL",
        "generated_at": now(),
        "server_root": str(server_root),
        "checks_passed": int((check_frame["status"] == "PASS").sum()),
        "checks_total": len(check_frame),
        "sections_passed": section_passed,
        "sections_total": len(section_frame),
        "key_files_hashed": len(manifest_rows),
    }
    (part / "validation_summary_v2.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    marker_name = "POSTRUN_VALIDATION_PASS_v2" if overall_pass else "POSTRUN_VALIDATION_FAIL_v2"
    (part / marker_name).write_text(
        f"status={summary['status']}\nchecks_passed={summary['checks_passed']}\nchecks_total={summary['checks_total']}\nsections_passed={section_passed}\n",
        encoding="utf-8",
    )
    os.replace(part, final)
    print(json.dumps(summary, sort_keys=True))
    return 0 if overall_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())

