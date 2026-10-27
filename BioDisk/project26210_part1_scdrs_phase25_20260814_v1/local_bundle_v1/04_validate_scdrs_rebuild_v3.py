#!/usr/bin/env python3
"""Independent terminal validation of Phase 25 rebuilt scDRS analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True, type=Path)
    parser.add_argument("--cnmf-usage", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--compute-dir", required=True, type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    parser.add_argument("--min-cells-per-patient", type=int, default=20)
    args = parser.parse_args()

    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(f"refusing stale partial output: {part}")
    part.mkdir(parents=True)
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    check("input_terminal_status", load_json(args.input_dir / "status.json")["status"] == "PASS_INPUT_COMPLETE", "input status")
    check("compute_terminal_status", load_json(args.compute_dir / "status.json")["status"] == "PASS_SCDRS_COMPUTE", "compute status")
    check("analysis_terminal_status", load_json(args.analysis_dir / "status.json")["status"] == "PASS_ANALYSIS_COMPLETE", "analysis status")

    h5ad_path = args.input_dir / "malignant_eligible_raw_counts_phase25_v1.h5ad"
    data = ad.read_h5ad(h5ad_path, backed="r")
    ids = pd.Index(data.obs_names.astype(str))
    check("h5ad_shape", data.n_obs == args.expected_cells and data.n_vars > 0, f"{data.n_obs}x{data.n_vars}")
    check("h5ad_unique_ids", ids.is_unique, f"unique={ids.is_unique}")

    registry = pd.read_csv(args.registry, sep="\t", low_memory=False)
    flag = registry["eligible_for_subclone_selection_v1"].astype(str).str.lower().isin({"true", "t", "1", "yes"})
    reg_ids = pd.Index(registry.loc[flag, "exact_cell_id"].astype(str))
    check("registry_identity_and_order", ids.equals(reg_ids), f"registry_n={len(reg_ids)}")

    usage = pd.read_csv(args.cnmf_usage, sep="\t").set_index("exact_cell_id")
    source_f1 = pd.to_numeric(usage.reindex(ids)["F1"], errors="coerce").to_numpy(float)
    h5ad_f1 = pd.to_numeric(data.obs["F1"], errors="coerce").to_numpy(float)
    check("cnmf_f1_exact", np.array_equal(source_f1, h5ad_f1), "exact F1 comparison")

    finite = nonnegative = integer = True
    library = []
    detected = []
    for start in range(0, data.n_obs, 512):
        block = sp.csr_matrix(data.X[start:start + 512, :])
        finite &= bool(np.isfinite(block.data).all())
        nonnegative &= bool((block.data >= 0).all())
        integer &= bool(np.equal(block.data, np.floor(block.data)).all())
        library.append(np.asarray(block.sum(axis=1)).ravel())
        detected.append(np.diff(block.indptr))
    library = np.concatenate(library).astype(float)
    detected = np.concatenate(detected)
    check("raw_count_safety", finite and nonnegative and integer, f"finite={finite}; nonnegative={nonnegative}; integer={integer}")
    check(
        "raw_count_metadata_conservation",
        np.allclose(library, pd.to_numeric(data.obs["raw_library_size_v1"]), rtol=0, atol=1e-8)
        and np.array_equal(detected, pd.to_numeric(data.obs["raw_detected_genes_v1"]).to_numpy(int)),
        "library size and detected genes",
    )

    cov = pd.read_csv(args.input_dir / "malignant_eligible_phase25_v1.cov", sep="\t", index_col=0)
    cov.index = cov.index.astype(str)
    check(
        "covariate_identity_and_values",
        cov.index.equals(ids)
        and np.array_equal(cov["n_genes"].to_numpy(int), detected)
        and np.allclose(cov["log_n_counts"].to_numpy(float), np.log1p(library), rtol=0, atol=1e-12),
        f"cov_n={len(cov)}",
    )

    coverage = pd.read_csv(args.input_dir / "gene_set_coverage_v1.tsv", sep="\t")
    n_present = int(coverage["present_in_matrix"].astype(bool).sum())
    check("gene_set_coverage", len(coverage) == 20 and n_present >= 10, f"{n_present}/20 present")

    score_paths = [
        path for path in (args.compute_dir / "scores").glob("*.score.gz")
        if not path.name.endswith(".full_score.gz")
    ]
    check("one_primary_score_file", len(score_paths) == 1, f"n={len(score_paths)}")
    if len(score_paths) != 1:
        raise ValueError("primary score file count invalid")
    score = pd.read_csv(score_paths[0], sep="\t", index_col=0)
    score.index = score.index.astype(str)
    check("score_cell_identity", score.index.is_unique and set(score.index) == set(ids), f"score_n={len(score)}")
    score = score.reindex(ids)
    norm = pd.to_numeric(score["norm_score"], errors="coerce").to_numpy(float)
    pval = pd.to_numeric(score["pval"], errors="coerce").to_numpy(float)
    check("score_numeric_safety", np.isfinite(norm).all() and np.isfinite(pval).all() and ((pval >= 0) & (pval <= 1)).all(), "norm finite; P in [0,1]")

    patients = pd.read_csv(args.analysis_dir / "patient_within_association_v1.tsv", sep="\t")
    recomputed = []
    obs = data.obs[["study", "curated_patient_id", "F1"]].copy()
    obs["norm"] = norm
    for (study, patient), frame in obs.groupby(["study", "curated_patient_id"], sort=True):
        f1 = pd.to_numeric(frame["F1"]).to_numpy(float)
        sc = frame["norm"].to_numpy(float)
        eligible = len(frame) >= args.min_cells_per_patient and np.ptp(f1) > 0 and np.ptp(sc) > 0
        rho = float(spearmanr(f1, sc).statistic) if eligible else np.nan
        recomputed.append((study, patient, len(frame), eligible, rho))
    rec = pd.DataFrame(recomputed, columns=["study", "curated_patient_id", "n_cells", "eligible_primary", "rho"])
    joined = patients.merge(rec, on=["study", "curated_patient_id"], validate="one_to_one")
    same_rho = np.allclose(joined["spearman_rho_f1_vs_scdrs"], joined["rho"], equal_nan=True, rtol=0, atol=1e-12)
    same_meta = np.array_equal(joined["n_cells_x"], joined["n_cells_y"]) and np.array_equal(joined["eligible_primary_x"], joined["eligible_primary_y"])
    check("patient_correlations_reproduced", same_rho and same_meta, f"patients={len(joined)}")

    overall = load_json(args.analysis_dir / "overall_patient_association_v1.json")
    n_eligible = int(rec["eligible_primary"].sum())
    tier_ok = overall["reporting_tier"] == ("full_statistics" if n_eligible >= 10 else "direction_only")
    inferential_fields = {"bootstrap_ci95_low", "bootstrap_ci95_high", "signflip_p_two_sided"}
    fields_ok = inferential_fields.issubset(overall) if n_eligible >= 10 else inferential_fields.isdisjoint(overall)
    check("patient_reporting_rule", tier_ok and fields_ok, f"eligible={n_eligible}; tier={overall['reporting_tier']}")

    descriptive = load_json(args.analysis_dir / "cell_level_descriptive_v1.json")
    check("no_pooled_cell_inference", descriptive.get("pooled_cell_p_value_calculated") is False, "descriptive pooled r only")
    check("five_studies_traceable", data.obs["study"].nunique() == 5, f"studies={data.obs['study'].nunique()}")
    data.file.close()

    check_df = pd.DataFrame(checks)
    check_df.to_csv(part / "INDEPENDENT_VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    passed = int(check_df["passed"].sum())
    total = int(len(check_df))
    terminal = "PASS_INDEPENDENT_VALIDATION" if passed == total else "FAIL_INDEPENDENT_VALIDATION"
    status = {"status": terminal, "checks_passed": passed, "checks_total": total}
    status_path = part / "status.json"
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8", newline="\n")
    manifest = pd.DataFrame(
        [
            {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted(part.iterdir()) if path.is_file()
        ]
    )
    manifest.to_csv(part / "VALIDATION_FILE_MANIFEST_v1.tsv", sep="\t", index=False)
    os.replace(part, args.output)
    if terminal != "PASS_INDEPENDENT_VALIDATION":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
