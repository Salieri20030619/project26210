#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


GENES = ["TYMP", "TK1", "SLC29A1", "SLC29A2", "ODC1", "SRM", "SMS", "MTAP",
         "CHKA", "PCYT1A", "PCYT2", "CEPT1"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def add(rows: list[dict], name: str, passed: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def main() -> None:
    root = Path(__file__).resolve().parent
    results = root / "results_local_v1"
    output = root / "validation_local_v1"
    part = root / "validation_local_v1.part"
    if output.exists() or part.exists():
        raise FileExistsError("refusing to overwrite M17 validation")
    part.mkdir(parents=True)
    checks: list[dict] = []

    required = ["input_manifest_v1.tsv", "m17_gene_registry_v1.tsv", "m17_gene_detection_qc_v1.tsv",
                "m17_cell_expression_v1.tsv.gz", "m17_patient_gene_f1_effects_v1.tsv",
                "m17_gene_f1_summary_v1.tsv", "m17_gene_integrated_evidence_v1.tsv",
                "m17_patient_gene_pair_correlations_v1.tsv", "m17_gene_pair_correlation_summary_v1.tsv",
                "status.json"]
    missing = [x for x in required if not (results / x).is_file()]
    add(checks, "required_outputs", not missing, f"missing={missing}")

    manifest = pd.read_csv(results / "input_manifest_v1.tsv", sep="\t")
    hashes_ok = len(manifest) == 4
    for row in manifest.itertuples(index=False):
        p = root / "input_snapshot_local_v1" / row.file
        hashes_ok &= p.is_file() and sha256(p) == row.sha256
    add(checks, "input_hashes", bool(hashes_ok), f"n_inputs={len(manifest)}")

    registry = pd.read_csv(results / "m17_gene_registry_v1.tsv", sep="\t")
    registry_ok = len(registry) == 12 and set(registry["gene"]) == set(GENES) and registry["prespecified"].all()
    add(checks, "gene_registry", bool(registry_ok), f"n_genes={len(registry)}; n_panels={registry['panel'].nunique()}")

    expression = pd.read_csv(results / "m17_cell_expression_v1.tsv.gz", sep="\t")
    identity_ok = len(expression) == 11950 and expression["exact_cell_id"].nunique() == 11950
    identity_ok &= expression["curated_patient_id"].nunique() == 26 and expression["study"].nunique() == 5
    add(checks, "cell_cohort_identity", bool(identity_ok),
        f"cells={len(expression)}; patients={expression['curated_patient_id'].nunique()}; studies={expression['study'].nunique()}")

    safety_ok = True
    formula_max = 0.0
    for gene in GENES:
        raw = expression[f"{gene}_raw_count"].to_numpy(dtype=float)
        logv = expression[f"{gene}_log1p_cp10k"].to_numpy(dtype=float)
        safety_ok &= np.isfinite(raw).all() and (raw >= 0).all() and np.equal(raw, np.floor(raw)).all()
        expected = np.log1p(raw / expression["raw_library_size"].to_numpy(dtype=float) * 10000.0)
        formula_max = max(formula_max, float(np.max(np.abs(expected - logv))))
    add(checks, "raw_count_safety", bool(safety_ok), "12 target raw-count columns")
    add(checks, "log_cp10k_formula", formula_max < 1e-12, f"max_abs_diff={formula_max:.3g}")

    qc = pd.read_csv(results / "m17_gene_detection_qc_v1.tsv", sep="\t")
    qc_ok = len(qc) == 12
    for row in qc.itertuples(index=False):
        raw = expression[f"{row.gene}_raw_count"]
        qc_ok &= row.detected_cells == int((raw > 0).sum()) and np.isclose(row.total_counts, raw.sum())
    add(checks, "detection_qc_recalculation", bool(qc_ok), f"rows={len(qc)}")

    effects = pd.read_csv(results / "m17_patient_gene_f1_effects_v1.tsv", sep="\t")
    endpoints = {"F1_continuous", "F1_dominant_group", "F1_patient_quartile"}
    coverage_ok = set(effects["gene"]) == set(GENES) and set(effects["endpoint"]) == endpoints
    add(checks, "endpoint_coverage", bool(coverage_ok), f"rows={len(effects)}")
    thresholds_ok = (effects.loc[effects["endpoint"] == "F1_continuous", "n_cells"] >= 20).all()
    group = effects[effects["endpoint"] == "F1_dominant_group"]
    thresholds_ok &= (group["n_low"] >= 10).all() and (group["n_high"] >= 10).all()
    quartile = effects[effects["endpoint"] == "F1_patient_quartile"]
    thresholds_ok &= (quartile["n_low"] >= 5).all() and (quartile["n_high"] >= 5).all()
    add(checks, "patient_cell_thresholds", bool(thresholds_ok), "20/10/5 cell gates")

    summary = pd.read_csv(results / "m17_gene_f1_summary_v1.tsv", sep="\t")
    recomputed = effects.groupby(["gene", "panel", "endpoint", "method"])["effect"].agg(n="count", median="median").reset_index()
    joined = summary.merge(recomputed, on=["gene", "panel", "endpoint", "method"], validate="one_to_one")
    summary_ok = len(summary) == 36 and len(joined) == 36
    summary_ok &= (joined["n_eligible_patients"] == joined["n"]).all()
    summary_ok &= np.allclose(joined["median_patient_effect"], joined["median"], atol=1e-12)
    add(checks, "gene_summary_recalculation", bool(summary_ok), f"rows={len(summary)}")

    two_sided_ok = ((summary["statistical_support"] & (summary["direction"] == "positive")) ==
                    (summary["supported_direction"] == "positive")).all()
    two_sided_ok &= ((summary["statistical_support"] & (summary["direction"] == "negative")) ==
                     (summary["supported_direction"] == "negative")).all()
    add(checks, "two_sided_direction_rule", bool(two_sided_ok),
        str(summary["supported_direction"].value_counts().to_dict()))

    integrated = pd.read_csv(results / "m17_gene_integrated_evidence_v1.tsv", sep="\t")
    tier_expected = np.select([integrated["n_supported_endpoints_of_3"] == 3,
                               integrated["n_supported_endpoints_of_3"] == 2,
                               integrated["n_supported_endpoints_of_3"] == 1],
                              ["STRONG_ASSOCIATION", "SECONDARY_ASSOCIATION", "EXPLORATORY_ASSOCIATION"],
                              default="NO_SUPPORTED_ASSOCIATION")
    integrated_ok = len(integrated) == 12 and set(integrated["gene"]) == set(GENES)
    integrated_ok &= (tier_expected == integrated["evidence_tier"]).all()
    add(checks, "integrated_evidence_rule", bool(integrated_ok),
        str(integrated["evidence_tier"].value_counts().to_dict()))

    pair_detail = pd.read_csv(results / "m17_patient_gene_pair_correlations_v1.tsv", sep="\t")
    pair_summary = pd.read_csv(results / "m17_gene_pair_correlation_summary_v1.tsv", sep="\t")
    expected_pairs = {tuple(sorted(x)) for x in combinations(GENES, 2)}
    observed_pairs = {tuple(sorted((r.gene_a, r.gene_b))) for r in pair_summary.itertuples(index=False)}
    add(checks, "pairwise_66_pairs", observed_pairs == expected_pairs and len(pair_summary) == 66,
        f"pairs={len(pair_summary)}")
    pair_recomputed = pair_detail.groupby(["gene_a", "gene_b", "panel_a", "panel_b"])["effect"].agg(n="count", median="median").reset_index()
    pair_joined = pair_summary.merge(pair_recomputed, on=["gene_a", "gene_b", "panel_a", "panel_b"], validate="one_to_one")
    pair_ok = len(pair_joined) == 66 and (pair_joined["n_eligible_patients"] == pair_joined["n"]).all()
    pair_ok &= np.allclose(pair_joined["median_patient_effect"], pair_joined["median"], atol=1e-12)
    add(checks, "pairwise_summary_recalculation", bool(pair_ok), f"patient_pair_rows={len(pair_detail)}")

    status = json.loads((results / "status.json").read_text(encoding="utf-8"))
    status_ok = status["status"] == "PASS_M17_ANALYSIS_AWAITING_INDEPENDENT_VALIDATION"
    status_ok &= status["direction_prespecified"] is False and status["two_sided_testing"] is True
    status_ok &= status["figures_created"] is False
    add(checks, "analysis_status", bool(status_ok), status["status"])
    figures = [p.name for p in results.rglob("*") if p.suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg"}]
    add(checks, "no_figures", not figures, f"figures={figures}")

    frame = pd.DataFrame(checks)
    frame.to_csv(part / "verification_checks_v1.tsv", sep="\t", index=False)
    passed = int((frame["status"] == "PASS").sum())
    validation_status = {"status": "PASS_INDEPENDENT_VALIDATION" if passed == len(frame) else "FAIL",
                         "checks_passed": passed, "checks_total": len(frame),
                         "analysis_file_sha256": {p.name: sha256(p) for p in sorted(results.iterdir()) if p.is_file()}}
    (part / "status.json").write_text(json.dumps(validation_status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, output)
    if passed != len(frame):
        raise RuntimeError(f"M17 validation failed: {passed}/{len(frame)}")


if __name__ == "__main__":
    main()

