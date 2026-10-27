#!/usr/bin/env python3
"""Independent structural and claim validation for Phase 16 TJRCC12 analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon


SAMPLES = ["H46_T", "J38_T", "R114_T", "R29_T", "R51_T", "S15_T", "X49_T", "X98_T", "Y12_T", "Y27_T", "Y7_T", "Z43_T"]
PRIMARY = ["F1polar_to_C1QC_expression", "F1polar_to_hypoxia_expression", "SPP1polar_to_Tdysfunction_expression"]
EXPECTED_SPOTS = 34609


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def bh(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty(len(values), float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def decision(row: pd.Series) -> str:
    direction = 1 if row.median_rho > 0 else -1 if row.median_rho < 0 else 0
    direction_count = row.positive_sections if direction > 0 else row.negative_sections
    if row.wilcoxon_fdr_bh < 0.05 and direction_count >= 8 and abs(row.median_rho) >= 0.10 and row.loo_direction_fraction == 1.0:
        return "AUXILIARY_REPRODUCED"
    if row.positive_sections >= 8 or row.negative_sections >= 8:
        return "DIRECTIONAL_NOT_FDR"
    if row.positive_sections >= 4 and row.negative_sections >= 4:
        return "HETEROGENEOUS"
    return "NOT_SUPPORTED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--signature-genes", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    checks: list[dict[str, str]] = []

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "result": "PASS" if passed else "FAIL", "detail": detail})

    required = [
        "PHASE16_TJRCC12_AUX_STATUS_v1.json", "OUTPUT_MANIFEST_v1.tsv",
        "figures/primary_section_rho_v1.png", "tables/input_provenance_v1.tsv",
        "tables/signature_gene_overlap_v1.tsv", "tables/spot_program_scores_v1.tsv.gz",
        "tables/section_correlations_all_v1.tsv", "tables/primary_section_statistics_v1.tsv",
        "tables/primary_across_section_inference_v1.tsv", "tables/primary_leave_one_section_out_v1.tsv",
    ]
    missing = [value for value in required if not (args.results / value).is_file() or (args.results / value).stat().st_size == 0]
    record("required_outputs_nonempty", not missing, f"required={len(required)} missing={';'.join(missing) or 'none'}")
    parts = list(args.results.rglob("*.part"))
    record("no_partial_files", not parts, f"part_files={len(parts)}")

    status = json.loads((args.results / required[0]).read_text(encoding="utf-8"))
    status_ok = status.get("status") == "PASS_ANALYSIS_COMPLETE" and status.get("biological_units") == 12 and status.get("spots") == EXPECTED_SPOTS and status.get("spot_as_replicate_inference") is False
    record("status_and_inference_scope", status_ok, f"status={status.get('status')} units={status.get('biological_units')} spots={status.get('spots')} spot_inference={status.get('spot_as_replicate_inference')}")
    record("frozen_input_hashes", status.get("registry_sha256") == sha256(args.registry) and status.get("signature_sha256") == sha256(args.signature_genes), f"registry={status.get('registry_sha256')} signature={status.get('signature_sha256')}")

    manifest = pd.read_csv(args.results / "OUTPUT_MANIFEST_v1.tsv", sep="\t")
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = args.results / row.file
        manifest_ok &= path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256
    record("output_manifest_rehash", bool(manifest_ok and len(manifest) == 8), f"rows={len(manifest)}")

    provenance = pd.read_csv(args.results / "tables/input_provenance_v1.tsv", sep="\t")
    source_ok = set(provenance["sample"]) == set(SAMPLES) and provenance["sample"].is_unique and provenance["n_spots"].sum() == EXPECTED_SPOTS
    for row in provenance.itertuples(index=False):
        path = Path(row.path)
        source_status = json.loads((args.input_root / row.sample / "status.json").read_text(encoding="utf-8"))
        source_ok &= path.is_file() and path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256
        source_ok &= source_status.get("status") == "completed" and int(source_status.get("epoch", -1)) == 10000 and source_status.get("output_sha256") == row.sha256
    record("tjrcc12_source_rehash_and_terminal_status", bool(source_ok), f"samples={len(provenance)} spots={provenance['n_spots'].sum()}")

    overlap = pd.read_csv(args.results / "tables/signature_gene_overlap_v1.tsv", sep="\t")
    minima = {"f1_top100_weighted": 74, "c1qc_state": 3, "spp1_state_no_spp1": 5, "spp1_state_full": 6, "hypoxia_angiogenesis": 5, "tcell_dysfunction": 7, "cytotoxicity": 5}
    overlap_ok = len(overlap) == 84 and overlap.groupby("signature_id")["sample"].nunique().eq(12).all()
    for signature_id, minimum in minima.items():
        overlap_ok &= int(overlap.loc[overlap.signature_id == signature_id, "present_genes"].min()) >= minimum
    record("frozen_signature_coverage", bool(overlap_ok), "84 rows; all prespecified minima checked")

    scores = pd.read_csv(args.results / "tables/spot_program_scores_v1.tsv.gz", sep="\t", usecols=["sample", "spot_id"])
    score_ok = len(scores) == EXPECTED_SPOTS and not scores.duplicated(["sample", "spot_id"]).any() and set(scores["sample"]) == set(SAMPLES)
    record("spot_identity_conservation", bool(score_ok), f"rows={len(scores)} duplicate_keys={scores.duplicated(['sample','spot_id']).sum()}")

    correlations = pd.read_csv(args.results / "tables/section_correlations_all_v1.tsv", sep="\t")
    p_columns = [column for column in correlations.columns if "p_value" in column.lower() or column.lower() in {"p", "pvalue", "fdr", "qvalue"}]
    primary = pd.read_csv(args.results / "tables/primary_section_statistics_v1.tsv", sep="\t")
    primary_ok = len(primary) == 36 and primary.groupby("relation")["sample"].nunique().eq(12).all() and set(primary.relation) == set(PRIMARY)
    record("section_descriptive_and_primary_grid", bool(primary_ok and not p_columns), f"primary_rows={len(primary)} spot_p_columns={';'.join(p_columns) or 'none'}")

    aggregate = pd.read_csv(args.results / "tables/primary_across_section_inference_v1.tsv", sep="\t")
    loo = pd.read_csv(args.results / "tables/primary_leave_one_section_out_v1.tsv", sep="\t")
    recomputed = []
    for relation in PRIMARY:
        values = primary.loc[primary.relation == relation, "spearman_rho"].to_numpy(float)
        nonzero = values[values != 0]
        sign_p = binomtest(int((nonzero > 0).sum()), len(nonzero), 0.5).pvalue if len(nonzero) else 1.0
        try:
            wilcoxon_p = wilcoxon(values, alternative="two-sided", zero_method="wilcox", method="auto").pvalue
        except ValueError:
            wilcoxon_p = 1.0
        loo_values = []
        for sample in SAMPLES:
            kept = primary[(primary.relation == relation) & (primary["sample"] != sample)].spearman_rho.to_numpy(float)
            loo_values.append(float(np.median(kept)))
        recomputed.append({"relation": relation, "median_rho": np.median(values), "positive_sections": (values > 0).sum(), "negative_sections": (values < 0).sum(), "sign_test_p_two_sided": sign_p, "wilcoxon_p_two_sided": wilcoxon_p, "loo_direction_fraction": np.mean(np.sign(loo_values) == np.sign(np.median(values)))})
    recomputed = pd.DataFrame(recomputed)
    recomputed["wilcoxon_fdr_bh"] = bh(recomputed.wilcoxon_p_two_sided.to_numpy(float))
    merged = aggregate.merge(recomputed, on="relation", suffixes=("_reported", "_audit"), validate="one_to_one")
    numeric = ["median_rho", "sign_test_p_two_sided", "wilcoxon_p_two_sided", "wilcoxon_fdr_bh", "loo_direction_fraction"]
    inference_ok = all(np.allclose(merged[f"{column}_reported"], merged[f"{column}_audit"], atol=1e-12, rtol=1e-10) for column in numeric)
    inference_ok &= (merged.positive_sections_reported == merged.positive_sections_audit).all() and (merged.negative_sections_reported == merged.negative_sections_audit).all()
    decisions = aggregate.apply(decision, axis=1)
    inference_ok &= (decisions.to_numpy() == aggregate.decision.to_numpy()).all()
    loo_ok = len(loo) == 36 and loo.groupby("relation")["left_out_sample"].nunique().eq(12).all()
    record("primary_section_level_inference_recomputed", bool(inference_ok and loo_ok), f"relations={len(aggregate)} loo_rows={len(loo)}")

    controls = correlations[(correlations.relation.isin(["C1QC_mapping_positive_control", "SPP1_mapping_positive_control"])) & (correlations.score_method == "primary") & (correlations.outcome == "lineage_depth_residual") & (correlations.analysis == "same_spot_q25")]
    control_summary = controls.groupby("relation").spearman_rho.agg(["median", lambda x: int((x > 0).sum())])
    controls_ok = len(controls) == 24 and (control_summary["median"] > 0).all() and (control_summary["<lambda_0>"] >= 8).all()
    record("mapping_positive_controls", bool(controls_ok), f"rows={len(controls)} medians={control_summary['median'].round(3).to_dict() if len(control_summary) else {}}")

    expected_decisions = {"F1polar_to_C1QC_expression": "HETEROGENEOUS", "F1polar_to_hypoxia_expression": "DIRECTIONAL_NOT_FDR", "SPP1polar_to_Tdysfunction_expression": "DIRECTIONAL_NOT_FDR"}
    claim_ok = aggregate.set_index("relation")["decision"].to_dict() == expected_decisions and (aggregate.wilcoxon_fdr_bh >= 0.05).all()
    record("hypothesis_claim_boundaries", bool(claim_ok), f"decisions={aggregate.set_index('relation')['decision'].to_dict()} no_fdr={(aggregate.wilcoxon_fdr_bh >= 0.05).all()}")

    overall = all(row["result"] == "PASS" for row in checks)
    checks.append({"check": "OVERALL", "result": "PASS" if overall else "FAIL", "detail": datetime.now().astimezone().isoformat(timespec="seconds")})
    atomic_tsv(pd.DataFrame(checks), args.output / "INDEPENDENT_VALIDATION_STATUS_v1.tsv")
    atomic_json({"status": "PASS" if overall else "FAIL", "checks_passed": sum(row["result"] == "PASS" for row in checks[:-1]), "checks_total": len(checks) - 1, "completed_at": checks[-1]["detail"]}, args.output / "INDEPENDENT_VALIDATION_SUMMARY_v1.json")
    print(pd.DataFrame(checks).to_string(index=False))
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
