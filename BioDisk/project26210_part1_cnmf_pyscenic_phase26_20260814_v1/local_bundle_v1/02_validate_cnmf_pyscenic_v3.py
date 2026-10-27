#!/usr/bin/env python3
"""Independent terminal validation for Phase 26 cNMF-pySCENIC analysis."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(series: pd.Series) -> pd.Series:
    valid = series.notna()
    output = pd.Series(np.nan, index=series.index, dtype=float)
    if not valid.any():
        return output
    x = series.loc[valid].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(ranked)
    result[order] = np.clip(ranked, 0, 1)
    output.loc[valid] = result
    return output


def tf_name(value: str) -> str:
    return re.sub(r"\([+-]\)$", "", str(value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--genie3-auc", required=True, type=Path)
    parser.add_argument("--grnboost2-auc", required=True, type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    args = parser.parse_args()

    if args.output.exists():
        raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(part)
    part.mkdir(parents=True)
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    status = json.loads((args.analysis_dir / "status.json").read_text(encoding="utf-8"))
    check("analysis_terminal_status", status.get("status") == "PASS_ANALYSIS_COMPLETE", str(status.get("status")))

    meta = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    ids = pd.Index(meta["exact_cell_id"].astype(str))
    check("metadata_identity", len(ids) == args.expected_cells and ids.is_unique, f"n={len(ids)}; unique={ids.is_unique}")
    check("cohort_traceability", meta["study"].nunique() == 5 and meta["curated_patient_id"].nunique() == 26, f"studies={meta['study'].nunique()}; patients={meta['curated_patient_id'].nunique()}")
    factors = meta[[f"F{i}" for i in range(1, 8)]].to_numpy(float)
    dominant = np.array([f"F{i}" for i in range(1, 8)])[np.argmax(factors, axis=1)]
    check("cnmf_factor_safety", np.isfinite(factors).all() and (factors >= 0).all(), "F1-F7 finite and nonnegative")
    check("dominant_factor_reproduced", np.array_equal(dominant, meta["dominant_factor_rank7_v1"].astype(str)), "argmax F1-F7")
    meta = meta.set_index(ids, drop=False)

    input_manifest = pd.read_csv(args.analysis_dir / "INPUT_MANIFEST_v1.tsv", sep="\t")
    expected_sources = {"metadata": args.metadata, "genie3_auc": args.genie3_auc, "grnboost2_auc": args.grnboost2_auc}
    manifest_ok = True
    for role, path in expected_sources.items():
        row = input_manifest.loc[input_manifest["role"] == role]
        manifest_ok &= len(row) == 1 and row.iloc[0]["sha256"] == sha256(path) and int(row.iloc[0]["bytes"]) == path.stat().st_size
    check("input_manifest_hashes", manifest_ok, "three source hashes and sizes")

    matrices = {}
    auc_ok = True
    for method, path, n_reg in (
        ("GENIE3", args.genie3_auc, 299),
        ("GRNBoost2", args.grnboost2_auc, 314),
    ):
        auc = pd.read_csv(path, index_col=0)
        auc.index = auc.index.astype(str)
        matrices[method] = auc
        values = auc.to_numpy(float)
        auc_ok &= (
            auc.shape == (args.expected_cells, n_reg)
            and auc.index.equals(ids)
            and np.isfinite(values).all()
            and (values >= 0).all()
            and (np.ptp(values, axis=0) > 0).all()
        )
    check("aucell_identity_shape_numeric", auc_ok, "GENIE3=11950x299; GRNBoost2=11950x314")
    hif1a_absent = all("HIF1A" not in {tf_name(x) for x in auc.columns} for auc in matrices.values())
    availability = pd.read_csv(args.analysis_dir / "prespecified_tf_availability_v1.tsv", sep="\t")
    check("hif1a_unavailable_reproduced", hif1a_absent and not availability["regulon_present"].astype(bool).any(), "absent in both methods")

    cont_detail = pd.read_csv(args.analysis_dir / "patient_continuous_regulon_effects_v1.tsv.gz", sep="\t")
    grp_detail = pd.read_csv(args.analysis_dir / "patient_cnmf_group_regulon_effects_v1.tsv.gz", sep="\t")
    cont = pd.read_csv(args.analysis_dir / "continuous_regulon_patient_summary_v1.tsv", sep="\t", na_values="NA")
    grp = pd.read_csv(args.analysis_dir / "cnmf_group_regulon_patient_summary_v1.tsv", sep="\t", na_values="NA")
    consensus = pd.read_csv(args.analysis_dir / "f1_driver_consensus_v1.tsv", sep="\t", na_values="NA")
    check("output_row_counts", len(cont_detail) == 26 * 613 and len(grp_detail) == 26 * 7 * 613 and len(cont) == 613 and len(grp) == 7 * 613, f"{len(cont_detail)},{len(grp_detail)},{len(cont)},{len(grp)}")
    check("summary_unique_keys", not cont.duplicated(["method", "regulon"]).any() and not grp.duplicated(["method", "factor", "regulon"]).any(), "continuous and grouped unique")

    tier_ok = True
    for frame in (cont, grp):
        full = frame["n_eligible_patients"] >= 10
        small_sig = (
            (frame["n_eligible_patients"] >= 3)
            & (frame["n_eligible_patients"] < 10)
            & (frame["wilcoxon_p_two_sided"] < 0.05)
        )
        direction = ~(full | small_sig)
        tier_ok &= (frame.loc[full, "reporting_tier"] == "full_statistics_fdr").all()
        tier_ok &= frame.loc[full, ["median_patient_effect", "wilcoxon_p_two_sided", "BH_FDR"]].notna().all().all()
        tier_ok &= (frame.loc[small_sig, "reporting_tier"] == "small_n_nominal_significant").all()
        tier_ok &= frame.loc[small_sig, ["median_patient_effect", "wilcoxon_p_two_sided"]].notna().all().all()
        tier_ok &= (frame.loc[direction, "reporting_tier"] == "direction_only").all()
        tier_ok &= frame.loc[direction, "median_patient_effect"].isna().all()
        has_effect = frame["n_eligible_patients"] > 0
        tier_ok &= frame.loc[has_effect, "median_patient_effect_internal"].notna().all()
        tier_ok &= frame.loc[~has_effect, "median_patient_effect_internal"].isna().all()
    check("patient_reporting_rule", tier_ok, "<10 nominal P<0.05 reportable; otherwise direction only; >=10 FDR tier")

    agg_ok = True
    for _, row in cont.iterrows():
        values = cont_detail.loc[
            (cont_detail["method"] == row["method"]) & (cont_detail["regulon"] == row["regulon"]),
            "spearman_rho",
        ].dropna().to_numpy(float)
        agg_ok &= len(values) == int(row["n_eligible_patients"])
        agg_ok &= int((values > 0).sum()) == int(row["n_positive"])
        agg_ok &= int((values < 0).sum()) == int(row["n_negative"])
        agg_ok &= np.isclose(np.median(values), row["median_patient_effect_internal"], rtol=0, atol=1e-12)
    check("continuous_patient_aggregation", agg_ok, "all 613 regulons")

    grp_agg_ok = True
    for _, row in grp.iterrows():
        values = grp_detail.loc[
            (grp_detail["method"] == row["method"])
            & (grp_detail["factor"] == row["factor"])
            & (grp_detail["regulon"] == row["regulon"]),
            "median_auc_difference",
        ].dropna().to_numpy(float)
        grp_agg_ok &= len(values) == int(row["n_eligible_patients"])
        grp_agg_ok &= int((values > 0).sum()) == int(row["n_positive"])
        grp_agg_ok &= int((values < 0).sum()) == int(row["n_negative"])
        if len(values):
            grp_agg_ok &= np.isclose(np.median(values), row["median_patient_effect_internal"], rtol=0, atol=1e-12)
        else:
            grp_agg_ok &= pd.isna(row["median_patient_effect_internal"])
    check("group_patient_aggregation", grp_agg_ok, "all 4,291 method-factor-regulon tests")

    fdr_ok = True
    for _, frame in cont.groupby("method"):
        fdr_ok &= np.allclose(bh(frame["wilcoxon_p_two_sided"]), frame["BH_FDR"], equal_nan=True, rtol=0, atol=1e-12)
    for _, frame in grp.groupby(["method", "factor"]):
        fdr_ok &= np.allclose(bh(frame["wilcoxon_p_two_sided"]), frame["BH_FDR"], equal_nan=True, rtol=0, atol=1e-12)
    check("bh_fdr_reproduced", fdr_ok, "2 continuous and 14 factor families")

    endpoint_ok = True
    for frame in (cont, grp):
        expected_pass = (
            (frame["median_patient_effect_internal"] > 0)
            & (
                ((frame["n_eligible_patients"] >= 10) & (frame["BH_FDR"] < 0.05))
                | ((frame["n_eligible_patients"] < 10) & (frame["wilcoxon_p_two_sided"] < 0.05))
            )
        )
        endpoint_ok &= np.array_equal(expected_pass.to_numpy(bool), frame["endpoint_pass"].astype(bool).to_numpy())
    check("endpoint_pass_rule", endpoint_ok, "FDR for >=10; nominal P for <10")

    sentinel_ok = True
    for method, auc in matrices.items():
        regulon = str(auc.columns[0])
        for (study, patient), idx in meta.groupby(["study", "curated_patient_id"], sort=True).groups.items():
            f1 = meta.loc[idx, "F1"].to_numpy(float)
            activity = auc.loc[pd.Index(idx).astype(str), regulon].to_numpy(float)
            expected = float(spearmanr(f1, activity).statistic) if len(idx) >= 20 and np.ptp(f1) > 0 and np.ptp(activity) > 0 else np.nan
            observed = cont_detail.loc[
                (cont_detail["method"] == method)
                & (cont_detail["regulon"] == regulon)
                & (cont_detail["study"] == study)
                & (cont_detail["curated_patient_id"] == patient),
                "spearman_rho",
            ].iloc[0]
            sentinel_ok &= np.isclose(expected, observed, equal_nan=True, rtol=0, atol=1e-12)
    check("raw_auc_sentinel_correlations", sentinel_ok, "first regulon per method across 26 patients")

    pass_cols = [x for x in consensus.columns if x.startswith("endpoint_pass__")]
    n_pass = consensus[pass_cols].astype(bool).sum(axis=1)
    high = n_pass == 4
    secondary = n_pass == 3
    rule_ok = (
        np.array_equal(n_pass.to_numpy(int), consensus["n_endpoints_passed"].to_numpy(int))
        and np.array_equal(high.to_numpy(bool), consensus["high_confidence_f1_driver"].astype(bool).to_numpy())
        and np.array_equal(secondary.to_numpy(bool), consensus["secondary_support_f1_driver"].astype(bool).to_numpy())
    )
    check("consensus_driver_rule", rule_ok, f"high={int(high.sum())}; secondary={int(secondary.sum())}")
    check("consensus_shared_count", len(consensus) == 197, f"shared={len(consensus)}")

    output_manifest = pd.read_csv(args.analysis_dir / "OUTPUT_MANIFEST_v1.tsv", sep="\t")
    manifest_ok = True
    for _, row in output_manifest.iterrows():
        path = args.analysis_dir / row["file"]
        manifest_ok &= path.is_file() and path.stat().st_size == int(row["bytes"]) and sha256(path) == row["sha256"]
    check("output_manifest_hashes", manifest_ok, f"files={len(output_manifest)}")
    check(
        "status_counts_reproduced",
        status["high_confidence_f1_drivers"] == int(high.sum())
        and status["secondary_support_f1_drivers"] == int(secondary.sum())
        and status["shared_tf_regulons"] == len(consensus),
        "status versus consensus",
    )

    table = pd.DataFrame(checks)
    table.to_csv(part / "INDEPENDENT_VALIDATION_CHECKS_v1.tsv", sep="\t", index=False)
    passed = int(table["passed"].sum())
    total = int(len(table))
    terminal = "PASS_INDEPENDENT_VALIDATION" if passed == total else "FAIL_INDEPENDENT_VALIDATION"
    result = {"status": terminal, "checks_passed": passed, "checks_total": total}
    (part / "status.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(
        [{"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(part.iterdir()) if p.is_file()]
    ).to_csv(part / "VALIDATION_MANIFEST_v1.tsv", sep="\t", index=False)
    os.replace(part, args.output)
    if terminal != "PASS_INDEPENDENT_VALIDATION":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
