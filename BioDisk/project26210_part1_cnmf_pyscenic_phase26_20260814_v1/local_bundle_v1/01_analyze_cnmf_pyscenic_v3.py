#!/usr/bin/env python3
"""Patient-level cNMF group and Factor-TYMP analysis of validated pySCENIC AUCell matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    valid = values.notna()
    result = pd.Series(np.nan, index=values.index, dtype=float)
    if not valid.any():
        return result
    x = values.loc[valid].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty_like(ranked)
    restored[order] = np.clip(ranked, 0, 1)
    result.loc[valid] = restored
    return result


def tf_name(regulon: str) -> str:
    return re.sub(r"\([+-]\)$", "", str(regulon))


def summarize(details: pd.DataFrame, keys: list[str], effect: str) -> pd.DataFrame:
    rows = []
    for key, frame in details.groupby(keys, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        values = frame[effect].dropna().to_numpy(float)
        row = dict(zip(keys, key))
        full_statistics = len(values) >= 10
        raw_median = float(np.median(values)) if len(values) else np.nan
        p_value = np.nan
        if len(values) >= 3:
            if np.all(values == 0):
                p_value = 1.0
            else:
                p_value = float(
                    wilcoxon(values, zero_method="wilcox", alternative="two-sided", method="auto").pvalue
                )
        small_n_significant = 3 <= len(values) < 10 and np.isfinite(p_value) and p_value < 0.05
        row.update(
            {
                "n_eligible_patients": int(len(values)),
                "n_positive": int((values > 0).sum()),
                "n_negative": int((values < 0).sum()),
                "n_zero": int((values == 0).sum()),
                "median_patient_effect_internal": raw_median,
                "median_patient_effect": raw_median if (full_statistics or small_n_significant) else np.nan,
                "reporting_tier": (
                    "full_statistics_fdr" if full_statistics
                    else "small_n_nominal_significant" if small_n_significant
                    else "direction_only"
                ),
                "wilcoxon_p_two_sided": p_value,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--genie3-auc", required=True, type=Path)
    parser.add_argument("--grnboost2-auc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    parser.add_argument("--min-cells-correlation", type=int, default=20)
    parser.add_argument("--min-cells-per-arm", type=int, default=10)
    args = parser.parse_args()

    for path in (args.metadata, args.genie3_auc, args.grnboost2_auc):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(f"refusing stale partial output: {part}")
    part.mkdir(parents=True)

    meta = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    factor_cols = [f"F{i}" for i in range(1, 8)]
    required = {"exact_cell_id", "study", "curated_patient_id", "dominant_factor_rank7_v1", *factor_cols}
    if not required.issubset(meta.columns):
        raise ValueError(f"metadata columns missing: {sorted(required - set(meta.columns))}")
    if len(meta) != args.expected_cells or meta["exact_cell_id"].duplicated().any():
        raise ValueError("metadata cell count/uniqueness gate failed")
    ids = pd.Index(meta["exact_cell_id"].astype(str))
    factors = meta[factor_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not np.isfinite(factors).all() or (factors < 0).any():
        raise ValueError("invalid cNMF usage values")
    reproduced = np.array(factor_cols)[np.argmax(factors, axis=1)]
    if not np.array_equal(reproduced, meta["dominant_factor_rank7_v1"].astype(str).to_numpy()):
        raise ValueError("dominant cNMF factor labels are not reproducible")
    if meta["study"].nunique() != 5 or meta["curated_patient_id"].nunique() != 26:
        raise ValueError("cohort traceability gate failed")
    meta = meta.set_index(ids, drop=False)

    input_manifest = []
    for role, path in (("metadata", args.metadata), ("genie3_auc", args.genie3_auc), ("grnboost2_auc", args.grnboost2_auc)):
        input_manifest.append({"role": role, "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(input_manifest).to_csv(part / "INPUT_MANIFEST_v1.tsv", sep="\t", index=False)

    continuous_details = []
    group_details = []
    availability = []
    method_summaries = []
    matrices = {"GENIE3": args.genie3_auc, "GRNBoost2": args.grnboost2_auc}
    for method, path in matrices.items():
        auc = pd.read_csv(path, index_col=0)
        auc.index = auc.index.astype(str)
        if len(auc) != args.expected_cells or not auc.index.is_unique or not auc.index.equals(ids):
            raise ValueError(f"{method} AUCell cell identity/order gate failed")
        values = auc.to_numpy(float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"{method} AUCell numeric safety gate failed")
        if (np.ptp(values, axis=0) == 0).any():
            raise ValueError(f"{method} contains constant regulons")
        regulons = auc.columns.astype(str).tolist()
        tfs = [tf_name(x) for x in regulons]
        availability.append(
            {"method": method, "target_tf": "HIF1A", "regulon_present": "HIF1A" in set(tfs)}
        )
        method_summaries.append({"method": method, "n_cells": len(auc), "n_regulons": len(regulons)})

        joined = meta[["study", "curated_patient_id", "F1", "dominant_factor_rank7_v1"]].copy()
        for (study, patient), idx in joined.groupby(["study", "curated_patient_id"], sort=True).groups.items():
            loc = ids.get_indexer(pd.Index(idx).astype(str))
            frame_auc = auc.iloc[loc]
            f1 = joined.loc[idx, "F1"].to_numpy(float)
            if len(idx) >= args.min_cells_correlation and np.ptp(f1) > 0:
                rho = frame_auc.corrwith(pd.Series(f1, index=frame_auc.index), method="spearman")
            else:
                rho = pd.Series(np.nan, index=regulons)
            for regulon, tf, value in zip(regulons, tfs, rho.to_numpy(float)):
                continuous_details.append(
                    {
                        "method": method,
                        "regulon": regulon,
                        "tf": tf,
                        "study": study,
                        "curated_patient_id": patient,
                        "n_cells": int(len(idx)),
                        "eligible": bool(np.isfinite(value)),
                        "spearman_rho": value,
                    }
                )

            labels = joined.loc[idx, "dominant_factor_rank7_v1"].astype(str).to_numpy()
            for factor in factor_cols:
                mask = labels == factor
                n_factor = int(mask.sum())
                n_other = int((~mask).sum())
                eligible = n_factor >= args.min_cells_per_arm and n_other >= args.min_cells_per_arm
                diff = (
                    np.median(frame_auc.to_numpy(float)[mask, :], axis=0)
                    - np.median(frame_auc.to_numpy(float)[~mask, :], axis=0)
                    if eligible else np.full(len(regulons), np.nan)
                )
                for regulon, tf, value in zip(regulons, tfs, diff):
                    group_details.append(
                        {
                            "method": method,
                            "factor": factor,
                            "regulon": regulon,
                            "tf": tf,
                            "study": study,
                            "curated_patient_id": patient,
                            "n_factor_cells": n_factor,
                            "n_other_cells": n_other,
                            "eligible": eligible,
                            "median_auc_difference": value,
                        }
                    )

    continuous_details = pd.DataFrame(continuous_details)
    group_details = pd.DataFrame(group_details)
    continuous = summarize(continuous_details, ["method", "regulon", "tf"], "spearman_rho")
    continuous["BH_FDR"] = continuous.groupby("method")["wilcoxon_p_two_sided"].transform(bh)
    continuous["endpoint_pass"] = (
        (continuous["median_patient_effect_internal"] > 0)
        & (
            ((continuous["n_eligible_patients"] >= 10) & (continuous["BH_FDR"] < 0.05))
            | ((continuous["n_eligible_patients"] < 10) & (continuous["wilcoxon_p_two_sided"] < 0.05))
        )
    )
    grouped = summarize(group_details, ["method", "factor", "regulon", "tf"], "median_auc_difference")
    grouped["BH_FDR"] = grouped.groupby(["method", "factor"])["wilcoxon_p_two_sided"].transform(bh)
    grouped["endpoint_pass"] = (
        (grouped["median_patient_effect_internal"] > 0)
        & (
            ((grouped["n_eligible_patients"] >= 10) & (grouped["BH_FDR"] < 0.05))
            | ((grouped["n_eligible_patients"] < 10) & (grouped["wilcoxon_p_two_sided"] < 0.05))
        )
    )

    pivot_values = ["median_patient_effect_internal", "median_patient_effect", "wilcoxon_p_two_sided", "BH_FDR", "n_eligible_patients", "reporting_tier", "endpoint_pass"]
    cont_wide = continuous.pivot(index="tf", columns="method", values=pivot_values)
    grp_f1 = grouped.loc[grouped["factor"] == "F1"].pivot(index="tf", columns="method", values=pivot_values)
    cont_wide.columns = [f"{metric}__{method}__continuous" for metric, method in cont_wide.columns]
    grp_f1.columns = [f"{metric}__{method}__f1_group" for metric, method in grp_f1.columns]
    shared = cont_wide.join(grp_f1, how="inner")
    presence_cols = [
        "n_eligible_patients__GENIE3__continuous",
        "n_eligible_patients__GRNBoost2__continuous",
        "n_eligible_patients__GENIE3__f1_group",
        "n_eligible_patients__GRNBoost2__f1_group",
    ]
    shared = shared.dropna(subset=presence_cols)
    shared = shared.reset_index()
    effect_cols = [
        "median_patient_effect_internal__GENIE3__continuous",
        "median_patient_effect_internal__GRNBoost2__continuous",
        "median_patient_effect_internal__GENIE3__f1_group",
        "median_patient_effect_internal__GRNBoost2__f1_group",
    ]
    fdr_cols = [
        "BH_FDR__GENIE3__continuous",
        "BH_FDR__GRNBoost2__continuous",
        "BH_FDR__GENIE3__f1_group",
        "BH_FDR__GRNBoost2__f1_group",
    ]
    for col in effect_cols:
        shared[col + "__percentile_rank"] = shared[col].rank(pct=True)
    shared["consensus_effect_percentile"] = shared[[c + "__percentile_rank" for c in effect_cols]].mean(axis=1)
    shared["all_four_effects_positive"] = (shared[effect_cols] > 0).all(axis=1)
    shared["all_four_fdr_below_0_05"] = (shared[fdr_cols] < 0.05).all(axis=1)
    pass_cols = [
        "endpoint_pass__GENIE3__continuous",
        "endpoint_pass__GRNBoost2__continuous",
        "endpoint_pass__GENIE3__f1_group",
        "endpoint_pass__GRNBoost2__f1_group",
    ]
    shared["n_endpoints_passed"] = shared[pass_cols].astype(bool).sum(axis=1)
    shared["high_confidence_f1_driver"] = shared["n_endpoints_passed"] == 4
    shared["secondary_support_f1_driver"] = shared["n_endpoints_passed"] == 3
    shared["evidence_tier"] = np.select(
        [
            shared["n_endpoints_passed"] == 4,
            shared["n_endpoints_passed"] == 3,
            shared["n_endpoints_passed"] == 2,
        ],
        ["high_confidence_4_of_4", "secondary_support_3_of_4", "exploratory_2_of_4"],
        default="insufficient_0_or_1_of_4",
    )
    shared = shared.sort_values(
        ["n_endpoints_passed", "consensus_effect_percentile", "tf"],
        ascending=[False, False, True],
    )
    shared.insert(0, "consensus_rank", np.arange(1, len(shared) + 1))

    continuous_details.to_csv(part / "patient_continuous_regulon_effects_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    group_details.to_csv(part / "patient_cnmf_group_regulon_effects_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    continuous.to_csv(part / "continuous_regulon_patient_summary_v1.tsv", sep="\t", index=False, na_rep="NA")
    grouped.to_csv(part / "cnmf_group_regulon_patient_summary_v1.tsv", sep="\t", index=False, na_rep="NA")
    shared.to_csv(part / "f1_driver_consensus_v1.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(availability).to_csv(part / "prespecified_tf_availability_v1.tsv", sep="\t", index=False)
    pd.DataFrame(method_summaries).to_csv(part / "method_input_summary_v1.tsv", sep="\t", index=False)

    status = {
        "status": "PASS_ANALYSIS_COMPLETE",
        "n_cells": int(len(meta)),
        "n_patients": int(meta["curated_patient_id"].nunique()),
        "n_studies": int(meta["study"].nunique()),
        "genie3_regulons": int(method_summaries[0]["n_regulons"]),
        "grnboost2_regulons": int(method_summaries[1]["n_regulons"]),
        "shared_tf_regulons": int(len(shared)),
        "high_confidence_f1_drivers": int(shared["high_confidence_f1_driver"].sum()),
        "secondary_support_f1_drivers": int(shared["secondary_support_f1_driver"].sum()),
        "exploratory_two_of_four_f1_drivers": int((shared["n_endpoints_passed"] == 2).sum()),
        "hif1a_regulon_available": False,
        "patient_unit_primary": True,
        "cell_level_inference": False,
        "python": platform.python_version(),
    }
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    outputs = []
    for path in sorted(p for p in part.iterdir() if p.is_file()):
        outputs.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(outputs).to_csv(part / "OUTPUT_MANIFEST_v1.tsv", sep="\t", index=False)
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
