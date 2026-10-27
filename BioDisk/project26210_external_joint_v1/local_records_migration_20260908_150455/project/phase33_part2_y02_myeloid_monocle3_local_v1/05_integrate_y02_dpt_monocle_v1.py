#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


ROOT = Path(__file__).resolve().parent
CELL = ROOT / "analysis_local_v1" / "y02_monocle3_dpt_cell_scores_v1.tsv.gz"
MONO = ROOT / "analysis_local_v1" / "y02_monocle3_primary_summary_v1.tsv"
OUT = ROOT / "integration_local_v1"
ENDPOINTS = ["monocyte_program", "spp1_state_no_spp1", "spp1_vs_c1qc", "SPP1_minus_Monocyte_pseudotime"]


def rho(x: pd.Series, y: pd.Series) -> float:
    a, b = x.to_numpy(float), y.to_numpy(float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 20 or np.ptp(a[keep]) == 0 or np.ptp(b[keep]) == 0:
        return math.nan
    return float(spearmanr(a[keep], b[keep]).statistic)


def bh(series: pd.Series) -> np.ndarray:
    p = series.to_numpy(float); out = np.full(len(p), np.nan); valid = np.isfinite(p)
    if not valid.any(): return out
    v = p[valid]; order = np.argsort(v); q = np.minimum.accumulate((v[order] * len(v) / np.arange(1, len(v)+1))[::-1])[::-1]
    restored = np.empty(len(v)); restored[order] = np.minimum(q, 1); out[valid] = restored; return out


def signed_rank(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if not len(values) or np.all(np.abs(values) < 1e-15): return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox", correction=False).pvalue)


def main() -> None:
    if OUT.exists(): raise FileExistsError(OUT)
    OUT.mkdir()
    cells = pd.read_csv(CELL, sep="\t")
    rows = []
    for (study, patient), frame in cells.groupby(["study_v1", "patient"], sort=True):
        base = len(frame) >= 50 and frame["state_v2"].nunique() >= 2 and frame["dpt_pseudotime"].nunique() > 1
        for endpoint in ENDPOINTS[:3]:
            rows.append({"study": study, "patient": patient, "endpoint": endpoint, "eligible": base and frame[endpoint].nunique() > 1, "patient_effect": rho(frame["dpt_pseudotime"], frame[endpoint]) if base else math.nan})
        left = frame.loc[frame["state_v2"].eq("Monocyte"), "dpt_pseudotime"]
        right = frame.loc[frame["state_v2"].eq("SPP1-state Macrophage"), "dpt_pseudotime"]
        eligible = len(left) >= 10 and len(right) >= 10
        rows.append({"study": study, "patient": patient, "endpoint": ENDPOINTS[3], "eligible": eligible, "patient_effect": float(right.median()-left.median()) if eligible else math.nan})
    effects = pd.DataFrame(rows)
    summary_rows = []
    for (study, endpoint), frame in effects.groupby(["study", "endpoint"], sort=True):
        values = frame.loc[frame["eligible"] & frame["patient_effect"].notna(), "patient_effect"].to_numpy(float)
        summary_rows.append({"study": study, "endpoint": endpoint, "n_eligible_patients": len(values), "n_positive": int((values>0).sum()), "n_negative": int((values<0).sum()), "median_patient_effect": float(np.median(values)) if len(values) else math.nan, "wilcoxon_p_two_sided": signed_rank(values) if len(values)>=3 else math.nan})
    dpt = pd.DataFrame(summary_rows)
    dpt["BH_FDR_within_study_primary4"] = dpt.groupby("study", sort=False)["wilcoxon_p_two_sided"].transform(bh)
    dpt["statistical_support"] = np.where(dpt["n_eligible_patients"]>=10, dpt["BH_FDR_within_study_primary4"]<.05, np.where(dpt["n_eligible_patients"]>=3, dpt["wilcoxon_p_two_sided"]<.05, False))
    dpt["direction"] = np.where(dpt["median_patient_effect"]>0, "positive", np.where(dpt["median_patient_effect"]<0, "negative", "zero"))
    dpt.to_csv(OUT / "y02_dpt_reclassified_summary_v1.tsv", sep="\t", index=False, lineterminator="\n")

    mono = pd.read_csv(MONO, sep="\t")
    compare = mono.merge(dpt, on=["study", "endpoint"], suffixes=("_monocle3", "_dpt"), validate="one_to_one")
    compare["same_direction"] = compare["direction_monocle3"].eq(compare["direction_dpt"])
    compare["both_supported_same_direction"] = compare["statistical_support_monocle3"].astype(bool) & compare["statistical_support_dpt"].astype(bool) & compare["same_direction"]
    compare["comparison"] = np.where(compare["same_direction"], np.where(compare["both_supported_same_direction"], "SAME_DIRECTION_BOTH_SUPPORTED", "SAME_DIRECTION_LIMITED_SUPPORT"), "OPPOSITE_DIRECTION")
    compare.to_csv(OUT / "y02_dpt_monocle3_endpoint_comparison_v1.tsv", sep="\t", index=False, lineterminator="\n")

    overall = []
    for endpoint, frame in compare.groupby("endpoint", sort=False):
        mpos = int(frame["median_patient_effect_monocle3"].gt(0).sum()); dpos = int(frame["median_patient_effect_dpt"].gt(0).sum()); same = int(frame["same_direction"].sum()); opposite = int((~frame["same_direction"]).sum())
        if same >= 4 and (mpos >= 4 and dpos >= 4): conclusion = "CROSS_METHOD_POSITIVE_DIRECTION"
        elif same >= 4 and (mpos <= 1 and dpos <= 1): conclusion = "CROSS_METHOD_NEGATIVE_DIRECTION"
        elif opposite >= 3: conclusion = "CROSS_METHOD_DIRECTION_CONFLICT"
        else: conclusion = "HETEROGENEOUS_OR_LIMITED"
        overall.append({"endpoint": endpoint, "n_studies": len(frame), "same_direction_studies": same, "opposite_direction_studies": opposite, "monocle3_positive_studies": mpos, "dpt_positive_studies": dpos, "both_supported_same_direction_studies": int(frame["both_supported_same_direction"].sum()), "cross_method_conclusion": conclusion})
    overall = pd.DataFrame(overall)
    overall.to_csv(OUT / "y02_cross_method_endpoint_conclusions_v1.tsv", sep="\t", index=False, lineterminator="\n")
    status = {"status":"PASS_Y02_CROSS_METHOD_INTEGRATION", "study_endpoint_rows":len(compare), "endpoints":len(overall), "direction_prespecified":False, "figures_created":False}
    (OUT / "status.json").write_text(json.dumps(status, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__": main()
