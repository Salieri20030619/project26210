#!/usr/bin/env python3
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
from scipy.stats import spearmanr, wilcoxon


PRIMARY = ["monocyte_program", "spp1_state_no_spp1", "spp1_vs_c1qc", "SPP1_minus_Monocyte_pseudotime"]
STATES = ["Monocyte", "Other Macrophage", "C1QC-state Macrophage", "SPP1-state Macrophage"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rho(x: pd.Series, y: pd.Series) -> float:
    a, b = x.to_numpy(float), y.to_numpy(float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 20 or np.ptp(a[keep]) == 0 or np.ptp(b[keep]) == 0:
        return math.nan
    return float(spearmanr(a[keep], b[keep]).statistic)


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    v = p[valid]
    order = np.argsort(v)
    q = np.minimum.accumulate((v[order] * len(v) / np.arange(1, len(v) + 1))[::-1])[::-1]
    restored = np.empty(len(v)); restored[order] = np.minimum(q, 1)
    out[valid] = restored
    return out


def signed_rank(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or np.all(np.abs(values) < 1e-15):
        return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox", correction=False).pvalue)


def atomic_tsv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, compression=compression, lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--monocle", type=Path, required=True)
    parser.add_argument("--dpt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    mono = pd.read_csv(args.monocle, sep="\t", dtype={"cell_id": str})
    dpt = pd.read_csv(args.dpt, sep="\t", dtype={"cell_id": str})
    if len(mono) != 50761 or mono["cell_id"].duplicated().any() or set(mono["cell_id"]) != set(dpt["cell_id"]):
        raise RuntimeError("Monocle3 and DPT cell identities differ")
    mono = mono.rename(columns={"study": "study_m3", "patient": "patient_m3", "state_v2": "state_m3"})
    cells = dpt.merge(mono, on="cell_id", how="left", validate="one_to_one")
    if not cells["study_v1"].eq(cells["study_m3"]).all() or not cells["patient"].eq(cells["patient_m3"]).all() or not cells["state_v2"].eq(cells["state_m3"]).all():
        raise RuntimeError("Monocle3 metadata differs from validated DPT metadata")
    if not np.isfinite(cells["monocle3_pseudotime_v1"]).all():
        raise RuntimeError("Non-finite Monocle3 pseudotime")

    patient_rows = []
    for (study, patient), frame in cells.groupby(["study_v1", "patient"], sort=True):
        base_eligible = len(frame) >= 50 and frame["state_v2"].nunique() >= 2 and frame["monocle3_pseudotime_v1"].nunique() > 1
        for endpoint in PRIMARY[:3]:
            patient_rows.append({
                "study": study, "patient": patient, "endpoint": endpoint,
                "n_cells": len(frame), "n_states": frame["state_v2"].nunique(),
                "eligible": base_eligible and frame[endpoint].nunique() > 1,
                "patient_effect": rho(frame["monocle3_pseudotime_v1"], frame[endpoint]) if base_eligible else math.nan,
                "effect_definition": "within-patient Spearman(Monocle3 pseudotime, score)",
            })
        left = frame.loc[frame["state_v2"].eq("Monocyte"), "monocle3_pseudotime_v1"]
        right = frame.loc[frame["state_v2"].eq("SPP1-state Macrophage"), "monocle3_pseudotime_v1"]
        eligible_state = len(left) >= 10 and len(right) >= 10
        patient_rows.append({
            "study": study, "patient": patient, "endpoint": "SPP1_minus_Monocyte_pseudotime",
            "n_cells": len(frame), "n_states": frame["state_v2"].nunique(),
            "eligible": eligible_state,
            "patient_effect": float(right.median() - left.median()) if eligible_state else math.nan,
            "effect_definition": "within-patient median pseudotime(SPP1-state) minus median pseudotime(Monocyte)",
        })
    patient_effects = pd.DataFrame(patient_rows)

    summary_rows = []
    for (study, endpoint), frame in patient_effects.groupby(["study", "endpoint"], sort=True):
        values = frame.loc[frame["eligible"] & frame["patient_effect"].notna(), "patient_effect"].to_numpy(float)
        n = len(values)
        p = signed_rank(values) if n >= 3 else math.nan
        summary_rows.append({
            "study": study, "endpoint": endpoint, "n_eligible_patients": n,
            "n_positive": int((values > 0).sum()), "n_negative": int((values < 0).sum()), "n_zero": int((values == 0).sum()),
            "median_patient_effect": float(np.median(values)) if n else math.nan,
            "wilcoxon_p_two_sided": p,
        })
    summary = pd.DataFrame(summary_rows)
    summary["BH_FDR_within_study_primary4"] = summary.groupby("study", sort=False)["wilcoxon_p_two_sided"].transform(bh)
    summary["reporting_tier"] = np.select(
        [summary["n_eligible_patients"] >= 10, summary["n_eligible_patients"].between(3, 9)],
        ["FDR_AT_LEAST_10", "NOMINAL_P_3_TO_9"], default="DIRECTION_ONLY_LT3",
    )
    summary["statistical_support"] = np.where(
        summary["n_eligible_patients"] >= 10,
        summary["BH_FDR_within_study_primary4"] < 0.05,
        np.where(summary["n_eligible_patients"] >= 3, summary["wilcoxon_p_two_sided"] < 0.05, False),
    )
    summary["direction"] = np.where(summary["median_patient_effect"] > 0, "positive", np.where(summary["median_patient_effect"] < 0, "negative", "zero"))

    concordance_rows = []
    for (study, patient), frame in cells.groupby(["study_v1", "patient"], sort=True):
        eligible = len(frame) >= 50 and frame["state_v2"].nunique() >= 2
        concordance_rows.append({
            "study": study, "patient": patient, "n_cells": len(frame), "n_states": frame["state_v2"].nunique(), "eligible": eligible,
            "monocle3_dpt_spearman": rho(frame["monocle3_pseudotime_v1"], frame["dpt_pseudotime"]) if eligible else math.nan,
        })
    concordance = pd.DataFrame(concordance_rows)
    concordance_summary = []
    for study, frame in concordance.groupby("study", sort=True):
        values = frame.loc[frame["eligible"] & frame["monocle3_dpt_spearman"].notna(), "monocle3_dpt_spearman"].to_numpy(float)
        n = len(values); p = signed_rank(values) if n >= 3 else math.nan
        concordance_summary.append({
            "study": study, "n_eligible_patients": n, "n_positive": int((values > 0).sum()), "n_negative": int((values < 0).sum()),
            "median_patient_rho": float(np.median(values)) if n else math.nan, "wilcoxon_p_two_sided": p,
            "reporting_tier": "FDR_AT_LEAST_10" if n >= 10 else ("NOMINAL_P_3_TO_9" if n >= 3 else "DIRECTION_ONLY_LT3"),
        })
    concordance_summary = pd.DataFrame(concordance_summary)
    concordance_summary["BH_FDR_across_studies"] = bh(concordance_summary["wilcoxon_p_two_sided"])
    concordance_summary["statistical_support"] = np.where(
        concordance_summary["n_eligible_patients"] >= 10,
        concordance_summary["BH_FDR_across_studies"] < 0.05,
        np.where(concordance_summary["n_eligible_patients"] >= 3, concordance_summary["wilcoxon_p_two_sided"] < 0.05, False),
    )

    state_summary = (
        cells.groupby(["study_v1", "state_v2"], observed=True)
        .agg(n_cells=("cell_id", "size"), n_patients=("patient", "nunique"), median_monocle3_pseudotime=("monocle3_pseudotime_v1", "median"), median_dpt_pseudotime=("dpt_pseudotime", "median"))
        .reset_index().rename(columns={"study_v1": "study"})
    )
    cross_study = (
        summary.groupby("endpoint", observed=True)
        .agg(studies=("study", "size"), positive_studies=("median_patient_effect", lambda x: int((x > 0).sum())), negative_studies=("median_patient_effect", lambda x: int((x < 0).sum())), supported_studies=("statistical_support", "sum"), median_study_effect=("median_patient_effect", "median"))
        .reset_index()
    )
    cross_study["interpretation"] = "descriptive synthesis of study-separated Monocle3 graphs"

    atomic_tsv(cells, args.output / "y02_monocle3_dpt_cell_scores_v1.tsv.gz", compression="gzip")
    atomic_tsv(patient_effects, args.output / "y02_monocle3_patient_effects_v1.tsv")
    atomic_tsv(summary, args.output / "y02_monocle3_primary_summary_v1.tsv")
    atomic_tsv(concordance, args.output / "y02_monocle3_dpt_patient_concordance_v1.tsv")
    atomic_tsv(concordance_summary, args.output / "y02_monocle3_dpt_concordance_summary_v1.tsv")
    atomic_tsv(state_summary, args.output / "y02_monocle3_state_summary_v1.tsv")
    atomic_tsv(cross_study, args.output / "y02_monocle3_cross_study_synthesis_v1.tsv")
    status = {
        "status": "PASS_Y02_MONOCLE3_ANALYSIS_COMPLETE_PRESENTATION_PENDING",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_cells": len(cells), "n_studies": int(cells["study_v1"].nunique()), "n_patients": int(cells["patient"].nunique()),
        "direction_prespecified": False, "two_sided_testing": True, "figures_created": False,
        "monocle_input_sha256": sha256(args.monocle), "dpt_input_sha256": sha256(args.dpt),
        "interpretation": "cross-sectional state ordering; not lineage tracing, temporal conversion, RNA velocity or causality",
    }
    part = args.output / "status.json.part"; part.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8"); os.replace(part, args.output / "status.json")


if __name__ == "__main__":
    main()
