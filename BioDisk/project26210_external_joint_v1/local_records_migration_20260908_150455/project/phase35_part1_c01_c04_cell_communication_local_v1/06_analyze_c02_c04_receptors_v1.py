#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


FEATURES = ["MST1R", "ITGAL", "ITGB2", "ITGAL_ITGB2", "CD209", "PRTG", "IGDCC4", "ITGAM_ITGB2_sensitivity"]
PAIRS = {
    "MST1_MST1R": ("MST1", ["MST1R"], "NICHENET_OMNIPATH_SUPPORTED"),
    "ICAM5_ITGAL_ITGB2": ("ICAM5", ["ITGAL", "ITGB2"], "NICHENET_OMNIPATH_SUPPORTED_PRIMARY"),
    "ICAM5_CD209": ("ICAM5", ["CD209"], "NICHENET_OMNIPATH_SUPPORTED"),
    "ICAM5_PRTG": ("ICAM5", ["PRTG"], "NICHENET_OMNIPATH_SUPPORTED"),
    "ICAM5_IGDCC4": ("ICAM5", ["IGDCC4"], "NICHENET_OMNIPATH_SUPPORTED"),
    "ICAM5_ITGAM_ITGB2_sensitivity": ("ICAM5", ["ITGAM", "ITGB2"], "CHECKLIST_RECEPTOR_SENSITIVITY_NOT_DATABASE_PAIR"),
}


def atomic_tsv(frame: pd.DataFrame, path: Path, gzip_output: bool = False) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, compression="gzip" if gzip_output else None, lineterminator="\n")
    os.replace(part, path)


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float); order = np.argsort(p); q = np.empty(len(p)); ranked = p[order] * len(p) / np.arange(1, len(p)+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]; q[order] = np.minimum(ranked, 1); return q


def signed_rank(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if len(values) < 3 or np.all(values == 0): return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox", correction=False).pvalue)


def receptor_features(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["MST1R"] = out["MST1R_log1p_CP10K"]
    out["ITGAL"] = out["ITGAL_log1p_CP10K"]
    out["ITGB2"] = out["ITGB2_log1p_CP10K"]
    out["ITGAL_ITGB2"] = np.sqrt(out["ITGAL_log1p_CP10K"] * out["ITGB2_log1p_CP10K"])
    out["CD209"] = out["CD209_log1p_CP10K"]
    out["PRTG"] = out["PRTG_log1p_CP10K"]
    out["IGDCC4"] = out["IGDCC4_log1p_CP10K"]
    out["ITGAM_ITGB2_sensitivity"] = np.sqrt(out["ITGAM_log1p_CP10K"] * out["ITGB2_log1p_CP10K"])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--input", type=Path, required=True); parser.add_argument("--output", type=Path, required=True); args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    group = pd.read_csv(args.input / "patient_group_targeted_expression_v1.tsv.gz", sep="\t")
    cells = pd.read_csv(args.input / "malignant_targeted_expression_umap_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str})
    cells = receptor_features(cells)

    gp = group.pivot_table(index=["study", "patient", "communication_group", "n_cells"], columns="gene", values="mean_log1p_CP10K").reset_index()
    receptor_rows = []
    pair_rows = []
    for row in gp.itertuples(index=False):
        data = row._asdict()
        features = {
            "MST1R": data["MST1R"], "ITGAL": data["ITGAL"], "ITGB2": data["ITGB2"],
            "ITGAL_ITGB2": math.sqrt(data["ITGAL"] * data["ITGB2"]), "CD209": data["CD209"],
            "PRTG": data["PRTG"], "IGDCC4": data["IGDCC4"],
            "ITGAM_ITGB2_sensitivity": math.sqrt(data["ITGAM"] * data["ITGB2"]),
        }
        for feature, value in features.items():
            receptor_rows.append({"study": data["study"], "patient": data["patient"], "communication_group": data["communication_group"],
                                  "n_cells": data["n_cells"], "receptor_feature": feature, "mean_receptor_availability": value})
    for (study, patient), frame in gp.groupby(["study", "patient"], sort=True):
        by_group = frame.set_index("communication_group")
        for sender in by_group.index:
            for receiver in by_group.index:
                for pair, (ligand, receptors, evidence) in PAIRS.items():
                    ligand_value = float(by_group.loc[sender, ligand])
                    receptor_values = np.array([float(by_group.loc[receiver, gene]) for gene in receptors])
                    receptor_value = float(np.prod(receptor_values) ** (1/len(receptor_values)))
                    pair_rows.append({"study": study, "patient": patient, "sender": sender, "receiver": receiver, "pair": pair,
                                      "evidence_status": evidence, "ligand_mean": ligand_value, "receptor_complex_mean": receptor_value,
                                      "expression_product_score": ligand_value * receptor_value})
    receptor_table = pd.DataFrame(receptor_rows)
    receptor_summary = receptor_table.groupby(["communication_group", "receptor_feature"], sort=True).agg(
        n_patients=("patient", "nunique"), median_receptor_availability=("mean_receptor_availability", "median"),
        n_positive_patients=("mean_receptor_availability", lambda x: int((x > 0).sum()))).reset_index()
    atomic_tsv(receptor_table, args.output / "c02_patient_state_receptor_availability_v1.tsv")
    atomic_tsv(receptor_summary, args.output / "c02_receptor_specificity_summary_v1.tsv")
    atomic_tsv(pd.DataFrame(pair_rows), args.output / "c02_targeted_lr_expression_scores_v1.tsv.gz", gzip_output=True)

    effects = []
    for (study, patient), frame in cells.groupby(["study_v1", "patient"], sort=True):
        for outcome, outcome_col in [("TYMP_continuous", "TYMP_log1p_CP10K"), ("F1_continuous", "F1")]:
            for feature in FEATURES:
                x = frame[feature].to_numpy(float); y = frame[outcome_col].to_numpy(float); keep = np.isfinite(x) & np.isfinite(y)
                eligible = keep.sum() >= 20 and np.ptp(x[keep]) > 0 and np.ptp(y[keep]) > 0
                value = float(spearmanr(x[keep], y[keep]).statistic) if eligible else math.nan
                effects.append({"study": study, "patient": patient, "outcome": outcome, "receptor_feature": feature,
                                "n_cells": len(frame), "eligible": eligible, "patient_spearman_rho": value})
    effects = pd.DataFrame(effects)
    summaries = []
    for (outcome, feature), frame in effects.groupby(["outcome", "receptor_feature"], sort=True):
        values = frame.loc[frame["eligible"], "patient_spearman_rho"].to_numpy(float); p = signed_rank(values)
        summaries.append({"outcome": outcome, "receptor_feature": feature, "n_eligible_patients": len(values),
                          "n_positive": int((values > 0).sum()), "n_negative": int((values < 0).sum()), "n_zero": int((values == 0).sum()),
                          "median_patient_rho": float(np.median(values)) if len(values) else math.nan, "wilcoxon_p_two_sided": p})
    summaries = pd.DataFrame(summaries)
    summaries["BH_FDR_within_outcome8features"] = summaries.groupby("outcome", sort=False)["wilcoxon_p_two_sided"].transform(lambda x: bh(x))
    summaries["statistical_support"] = np.where(summaries["n_eligible_patients"] >= 10, summaries["BH_FDR_within_outcome8features"] < 0.05,
        np.where(summaries["n_eligible_patients"] >= 3, summaries["wilcoxon_p_two_sided"] < 0.05, False))
    summaries["direction"] = np.where(summaries["median_patient_rho"] > 0, "positive", np.where(summaries["median_patient_rho"] < 0, "negative", "zero"))
    atomic_tsv(effects, args.output / "c04_patient_receptor_tymp_f1_effects_v1.tsv")
    atomic_tsv(summaries, args.output / "c04_receptor_tymp_f1_summary_v1.tsv")
    umap_cols = ["exact_cell_id", "study_v1", "patient", "UMAP1", "UMAP2", "TYMP_log1p_CP10K", "F1"] + FEATURES
    atomic_tsv(cells[umap_cols], args.output / "c04_receptor_tymp_f1_umap_source_v1.tsv.gz", gzip_output=True)
    status = {"status":"PASS_C02_C04_ANALYSIS_COMPLETE_PRESENTATION_PENDING", "completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),
              "n_malignant_cells":len(cells), "n_patients":int(cells["patient"].nunique()), "two_sided_testing":True,
              "direction_prespecified":False, "figures_created":False,
              "interpretation":"RNA receptor availability and expression-product scores; not measured protein binding, circulating ligand level or causality"}
    part=args.output/"status.json.part";part.write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8");os.replace(part,args.output/"status.json")


if __name__ == "__main__": main()
