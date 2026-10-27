#!/usr/bin/env python3
"""Method-specific validation of supportive malignant programs for key relations."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
PROGRAMS = HERE / "posthoc_malignancy_evidence_tier_programs_v4" / "cross_sectional_patient_evidence_tier_programs_v1.tsv"
EXPANDED_RUN = HERE / "runs" / "20260812_expanded_malignancy_v2"
OUT = HERE / "posthoc_method_specific_supportive_validation_v1"
DATASETS = ["GSE178481", "GSE202374", "GSE207493", "GSE264586", "GSE294109"]
METHODS = [("primary_gene_z_mean", "mean_primary"), ("sensitivity_mean_log_normalized", "mean_sensitivity")]
SUPPORT_GROUPS = [
    "supportive_part1_infercnv_only", "supportive_part1_copykat_only",
    "supportive_external_infercnv_only", "supportive_external_copykat_only",
]
RELATIONS = [
    ("H2_F1_to_C1QC", "f1_top100_weighted", "c1qc_state", "myeloid"),
    ("H2_F1_to_APOE_C1QC", "f1_top100_weighted", "apoe_tam_expression_proxy", "myeloid"),
    ("H3_F1_to_hypoxia", "f1_top100_weighted", "hypoxia_angiogenesis", "malignant"),
    ("H1_TYMP_to_SPP1", "tymp_expression", "spp1_state_no_spp1", "myeloid"),
]


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def nonmalignant(dataset: str) -> pd.DataFrame:
    frame = pd.read_csv(EXPANDED_RUN / dataset / "tables" / "sample_program_summary_v1.tsv", sep="\t")
    frame = frame[frame["analysis_compartment"] != "malignant_provisional"].copy()
    if dataset == "GSE178481":
        frame = frame[~frame["sample_id"].astype(str).str.contains("Normal", case=False, regex=False)]
    if dataset == "GSE264586":
        frame = frame[frame["timepoint"] == "Baseline"]
    rows = []
    keys = ["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"]
    for group_keys, group in frame.groupby(keys, sort=True, dropna=False):
        w = group["n_cells"].to_numpy(float)
        rows.append(dict(zip(keys, group_keys), n_cells=int(w.sum()),
                         mean_primary=float(np.average(group["mean_primary"], weights=w)),
                         mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=w))))
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir()
    programs = pd.read_csv(PROGRAMS, sep="\t")
    rows = []
    for dataset in DATASETS:
        non = nonmalignant(dataset)
        for method, value in METHODS:
            n = non[non["n_cells"] >= 20].copy()
            n["key"] = n["signature_id"] + "__" + n["analysis_compartment"]
            nw = n.pivot(index=["dataset", "patient_key", "patient_id"], columns="key", values=value).reset_index()
            for group_name in SUPPORT_GROUPS:
                m = programs[(programs["dataset"] == dataset) & (programs["evidence_tier"] == group_name) & (programs["n_cells"] >= 20)]
                mw = m.pivot(index=["dataset", "patient_key", "patient_id"], columns="signature_id", values=value).reset_index()
                wide = mw.merge(nw, on=["dataset", "patient_key", "patient_id"], how="outer")
                for relation, x_sig, y_sig, y_comp in RELATIONS:
                    y_col = y_sig if y_comp == "malignant" else f"{y_sig}__{y_comp}"
                    cols = [c for c in [x_sig, y_col] if c in wide.columns]
                    sub = wide[["patient_key", *cols]].dropna() if len(cols) == 2 else pd.DataFrame()
                    count = len(sub)
                    if count >= 5 and sub[x_sig].nunique() > 1 and sub[y_col].nunique() > 1:
                        test = stats.spearmanr(sub[x_sig], sub[y_col])
                        rho, p = float(test.statistic), float(test.pvalue)
                        loo = []
                        x, y = sub[x_sig].to_numpy(float), sub[y_col].to_numpy(float)
                        for i in range(count):
                            loo.append(float(stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic))
                        loo_same = float(np.mean(np.sign(loo) == np.sign(rho))) if rho != 0 else np.nan
                    else:
                        rho = p = loo_same = np.nan
                    rows.append({"dataset": dataset, "method": method, "support_group": group_name,
                                 "relation": relation, "n_patients": count, "spearman_rho": rho,
                                 "p_value": p, "loo_sign_stability_fraction": loo_same})
    result = pd.DataFrame(rows)
    atomic_tsv(result, OUT / "method_specific_supportive_key_relations_v1.tsv")

    g207 = result[(result["dataset"] == "GSE207493") & (result["relation"] == "H2_F1_to_C1QC")]
    status = {
        "schema_version": "phase14_method_specific_supportive_validation_v1",
        "status": "PASS_POSTHOC_METHOD_SPECIFIC_SUPPORTIVE_VALIDATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "rows": int(len(result)),
        "gse207493_f1_c1qc_groups_testable": int(g207["spearman_rho"].notna().sum()),
        "gse207493_f1_c1qc_groups_positive": int((g207["spearman_rho"] > 0).sum()),
        "formal_strict_v6_replaced": False,
    }
    part = OUT / "METHOD_SPECIFIC_SUPPORTIVE_VALIDATION_STATUS_v1.json.part"
    part.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, OUT / "METHOD_SPECIFIC_SUPPORTIVE_VALIDATION_STATUS_v1.json")
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
