#!/usr/bin/env python3
"""Post-hoc association sensitivity using fresh external strict consensus for all cells."""

from __future__ import annotations

import importlib.util
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
FRESH = HERE / "posthoc_fresh_external_consensus_inputs_v1"
OUT = HERE / "posthoc_fresh_external_consensus_analysis_v1"
BASE_SCRIPT = HERE / "10_analyze_malignancy_gate_tiers_v1.py"

spec = importlib.util.spec_from_file_location("gate_tiers", BASE_SCRIPT)
gate = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(gate)


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", na_rep="NA", lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    if path.exists() or part.exists():
        raise FileExistsError(path)
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    rows = []
    for dataset_index, dataset in enumerate(gate.DATASETS):
        malignant = pd.read_csv(FRESH / dataset / "malignant_patient_program_summary_fresh_external_v1.tsv", sep="\t")
        nonmal = pd.read_csv(HERE / "runs" / "20260812_final_cnv_gated_v6" / dataset / "tables" / "patient_program_summary_v1.tsv", sep="\t")
        nonmal = nonmal[nonmal["analysis_compartment"] != "malignant_provisional"]
        for method_index, method in enumerate(gate.METHODS):
            value = "mean_primary" if method == "primary_gene_z_mean" else "mean_sensitivity"
            mw = malignant[malignant["n_cells"] >= gate.MIN_CELLS].pivot(
                index=["dataset", "patient_key", "patient_id"], columns="signature_id", values=value
            ).reset_index()
            nw = nonmal[nonmal["n_cells"] >= gate.MIN_CELLS].pivot(
                index=["dataset", "patient_key", "patient_id"], columns=["signature_id", "analysis_compartment"], values=value
            )
            nw.columns = [f"{a}__{b}" for a, b in nw.columns]
            wide = mw.merge(nw.reset_index(), on=["dataset", "patient_key", "patient_id"], how="outer")
            for relation_index, (relation, x_sig, y_sig, y_comp, expected) in enumerate(gate.RELATIONS):
                y_col = y_sig if y_comp == "malignant" else f"{y_sig}__{y_comp}"
                cols = [c for c in [x_sig, y_col] if c in wide]
                sub = wide[["dataset", "patient_key", "patient_id", *cols]].dropna()
                n = len(sub) if len(cols) == 2 else 0
                if n >= gate.MIN_PATIENTS and sub[x_sig].nunique() > 1 and sub[y_col].nunique() > 1:
                    x = sub[x_sig].to_numpy(float)
                    y = sub[y_col].to_numpy(float)
                    test = stats.spearmanr(x, y)
                    rho, p = float(test.statistic), float(test.pvalue)
                    low, high, valid = gate.bootstrap_ci(x, y, gate.SEED + dataset_index * 1000 + method_index * 100 + relation_index)
                    loo = [float(stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic) for i in range(n)]
                    status = "TESTED_POSTHOC_FRESH_EXTERNAL_CONSENSUS"
                else:
                    rho = p = low = high = valid = np.nan
                    loo = []
                    status = "NOT_TESTED_LT5_ELIGIBLE_PATIENTS"
                rows.append({
                    "gate_tier": "fresh_external_strict_consensus",
                    "dataset": dataset,
                    "method": method,
                    "relation": relation,
                    "n_patients": n,
                    "spearman_rho": rho,
                    "bootstrap_ci_low": low,
                    "bootstrap_ci_high": high,
                    "bootstrap_valid_fraction": valid,
                    "p_value": p,
                    "direction_matches_expected": bool(np.sign(rho) == expected) if expected and np.isfinite(rho) else np.nan,
                    "loo_min_rho": min(loo) if loo else np.nan,
                    "loo_max_rho": max(loo) if loo else np.nan,
                    "loo_sign_stability_fraction": float(np.mean(np.sign(loo) == np.sign(rho))) if loo and rho != 0 else np.nan,
                    "status": status,
                })
    frame = pd.DataFrame(rows)
    frame["p_fdr_bh_within_dataset_method_9_malignant"] = frame.groupby(["dataset", "method"], group_keys=False)["p_value"].apply(gate.bh)
    atomic_tsv(frame, OUT / "fresh_external_strict_associations_v1.tsv")

    existing = pd.read_csv(HERE / "posthoc_malignancy_gate_tiers_v3" / "malignant_relations_by_gate_tier_v1.tsv", sep="\t")
    common = ["gate_tier", "dataset", "method", "relation", "n_patients", "spearman_rho", "bootstrap_ci_low", "bootstrap_ci_high", "p_value", "loo_sign_stability_fraction", "direction_matches_expected", "status"]
    combined = pd.concat([existing[common], frame[common]], ignore_index=True)
    atomic_tsv(combined.sort_values(["dataset", "method", "relation", "gate_tier"]), OUT / "four_gate_tier_key_comparison_source_v1.tsv")
    status = {
        "schema_version": "phase14_fresh_external_consensus_posthoc_v1",
        "status": "PASS_POSTHOC_SENSITIVITY",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "rows": int(len(frame)),
        "rows_tested": int(frame["spearman_rho"].notna().sum()),
        "nominal_p_lt_0_05": int((frame["p_value"] < 0.05).sum()),
        "fdr_q_lt_0_05": int((frame["p_fdr_bh_within_dataset_method_9_malignant"] < 0.05).sum()),
        "warning": "Fresh external consensus tier is post-hoc and tests source-priority sensitivity; it does not replace Part1-priority v6.",
    }
    atomic_json(status, OUT / "FRESH_EXTERNAL_ANALYSIS_STATUS_v1.json")


if __name__ == "__main__":
    main()

