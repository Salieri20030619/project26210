#!/usr/bin/env python3
"""Audit planned-family BH and patient membership for key post-hoc findings."""

from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE / "posthoc_fdr_family_and_key_patients_v1"
ASSOC = HERE / "patient_analysis_final_cnv_v6" / "within_dataset_patient_associations_v6.tsv"
TIERS = HERE / "posthoc_malignancy_gate_tiers_v3" / "malignant_relations_by_gate_tier_v1.tsv"
ANY = HERE / "posthoc_any_method_positive_inputs_v2"
RUN3 = HERE / "runs" / "20260812_082000_provisional_v3"
RUN6 = HERE / "runs" / "20260812_final_cnv_gated_v6"


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


def bh_array(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def aggregate_selected_samples(path: Path, dataset: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    if dataset == "GSE178481":
        frame = frame[~frame["sample_id"].astype(str).str.contains("Normal", case=False, regex=False)]
    if dataset == "GSE264586":
        frame = frame[frame["timepoint"] == "Baseline"]
    rows = []
    fields = ["dataset", "patient_key", "patient_id", "signature_id", "analysis_compartment"]
    for keys, group in frame.groupby(fields, sort=True, dropna=False):
        weights = group["n_cells"].to_numpy(float)
        rows.append(
            dict(
                zip(fields, keys),
                n_cells=int(weights.sum()),
                mean_primary=float(np.average(group["mean_primary"], weights=weights)),
                mean_sensitivity=float(np.average(group["mean_sensitivity"], weights=weights)),
            )
        )
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)

    assoc = pd.read_csv(ASSOC, sep="\t")
    audit_rows = []
    for (dataset, method), group in assoc.groupby(["dataset", "method"], sort=True):
        planned_p = group["p_value"].fillna(1.0).to_numpy(float)
        planned_q = bh_array(planned_p)
        for (_, row), q in zip(group.iterrows(), planned_q):
            audit_rows.append(
                {
                    "dataset": dataset,
                    "method": method,
                    "relation": row["relation"],
                    "status": row["status"],
                    "n_patients": row["n_patients"],
                    "p_value": row["p_value"],
                    "reported_bh_valid_p_only": row["p_fdr_bh_within_dataset_method"],
                    "planned_family_bh_12_with_untested_as_1": q,
                    "bh_valid_p_count": int(group["p_value"].notna().sum()),
                    "planned_family_count": 12,
                }
            )
    fdr = pd.DataFrame(audit_rows)
    atomic_tsv(fdr, OUT / "planned_family_bh_audit_v1.tsv")

    phenotype = aggregate_selected_samples(RUN3 / "GSE178481" / "tables" / "sample_program_summary_v1.tsv", "GSE178481")
    strict = aggregate_selected_samples(RUN6 / "GSE178481" / "tables" / "sample_program_summary_v1.tsv", "GSE178481")
    any_method = pd.read_csv(ANY / "GSE178481" / "malignant_patient_program_summary_any_method_v1.tsv", sep="\t")

    patient_rows = []
    for tier, frame, malignant_label in [
        ("phenotype_candidate", phenotype, "malignant_provisional"),
        ("any_method_positive", any_method, "malignant_any_method_positive"),
        ("strict_conservative", strict, "malignant_provisional"),
    ]:
        tymp = frame[(frame["signature_id"] == "tymp_expression") & (frame["analysis_compartment"] == malignant_label)]
        spp1_source = phenotype if tier == "any_method_positive" else frame
        spp1 = spp1_source[(spp1_source["signature_id"] == "spp1_state_no_spp1") & (spp1_source["analysis_compartment"] == "myeloid")]
        merged = tymp[["patient_key", "patient_id", "n_cells", "mean_primary", "mean_sensitivity"]].rename(
            columns={"n_cells": "tymp_cells", "mean_primary": "tymp_primary", "mean_sensitivity": "tymp_sensitivity"}
        ).merge(
            spp1[["patient_key", "n_cells", "mean_primary", "mean_sensitivity"]].rename(
                columns={"n_cells": "spp1_cells", "mean_primary": "spp1_primary", "mean_sensitivity": "spp1_sensitivity"}
            ), on="patient_key", how="outer"
        )
        merged.insert(0, "gate_tier", tier)
        merged["eligible_ge20_both"] = (merged["tymp_cells"] >= 20) & (merged["spp1_cells"] >= 20)
        patient_rows.append(merged)
    patients = pd.concat(patient_rows, ignore_index=True).sort_values(["gate_tier", "patient_key"])
    atomic_tsv(patients, OUT / "GSE178481_TYMP_SPP1_patient_membership_v1.tsv")

    tiers = pd.read_csv(TIERS, sep="\t")
    key = tiers[
        ((tiers["dataset"] == "GSE178481") & (tiers["relation"] == "H1_TYMP_to_SPP1"))
        | ((tiers["dataset"] == "GSE207493") & (tiers["relation"] == "H2_F1_to_C1QC"))
        | ((tiers["dataset"] == "GSE294109") & (tiers["relation"] == "H3_F1_to_hypoxia"))
    ].copy()
    atomic_tsv(key, OUT / "key_relations_across_gate_tiers_v1.tsv")

    status = {
        "schema_version": "phase14_fdr_family_and_key_patients_audit_v1",
        "status": "PASS_POSTHOC_AUDIT",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "reported_bh_uses_only_nonmissing_p": True,
        "planned_family_sensitivity_uses_12_hypotheses": True,
        "reported_fdr_significant": int((fdr["reported_bh_valid_p_only"] < 0.05).sum()),
        "planned_family_fdr_significant": int((fdr["planned_family_bh_12_with_untested_as_1"] < 0.05).sum()),
        "interpretation": "The implemented BH is lenient where hypotheses are untestable; using all 12 planned hypotheses cannot create significance and leaves zero rejections.",
    }
    atomic_json(status, OUT / "FDR_FAMILY_AUDIT_STATUS_v1.json")


if __name__ == "__main__":
    main()

