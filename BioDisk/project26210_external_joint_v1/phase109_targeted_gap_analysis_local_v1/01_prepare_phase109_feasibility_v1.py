from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("--assignments", type=Path, required=True)
parser.add_argument("--axis-scores", type=Path, required=True)
parser.add_argument("--cellchat-qc", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)

assign = pd.read_csv(args.assignments, sep="\t")
axis = pd.read_csv(args.axis_scores, sep="\t")
chat_qc = pd.read_csv(args.cellchat_qc, sep="\t")

state_names = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
state_counts = (
    assign.loc[assign.state_v2.isin(state_names)]
    .groupby(["study_v1", "patient", "state_v2"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)
for state in state_names:
    if state not in state_counts:
        state_counts[state] = 0
state_counts["macrophage_total"] = state_counts[state_names].sum(axis=1)

axis_counts = (
    axis.groupby(["study", "curated_patient_id"])
    .agg(
        malignant_axis_cells=("exact_cell_id", "size"),
        f1_nonmissing=("F1", lambda x: int(x.notna().sum())),
        strict_tymp_axis_nonmissing=("TYMP_axis_primary", lambda x: int(x.notna().sum())),
        broad_tymp_axis_nonmissing=("TYMP_axis_broad_sensitivity", lambda x: int(x.notna().sum())),
    )
    .reset_index()
    .rename(columns={"study": "study_v1", "curated_patient_id": "patient"})
)

table = state_counts.merge(axis_counts, on=["study_v1", "patient"], how="outer").fillna(0)
integer_columns = state_names + [
    "macrophage_total",
    "malignant_axis_cells",
    "f1_nonmissing",
    "strict_tymp_axis_nonmissing",
    "broad_tymp_axis_nonmissing",
]
for column in integer_columns:
    table[column] = table[column].astype(int)

table["eligible_abundance_primary"] = (
    (table.strict_tymp_axis_nonmissing >= 20)
    & (table.macrophage_total >= 50)
    & (table["SPP1-state Macrophage"] >= 20)
    & (table["C1QC-state Macrophage"] >= 20)
    & (table["Other Macrophage"] >= 20)
)
table["eligible_abundance_sensitivity"] = (
    (table.f1_nonmissing >= 20)
    & (table.macrophage_total >= 50)
    & (table["SPP1-state Macrophage"] >= 10)
    & (table["C1QC-state Macrophage"] >= 10)
)

chat_patients = set(chat_qc.patient.astype(str))
table["phase35_cellchat_run_available"] = table.patient.astype(str).isin(chat_patients)
table.to_csv(args.output / "phase109_patient_feasibility_v1.tsv", sep="\t", index=False)

summary = (
    table.groupby("study_v1")
    .agg(
        patients=("patient", "nunique"),
        primary_eligible=("eligible_abundance_primary", "sum"),
        sensitivity_eligible=("eligible_abundance_sensitivity", "sum"),
        phase35_cellchat_patients=("phase35_cellchat_run_available", "sum"),
    )
    .reset_index()
)
summary.to_csv(args.output / "phase109_study_feasibility_v1.tsv", sep="\t", index=False)

payload = {
    "status": "PASS_FEASIBILITY_COMPLETE",
    "assignment_rows": int(len(assign)),
    "axis_rows": int(len(axis)),
    "phase35_cellchat_patients": int(chat_qc.patient.nunique()),
    "patients_in_bridge": int(table.patient.nunique()),
    "primary_eligible_patients": int(table.eligible_abundance_primary.sum()),
    "sensitivity_eligible_patients": int(table.eligible_abundance_sensitivity.sum()),
    "primary_rule": "strict axis>=20; total macrophage>=50; each SPP1/C1QC/Other>=20",
    "sensitivity_rule": "F1>=20; total macrophage>=50; SPP1/C1QC>=10",
}
(args.output / "PHASE109_FEASIBILITY_STATUS_v1.json").write_text(
    json.dumps(payload, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(payload, indent=2))
