from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, wilcoxon


PRIMARY_EXPOSURES = ["TYMP_axis_primary", "z_TYMP_RNA"]
SENSITIVITY_EXPOSURES = ["TYMP_axis_broad_sensitivity", "F1"]
OUTCOMES = ["spp1_fraction", "spp1_vs_c1qc_logratio"]
STATES = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
SEED = 262109
PERMUTATIONS = 100_000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def within_study_rank(values: np.ndarray, studies: np.ndarray) -> np.ndarray:
    out = np.zeros(len(values), dtype=float)
    for study in sorted(set(studies)):
        idx = np.flatnonzero(studies == study)
        ranks = rankdata(values[idx], method="average")
        sd = ranks.std(ddof=0)
        out[idx] = (ranks - ranks.mean()) / sd if sd > 0 else 0.0
    return out


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def fixed_study_test(frame: pd.DataFrame, exposure: str, outcome: str, count_adjust: bool = False, permutations: int = PERMUTATIONS) -> dict:
    studies = frame.study_v1.to_numpy(str)
    x = within_study_rank(frame[exposure].to_numpy(float), studies)
    y = within_study_rank(frame[outcome].to_numpy(float), studies)
    levels = sorted(set(studies))
    design_parts = [np.ones(len(frame))]
    for level in levels[1:]:
        design_parts.append((studies == level).astype(float))
    if count_adjust:
        design_parts.append(within_study_rank(np.log1p(frame.malignant_axis_cells.to_numpy(float)), studies))
        design_parts.append(within_study_rank(np.log1p(frame.macrophage_total.to_numpy(float)), studies))
    design = np.column_stack(design_parts)
    xr = residualize(x, design)
    yr = residualize(y, design)
    observed = correlation(xr, yr)
    p = float("nan")
    if permutations > 0:
        rng = np.random.default_rng(SEED)
        exceed = 0
        for _ in range(permutations):
            perm = yr.copy()
            for level in levels:
                idx = np.flatnonzero(studies == level)
                perm[idx] = rng.permutation(perm[idx])
            value = correlation(xr, perm)
            if np.isfinite(value) and abs(value) >= abs(observed) - 1e-15:
                exceed += 1
        p = (exceed + 1) / (permutations + 1)
    return {"estimate": observed, "permutation_p": p, "n_patients": len(frame), "count_adjusted": count_adjust}


parser = argparse.ArgumentParser()
parser.add_argument("--assignments", type=Path, required=True)
parser.add_argument("--axis-scores", type=Path, required=True)
parser.add_argument("--cellchat-interactions", type=Path, required=True)
parser.add_argument("--cellchat-qc", type=Path, required=True)
parser.add_argument("--spatial-decisions", type=Path, required=True)
parser.add_argument("--spatial-meta", type=Path, required=True)
parser.add_argument("--feasibility", type=Path, required=True)
parser.add_argument("--preregistration", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
if args.output.exists():
    raise FileExistsError(args.output)
args.output.mkdir(parents=True)

input_paths = {
    "assignments": args.assignments,
    "axis_scores": args.axis_scores,
    "cellchat_interactions": args.cellchat_interactions,
    "cellchat_qc": args.cellchat_qc,
    "spatial_decisions": args.spatial_decisions,
    "spatial_meta": args.spatial_meta,
    "feasibility": args.feasibility,
    "preregistration": args.preregistration,
}
input_hashes = {key: sha256(path) for key, path in input_paths.items()}

assign = pd.read_csv(args.assignments, sep="\t")
axis = pd.read_csv(args.axis_scores, sep="\t")
feas = pd.read_csv(args.feasibility, sep="\t")

state_counts = (
    assign.loc[assign.state_v2.isin(STATES)]
    .groupby(["study_v1", "patient", "state_v2"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)
for state in STATES:
    if state not in state_counts:
        state_counts[state] = 0
state_counts["macrophage_total"] = state_counts[STATES].sum(axis=1)
state_counts["spp1_fraction"] = state_counts["SPP1-state Macrophage"] / state_counts.macrophage_total
state_counts["spp1_vs_c1qc_logratio"] = np.log(
    (state_counts["SPP1-state Macrophage"] + 0.5) / (state_counts["C1QC-state Macrophage"] + 0.5)
)

axis_patient = (
    axis.groupby(["study", "curated_patient_id"])[PRIMARY_EXPOSURES + SENSITIVITY_EXPOSURES]
    .median()
    .reset_index()
    .rename(columns={"study": "study_v1", "curated_patient_id": "patient"})
)
patient = state_counts.merge(axis_patient, on=["study_v1", "patient"], how="inner")
patient = patient.merge(
    feas[["study_v1", "patient", "malignant_axis_cells", "eligible_abundance_primary"]],
    on=["study_v1", "patient"],
    how="left",
)
eligible = patient.loc[patient.eligible_abundance_primary.astype(bool)].copy()
if len(eligible) != 20:
    raise RuntimeError(f"Frozen primary eligibility changed: {len(eligible)}")
patient.to_csv(args.output / "phase109_patient_summary_v1.tsv", sep="\t", index=False)

primary_rows = []
adjusted_rows = []
loo_rows = []
study_rows = []
for exposure in PRIMARY_EXPOSURES:
    for outcome in OUTCOMES:
        result = fixed_study_test(eligible, exposure, outcome, count_adjust=False)
        primary_rows.append({"exposure": exposure, "outcome": outcome, **result})
        adjusted = fixed_study_test(eligible, exposure, outcome, count_adjust=True)
        adjusted_rows.append({"exposure": exposure, "outcome": outcome, **adjusted})
        for omitted in sorted(eligible.study_v1.unique()):
            subset = eligible.loc[eligible.study_v1 != omitted]
            loo = fixed_study_test(subset, exposure, outcome, count_adjust=False, permutations=0)
            loo_rows.append({"exposure": exposure, "outcome": outcome, "omitted_study": omitted, **loo})
        for study, group in eligible.groupby("study_v1"):
            rho, p = spearmanr(group[exposure], group[outcome]) if len(group) >= 3 else (np.nan, np.nan)
            study_rows.append({"exposure": exposure, "outcome": outcome, "study": study, "n_patients": len(group), "spearman_rho": rho, "spearman_p": p})

primary = pd.DataFrame(primary_rows)
primary["bh_q_4"] = bh(primary.permutation_p)
loo = pd.DataFrame(loo_rows)
for index, row in primary.iterrows():
    subset = loo.loc[(loo.exposure == row.exposure) & (loo.outcome == row.outcome)]
    positive_all_loo = bool((subset.estimate > 0).all())
    no_gse207493 = float(subset.loc[subset.omitted_study == "GSE207493", "estimate"].iloc[0])
    support = bool(row.estimate > 0 and row.bh_q_4 < 0.05 and positive_all_loo and no_gse207493 > 0)
    primary.loc[index, "all_leave_one_study_positive"] = positive_all_loo
    primary.loc[index, "estimate_without_GSE207493"] = no_gse207493
    primary.loc[index, "support"] = support
    primary.loc[index, "decision"] = "SUPPORTED" if support else ("CONTEXT_DEPENDENT" if row.estimate > 0 else "NOT_SUPPORTED")

primary.to_csv(args.output / "abundance_primary_results_v1.tsv", sep="\t", index=False)
pd.DataFrame(adjusted_rows).to_csv(args.output / "abundance_count_adjusted_sensitivity_v1.tsv", sep="\t", index=False)
loo.to_csv(args.output / "abundance_leave_one_study_out_v1.tsv", sep="\t", index=False)
pd.DataFrame(study_rows).to_csv(args.output / "abundance_study_specific_v1.tsv", sep="\t", index=False)

sensitivity_rows = []
for exposure in SENSITIVITY_EXPOSURES:
    for outcome in OUTCOMES:
        sensitivity_rows.append({"exposure": exposure, "outcome": outcome, **fixed_study_test(eligible, exposure, outcome, count_adjust=False)})
sensitivity = pd.DataFrame(sensitivity_rows)
sensitivity["bh_q_4"] = bh(sensitivity.permutation_p)
sensitivity.to_csv(args.output / "abundance_exposure_sensitivity_v1.tsv", sep="\t", index=False)

interactions = pd.read_csv(args.cellchat_interactions, sep="\t")
qc = pd.read_csv(args.cellchat_qc, sep="\t")
groups = ["TYMP-high Malignant", "TYMP-low Malignant", *STATES, "Endothelial"]
pair_sum = interactions.groupby(["study", "patient", "source", "target"], as_index=False).prob.sum()

def value(frame: pd.DataFrame, source: str, target: str) -> float:
    hit = frame.loc[(frame.source == source) & (frame.target == target), "prob"]
    return float(hit.iloc[0]) if len(hit) else 0.0

chat_rows = []
for _, qrow in qc.iterrows():
    z = pair_sum.loc[pair_sum.patient == qrow.patient]
    high_to = {state: value(z, "TYMP-high Malignant", state) for state in STATES}
    low_to = {state: value(z, "TYMP-low Malignant", state) for state in STATES}
    to_high = {state: value(z, state, "TYMP-high Malignant") for state in STATES}
    to_low = {state: value(z, state, "TYMP-low Malignant") for state in STATES}
    malignant_to = (high_to[STATES[0]] - low_to[STATES[0]]) - np.mean([high_to[s] - low_to[s] for s in STATES[1:]])
    state_to = (to_high[STATES[0]] - to_low[STATES[0]]) - np.mean([to_high[s] - to_low[s] for s in STATES[1:]])
    spp1_endo = value(z, STATES[0], "Endothelial") - np.mean([value(z, s, "Endothelial") for s in STATES[1:]])
    endo_spp1 = value(z, "Endothelial", STATES[0]) - np.mean([value(z, "Endothelial", s) for s in STATES[1:]])
    chat_rows.append({"study": qrow.study, "patient": qrow.patient, "malignant_to_spp1_specificity": malignant_to, "spp1_to_malignant_specificity": state_to, "spp1_to_endothelial_specificity": spp1_endo, "endothelial_to_spp1_specificity": endo_spp1})
chat = pd.DataFrame(chat_rows)
chat.to_csv(args.output / "cellchat_patient_fine_state_specificity_v1.tsv", sep="\t", index=False)

chat_summary = []
chat_loo = []
for endpoint in chat.columns[2:]:
    values = chat[endpoint].to_numpy(float)
    test = wilcoxon(values, zero_method="wilcox", alternative="two-sided", method="approx") if np.any(values != 0) else None
    row = {"endpoint": endpoint, "n_patients": len(values), "n_positive": int((values > 0).sum()), "n_negative": int((values < 0).sum()), "n_zero": int((values == 0).sum()), "positive_fraction": float((values > 0).mean()), "median": float(np.median(values)), "wilcoxon_p": float(test.pvalue) if test else 1.0}
    chat_summary.append(row)
    for patient_name in chat.patient:
        kept = chat.loc[chat.patient != patient_name, endpoint]
        chat_loo.append({"endpoint": endpoint, "omitted_patient": patient_name, "median": float(kept.median())})
chat_summary = pd.DataFrame(chat_summary)
chat_summary["bh_q_4"] = bh(chat_summary.wilcoxon_p)
chat_loo = pd.DataFrame(chat_loo)
for index, row in chat_summary.iterrows():
    all_loo = bool((chat_loo.loc[chat_loo.endpoint == row.endpoint, "median"] > 0).all())
    support = bool(row.n_patients >= 10 and row["median"] > 0 and row.positive_fraction >= 2 / 3 and row.bh_q_4 < 0.05 and all_loo)
    chat_summary.loc[index, "all_leave_one_patient_positive"] = all_loo
    chat_summary.loc[index, "support"] = support
    chat_summary.loc[index, "decision"] = "SUPPORTED" if support else ("CONTEXT_DEPENDENT" if row["median"] > 0 else "NOT_SUPPORTED")
chat_summary.to_csv(args.output / "cellchat_fine_state_specificity_summary_v1.tsv", sep="\t", index=False)
chat_loo.to_csv(args.output / "cellchat_leave_one_patient_out_v1.tsv", sep="\t", index=False)

decisions = pd.read_csv(args.spatial_decisions, sep="\t")
meta = pd.read_csv(args.spatial_meta, sep="\t")
selected_ids = ["B01_TYMP_to_SPP1", "B08_SPP1_to_Tdysfunction", "B09_SPP1_to_Vascular"]
spatial = meta.loc[meta.relation_id.isin(selected_ids)].merge(decisions, on="relation_id", how="left")
spatial.to_csv(args.output / "spatial_frozen_edge_synthesis_v1.tsv", sep="\t", index=False)

summary_lines = [
    "# Phase 109 targeted gap analysis",
    "",
    "## Abundance bridge",
]
for _, row in primary.iterrows():
    summary_lines.append(f"- {row.exposure} -> {row.outcome}: estimate={row.estimate:.4f}, P={row.permutation_p:.5g}, q={row.bh_q_4:.5g}, decision={row.decision}.")
summary_lines += ["", "## Fine-state CellChat"]
for _, row in chat_summary.iterrows():
    summary_lines.append(f"- {row.endpoint}: median={row['median']:.6g}, positive={int(row.n_positive)}/{int(row.n_patients)}, P={row.wilcoxon_p:.5g}, q={row.bh_q_4:.5g}, decision={row.decision}.")
summary_lines += [
    "",
    "## Spatial frozen synthesis",
    "- B01 TYMP-to-SPP1: unsupported in FFPE and frozen.",
    "- B08 SPP1-to-T-cell dysfunction: FFPE-only support.",
    "- B09 SPP1-to-vascular: supported in FFPE and frozen.",
    "",
    "## Boundary",
    "These analyses do not establish direct TYMP-to-MTA or MTA-to-SPP1 causality. CellChat is computational compatibility and the spatial synthesis introduces no new test.",
]
(args.output / "PHASE109_TARGETED_ANALYSIS_SUMMARY_v1.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

status = {
    "status": "PASS_ANALYSIS_COMPLETE",
    "primary_abundance_supported": int(primary.support.sum()),
    "primary_abundance_total": len(primary),
    "fine_state_cellchat_supported": int(chat_summary.support.sum()),
    "fine_state_cellchat_total": len(chat_summary),
    "eligible_abundance_patients": len(eligible),
    "cellchat_patients": len(chat),
    "permutations": PERMUTATIONS,
    "seed": SEED,
    "input_hashes": input_hashes,
    "direct_mta_tested": False,
    "causal_mta_to_spp1_tested": False,
}
(args.output / "PHASE109_ANALYSIS_STATUS_v1.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
print(json.dumps(status, indent=2))
