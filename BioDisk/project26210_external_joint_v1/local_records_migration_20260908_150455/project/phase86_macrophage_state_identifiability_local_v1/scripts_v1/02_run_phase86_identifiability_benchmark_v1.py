#!/usr/bin/env python

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import rankdata, spearmanr


MACRO = ["SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage"]
STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]
ANCHORS = {
    "SPP1", "GPNMB", "LGALS3", "FABP5", "CD9", "LPL",
    "C1QA", "C1QB", "C1QC", "APOE", "TYMP",
}
BASE_SEED = 2621086
READS_PER_MIXTURE = 500_000
MIXTURES_PER_TIER = 24


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def stable_seed(text: str) -> int:
    value = int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)
    return (BASE_SEED + value) % (2**32 - 1)


def normalize_reference(frame: pd.DataFrame) -> pd.DataFrame:
    norms = np.sqrt(np.square(frame.to_numpy(dtype=float)).sum(axis=0))
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0):
        raise RuntimeError("Zero or invalid reference column")
    return frame / norms


def rebuild_reference(
    profiles: pd.DataFrame,
    genes: list[str],
    held_out: str,
    states: list[str],
    anchor_excluded: bool,
) -> pd.DataFrame:
    train = profiles.loc[profiles["dataset"] != held_out].copy()
    if anchor_excluded:
        genes = [g for g in genes if g not in ANCHORS]
    if not genes:
        raise RuntimeError("No genes remain for rebuilt reference")

    values = train[genes].to_numpy(dtype=float) / train["n_cells"].to_numpy(dtype=float)[:, None]
    mean_profiles = pd.DataFrame(values, columns=genes)
    mean_profiles.insert(0, "state", train["state"].to_numpy())
    mean_profiles.insert(0, "dataset", train["dataset"].to_numpy())
    dataset_state = mean_profiles.groupby(["dataset", "state"], sort=True)[genes].mean().reset_index()
    state_mean = dataset_state.groupby("state", sort=True)[genes].mean()
    missing = [state for state in states if state not in state_mean.index]
    if missing:
        raise RuntimeError(f"Missing rebuilt-reference states for {held_out}: {missing}")
    state_mean = state_mean.loc[states].T

    logged = np.log1p(state_mean)
    specificity = pd.DataFrame(index=logged.index, columns=states, dtype=float)
    for state in states:
        other = logged.drop(columns=state).mean(axis=1)
        specificity[state] = (logged[state] - other).clip(lower=0)

    members: set[str] = set()
    for state in states:
        positive = specificity.loc[specificity[state] > 0, state].sort_values(ascending=False)
        if len(positive) < 20:
            raise RuntimeError(f"Too few positive markers for {held_out}/{state}: {len(positive)}")
        members.update(positive.head(60).index.tolist())
    reference = specificity.loc[sorted(members), states]
    return normalize_reference(reference)


def project(counts: np.ndarray, reference: pd.DataFrame) -> dict[str, float]:
    ranks = rankdata(counts, method="average") / len(counts)
    coefficients, _ = nnls(reference.to_numpy(dtype=float), ranks)
    total = float(coefficients.sum())
    if not math.isfinite(total) or total <= 0:
        coefficients = np.zeros_like(coefficients)
    else:
        coefficients = coefficients / total
    return {state: float(value) for state, value in zip(reference.columns, coefficients)}


def composition_schedule(rng: np.random.Generator) -> np.ndarray:
    fixed = np.asarray([
        [0.80, 0.10, 0.10],
        [0.10, 0.80, 0.10],
        [0.10, 0.10, 0.80],
        [1 / 3, 1 / 3, 1 / 3],
    ])
    random = rng.dirichlet(np.ones(3), size=20)
    return np.vstack([fixed, random])


def patient_profiles(frame: pd.DataFrame, genes: list[str]) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    profiles: dict[str, np.ndarray] = {}
    cells: dict[str, int] = {}
    for _, row in frame.iterrows():
        state = str(row["state"])
        n_cells = int(row["n_cells"])
        values = row[genes].to_numpy(dtype=float) / n_cells
        profiles[state] = values
        cells[state] = n_cells
    return profiles, cells


def simulate_mixtures(profiles: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    count_rows: list[dict] = []
    truth_rows: list[dict] = []
    grouped = profiles.groupby(["dataset", "patient_id"], sort=True)

    for (dataset, patient), frame in grouped:
        state_profiles, state_cells = patient_profiles(frame, genes)
        if any(state not in state_profiles or state_cells[state] < 25 for state in MACRO):
            continue
        if sum(state_cells[state] for state in MACRO) < 150:
            continue
        rng = np.random.default_rng(stable_seed(f"{dataset}|{patient}"))
        macro_compositions = composition_schedule(rng)

        for tier in ["macrophage_only", "tissue_context"]:
            nonmacro = [state for state in STATES if state not in MACRO and state in state_profiles]
            if tier == "tissue_context":
                if len(nonmacro) < 3 or sum(state_cells[state] for state in nonmacro) < 200:
                    continue
                nonmacro_weights = np.asarray([state_cells[state] for state in nonmacro], dtype=float)
                nonmacro_weights /= nonmacro_weights.sum()
            else:
                nonmacro_weights = np.asarray([], dtype=float)

            for index, macro_within in enumerate(macro_compositions):
                macro_total = 1.0 if tier == "macrophage_only" else [0.05, 0.10, 0.20, 0.40][index % 4]
                truth = {state: 0.0 for state in STATES}
                expected = np.zeros(len(genes), dtype=float)
                for state, weight in zip(MACRO, macro_within):
                    full_weight = macro_total * float(weight)
                    truth[state] = full_weight
                    expected += full_weight * state_profiles[state]
                if tier == "tissue_context":
                    for state, weight in zip(nonmacro, nonmacro_weights):
                        full_weight = (1.0 - macro_total) * float(weight)
                        truth[state] = full_weight
                        expected += full_weight * state_profiles[state]

                expected = np.clip(expected, 0, None)
                if not np.isfinite(expected).all() or expected.sum() <= 0:
                    raise RuntimeError(f"Invalid expected profile: {dataset}/{patient}/{tier}/{index}")
                probability = expected / expected.sum()
                mixture_id = f"{dataset}__{patient}__{tier}__{index + 1:02d}"
                mixture_rng = np.random.default_rng(stable_seed(mixture_id))
                counts = mixture_rng.multinomial(READS_PER_MIXTURE, probability)

                count_row = {
                    "mixture_id": mixture_id,
                    "dataset": dataset,
                    "patient_id": patient,
                    "tier": tier,
                    "mixture_index": index + 1,
                }
                count_row.update(dict(zip(genes, counts.astype(int))))
                count_rows.append(count_row)

                truth_row = {
                    "mixture_id": mixture_id,
                    "dataset": dataset,
                    "patient_id": patient,
                    "tier": tier,
                    "true_total_macrophage": macro_total,
                }
                truth_row.update({f"true_{state}": truth[state] for state in STATES})
                truth_row.update({f"true_within_{state}": float(weight) for state, weight in zip(MACRO, macro_within)})
                truth_rows.append(truth_row)

    if not count_rows:
        raise RuntimeError("No eligible synthetic mixtures")
    return pd.DataFrame(count_rows), pd.DataFrame(truth_rows)


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    if x.nunique() < 2 or y.nunique() < 2:
        return float("nan")
    return float(spearmanr(x.to_numpy(dtype=float), y.to_numpy(dtype=float)).statistic)


def build_metrics(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict] = []
    accuracy_rows: list[dict] = []
    total_rows: list[dict] = []

    keys = ["dataset", "tier", "method"]
    for key, frame in predictions.groupby(keys, sort=True):
        base = dict(zip(keys, key))
        true_matrix = frame[[f"true_within_{state}" for state in MACRO]].to_numpy(dtype=float)
        est_matrix = frame[[f"estimated_within_{state}" for state in MACRO]].to_numpy(dtype=float)
        accuracy_rows.append({
            **base,
            "n_mixtures": len(frame),
            "dominant_state_accuracy": float(np.mean(np.argmax(true_matrix, axis=1) == np.argmax(est_matrix, axis=1))),
        })
        for state in MACRO:
            truth = frame[f"true_within_{state}"]
            estimate = frame[f"estimated_within_{state}"]
            difference = estimate.to_numpy(dtype=float) - truth.to_numpy(dtype=float)
            metric_rows.append({
                **base,
                "state": state,
                "n_mixtures": len(frame),
                "spearman_rho": safe_spearman(truth, estimate),
                "mae": float(np.mean(np.abs(difference))),
                "rmse": float(np.sqrt(np.mean(np.square(difference)))),
                "estimated_zero_fraction": float(np.mean(estimate.to_numpy(dtype=float) <= 1e-12)),
            })
        if frame["method_states"].iloc[0] == 13:
            truth = frame["true_total_macrophage"]
            estimate = frame["estimated_total_macrophage"]
            total_rows.append({
                **base,
                "n_mixtures": len(frame),
                "spearman_rho": safe_spearman(truth, estimate),
                "mae": float(np.mean(np.abs(estimate.to_numpy(dtype=float) - truth.to_numpy(dtype=float)))),
            })
    return pd.DataFrame(metric_rows), pd.DataFrame(accuracy_rows), pd.DataFrame(total_rows)


def evaluate_gates(metrics: pd.DataFrame, accuracy: pd.DataFrame, total: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    methods = {
        "frozen13": dict(rho=0.70, mae=0.15, zero=0.20, accuracy=0.80, total_rho=0.70, total_mae=0.10),
        "lodo13": dict(rho=0.60, mae=0.18, zero=1.00, accuracy=0.75, total_rho=0.60, total_mae=0.12),
        "frozen13_anchor_excluded": dict(rho=0.50, mae=0.20, zero=0.50, accuracy=0.70),
        "lodo3_anchor_excluded": dict(rho=0.50, mae=0.20, zero=0.50, accuracy=0.70),
    }
    for method, threshold in methods.items():
        subset = metrics.loc[metrics["method"] == method]
        for (dataset, tier), frame in subset.groupby(["dataset", "tier"], sort=True):
            acc = accuracy.loc[
                (accuracy["method"] == method) & (accuracy["dataset"] == dataset) & (accuracy["tier"] == tier),
                "dominant_state_accuracy",
            ]
            if len(frame) != 3 or len(acc) != 1:
                raise RuntimeError(f"Incomplete metric group: {method}/{dataset}/{tier}")
            state_pass = bool(
                (frame["spearman_rho"] >= threshold["rho"]).all()
                and (frame["mae"] <= threshold["mae"]).all()
                and (frame["estimated_zero_fraction"] <= threshold["zero"]).all()
            )
            accuracy_pass = bool(float(acc.iloc[0]) >= threshold["accuracy"])
            total_pass = True
            total_rho = float("nan")
            total_mae = float("nan")
            if method in {"frozen13", "lodo13"} and tier == "tissue_context":
                current = total.loc[
                    (total["method"] == method) & (total["dataset"] == dataset) & (total["tier"] == tier)
                ]
                if len(current) != 1:
                    raise RuntimeError(f"Missing total-macrophage metric: {method}/{dataset}/{tier}")
                total_rho = float(current["spearman_rho"].iloc[0])
                total_mae = float(current["mae"].iloc[0])
                total_pass = bool(total_rho >= threshold["total_rho"] and total_mae <= threshold["total_mae"])
            rows.append({
                "method": method,
                "dataset": dataset,
                "tier": tier,
                "minimum_state_rho": float(frame["spearman_rho"].min()),
                "maximum_state_mae": float(frame["mae"].max()),
                "maximum_state_zero_fraction": float(frame["estimated_zero_fraction"].max()),
                "dominant_state_accuracy": float(acc.iloc[0]),
                "total_macrophage_rho": total_rho,
                "total_macrophage_mae": total_mae,
                "state_gate_pass": state_pass,
                "accuracy_gate_pass": accuracy_pass,
                "total_gate_pass": total_pass,
                "group_gate_pass": bool(state_pass and accuracy_pass and total_pass),
            })
    return pd.DataFrame(rows)


def write_summary(path: Path, decision: dict, gates: pd.DataFrame, metrics: pd.DataFrame) -> None:
    lines = [
        "# Phase 86 macrophage-state identifiability benchmark",
        "",
        f"- Computational status: **{decision['status']}**.",
        f"- Scientific decision: **{decision['decision']}**.",
        f"- Eligible patients: {decision['eligible_patients']}.",
        f"- Synthetic mixtures: {decision['mixtures']}.",
        f"- Operational frozen-reference gate: {decision['operational_gate_pass']}.",
        f"- Cross-dataset generalization gate: {decision['generalization_gate_pass']}.",
        f"- Anchor-independent sensitivity gate: {decision['anchor_gate_pass']}.",
        f"- Failure interpretation: `{decision['failure_interpretation']}`.",
        "",
        "## Gate table",
        "",
        "| method | dataset | tier | min rho | max MAE | max zero | dominant accuracy | total rho | total MAE | pass |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in gates.sort_values(["method", "dataset", "tier"]).itertuples(index=False):
        lines.append(
            f"| {row.method} | {row.dataset} | {row.tier} | {row.minimum_state_rho:.3f} | "
            f"{row.maximum_state_mae:.3f} | {row.maximum_state_zero_fraction:.3f} | "
            f"{row.dominant_state_accuracy:.3f} | "
            f"{row.total_macrophage_rho if math.isfinite(row.total_macrophage_rho) else float('nan'):.3f} | "
            f"{row.total_macrophage_mae if math.isfinite(row.total_macrophage_mae) else float('nan'):.3f} | "
            f"{row.group_gate_pass} |"
        )
    lines.extend([
        "",
        "## Boundary",
        "",
        "The benchmark tests recoverability of frozen Phase 8 expression-program labels. It does not validate lineage, causal polarization, 5′-MTA mediation or clinical prediction.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-root", required=True, type=Path)
    args = parser.parse_args()
    phase_root = args.phase_root.resolve()
    output = phase_root / "results_v1"
    part = phase_root / "results_v1.part"
    if output.exists() or part.exists():
        raise RuntimeError("Refusing existing benchmark output")
    part.mkdir(parents=True)

    prereg = phase_root / "PHASE86_PREREGISTRATION_v1.md"
    profile_path = phase_root / "profiles_v1" / "patient_state_profiles_v1.tsv.gz"
    reference_root = phase_root / "inputs_v1" / "phase83_reference_v1"
    reference_path = reference_root / "phase83_reference_specificity_v1.tsv.gz"
    coverage_path = reference_root / "phase83_marker_coverage_v1.tsv"
    for path in [prereg, profile_path, reference_path, coverage_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    profiles = pd.read_csv(profile_path, sep="\t")
    coverage = pd.read_csv(coverage_path, sep="\t")
    coverage["present_bool"] = coverage["present"].astype(str).str.lower().eq("true")
    common = sorted(
        gene.upper()
        for gene, frame in coverage.groupby("gene")
        if frame["dataset"].nunique() == 2 and frame["present_bool"].all()
    )
    if len(common) != 478:
        raise RuntimeError(f"Expected 478 common genes, found {len(common)}")
    missing_profile = sorted(set(common) - set(profiles.columns))
    if missing_profile:
        raise RuntimeError(f"Profiles missing common genes: {missing_profile[:10]}")
    if profiles.duplicated(["dataset", "patient_id", "state"]).any():
        raise RuntimeError("Duplicate patient-state profile")

    frozen = pd.read_csv(reference_path, sep="\t", compression="gzip")
    frozen.iloc[:, 0] = frozen.iloc[:, 0].astype(str).str.upper()
    frozen = frozen.set_index(frozen.columns[0])
    if frozen.columns.tolist() != STATES:
        raise RuntimeError("Unexpected frozen reference states")
    frozen = normalize_reference(frozen.loc[common, STATES])
    frozen_anchor = normalize_reference(frozen.loc[[g for g in common if g not in ANCHORS], STATES])

    counts, truth = simulate_mixtures(profiles, common)
    eligible_patients = truth[["dataset", "patient_id"]].drop_duplicates().shape[0]
    count_gene_columns = common
    predictions: list[dict] = []

    references: dict[tuple[str, str], pd.DataFrame] = {}
    for dataset in sorted(profiles["dataset"].unique()):
        references[(dataset, "frozen13")] = frozen
        references[(dataset, "frozen13_anchor_excluded")] = frozen_anchor
        references[(dataset, "lodo13")] = rebuild_reference(profiles, common, dataset, STATES, False)
        references[(dataset, "lodo3_anchor_excluded")] = rebuild_reference(profiles, common, dataset, MACRO, True)

    truth_index = truth.set_index("mixture_id")
    for row_index, row in counts.iterrows():
        mixture_id = str(row["mixture_id"])
        dataset = str(row["dataset"])
        truth_row = truth_index.loc[mixture_id]
        full_counts = dict(zip(count_gene_columns, row[count_gene_columns].to_numpy(dtype=float)))
        for method in ["frozen13", "frozen13_anchor_excluded", "lodo13", "lodo3_anchor_excluded"]:
            reference = references[(dataset, method)]
            vector = np.asarray([full_counts[gene] for gene in reference.index], dtype=float)
            estimate = project(vector, reference)
            macro_values = np.asarray([estimate.get(state, 0.0) for state in MACRO], dtype=float)
            macro_total = float(macro_values.sum())
            within = macro_values / macro_total if macro_total > 0 else np.zeros(3)
            prediction = {
                "mixture_id": mixture_id,
                "dataset": dataset,
                "patient_id": row["patient_id"],
                "tier": row["tier"],
                "method": method,
                "method_states": len(reference.columns),
                "true_total_macrophage": float(truth_row["true_total_macrophage"]),
                "estimated_total_macrophage": macro_total,
            }
            for state, value in zip(MACRO, within):
                prediction[f"true_within_{state}"] = float(truth_row[f"true_within_{state}"])
                prediction[f"estimated_within_{state}"] = float(value)
            predictions.append(prediction)

    prediction_frame = pd.DataFrame(predictions)
    metrics, accuracy, total = build_metrics(prediction_frame)
    gates = evaluate_gates(metrics, accuracy, total)

    operational = bool(gates.loc[gates["method"] == "frozen13", "group_gate_pass"].all())
    generalization = bool(gates.loc[gates["method"] == "lodo13", "group_gate_pass"].all())
    anchor = bool(gates.loc[gates["method"].isin(["frozen13_anchor_excluded", "lodo3_anchor_excluded"]), "group_gate_pass"].all())
    overall = operational and generalization and anchor
    frozen_macro = gates.loc[(gates["method"] == "frozen13") & (gates["tier"] == "macrophage_only"), "group_gate_pass"].all()
    frozen_tissue = gates.loc[(gates["method"] == "frozen13") & (gates["tier"] == "tissue_context"), "group_gate_pass"].all()
    if not frozen_macro:
        failure = "INTRINSIC_MACROPHAGE_STATE_NONIDENTIFIABILITY"
    elif not frozen_tissue:
        failure = "TISSUE_CONTEXT_COMPARTMENT_CONFOUNDING"
    elif not generalization:
        failure = "CROSS_DATASET_GENERALIZATION_FAILURE"
    elif not anchor:
        failure = "LABEL_ANCHOR_DEPENDENCE"
    else:
        failure = "NONE"

    decision = {
        "status": "PASS_COMPUTATION_COMPLETE",
        "decision": "IDENTIFIABLE_FOR_FINE_STATE_BULK_USE" if overall else "NOT_IDENTIFIABLE_FOR_FINE_STATE_BULK_USE",
        "operational_gate_pass": operational,
        "generalization_gate_pass": generalization,
        "anchor_gate_pass": anchor,
        "failure_interpretation": failure,
        "eligible_patients": int(eligible_patients),
        "mixtures": int(len(counts)),
        "prediction_rows": int(len(prediction_frame)),
        "preregistration_sha256": sha256(prereg),
        "profile_sha256": sha256(profile_path),
        "frozen_reference_sha256": sha256(reference_path),
        "marker_coverage_sha256": sha256(coverage_path),
        "seed": BASE_SEED,
        "reads_per_mixture": READS_PER_MIXTURE,
    }

    counts.to_csv(part / "phase86_mixture_marker_counts_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    truth.to_csv(part / "phase86_mixture_truth_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    prediction_frame.to_csv(part / "phase86_mixture_predictions_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    metrics.to_csv(part / "phase86_state_metrics_v1.tsv", sep="\t", index=False)
    accuracy.to_csv(part / "phase86_dominant_accuracy_v1.tsv", sep="\t", index=False)
    total.to_csv(part / "phase86_total_macrophage_metrics_v1.tsv", sep="\t", index=False)
    gates.to_csv(part / "phase86_gate_decisions_v1.tsv", sep="\t", index=False)
    with (part / "PHASE86_DECISION_v1.json").open("w", encoding="utf-8") as handle:
        json.dump(decision, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    write_summary(part / "PHASE86_IDENTIFIABILITY_SUMMARY_v1.md", decision, gates, metrics)
    (part / "PHASE86_BENCHMARK_COMPUTATION_PASS_v1").write_text("PASS_COMPUTATION_COMPLETE\n", encoding="utf-8")

    manifest_rows = []
    for path in sorted(part.iterdir()):
        if path.name == "FILE_MANIFEST_SHA256_v1.tsv":
            continue
        manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(part / "FILE_MANIFEST_SHA256_v1.tsv", sep="\t", index=False)
    os.replace(part, output)
    print(json.dumps(decision, ensure_ascii=False))


if __name__ == "__main__":
    main()
