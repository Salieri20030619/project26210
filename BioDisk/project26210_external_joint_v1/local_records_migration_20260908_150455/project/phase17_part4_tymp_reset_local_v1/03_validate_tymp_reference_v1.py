#!/usr/bin/env python3
"""Independent validation for the Phase 17 TYMP-focused reference."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr


EXPECTED_STATES = sorted([
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
    "TYMP-high Malignant", "TYMP-low Malignant",
])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--old-signature", type=Path, required=True)
    p.add_argument("--bridge", type=Path, required=True)
    p.add_argument("--thresholds", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, sep="\t", index=False, encoding="utf-8", lineterminator="\n")
    os.replace(part, path)


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output}")
    work = args.output.with_name(args.output.name + ".part")
    if work.exists():
        raise FileExistsError(f"Stale validation work directory: {work}")
    work.mkdir(parents=True)

    files = {
        "signature": args.reference / "reference_tymp_v1_signatures_study_balanced_mean_counts.tsv.gz",
        "state_qc": args.reference / "reference_tymp_v1_state_qc.tsv",
        "patient_qc": args.reference / "reference_tymp_v1_patient_state_qc.tsv",
        "assignments": args.reference / "reference_tymp_v1_cell_assignments.tsv.gz",
        "markers": args.reference / "reference_tymp_v1_marker_ranking.tsv",
        "internal": args.reference / "reference_tymp_v1_internal_validation.tsv",
        "excluded": args.reference / "reference_tymp_v1_excluded_malignant_patients.tsv",
        "metadata": args.reference / "run_metadata_tymp_v1.json",
        "manifest": args.reference / "file_manifest_tymp_v1.tsv",
    }
    missing = [str(p) for p in files.values() if not p.is_file()]
    if missing:
        raise RuntimeError(f"Missing reference deliverables: {missing}")

    checks: list[dict[str, object]] = []
    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({"check": name, "pass": bool(passed), "observed": observed, "expected": expected})

    manifest = pd.read_csv(files["manifest"], sep="\t")
    manifest_ok = True
    for row in manifest.itertuples(index=False):
        path = args.reference / row.file
        ok = path.is_file() and int(path.stat().st_size) == int(row.size_bytes) and sha256(path) == row.sha256
        manifest_ok &= ok
    check("manifest_hashes", manifest_ok, int(manifest_ok), 1)

    signature = pd.read_csv(files["signature"], sep="\t", index_col=0)
    signature.index = signature.index.astype(str)
    check("signature_shape", signature.shape == (12758, 14), f"{signature.shape[0]}x{signature.shape[1]}", "12758x14")
    check("state_names", sorted(signature.columns) == EXPECTED_STATES, "|".join(sorted(signature.columns)), "|".join(EXPECTED_STATES))
    values = signature.to_numpy()
    check("signature_finite", np.isfinite(values).all(), int(np.isfinite(values).all()), 1)
    check("signature_nonnegative", (values >= 0).all(), float(values.min()), ">=0")

    old = pd.read_csv(args.old_signature, sep="\t", index_col=0)
    old.index = old.index.astype(str)
    nonmalignant = sorted(set(EXPECTED_STATES).difference({"TYMP-high Malignant", "TYMP-low Malignant"}))
    delta = np.abs(signature[nonmalignant].to_numpy() - old[nonmalignant].to_numpy())
    max_delta = float(delta.max())
    check("nonmalignant_signatures_exact", np.allclose(delta, 0, rtol=1e-10, atol=1e-12), max_delta, "<=1e-12")

    assignments = pd.read_csv(files["assignments"], sep="\t", index_col="exact_cell_id")
    assignments.index = assignments.index.astype(str)
    check("assignment_ids_unique", not assignments.index.has_duplicates, int(assignments.index.duplicated().sum()), 0)
    check("assignment_state_names", sorted(assignments["state_tymp_v1"].unique()) == EXPECTED_STATES,
          assignments["state_tymp_v1"].nunique(), 14)
    qc = pd.read_csv(files["state_qc"], sep="\t").set_index("state_tymp_v1")
    observed_counts = assignments["state_tymp_v1"].value_counts().sort_index()
    qc_counts = qc["n_cells"].astype(int).sort_index()
    check("state_cell_counts_conserved", observed_counts.equals(qc_counts), int((observed_counts - qc_counts).abs().sum()), 0)

    bridge = pd.read_csv(args.bridge, sep="\t", index_col="exact_cell_id")
    bridge.index = bridge.index.astype(str)
    thresholds = pd.read_csv(args.thresholds, sep="\t")
    positive_patients = set(thresholds.loc[thresholds["tymp_q75"] > 0, "patient"].astype(str))
    malignant = assignments[assignments["state_tymp_v1"].isin(["TYMP-high Malignant", "TYMP-low Malignant"])]
    check("positive_q75_patient_count", len(positive_patients) == 17, len(positive_patients), 17)
    check("malignant_only_positive_q75_patients", set(malignant["patient"]) == positive_patients,
          len(set(malignant["patient"])), 17)
    check("malignant_cells_match_frozen_bridge", set(malignant.index) == set(bridge.index[bridge["patient"].isin(positive_patients)]),
          len(malignant), int(bridge["patient"].isin(positive_patients).sum()))
    expected_high = malignant["TYMP_log1p_CP10K"] >= malignant["TYMP_q75_patient"]
    observed_high = malignant["state_tymp_v1"].eq("TYMP-high Malignant")
    check("malignant_labels_match_patient_q75", np.array_equal(expected_high, observed_high),
          int((expected_high != observed_high).sum()), 0)
    check("no_zero_expression_in_tymp_high", (malignant.loc[observed_high, "TYMP_raw_count"] > 0).all(),
          int((malignant.loc[observed_high, "TYMP_raw_count"] <= 0).sum()), 0)
    state_counts = malignant["state_tymp_v1"].value_counts()
    check("tymp_high_cells_ge_500", int(state_counts["TYMP-high Malignant"]) >= 500,
          int(state_counts["TYMP-high Malignant"]), ">=500")
    check("tymp_low_cells_ge_500", int(state_counts["TYMP-low Malignant"]) >= 500,
          int(state_counts["TYMP-low Malignant"]), ">=500")
    per_state_patients = malignant.groupby("state_tymp_v1", observed=True)["patient"].nunique()
    check("both_malignant_states_17_patients", (per_state_patients == 17).all(),
          "|".join(f"{k}:{v}" for k, v in per_state_patients.items()), "both=17")

    high = signature["TYMP-high Malignant"]
    low = signature["TYMP-low Malignant"]
    check("tymp_signature_high_gt_low", high.loc["TYMP"] > low.loc["TYMP"],
          f"{high.loc['TYMP']:.9g}>{low.loc['TYMP']:.9g}", "high>low")
    rho = float(spearmanr(high.to_numpy(), low.to_numpy()).statistic)
    check("malignant_signature_rho_lt_0.98", rho < 0.98, rho, "<0.98")
    excluded_genes = signature.index.str.upper().str.startswith(("MT-", "RPL", "RPS", "HBA", "HBB"))
    log2fc = np.log2((high + 0.01) / (low + 0.01))
    n_markers = int(((np.abs(log2fc) >= 0.5) & ~excluded_genes).sum())
    check("malignant_markers_ge_20", n_markers >= 20, n_markers, ">=20")

    # Explicitly document that F1 and TYMP labels are not interchangeable.
    f1_q75 = bridge.loc[malignant.index, "F1_q75_patient"]
    f1_high = malignant["F1"] >= f1_q75
    concordance = float((f1_high.to_numpy() == observed_high.to_numpy()).mean())
    check("f1_not_used_as_tymp_label", concordance < 0.90, concordance, "<0.90 concordance")

    result = pd.DataFrame(checks)
    passed = bool(result["pass"].all())
    write_tsv(result, work / "reference_tymp_v1_independent_checks.tsv")
    status = {
        "schema_version": "phase17_reference_tymp_independent_validation_v1",
        "status": "PASS" if passed else "FAIL",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_checks": int(len(result)),
        "n_pass": int(result["pass"].sum()),
        "n_fail": int((~result["pass"]).sum()),
        "reference": str(args.reference),
        "reference_manifest_sha256": sha256(files["manifest"]),
        "reference_signature_sha256": sha256(files["signature"]),
        "validation_metrics": {
            "tymp_high_cells": int(state_counts["TYMP-high Malignant"]),
            "tymp_low_cells": int(state_counts["TYMP-low Malignant"]),
            "tymp_high_mean_count": float(high.loc["TYMP"]),
            "tymp_low_mean_count": float(low.loc["TYMP"]),
            "malignant_signature_spearman_rho": rho,
            "n_abs_log2fc_ge_0.5_markers": n_markers,
            "f1_tymp_binary_concordance": concordance,
            "max_nonmalignant_signature_abs_delta": max_delta,
        },
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    part = work / "PHASE17_REFERENCE_TYMP_INDEPENDENT_VALIDATION_v1.json.part"
    part.write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, work / "PHASE17_REFERENCE_TYMP_INDEPENDENT_VALIDATION_v1.json")
    os.replace(work, args.output)
    print(json.dumps(status, ensure_ascii=False, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
