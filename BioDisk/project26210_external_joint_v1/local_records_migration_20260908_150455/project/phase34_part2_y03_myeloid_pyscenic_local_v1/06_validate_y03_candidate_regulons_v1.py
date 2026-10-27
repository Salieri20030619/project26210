#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


CANDIDATES = {"JUN", "JUNB", "CEBPD"}
ENDPOINTS = {
    "Monocle3_continuous",
    "DPT_continuous",
    "SPP1_minus_Monocyte",
    "SPP1_minus_C1QC",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    part = Path(str(path) + ".part")
    frame.to_csv(part, sep="\t", index=False, lineterminator="\n")
    os.replace(part, path)


def atomic_json(payload: dict, path: Path) -> None:
    part = Path(str(path) + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(part, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--run-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, observed: object, expected: object) -> None:
        checks.append({
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": str(observed),
            "expected": str(expected),
        })

    core_files = {
        "adjacency": args.core / "adjacencies_grnboost2_v1.tsv",
        "auc": args.core / "auc_grnboost2_all_myeloid_v1.csv",
        "regulons": args.core / "regulons_grnboost2_v1.dat",
        "registry": args.core / "regulon_registry_v1.tsv",
        "marker": args.core / "RUN_COMPLETE",
    }
    for label, path in core_files.items():
        check(f"core_file_{label}", path.is_file() and path.stat().st_size > 0,
              path.stat().st_size if path.is_file() else "missing", ">0 bytes")
    marker_text = core_files["marker"].read_text(encoding="utf-8").strip() if core_files["marker"].is_file() else ""
    try:
        marker_time_valid = datetime.fromisoformat(marker_text).tzinfo is not None
    except ValueError:
        marker_time_valid = False
    check("core_completion_marker", marker_time_valid, marker_text, "timezone-aware ISO completion timestamp")
    log = args.run_metadata / "logs" / "pyscenic_core_v1.log"
    log_text = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
    check("core_log_terminal_pass", "Y03_PYSCENIC_CORE_PASS" in log_text,
          "Y03_PYSCENIC_CORE_PASS" in log_text, True)

    input_status = json.loads((args.run_metadata / "input_status_v1.json").read_text(encoding="utf-8"))
    loom_status = json.loads((args.run_metadata / "loom_status_v1.json").read_text(encoding="utf-8"))
    check("input_status", input_status.get("status") == "PASS_Y03_PYSCENIC_INPUT_PREPARATION",
          input_status.get("status"), "PASS_Y03_PYSCENIC_INPUT_PREPARATION")
    check("input_dimensions", input_status.get("n_cells") == 50761 and input_status.get("n_patients") == 40,
          f"cells={input_status.get('n_cells')};patients={input_status.get('n_patients')}", "cells=50761;patients=40")
    check("loom_status", loom_status.get("status") == "PASS_Y03_PYSCENIC_LOOM_BUILD",
          loom_status.get("status"), "PASS_Y03_PYSCENIC_LOOM_BUILD")

    registry = pd.read_csv(args.analysis / "y03_candidate_regulon_registry_v1.tsv", sep="\t")
    check("candidate_registry_exact", set(registry["tf"]) == CANDIDATES and len(registry) == 3,
          sorted(registry["tf"].tolist()), sorted(CANDIDATES))
    check("all_candidates_available", registry["available"].astype(bool).all(),
          int(registry["available"].astype(bool).sum()), 3)
    check("candidate_targets_positive", (registry["n_targets"] > 0).all(),
          registry.set_index("tf")["n_targets"].to_dict(), "all > 0")

    cells = pd.read_csv(args.analysis / "y03_candidate_regulon_cell_scores_v1.tsv.gz", sep="\t", dtype={"cell_id": str})
    score_cols = {f"{tf}_AUCell_v1" for tf in CANDIDATES}
    check("cell_score_dimensions", len(cells) == 50761 and not cells["cell_id"].duplicated().any(),
          f"rows={len(cells)};unique={cells['cell_id'].nunique()}", "rows=unique=50761")
    check("cell_score_columns", score_cols.issubset(cells.columns), sorted(set(cells.columns) & score_cols), sorted(score_cols))
    check("cell_scores_finite", cells[list(score_cols)].notna().all().all(),
          int(cells[list(score_cols)].isna().sum().sum()), 0)

    effects = pd.read_csv(args.analysis / "y03_candidate_patient_effects_v1.tsv", sep="\t")
    check("patient_effect_grid", len(effects) == 40 * 3 * 4,
          len(effects), 40 * 3 * 4)
    check("patient_effect_unique", not effects.duplicated(["study", "patient", "tf", "endpoint"]).any(),
          int(effects.duplicated(["study", "patient", "tf", "endpoint"]).sum()), 0)
    check("patient_effect_candidates", set(effects["tf"]) == CANDIDATES, sorted(effects["tf"].unique()), sorted(CANDIDATES))
    check("patient_effect_endpoints", set(effects["endpoint"]) == ENDPOINTS, sorted(effects["endpoint"].unique()), sorted(ENDPOINTS))

    summary = pd.read_csv(args.analysis / "y03_candidate_endpoint_summary_v1.tsv", sep="\t")
    check("endpoint_summary_grid", len(summary) == 12 and not summary.duplicated(["endpoint", "tf"]).any(),
          len(summary), 12)
    check("endpoint_summary_patient_unit", summary["n_eligible_patients"].between(3, 40).all(),
          sorted(summary["n_eligible_patients"].unique().tolist()), "3..40 patients")
    check("reporting_rule", (summary["reporting_tier"] == "FDR_AT_LEAST_10").all(),
          sorted(summary["reporting_tier"].unique().tolist()), "FDR_AT_LEAST_10")
    recomputed_support = summary["BH_FDR_within_endpoint3TF"].lt(0.05)
    check("support_matches_fdr", (summary["statistical_support"].astype(bool) == recomputed_support).all(),
          int((summary["statistical_support"].astype(bool) != recomputed_support).sum()), 0)
    recomputed_direction = summary["median_patient_effect"].map(lambda value: "positive" if value > 0 else ("negative" if value < 0 else "zero"))
    check("direction_matches_effect", (summary["direction"] == recomputed_direction).all(),
          int((summary["direction"] != recomputed_direction).sum()), 0)

    integrated = pd.read_csv(args.analysis / "y03_candidate_integrated_evidence_v1.tsv", sep="\t")
    supported = summary.loc[summary["statistical_support"].astype(bool)]
    by_tf = supported.groupby("tf").agg(
        n=("endpoint", "size"),
        pos=("direction", lambda x: int((x == "positive").sum())),
        neg=("direction", lambda x: int((x == "negative").sum())),
    ).reindex(sorted(CANDIDATES), fill_value=0)
    observed = integrated.set_index("tf").sort_index()
    integration_ok = (
        observed["n_supported_endpoints_of_4"].equals(by_tf["n"])
        and observed["n_positive_supported"].equals(by_tf["pos"])
        and observed["n_negative_supported"].equals(by_tf["neg"])
    )
    check("integrated_evidence_recomputed", integration_ok,
          observed[["n_supported_endpoints_of_4", "n_positive_supported", "n_negative_supported"]].to_dict("index"),
          by_tf.to_dict("index"))

    analysis_status = json.loads((args.analysis / "status.json").read_text(encoding="utf-8"))
    check("two_sided_no_prespecified_direction",
          analysis_status.get("two_sided_testing") is True and analysis_status.get("direction_prespecified") is False,
          f"two_sided={analysis_status.get('two_sided_testing')};prespecified={analysis_status.get('direction_prespecified')}",
          "two_sided=True;prespecified=False")
    figure_files = [p.name for p in args.analysis.rglob("*") if p.suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg"}]
    check("no_figures_before_review", not figure_files and analysis_status.get("figures_created") is False,
          figure_files, [])
    interpretation = str(analysis_status.get("interpretation", "")).lower()
    check("interpretation_guard", "not direct tf binding" in interpretation and "causality" in interpretation,
          analysis_status.get("interpretation"), "explicit non-binding/non-causal guard")

    check_frame = pd.DataFrame(checks)
    atomic_tsv(check_frame, args.output / "y03_independent_validation_checks_v1.tsv")
    failures = check_frame.loc[check_frame["status"].eq("FAIL"), "check"].tolist()
    payload = {
        "status": "PASS_Y03_INDEPENDENT_VALIDATION" if not failures else "FAIL_Y03_INDEPENDENT_VALIDATION",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_checks": len(check_frame),
        "n_pass": int(check_frame["status"].eq("PASS").sum()),
        "n_fail": len(failures),
        "failures": failures,
        "core_sha256": {label: sha256(path) for label, path in core_files.items()},
        "analysis_status_sha256": sha256(args.analysis / "status.json"),
    }
    atomic_json(payload, args.output / "status.json")
    if failures:
        raise RuntimeError(f"Independent validation failed: {failures}")


if __name__ == "__main__":
    main()
