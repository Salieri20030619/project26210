#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--score", required=True, type=Path)
    parser.add_argument("--ctx", required=True, type=Path)
    parser.add_argument("--analysis", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(part)
    part.mkdir(parents=True)
    checks: list[dict[str, str]] = []

    def add(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "status": "PASS" if bool(passed) else "FAIL", "detail": str(detail)})

    prep = pd.read_csv(args.prepared / "hif1a_adjacency_summary_v1.tsv", sep="\t")
    add("prepared_two_methods", len(prep) == 2 and set(prep["method"]) == {"GENIE3", "GRNBoost2"}, len(prep))
    expected_edges = {"GENIE3": 2632, "GRNBoost2": 5298}
    add("prepared_edge_counts", all(int(r.n_hif1a_edges) == expected_edges[r.method] for r in prep.itertuples()), prep[["method", "n_hif1a_edges"]].to_dict("records"))
    add("prepared_unique_positive_edges", bool((prep["n_hif1a_edges"] == prep["n_unique_targets"]).all() and (prep["min_importance"] > 0).all()), "unique targets and positive importance")

    score_status = json.loads((args.score / "status.json").read_text(encoding="utf-8"))
    add("score_status", score_status.get("status") == "PASS", score_status.get("status"))
    add("score_scope", score_status.get("n_cells") == 11950 and score_status.get("signature_genes") == ["CA9", "VEGFA", "ADM", "ANGPTL4", "LDHA"], score_status)
    add("hif1a_detection_reproduced", score_status.get("HIF1A_detected_cells") == 2539, score_status.get("HIF1A_detected_cells"))
    scores = pd.read_csv(args.score / "cell_hif1a_downstream_score_v1.tsv.gz", sep="\t")
    add("cell_score_keys", len(scores) == 11950 and scores["exact_cell_id"].is_unique, len(scores))
    numeric = scores[["hif1a_downstream_score", "HIF1A_raw_count", "HIF1A_log1p_cp10k", "library_size"]].to_numpy(float)
    add("cell_score_numeric_safety", np.isfinite(numeric).all() and (scores["HIF1A_raw_count"] >= 0).all() and (scores["library_size"] > 0).all(), numeric.shape)

    ctx = pd.read_csv(args.ctx / "run_summary_v1.tsv", sep="\t")
    add("ctx_eight_runs", len(ctx) == 8 and set(ctx["method"]) == {"GENIE3", "GRNBoost2"} and ctx["variant"].nunique() == 4, len(ctx))
    accepted_ctx = {"PASS_REGULON", "PASS_NO_REGULON"}
    add("ctx_all_completed", bool(ctx["ctx_status"].isin(accepted_ctx).all()), ctx["ctx_status"].value_counts().to_dict())
    present_files = list((args.ctx / "regulons").glob("*_HIF1A_PRESENT.json"))
    add("ctx_presence_markers_match", int(ctx["hif1a_positive_present"].astype(str).str.lower().eq("true").sum()) == len(present_files), len(present_files))
    auc_files = list((args.ctx / "auc").glob("*_auc.csv"))
    add("ctx_auc_matches_rescues", len(auc_files) == len(present_files) and int(ctx["auc_status"].eq("PASS").sum()) == len(auc_files), len(auc_files))
    auc_shapes_ok = True
    for path in auc_files:
        x = pd.read_csv(path)
        auc_shapes_ok &= len(x) == 11950 and any(str(c) == "HIF1A(+)" for c in x.columns)
    add("targeted_auc_shape", auc_shapes_ok, len(auc_files))

    analysis_status = json.loads((args.analysis / "status.json").read_text(encoding="utf-8"))
    add("analysis_status", analysis_status.get("status") == "PASS_ANALYSIS_COMPLETE", analysis_status.get("status"))
    add("analysis_common_scope", analysis_status.get("n_common_cells") == 11921 and analysis_status.get("n_patients") == 20, analysis_status)
    add("analysis_auc_count", analysis_status.get("n_targeted_auc_profiles") == len(auc_files), analysis_status.get("n_targeted_auc_profiles"))
    summary = pd.read_csv(args.analysis / "hif1a_sensitivity_patient_summary_v1.tsv", sep="\t")
    effects = pd.read_csv(args.analysis / "patient_hif1a_sensitivity_effects_v1.tsv", sep="\t")
    add("independent_modalities_present", {"metaVIPER_HIF1A", "frozen_downstream_5gene"}.issubset(set(summary["modality"])), sorted(summary["modality"].unique()))
    expected_rows = int(summary["modality"].nunique()) * 20 * 2
    add("patient_effect_row_count", len(effects) == expected_rows, f"{len(effects)}/{expected_rows}")
    recompute_ok = True
    for row in summary.itertuples(index=False):
        vals = effects.loc[(effects["modality"] == row.modality) & (effects["endpoint"] == row.endpoint), "effect"].dropna().to_numpy(float)
        if len(vals) != int(row.n_eligible_patients) or not np.isclose(np.median(vals), float(row.median_patient_effect), atol=1e-12, rtol=1e-10):
            recompute_ok = False
        if len(vals) >= 3:
            p = 1.0 if np.all(vals == 0) else float(wilcoxon(vals, zero_method="wilcox", alternative="two-sided", method="auto").pvalue)
            if not np.isclose(p, float(row.wilcoxon_p_two_sided), atol=1e-8, rtol=1e-6):
                recompute_ok = False
    add("patient_summary_recomputed", recompute_ok, f"rows={len(summary)}")
    conc = pd.read_csv(args.analysis / "metaviper_downstream_patient_concordance_v1.tsv", sep="\t")
    add("independent_concordance_scope", len(conc) == 20 and conc["curated_patient_id"].is_unique and np.isfinite(conc["spearman_rho"]).all(), len(conc))
    overlap = pd.read_csv(args.analysis / "cross_method_target_overlap_v1.tsv", sep="\t")
    add("target_overlap_complete", len(overlap) == 4 and set(overlap["variant"]) == {"default", "min10_only", "nes2_5_only", "combined_light"}, len(overlap))

    check_frame = pd.DataFrame(checks)
    check_frame.to_csv(part / "verification_checks_v1.tsv", sep="\t", index=False)
    failed = check_frame.loc[check_frame["status"] != "PASS", "check"].tolist()
    status = {"status": "PASS_INDEPENDENT_VALIDATION" if not failed else "FAIL_INDEPENDENT_VALIDATION", "checks_passed": int((check_frame["status"] == "PASS").sum()), "checks_total": len(check_frame), "failed_checks": failed}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    manifest_rows = []
    for root in (args.prepared, args.score, args.ctx, args.analysis):
        for path in sorted(root.rglob("*")):
            if path.is_file():
                manifest_rows.append({"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(part / "VALIDATION_MANIFEST_v1.tsv", sep="\t", index=False)
    os.replace(part, args.output)
    if failed:
        raise SystemExit("validation failures: " + ", ".join(failed))


if __name__ == "__main__":
    main()
