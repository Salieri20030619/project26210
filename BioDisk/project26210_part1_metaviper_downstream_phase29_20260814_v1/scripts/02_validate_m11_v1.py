#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, wilcoxon


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists(): raise FileExistsError(part)
    part.mkdir(parents=True)
    checks = []
    def add(name, ok, detail): checks.append({"check": name, "status": "PASS" if bool(ok) else "FAIL", "detail": str(detail)})

    status = json.loads((args.analysis / "status.json").read_text(encoding="utf-8"))
    add("analysis_status", status.get("status") == "PASS_M11_ANALYSIS_COMPLETE", status.get("status"))
    add("scope", status.get("n_cells") == 11921 and status.get("n_patients") == 20, status)
    targets = pd.read_csv(args.analysis / "m11_final_mr_target_inference_status_v1.tsv", sep="\t")
    add("eight_final_targets", len(targets) == 8 and targets["protein"].is_unique, len(targets))
    infer = targets.loc[targets["activity_row_present"].astype(bool)]
    add("four_inferable_targets", set(infer["protein"]) == {"CASP9", "CHMP1A", "PRSS8", "PXN"}, sorted(infer["protein"]))
    add("beta_reproduced", np.allclose(targets["mr_beta_log_or"], np.log(targets["odds_ratio"]), atol=1e-14), "log(OR)")
    means = pd.read_csv(args.analysis / "m11_patient_mean_activity_v1.tsv", sep="\t")
    add("patient_mean_rows", len(means) == 80 and means.groupby("protein")["curated_patient_id"].nunique().eq(20).all(), len(means))
    points = pd.read_csv(args.analysis / "m11_mr_effect_activity_points_v1.tsv", sep="\t")
    reproduced = means.groupby("protein", as_index=False)["mean_activity"].mean().rename(columns={"mean_activity":"reproduced"}).merge(points, on="protein")
    add("patient_equal_mean_reproduced", np.allclose(reproduced["reproduced"], reproduced["patient_equal_mean_activity"], atol=1e-12), len(reproduced))
    glob = pd.read_csv(args.analysis / "m11_global_correlation_v1.tsv", sep="\t").iloc[0]
    pr = pearsonr(points["mr_beta_log_or"], points["patient_equal_mean_activity"])
    sr = spearmanr(points["mr_beta_log_or"], points["patient_equal_mean_activity"])
    add("global_correlations_reproduced", np.isclose(glob.pearson_r, pr.statistic, atol=1e-12) and np.isclose(glob.spearman_rho, sr.statistic, atol=1e-12), (pr.statistic, sr.statistic))
    pc = pd.read_csv(args.analysis / "m11_patient_cross_protein_correlations_v1.tsv", sep="\t")
    add("patient_correlation_scope", len(pc) == 20 and pc["curated_patient_id"].is_unique and np.isfinite(pc[["pearson_r","spearman_rho"]]).all().all(), len(pc))
    ps = pd.read_csv(args.analysis / "m11_patient_correlation_summary_v1.tsv", sep="\t")
    summary_ok = True
    for row in ps.itertuples(index=False):
        vals = pc[row.metric].to_numpy(float)
        nonzero = vals[vals != 0]
        has_tied_absolute_ranks = len(np.unique(np.abs(nonzero))) < len(nonzero)
        if (vals == 0).any() or has_tied_absolute_ranks:
            pv = wilcoxon(vals, zero_method="wilcox", alternative="two-sided", method="approx", correction=True).pvalue
        else:
            pv = wilcoxon(vals, zero_method="wilcox", alternative="two-sided", method="exact").pvalue
        summary_ok &= row.n_patients == 20 and np.isclose(row.median_patient_correlation, np.median(vals), atol=1e-12) and np.isclose(row.wilcoxon_p_two_sided, pv, atol=1e-8)
    add("patient_summary_reproduced", summary_ok, len(ps))
    add("no_figures", not any(args.analysis.rglob("*.png")) and not any(args.analysis.rglob("*.pdf")), "tables only")
    frame = pd.DataFrame(checks)
    frame.to_csv(part / "verification_checks_v1.tsv", sep="\t", index=False)
    failed = frame.loc[frame.status != "PASS", "check"].tolist()
    out = {"status": "PASS_INDEPENDENT_VALIDATION" if not failed else "FAIL_INDEPENDENT_VALIDATION", "checks_passed": int((frame.status == "PASS").sum()), "checks_total": len(frame), "failed_checks": failed}
    (part / "status.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output)
    if failed: raise SystemExit("validation failures: " + ",".join(failed))


if __name__ == "__main__": main()
