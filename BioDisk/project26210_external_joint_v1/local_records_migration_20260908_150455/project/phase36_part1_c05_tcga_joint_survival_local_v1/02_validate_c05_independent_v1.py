#!/usr/bin/env python3
"""Independent output-only validation for Phase 36 C05 TCGA survival."""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "server_results_sync_v1" / "analysis_v1"
OUT = ROOT / "independent_validation_local_v1"
AXES = ["MST1R", "ICAM5_ITGAL_ITGB2_primary", "CD209", "PRTG", "IGDCC4", "ICAM5_ITGAM_ITGB2_sensitivity"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values):
    p = np.asarray(values, float); order = np.argsort(p)
    z = p[order] * len(p) / np.arange(1, len(p) + 1)
    z = np.minimum.accumulate(z[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(z, 1); return out


def logrank_k(time, event, group):
    """Mantel-Haenszel k-sample log-rank, matching survival::survdiff rho=0."""
    time = np.asarray(time, float); event = np.asarray(event, int); group = np.asarray(group)
    levels = list(pd.unique(group)); k = len(levels)
    obs = np.zeros(k); exp = np.zeros(k); var = np.zeros((k, k))
    for t in np.sort(np.unique(time[event == 1])):
        risk = np.array([np.sum((group == g) & (time >= t)) for g in levels], float)
        deaths = np.array([np.sum((group == g) & (time == t) & (event == 1)) for g in levels], float)
        n = risk.sum(); d = deaths.sum()
        if n <= 0:
            continue
        obs += deaths; exp += d * risk / n
        if n > 1:
            factor = d * (n - d) / (n - 1)
            frac = risk / n
            var += factor * (np.diag(frac) - np.outer(frac, frac))
    oe = (obs - exp)[:-1]; vv = var[:-1, :-1]
    stat = float(oe @ np.linalg.pinv(vv) @ oe)
    return stat, k - 1, float(chi2.sf(stat, k - 1))


checks = []


def check(name, passed, detail):
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


def close(a, b, atol=1e-8):
    return np.allclose(np.asarray(a, float), np.asarray(b, float), rtol=1e-7, atol=atol, equal_nan=True)


def main():
    status = json.loads((SRC / "status.json").read_text(encoding="utf-8"))
    cohort = pd.read_csv(SRC / "c05_tcga_analysis_cohort_v1.tsv", sep="\t")
    groups = pd.read_csv(SRC / "c05_joint_group_survival_summary_v1.tsv", sep="\t")
    lr = pd.read_csv(SRC / "c05_joint_group_logrank_summary_v1.tsv", sep="\t")
    cox = pd.read_csv(SRC / "c05_joint_group_cox_models_v1.tsv", sep="\t")
    inter = pd.read_csv(SRC / "c05_continuous_interaction_cox_v1.tsv", sep="\t")
    check("analysis status PASS", status["status"].startswith("PASS"), status["status"])
    check("median cutoffs only", status.get("median_cutoffs_only") is True, status.get("median_cutoffs_only"))
    check("no prespecified direction", status.get("direction_prespecified") is False, status.get("direction_prespecified"))
    check("no figures generated", status.get("figures_created") is False, status.get("figures_created"))
    check("510 unique merged patients", len(cohort) == 510 and cohort.patient.nunique() == 510, cohort.shape)
    check("all six receptor axes present and finite", all(a in cohort.columns for a in AXES) and
          np.isfinite(cohort[["TYMP"] + AXES].to_numpy(float)).all(), AXES)
    check("OS analysis population and events", (cohort.OS_time > 0).sum() == 508 and
          cohort.loc[cohort.OS_time > 0, "OS_event"].sum() == 168, "508 patients, 168 events")
    check("PFS analysis population and events", (cohort.PFS_time > 0).sum() == 506 and
          cohort.loc[cohort.PFS_time > 0, "PFS_event"].sum() == 156, "506 patients, 156 events")
    check("12 axis-endpoint log-rank tests", len(lr) == 12 and not lr.duplicated(["receptor_axis", "endpoint"]).any(), lr.shape)
    check("36 axis-endpoint-group summaries", len(groups) == 36 and not groups.duplicated(["receptor_axis", "endpoint", "joint_group"]).any(), groups.shape)

    rebuilt = []
    rebuilt_groups = []
    for endpoint in ["OS", "PFS"]:
        tcol, ecol = f"{endpoint}_time", f"{endpoint}_event"
        for axis in AXES:
            d = cohort[np.isfinite(cohort[tcol]) & (cohort[tcol] > 0) & np.isfinite(cohort[ecol]) &
                       np.isfinite(cohort.TYMP) & np.isfinite(cohort[axis])].copy()
            th = d.TYMP > d.TYMP.median(); rh = d[axis] > d[axis].median()
            d["joint_group"] = np.where(~th & ~rh, "Both_low", np.where(th & rh, "Both_high", "Single_high"))
            stat, df, p = logrank_k(d[tcol], d[ecol], d.joint_group)
            rebuilt.append([axis, endpoint, len(d), int(d[ecol].sum()), stat, df, p])
            for g in ["Both_low", "Single_high", "Both_high"]:
                z = d[d.joint_group == g]
                rebuilt_groups.append([axis, endpoint, g, len(z), int(z[ecol].sum()), z.TYMP.median(), z[axis].median()])
    rb = pd.DataFrame(rebuilt, columns=["receptor_axis", "endpoint", "n_patients", "n_events", "chisq", "df", "p"])
    rb["q"] = rb.groupby("endpoint", sort=False).p.transform(bh)
    key = ["receptor_axis", "endpoint"]
    lr2 = lr.sort_values(key).reset_index(drop=True); rb2 = rb.sort_values(key).reset_index(drop=True)
    check("log-rank population counts reproduced", np.array_equal(lr2.n_patients, rb2.n_patients) and
          np.array_equal(lr2.n_events, rb2.n_events), "12/12")
    check("log-rank chi-square and P reproduced", close(lr2.chisq, rb2.chisq, 1e-7) and
          close(lr2.logrank_p_two_sided, rb2.p, 1e-9), "12/12")
    check("endpoint-wise six-axis BH reproduced", close(lr2.BH_FDR_within_endpoint6axes, rb2.q, 1e-9), "12/12")
    check("support flags follow BH FDR", np.array_equal(lr2.statistical_support.astype(bool), rb2.q < .05), "12/12")
    rg = pd.DataFrame(rebuilt_groups, columns=["receptor_axis", "endpoint", "joint_group", "n_patients", "n_events", "median_TYMP", "median_receptor"])
    gkey = ["receptor_axis", "endpoint", "joint_group"]
    g1 = groups.sort_values(gkey).reset_index(drop=True); g2 = rg.sort_values(gkey).reset_index(drop=True)
    check("joint-group patient and event counts reproduced", np.array_equal(g1.n_patients, g2.n_patients) and
          np.array_equal(g1.n_events, g2.n_events), "36/36")
    check("joint-group expression medians reproduced", close(g1.median_TYMP, g2.median_TYMP) and
          close(g1.median_receptor, g2.median_receptor), "36/36")
    check("Cox outputs cover six axes and two endpoints", set(cox.receptor_axis) == set(AXES) and
          set(inter.receptor_axis) == set(AXES) and set(cox.endpoint) == {"OS", "PFS"} and set(inter.endpoint) == {"OS", "PFS"},
          f"group={cox.shape}, interaction={inter.shape}")
    for name, frame in [("group Cox", cox), ("interaction Cox", inter)]:
        finite = frame[["beta", "hazard_ratio", "standard_error", "z", "p_two_sided", "ci95_low", "ci95_high"]].apply(pd.to_numeric, errors="coerce")
        ok = (np.isfinite(finite.to_numpy()).all() and
              close(frame.hazard_ratio, np.exp(frame.beta)) and
              close(frame.z, frame.beta / frame.standard_error) and
              close(frame.ci95_low, np.exp(frame.beta - 1.96 * frame.standard_error), 1e-6) and
              close(frame.ci95_high, np.exp(frame.beta + 1.96 * frame.standard_error), 1e-6))
        check(f"{name} coefficient algebra", ok, f"rows={len(frame)}")
    check("primary and sensitivity ICAM5 axes remain distinct", "ICAM5_ITGAL_ITGB2_primary" in set(lr.receptor_axis) and
          "ICAM5_ITGAM_ITGB2_sensitivity" in set(lr.receptor_axis), "both explicitly labeled")
    figure_ext = {".png", ".jpg", ".jpeg", ".pdf", ".svg", ".tif", ".tiff"}
    figures = [str(p) for p in SRC.rglob("*") if p.suffix.lower() in figure_ext]
    check("filesystem contains no figures", not figures, figures)

    OUT.mkdir(parents=False, exist_ok=False)
    frame = pd.DataFrame(checks)
    frame.to_csv(OUT / "phase36_c05_independent_validation_checks_v1.tsv", sep="\t", index=False)
    passed = int(frame.passed.sum())
    payload = {"status": "PASS_PHASE36_C05_INDEPENDENT_VALIDATION" if passed == len(frame) else "FAIL_PHASE36_C05_INDEPENDENT_VALIDATION",
               "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
               "checks_passed": passed, "checks_total": len(frame),
               "median_cutoffs_only": True, "direction_prespecified": False, "figures_created": False,
               "validated_files": {str(p.relative_to(ROOT)): sha256(p) for p in sorted(SRC.iterdir()) if p.is_file()}}
    (OUT / "status.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if passed == len(frame) else 1


if __name__ == "__main__":
    sys.exit(main())
