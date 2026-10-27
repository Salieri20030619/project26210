#!/usr/bin/env python3
"""Independent, output-only validation for Phase 35 C01--C04."""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parent
SYNC = ROOT / "server_results_sync_v1"
OUT = ROOT / "independent_validation_local_v1"


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def signed_rank(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 3 or np.all(values == 0):
        return 1.0
    return float(wilcoxon(values, alternative="two-sided", zero_method="wilcox",
                          correction=False, method="auto").pvalue)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


checks = []


def check(name, passed, detail):
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})


def close(a, b, atol=1e-10):
    return np.allclose(np.asarray(a, float), np.asarray(b, float), rtol=1e-8,
                       atol=atol, equal_nan=True)


def main():
    expected_dirs = ["c01_cellchat_integrated_v1", "c02_c04_analysis_v1", "c03_nichenet_v1"]
    check("synced result directories present", all((SYNC / d).is_dir() for d in expected_dirs), expected_dirs)
    statuses = {d: json.loads((SYNC / d / "status.json").read_text(encoding="utf-8")) for d in expected_dirs}
    check("all component status markers PASS", all(x["status"].startswith("PASS") for x in statuses.values()),
          {k: v["status"] for k, v in statuses.items()})
    check("two-sided and no prespecified direction",
          all(x.get("two_sided_testing") is True and x.get("direction_prespecified") is False for x in statuses.values()),
          "all three components")
    check("no figures generated", all(x.get("figures_created") is False for x in statuses.values()), "all three components")

    c01 = SYNC / "c01_cellchat_integrated_v1"
    effects = pd.read_csv(c01 / "c01_cellchat_patient_endpoint_effects_v1.tsv", sep="\t")
    summary = pd.read_csv(c01 / "c01_cellchat_endpoint_summary_v1.tsv", sep="\t")
    qc = pd.read_csv(c01 / "c01_cellchat_patient_run_qc_v1.tsv", sep="\t")
    interactions = pd.read_csv(c01 / "c01_cellchat_patient_interactions_v1.tsv.gz", sep="\t")
    check("C01 has 14 unique patients", effects.patient.nunique() == 14 and qc.patient.nunique() == 14, effects.patient.nunique())
    check("C01 has exactly four endpoints per patient", len(effects) == 56 and effects.groupby("patient").endpoint.nunique().eq(4).all(), len(effects))
    check("C01 patient runs all PASS", len(qc) == 14 and qc.status.astype(str).str.startswith("PASS").all(), qc.status.value_counts().to_dict())
    check("C01 interaction conservation", len(interactions) == 45266 and interactions.patient.nunique() == 14,
          f"rows={len(interactions)}, patients={interactions.patient.nunique()}")
    rebuilt = []
    for endpoint, frame in effects.groupby("endpoint", sort=True):
        v = frame.loc[frame.eligible.astype(bool), "high_minus_low"].to_numpy(float)
        rebuilt.append([endpoint, len(v), (v > 0).sum(), (v < 0).sum(), (v == 0).sum(), np.median(v), signed_rank(v)])
    rb = pd.DataFrame(rebuilt, columns=["endpoint", "n", "pos", "neg", "zero", "median", "p"])
    rb["q"] = bh(rb.p)
    s = summary.sort_values("endpoint").reset_index(drop=True)
    check("C01 endpoint identities", list(s.endpoint) == list(rb.endpoint), list(s.endpoint))
    check("C01 counts independently reproduced",
          np.array_equal(s.n_eligible_patients, rb.n) and np.array_equal(s.n_positive, rb.pos) and
          np.array_equal(s.n_negative, rb.neg) and np.array_equal(s.n_zero, rb.zero), "4/4 endpoints")
    check("C01 medians independently reproduced", close(s.median_high_minus_low, rb["median"]), "4/4 endpoints")
    check("C01 two-sided Wilcoxon independently reproduced", close(s.wilcoxon_p_two_sided, rb.p), "4/4 endpoints")
    check("C01 four-endpoint BH independently reproduced", close(s.BH_FDR_4_endpoints, rb.q), "4/4 endpoints")
    check("C01 support flags follow FDR", np.array_equal(s.statistical_support.astype(bool), rb.q < .05), "4/4 endpoints")

    c03 = SYNC / "c03_nichenet_v1"
    p03 = pd.read_csv(c03 / "c03_nichenet_patient_ligand_activity_v1.tsv", sep="\t")
    s03 = pd.read_csv(c03 / "c03_nichenet_ligand_activity_summary_v1.tsv", sep="\t")
    weights = pd.read_csv(c03 / "c03_nichenet_direct_target_weights_v1.tsv", sep="\t")
    check("C03 has 14 patients and two ligands", p03.patient.nunique() == 14 and set(p03.ligand) == {"ICAM5", "MST1"}, p03.shape)
    check("C03 one row per patient-ligand", len(p03) == 28 and not p03.duplicated(["patient", "ligand"]).any(), len(p03))
    check("C03 direct weights cover two ligands and seven targets", len(weights) == 14 and weights.ligand.nunique() == 2 and weights.target.nunique() == 7, weights.shape)
    endpoint_map = {"complete_target_pearson": "complete_target_pearson",
                    "complete_target_spearman": "complete_target_spearman",
                    "hif1a_axis_weighted_response": "hif1a_axis_weighted_response"}
    c03_ok = []
    for endpoint, col in endpoint_map.items():
        sub = s03[s03.endpoint == endpoint].sort_values("ligand")
        vals, ps = [], []
        for ligand in sorted(p03.ligand.unique()):
            v = p03.loc[p03.ligand == ligand, col].dropna().to_numpy(float)
            vals.append(v); ps.append(signed_rank(v))
        q = bh(ps)
        c03_ok.append(len(sub) == 2 and np.array_equal(sub.n_eligible_patients, [len(x) for x in vals]) and
                      close(sub.median_patient_effect, [np.median(x) for x in vals]) and
                      close(sub.wilcoxon_p_two_sided, ps) and close(sub.BH_FDR_2_ligands, q))
    check("C03 patient-level summaries independently reproduced", all(c03_ok), "6/6 ligand-endpoints")
    expected_support = np.where(s03.n_eligible_patients >= 10, s03.BH_FDR_2_ligands < .05,
                                np.where(s03.n_eligible_patients >= 3, s03.wilcoxon_p_two_sided < .05, False))
    check("C03 support policy reproduced", np.array_equal(s03.statistical_support.astype(bool), expected_support), "6/6")

    c24 = SYNC / "c02_c04_analysis_v1"
    avail = pd.read_csv(c24 / "c02_patient_state_receptor_availability_v1.tsv", sep="\t")
    a_sum = pd.read_csv(c24 / "c02_receptor_specificity_summary_v1.tsv", sep="\t")
    check("C02 patient-state-feature grid complete", len(avail) == 960 and avail.patient.nunique() == 20 and
          avail.communication_group.nunique() == 6 and avail.receptor_feature.nunique() == 8, avail.shape)
    check("C02 summaries cover 48 group-feature combinations", len(a_sum) == 48 and not a_sum.duplicated(["communication_group", "receptor_feature"]).any(), a_sum.shape)
    check("C02 sensitivity mapping is explicitly labeled",
          "ITGAM_ITGB2_sensitivity" in set(avail.receptor_feature) and
          not any((avail.receptor_feature == "ITGAM_ITGB2_sensitivity") & avail.receptor_feature.str.contains("primary")),
          "ITGAM/ITGB2 is not represented as the primary ICAM5 receptor")
    eff = pd.read_csv(c24 / "c04_patient_receptor_tymp_f1_effects_v1.tsv", sep="\t")
    s04 = pd.read_csv(c24 / "c04_receptor_tymp_f1_summary_v1.tsv", sep="\t")
    check("C04 covers 20 patients, two outcomes and eight features", eff.patient.nunique() == 20 and
          eff.outcome.nunique() == 2 and eff.receptor_feature.nunique() == 8 and len(eff) == 320, eff.shape)
    rebuilt04 = []
    for (outcome, feature), frame in eff.groupby(["outcome", "receptor_feature"], sort=True):
        v = frame.loc[frame.eligible.astype(bool), "patient_spearman_rho"].dropna().to_numpy(float)
        rebuilt04.append([outcome, feature, len(v), (v > 0).sum(), (v < 0).sum(), (v == 0).sum(),
                          np.median(v) if len(v) else np.nan, signed_rank(v)])
    r04 = pd.DataFrame(rebuilt04, columns=["outcome", "receptor_feature", "n", "pos", "neg", "zero", "median", "p"])
    r04["q"] = r04.groupby("outcome", sort=False).p.transform(bh)
    s04 = s04.sort_values(["outcome", "receptor_feature"]).reset_index(drop=True)
    check("C04 row identities", list(zip(s04.outcome, s04.receptor_feature)) == list(zip(r04.outcome, r04.receptor_feature)), "16/16")
    check("C04 counts and medians independently reproduced",
          np.array_equal(s04.n_eligible_patients, r04.n) and close(s04.median_patient_rho, r04["median"]), "16/16")
    check("C04 Wilcoxon and BH independently reproduced", close(s04.wilcoxon_p_two_sided, r04.p) and
          close(s04.BH_FDR_within_outcome8features, r04.q), "16/16")
    expected04 = np.where(s04.n_eligible_patients >= 10, s04.BH_FDR_within_outcome8features < .05,
                          np.where(s04.n_eligible_patients >= 3, s04.wilcoxon_p_two_sided < .05, False))
    check("C04 small-n reporting policy reproduced", np.array_equal(s04.statistical_support.astype(bool), expected04), "16/16")

    figure_ext = {".png", ".jpg", ".jpeg", ".pdf", ".svg", ".tif", ".tiff"}
    figures = [str(p) for p in SYNC.rglob("*") if p.suffix.lower() in figure_ext]
    check("filesystem contains no figures", not figures, figures)

    OUT.mkdir(parents=False, exist_ok=False)
    frame = pd.DataFrame(checks)
    frame.to_csv(OUT / "phase35_independent_validation_checks_v1.tsv", sep="\t", index=False)
    passed = int(frame.passed.sum())
    payload = {"status": "PASS_PHASE35_C01_C04_INDEPENDENT_VALIDATION" if passed == len(frame) else "FAIL_PHASE35_C01_C04_INDEPENDENT_VALIDATION",
               "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
               "checks_passed": passed, "checks_total": len(frame),
               "patient_is_inferential_unit": True, "two_sided_testing": True,
               "direction_prespecified": False, "figures_created": False,
               "validated_files": {str(p.relative_to(ROOT)): sha256(p) for p in sorted(SYNC.rglob("*")) if p.is_file()}}
    (OUT / "status.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if passed == len(frame) else 1


if __name__ == "__main__":
    sys.exit(main())
