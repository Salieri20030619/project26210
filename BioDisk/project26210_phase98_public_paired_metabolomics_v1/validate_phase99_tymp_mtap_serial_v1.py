#!/usr/bin/env python3
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MTA = "5-methylthioadenosine (MTA)"
PROGRAM = ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
DATASETS = {
    "RC18": ("matched_Harmonized_RawData_RC18.csv", "matched_tpm_RC18.csv"),
    "RC20": ("matched_Harmonized_RawData_RC20.csv", "matched_tpm_RC20.csv"),
    "CPTAC": ("matched_Harmonized_Met_CPTAC.csv", "matched_tpm_CPTAC.csv"),
    "CPTAC_val": ("matched_Harmonized_Met_CPTAC_val.csv", "matched_tpm_CPTAC_val.csv"),
}
VALIDATION_SEED = 9926210
N_BOOT = 5000


def ranks(x):
    x = stats.rankdata(np.asarray(x, dtype=float), method="average")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if sd else np.zeros_like(x)


def ols(y, *x):
    design = np.column_stack((np.ones(len(y)),) + x)
    return np.linalg.pinv(design) @ y


def fixed_meta(rhos, ns):
    z = np.arctanh(np.clip(np.asarray(rhos), -0.999999, 0.999999))
    w = np.asarray(ns, dtype=float) - 3.0
    estimate = np.sum(w * z) / np.sum(w)
    se = np.sqrt(1.0 / np.sum(w))
    return math.tanh(estimate), float(2 * stats.norm.sf(abs(estimate / se)))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: validate_phase99_tymp_mtap_serial_v1.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    raw = root / "downloads_v1"
    analysis = root / "analysis_phase99_v1"
    out = root / "analysis_phase99_validation_v1"
    out.mkdir(parents=True, exist_ok=False)
    mapping = pd.read_csv(root / "MasterMapping_ccRCC.csv", dtype=str).fillna("")
    reported_assoc = pd.read_csv(analysis / "cohort_tymp_mtap_associations_v1.tsv", sep="\t")
    reported_meta = pd.read_csv(analysis / "meta_tymp_mtap_associations_v1.tsv", sep="\t")
    reported_cohort = pd.read_csv(analysis / "cohort_indirect_effects_v1.tsv", sep="\t")
    reported_pool = pd.read_csv(analysis / "pooled_indirect_effects_v1.tsv", sep="\t")
    status = json.loads((analysis / "analysis_status_phase99_v1.json").read_text())
    checks = []

    def check(name, passed, detail=""):
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    expected = {
        "cohort_tymp_mtap_associations_v1.tsv", "meta_tymp_mtap_associations_v1.tsv",
        "cohort_indirect_effects_v1.tsv", "pooled_indirect_effects_v1.tsv",
        "sample_counts_v1.tsv", "analysis_status_phase99_v1.json",
    }
    check("six_expected_outputs", expected == {p.name for p in analysis.iterdir() if p.is_file()})
    frames = {}
    primary_rho, primary_n = [], []
    pooled_blocks = []

    for dataset, (mf, rf) in DATASETS.items():
        mapd = mapping[mapping.Dataset == dataset]
        met = pd.read_csv(raw / mf, index_col=0).loc[mapd.MetabID]
        rna = pd.read_csv(raw / rf, index_col=0).loc[mapd.RNAID]
        frame = pd.DataFrame({
            "TN": mapd.TN.to_numpy(),
            "MTA": pd.to_numeric(met[MTA], errors="coerce").to_numpy(),
            "TYMP": np.log2(pd.to_numeric(rna.TYMP, errors="coerce").to_numpy() + 1),
            "MTAP": np.log2(pd.to_numeric(rna.MTAP, errors="coerce").to_numpy() + 1),
            "SPP1": np.log2(pd.to_numeric(rna.SPP1, errors="coerce").to_numpy() + 1),
        })
        for gene in PROGRAM:
            frame[gene] = np.log2(pd.to_numeric(rna[gene], errors="coerce").to_numpy() + 1)
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
        tumor = frame[frame.TN == "Tumor"].copy()
        pm = tumor[PROGRAM]
        tumor["SPP1_PROGRAM"] = ((pm - pm.mean()) / pm.std(ddof=0)).mean(axis=1)
        frames[dataset] = tumor
        rho, p = stats.spearmanr(tumor.TYMP, tumor.MTAP)
        row = reported_assoc[(reported_assoc.dataset == dataset) & (reported_assoc.endpoint == "TYMP_MTAP_MARGINAL")].iloc[0]
        check(f"{dataset}_sample_count", len(tumor) == int(row.n), f"n={len(tumor)}")
        check(f"{dataset}_rho", abs(rho - row.rho) < 1e-12, f"recomputed={rho:.16g}")
        check(f"{dataset}_p", abs(p - row.p) < 1e-12, f"recomputed={p:.16g}")
        primary_rho.append(rho)
        primary_n.append(len(tumor))
        pooled_blocks.append(pd.DataFrame({
            "dataset": dataset, "A": ranks(tumor.TYMP), "D": ranks(tumor.MTAP),
            "B": ranks(tumor.MTA), "SPP1": ranks(tumor.SPP1), "PROGRAM": ranks(tumor.SPP1_PROGRAM),
        }))

    meta_rho, meta_p = fixed_meta(primary_rho, primary_n)
    meta_row = reported_meta[reported_meta.endpoint == "TYMP_MTAP_MARGINAL"].iloc[0]
    check("fixed_meta_rho", abs(meta_rho - meta_row.fixed_rho) < 1e-12, f"recomputed={meta_rho:.16g}")
    check("fixed_meta_p", abs(meta_p - meta_row.fixed_p) < 1e-12, f"recomputed={meta_p:.16g}")
    check("primary_directions_2_plus_2_minus", sum(x > 0 for x in primary_rho) == 2 and sum(x < 0 for x in primary_rho) == 2)

    data = pd.concat(pooled_blocks, ignore_index=True)
    dummies = pd.get_dummies(data.dataset, drop_first=True, dtype=float).to_numpy()

    def estimate(sample, outcome=None):
        aa, dd, bb = sample.A.to_numpy(), sample.D.to_numpy(), sample.B.to_numpy()
        fixed = pd.get_dummies(sample.dataset, dtype=float).reindex(columns=list(DATASETS), fill_value=0).iloc[:, 1:].to_numpy()
        a1 = ols(dd, aa, *fixed.T)[1]
        a2 = ols(bb, dd, aa, *fixed.T)[1]
        if outcome is None:
            return a1, a2, np.nan, a1 * a2
        cc = sample[outcome].to_numpy()
        b1 = ols(cc, bb, dd, aa, *fixed.T)[1]
        return a1, a2, b1, a1 * a2 * b1

    names = [("SHORT_TYMP_MTAP_MTA", None), ("SERIAL_TO_SPP1", "SPP1"), ("SERIAL_TO_PROGRAM", "PROGRAM")]
    rng = np.random.default_rng(VALIDATION_SEED)
    indices = {ds: np.flatnonzero(data.dataset.to_numpy() == ds) for ds in DATASETS}
    validation_intervals = []
    for effect_type, outcome in names:
        point = estimate(data, outcome)
        row = reported_pool[reported_pool.effect_type == effect_type].iloc[0]
        check(f"{effect_type}_point", abs(point[3] - row.indirect_effect) < 1e-12, f"recomputed={point[3]:.16g}")
        draws = np.empty(N_BOOT)
        for i in range(N_BOOT):
            take = np.concatenate([rng.choice(idx, len(idx), replace=True) for idx in indices.values()])
            draws[i] = estimate(data.iloc[take], outcome)[3]
        low, high = np.quantile(draws, [0.025, 0.975])
        check(f"{effect_type}_independent_bootstrap_crosses_zero", low <= 0 <= high, f"CI=[{low:.8g},{high:.8g}]")
        validation_intervals.append({"effect_type": effect_type, "independent_ci_low": low, "independent_ci_high": high, "replicates": N_BOOT, "seed": VALIDATION_SEED})

    check("reported_all_cohort_indirect_cis_cross_zero", bool(((reported_cohort.bootstrap_ci_low <= 0) & (reported_cohort.bootstrap_ci_high >= 0)).all()))
    check("reported_all_pooled_indirect_cis_cross_zero", bool(((reported_pool.bootstrap_ci_low <= 0) & (reported_pool.bootstrap_ci_high >= 0)).all()))
    check("status_no_relation", status["reproducible_tymp_mtap_association"] is False)
    check("status_no_short_bridge", status["short_tymp_mtap_mta_support"] is False)
    check("status_no_serial_spp1", status["serial_support"]["SERIAL_TO_SPP1"] is False)
    check("status_no_serial_program", status["serial_support"]["SERIAL_TO_PROGRAM"] is False)
    check("status_decision", status["decision"] == "NO_SERIAL_MODEL_SUPPORT")

    pd.DataFrame(checks).to_csv(out / "validation_checks_phase99_v1.tsv", sep="\t", index=False)
    pd.DataFrame(validation_intervals).to_csv(out / "independent_bootstrap_intervals_v1.tsv", sep="\t", index=False)
    hashes = [{"file": p.name, "sha256": sha256(p)} for p in sorted(analysis.iterdir()) if p.is_file()]
    pd.DataFrame(hashes).to_csv(out / "analysis_output_hashes_v1.tsv", sep="\t", index=False)
    passed = sum(x["passed"] for x in checks)
    result = {"status": "PASS" if passed == len(checks) else "FAIL", "checks": len(checks), "passed": passed, "failed": len(checks) - passed, "validation_seed": VALIDATION_SEED, "bootstrap_replicates": N_BOOT}
    (out / "validation_status_phase99_v1.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
