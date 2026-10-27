#!/usr/bin/env python3
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MTA = "5-methylthioadenosine (MTA)"
PROGRAM_GENES = ["GPNMB", "LGALS3", "FABP5", "CD9", "LPL"]
DATASETS = {
    "RC18": ("matched_Harmonized_RawData_RC18.csv", "matched_tpm_RC18.csv"),
    "RC20": ("matched_Harmonized_RawData_RC20.csv", "matched_tpm_RC20.csv"),
    "CPTAC": ("matched_Harmonized_Met_CPTAC.csv", "matched_tpm_CPTAC.csv"),
    "CPTAC_val": ("matched_Harmonized_Met_CPTAC_val.csv", "matched_tpm_CPTAC_val.csv"),
}
SEED = 2621098
N_BOOT = 10000


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def zrank(values):
    ranks = stats.rankdata(np.asarray(values, dtype=float), method="average")
    sd = ranks.std(ddof=0)
    return (ranks - ranks.mean()) / sd if sd > 0 else np.zeros_like(ranks)


def partial_rank_corr(x, y, covariates):
    x = zrank(x)
    y = zrank(y)
    if covariates:
        design = np.column_stack([np.ones(len(x))] + [zrank(value) for value in covariates])
        x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    rho, p = stats.pearsonr(x, y)
    return float(rho), float(p)


def coefficient(y, columns, target_index):
    design = np.column_stack([np.ones(len(y))] + columns)
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(beta[target_index + 1])


def bootstrap_indirect(a, b, mtap, c, rng):
    a = zrank(a)
    b = zrank(b)
    mtap = zrank(mtap)
    c = zrank(c)
    point_a = coefficient(b, [a, mtap], 0)
    point_b = coefficient(c, [b, a, mtap], 0)
    draws = np.empty(N_BOOT, dtype=float)
    n = len(a)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        ai, bi, mi, ci = a[idx], b[idx], mtap[idx], c[idx]
        draws[i] = coefficient(bi, [ai, mi], 0) * coefficient(ci, [bi, ai, mi], 0)
    low, high = np.quantile(draws, [0.025, 0.975])
    p_boot = 2 * min(float(np.mean(draws <= 0)), float(np.mean(draws >= 0)))
    return point_a, point_b, point_a * point_b, float(low), float(high), min(p_boot, 1.0)


def meta_correlation(rows):
    rho = np.clip(np.array([row["rho"] for row in rows], dtype=float), -0.999999, 0.999999)
    n = np.array([row["n"] for row in rows], dtype=float)
    z = np.arctanh(rho)
    var = 1.0 / (n - 3.0)
    w = 1.0 / var
    fixed_z = float(np.sum(w * z) / np.sum(w))
    fixed_se = math.sqrt(1.0 / float(np.sum(w)))
    fixed_p = float(2 * stats.norm.sf(abs(fixed_z / fixed_se)))
    q_stat = float(np.sum(w * (z - fixed_z) ** 2))
    q_p = float(stats.chi2.sf(q_stat, len(rows) - 1))
    c_term = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q_stat - (len(rows) - 1)) / c_term) if c_term > 0 else 0.0
    wr = 1.0 / (var + tau2)
    random_z = float(np.sum(wr * z) / np.sum(wr))
    random_se = math.sqrt(1.0 / float(np.sum(wr)))
    random_p = float(2 * stats.norm.sf(abs(random_z / random_se)))
    return {
        "k": len(rows),
        "total_n": int(np.sum(n)),
        "fixed_rho": math.tanh(fixed_z),
        "fixed_ci_low": math.tanh(fixed_z - 1.96 * fixed_se),
        "fixed_ci_high": math.tanh(fixed_z + 1.96 * fixed_se),
        "fixed_p": fixed_p,
        "q_stat": q_stat,
        "q_p": q_p,
        "tau2": tau2,
        "random_rho": math.tanh(random_z),
        "random_ci_low": math.tanh(random_z - 1.96 * random_se),
        "random_ci_high": math.tanh(random_z + 1.96 * random_se),
        "random_p": random_p,
        "directions": ";".join("positive" if value > 0 else "negative" if value < 0 else "zero" for value in rho),
    }


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_phase98_paired_mta_analysis_v1.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    inputs = root / "downloads_v1"
    out = root / "analysis_v1"
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    associations = []
    mediation = []

    for dataset, (met_file, rna_file) in DATASETS.items():
        met = pd.read_csv(inputs / met_file, index_col=0)
        rna = pd.read_csv(inputs / rna_file, index_col=0)
        if list(met.index.astype(str)) != list(rna.index.astype(str)):
            raise RuntimeError(f"{dataset}: ordered sample IDs do not match")
        required = ["TYMP", "MTAP", "SPP1"] + PROGRAM_GENES
        missing = [gene for gene in required if gene not in rna.columns]
        if missing:
            raise RuntimeError(f"{dataset}: missing RNA genes {missing}")
        frame = pd.DataFrame(index=met.index)
        frame["MTA"] = pd.to_numeric(met[MTA], errors="coerce")
        for gene in required:
            frame[gene] = np.log2(pd.to_numeric(rna[gene], errors="coerce") + 1.0)
        program_matrix = frame[PROGRAM_GENES]
        frame["SPP1_PROGRAM"] = ((program_matrix - program_matrix.mean()) / program_matrix.std(ddof=0)).mean(axis=1)
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna()

        tests = [
            ("A_B_TYMP_MTA_MARGINAL", "primary_ab_marginal", "TYMP", "MTA", []),
            ("A_B_TYMP_MTA_PARTIAL_MTAP", "primary_ab_partial", "TYMP", "MTA", ["MTAP"]),
            ("B_C_MTA_SPP1_MARGINAL", "exploratory_bc_spp1", "MTA", "SPP1", []),
            ("B_C_MTA_SPP1_PARTIAL_TYMP_MTAP", "exploratory_bc_spp1_partial", "MTA", "SPP1", ["TYMP", "MTAP"]),
            ("B_C_MTA_SPP1_PROGRAM_MARGINAL", "exploratory_bc_program", "MTA", "SPP1_PROGRAM", []),
            ("B_C_MTA_SPP1_PROGRAM_PARTIAL_TYMP_MTAP", "exploratory_bc_program_partial", "MTA", "SPP1_PROGRAM", ["TYMP", "MTAP"]),
            ("A_C_TYMP_SPP1_DESCRIPTIVE", "descriptive_ac", "TYMP", "SPP1", []),
            ("A_C_TYMP_SPP1_PROGRAM_DESCRIPTIVE", "descriptive_ac", "TYMP", "SPP1_PROGRAM", []),
            ("MTAP_B_MTAP_MTA_CONTEXT", "context_mtap_mta", "MTAP", "MTA", []),
        ]
        for endpoint, family, x_name, y_name, cov_names in tests:
            values = frame[[x_name, y_name] + cov_names].dropna()
            if cov_names:
                rho, p = partial_rank_corr(values[x_name], values[y_name], [values[name] for name in cov_names])
            else:
                rho, p = stats.spearmanr(values[x_name], values[y_name])
            associations.append({
                "dataset": dataset,
                "endpoint": endpoint,
                "family": family,
                "x": x_name,
                "y": y_name,
                "covariates": ";".join(cov_names),
                "n": len(values),
                "rho": float(rho),
                "p": float(p),
            })

        for outcome in ["SPP1", "SPP1_PROGRAM"]:
            values = frame[["TYMP", "MTA", "MTAP", outcome]].dropna()
            pa, pb, indirect, low, high, p_boot = bootstrap_indirect(
                values["TYMP"].to_numpy(), values["MTA"].to_numpy(),
                values["MTAP"].to_numpy(), values[outcome].to_numpy(), rng
            )
            mediation.append({
                "dataset": dataset,
                "outcome": outcome,
                "n": len(values),
                "a_TYMP_to_MTA_adjusted_MTAP": pa,
                "b_MTA_to_C_adjusted_TYMP_MTAP": pb,
                "indirect_a_times_b": indirect,
                "bootstrap_ci_low": low,
                "bootstrap_ci_high": high,
                "bootstrap_p_two_sided": p_boot,
                "bootstrap_replicates": N_BOOT,
                "seed": SEED,
            })

    assoc_df = pd.DataFrame(associations)
    assoc_df["q_within_family"] = np.nan
    for family, idx in assoc_df.groupby("family").groups.items():
        assoc_df.loc[idx, "q_within_family"] = bh_adjust(assoc_df.loc[idx, "p"].to_numpy())
    assoc_df.to_csv(out / "cohort_associations_v1.tsv", sep="\t", index=False)
    pd.DataFrame(mediation).to_csv(out / "mediation_bootstrap_v1.tsv", sep="\t", index=False)

    meta_rows = []
    for endpoint, group in assoc_df.groupby("endpoint", sort=False):
        values = group[["dataset", "n", "rho"]].to_dict("records")
        meta_rows.append({"endpoint": endpoint, **meta_correlation(values)})
    meta_df = pd.DataFrame(meta_rows)
    meta_df["fixed_q_across_endpoints"] = bh_adjust(meta_df["fixed_p"].to_numpy())
    meta_df.to_csv(out / "cross_cohort_meta_analysis_v1.tsv", sep="\t", index=False)

    primary_marginal = assoc_df[assoc_df.endpoint == "A_B_TYMP_MTA_MARGINAL"]
    primary_partial = assoc_df[assoc_df.endpoint == "A_B_TYMP_MTA_PARTIAL_MTAP"]
    meta_marginal = meta_df[meta_df.endpoint == "A_B_TYMP_MTA_MARGINAL"].iloc[0]
    same_direction = primary_marginal.rho.gt(0).all() or primary_marginal.rho.lt(0).all()
    strong = same_direction and meta_marginal.fixed_p < 0.05 and int((primary_marginal.q_within_family < 0.05).sum()) >= 2
    partial_meta = meta_df[meta_df.endpoint == "A_B_TYMP_MTA_PARTIAL_MTAP"].iloc[0]
    if strong:
        decision = "STRONG_OBSERVATIONAL_A_TO_B_SUPPORT"
    elif meta_marginal.fixed_p < 0.05 or partial_meta.fixed_p < 0.05:
        decision = "PARTIAL_OBSERVATIONAL_A_TO_B_SUPPORT"
    else:
        decision = "NO_OBSERVATIONAL_A_TO_B_SUPPORT"
    status = {
        "status": "PASS_COMPUTATION_VALIDATION_PENDING",
        "decision": decision,
        "datasets": list(DATASETS),
        "total_samples": int(primary_marginal.n.sum()),
        "primary_same_direction": bool(same_direction),
        "primary_cohort_q_lt_0_05": int((primary_marginal.q_within_family < 0.05).sum()),
        "primary_fixed_meta_rho": float(meta_marginal.fixed_rho),
        "primary_fixed_meta_p": float(meta_marginal.fixed_p),
        "primary_partial_fixed_meta_rho": float(partial_meta.fixed_rho),
        "primary_partial_fixed_meta_p": float(partial_meta.fixed_p),
        "scope": "Observational paired bulk RNA and measured MTA. B-to-C and mediation are compatibility analyses, not macrophage fine-state or causal proof.",
    }
    (out / "analysis_status_v1.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    summary = [
        "# Phase 98 paired ccRCC RNA–MTA analysis", "",
        f"- Status: `{status['status']}`.",
        f"- Decision: `{decision}`.",
        f"- Four cohorts, total paired samples: {status['total_samples']}.",
        f"- TYMP–MTA fixed-effect meta rho={status['primary_fixed_meta_rho']:.4f}, P={status['primary_fixed_meta_p']:.4g}.",
        f"- MTAP-adjusted TYMP–MTA meta rho={status['primary_partial_fixed_meta_rho']:.4f}, P={status['primary_partial_fixed_meta_p']:.4g}.",
        "- MTA–SPP1, anchor-excluded bulk program and mediation-compatible results are exploratory; bulk data do not identify an SPP1 macrophage fine state.",
        "- No causal statement is allowed without perturbation, direct MTA measurement and phenocopy/depletion-rescue.",
    ]
    (out / "PHASE98_PAIRED_MTA_SUMMARY_v1.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(json.dumps(status))


if __name__ == "__main__":
    main()
