#!/usr/bin/env python3
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
SEED = 2621099
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


def coefficient(y, columns, target_index):
    design = np.column_stack([np.ones(len(y))] + columns)
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(beta[target_index + 1])


def partial_rank_corr(x, y, covariates):
    xz, yz = zrank(x), zrank(y)
    design = np.column_stack([np.ones(len(xz))] + [zrank(v) for v in covariates])
    xres = xz - design @ np.linalg.lstsq(design, xz, rcond=None)[0]
    yres = yz - design @ np.linalg.lstsq(design, yz, rcond=None)[0]
    rho, p = stats.pearsonr(xres, yres)
    return float(rho), float(p)


def bootstrap_p(draws):
    n = len(draws)
    lo = (np.count_nonzero(draws <= 0) + 1) / (n + 1)
    hi = (np.count_nonzero(draws >= 0) + 1) / (n + 1)
    return min(1.0, 2.0 * min(lo, hi))


def path_effects(a, d, b, c=None):
    a1 = coefficient(d, [a], 0)
    a2 = coefficient(b, [d, a], 0)
    if c is None:
        return a1, a2, np.nan, a1 * a2
    b1 = coefficient(c, [b, d, a], 0)
    return a1, a2, b1, a1 * a2 * b1


def bootstrap_path(a, d, b, c, rng):
    a, d, b = zrank(a), zrank(d), zrank(b)
    c = None if c is None else zrank(c)
    point = path_effects(a, d, b, c)
    draws = np.empty(N_BOOT, dtype=float)
    for i in range(N_BOOT):
        idx = rng.integers(0, len(a), size=len(a))
        ci = None if c is None else c[idx]
        draws[i] = path_effects(a[idx], d[idx], b[idx], ci)[3]
    low, high = np.quantile(draws, [0.025, 0.975])
    return point, float(low), float(high), bootstrap_p(draws)


def meta_correlation(rows):
    rho = np.clip(np.array([r["rho"] for r in rows]), -0.999999, 0.999999)
    n = np.array([r["n"] for r in rows], dtype=float)
    z = np.arctanh(rho)
    var = 1.0 / (n - 3.0)
    w = 1.0 / var
    fixed_z = float(np.sum(w * z) / np.sum(w))
    fixed_se = math.sqrt(1.0 / float(np.sum(w)))
    q_stat = float(np.sum(w * (z - fixed_z) ** 2))
    q_p = float(stats.chi2.sf(q_stat, len(rows) - 1))
    c_term = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q_stat - (len(rows) - 1)) / c_term) if c_term > 0 else 0.0
    wr = 1.0 / (var + tau2)
    random_z = float(np.sum(wr * z) / np.sum(wr))
    random_se = math.sqrt(1.0 / float(np.sum(wr)))
    return {
        "k": len(rows), "total_n": int(n.sum()),
        "fixed_rho": math.tanh(fixed_z),
        "fixed_ci_low": math.tanh(fixed_z - 1.96 * fixed_se),
        "fixed_ci_high": math.tanh(fixed_z + 1.96 * fixed_se),
        "fixed_p": float(2 * stats.norm.sf(abs(fixed_z / fixed_se))),
        "q_stat": q_stat, "q_p": q_p, "tau2": tau2,
        "random_rho": math.tanh(random_z),
        "random_ci_low": math.tanh(random_z - 1.96 * random_se),
        "random_ci_high": math.tanh(random_z + 1.96 * random_se),
        "random_p": float(2 * stats.norm.sf(abs(random_z / random_se))),
        "directions": ";".join("positive" if x > 0 else "negative" if x < 0 else "zero" for x in rho),
    }


def cohort_dummies(labels):
    levels = list(DATASETS)
    return [np.asarray(labels == level, dtype=float) for level in levels[1:]]


def pooled_effect(frames, outcome, rng):
    blocks = []
    for dataset, frame in frames.items():
        block = pd.DataFrame({
            "dataset": dataset,
            "A": zrank(frame["TYMP"]), "D": zrank(frame["MTAP"]),
            "B": zrank(frame["MTA"]), "C": zrank(frame[outcome]) if outcome else np.nan,
        })
        blocks.append(block)
    data = pd.concat(blocks, ignore_index=True)

    def estimate(sample):
        a, d, b = sample["A"].to_numpy(), sample["D"].to_numpy(), sample["B"].to_numpy()
        fixed = cohort_dummies(sample["dataset"].to_numpy())
        a1 = coefficient(d, [a] + fixed, 0)
        a2 = coefficient(b, [d, a] + fixed, 0)
        if not outcome:
            return a1, a2, np.nan, a1 * a2
        c = sample["C"].to_numpy()
        b1 = coefficient(c, [b, d, a] + fixed, 0)
        return a1, a2, b1, a1 * a2 * b1

    point = estimate(data)
    draws = np.empty(N_BOOT, dtype=float)
    indices = {ds: np.flatnonzero(data["dataset"].to_numpy() == ds) for ds in DATASETS}
    for i in range(N_BOOT):
        take = np.concatenate([rng.choice(idx, size=len(idx), replace=True) for idx in indices.values()])
        draws[i] = estimate(data.iloc[take])[3]
    low, high = np.quantile(draws, [0.025, 0.975])
    return point, float(low), float(high), bootstrap_p(draws), len(data)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_phase99_tymp_mtap_serial_v1.py PHASE98_ROOT")
    root = Path(sys.argv[1]).resolve()
    inputs = root / "downloads_v1"
    out = root / "analysis_phase99_v1"
    out.mkdir(parents=True, exist_ok=False)
    mapping = pd.read_csv(root / "MasterMapping_ccRCC.csv", dtype=str).fillna("")
    rng = np.random.default_rng(SEED)
    associations, indirect_rows, counts = [], [], []
    tumor_frames = {}

    for dataset, (met_file, rna_file) in DATASETS.items():
        met = pd.read_csv(inputs / met_file, index_col=0)
        rna = pd.read_csv(inputs / rna_file, index_col=0)
        mapd = mapping[mapping["Dataset"] == dataset].copy()
        if set(mapd["MetabID"]) != set(met.index.astype(str)) or set(mapd["RNAID"]) != set(rna.index.astype(str)):
            raise RuntimeError(f"{dataset}: identifier mapping mismatch")
        met = met.loc[mapd["MetabID"]].copy()
        rna = rna.loc[mapd["RNAID"]].copy()
        ids = mapd["CommonID"].where(mapd["CommonID"] != "", mapd["MetabID"]).tolist()
        met.index = ids
        rna.index = ids
        needed = ["TYMP", "MTAP", "SPP1"] + PROGRAM_GENES
        if any(g not in rna.columns for g in needed):
            raise RuntimeError(f"{dataset}: required gene missing")
        frame = pd.DataFrame(index=ids)
        frame["TN"] = mapd["TN"].to_numpy()
        frame["MTA"] = pd.to_numeric(met[MTA], errors="coerce")
        for gene in needed:
            frame[gene] = np.log2(pd.to_numeric(rna[gene], errors="coerce") + 1.0)
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
        counts.append({"dataset": dataset, "total": len(frame), "tumor": int((frame.TN == "Tumor").sum()), "normal": int((frame.TN == "Normal").sum())})
        tumor = frame[frame.TN == "Tumor"].copy()
        pm = tumor[PROGRAM_GENES]
        tumor["SPP1_PROGRAM"] = ((pm - pm.mean()) / pm.std(ddof=0)).mean(axis=1)
        tumor_frames[dataset] = tumor

        tests = [
            ("TYMP_MTAP_MARGINAL", tumor, []),
            ("TYMP_MTAP_PARTIAL_SPP1", tumor, ["SPP1"]),
            ("TYMP_MTAP_PARTIAL_PROGRAM", tumor, ["SPP1_PROGRAM"]),
        ]
        all_frame = frame.copy()
        all_frame["TN_BINARY"] = (all_frame.TN == "Tumor").astype(float)
        tests.append(("TYMP_MTAP_ALL_PARTIAL_TN", all_frame, ["TN_BINARY"]))
        for endpoint, values, covs in tests:
            use = values[["TYMP", "MTAP"] + covs].dropna()
            if covs:
                rho, p = partial_rank_corr(use.TYMP, use.MTAP, [use[c] for c in covs])
            else:
                rho, p = stats.spearmanr(use.TYMP, use.MTAP)
            associations.append({"dataset": dataset, "endpoint": endpoint, "n": len(use), "rho": float(rho), "p": float(p), "covariates": ";".join(covs)})

        for effect_type, outcome in [("SHORT_TYMP_MTAP_MTA", None), ("SERIAL_TO_SPP1", "SPP1"), ("SERIAL_TO_PROGRAM", "SPP1_PROGRAM")]:
            usecols = ["TYMP", "MTAP", "MTA"] + ([] if outcome is None else [outcome])
            use = tumor[usecols].dropna()
            point, low, high, p_boot = bootstrap_path(use.TYMP, use.MTAP, use.MTA, None if outcome is None else use[outcome], rng)
            indirect_rows.append({
                "scope": "cohort", "dataset": dataset, "effect_type": effect_type, "outcome": outcome or "MTA",
                "n": len(use), "a1_TYMP_to_MTAP": point[0], "a2_MTAP_to_MTA_adjusted_TYMP": point[1],
                "b_MTA_to_C_adjusted_TYMP_MTAP": point[2], "indirect_effect": point[3],
                "bootstrap_ci_low": low, "bootstrap_ci_high": high, "bootstrap_p_two_sided": p_boot,
                "bootstrap_replicates": N_BOOT, "seed": SEED,
            })

    assoc = pd.DataFrame(associations)
    assoc["q_within_endpoint"] = np.nan
    for endpoint, idx in assoc.groupby("endpoint").groups.items():
        assoc.loc[idx, "q_within_endpoint"] = bh_adjust(assoc.loc[idx, "p"].to_numpy())
    assoc.to_csv(out / "cohort_tymp_mtap_associations_v1.tsv", sep="\t", index=False)
    pd.DataFrame(counts).to_csv(out / "sample_counts_v1.tsv", sep="\t", index=False)

    meta_rows = []
    for endpoint, group in assoc.groupby("endpoint", sort=False):
        meta_rows.append({"endpoint": endpoint, **meta_correlation(group[["n", "rho"]].to_dict("records"))})
    meta = pd.DataFrame(meta_rows)
    meta["fixed_q_across_endpoints"] = bh_adjust(meta.fixed_p.to_numpy())
    meta.to_csv(out / "meta_tymp_mtap_associations_v1.tsv", sep="\t", index=False)

    pooled_rows = []
    for effect_type, outcome in [("SHORT_TYMP_MTAP_MTA", None), ("SERIAL_TO_SPP1", "SPP1"), ("SERIAL_TO_PROGRAM", "SPP1_PROGRAM")]:
        point, low, high, p_boot, n = pooled_effect(tumor_frames, outcome, rng)
        pooled_rows.append({
            "scope": "stratified_pooled", "dataset": "ALL", "effect_type": effect_type, "outcome": outcome or "MTA",
            "n": n, "a1_TYMP_to_MTAP": point[0], "a2_MTAP_to_MTA_adjusted_TYMP": point[1],
            "b_MTA_to_C_adjusted_TYMP_MTAP": point[2], "indirect_effect": point[3],
            "bootstrap_ci_low": low, "bootstrap_ci_high": high, "bootstrap_p_two_sided": p_boot,
            "bootstrap_replicates": N_BOOT, "seed": SEED,
        })
    pooled = pd.DataFrame(pooled_rows)
    serial_idx = pooled.effect_type.str.startswith("SERIAL")
    pooled["q_serial_outcomes"] = np.nan
    pooled.loc[serial_idx, "q_serial_outcomes"] = bh_adjust(pooled.loc[serial_idx, "bootstrap_p_two_sided"].to_numpy())
    pd.DataFrame(indirect_rows).to_csv(out / "cohort_indirect_effects_v1.tsv", sep="\t", index=False)
    pooled.to_csv(out / "pooled_indirect_effects_v1.tsv", sep="\t", index=False)

    primary = assoc[assoc.endpoint == "TYMP_MTAP_MARGINAL"].copy()
    primary_meta = meta[meta.endpoint == "TYMP_MTAP_MARGINAL"].iloc[0]
    signs = np.sign(primary.rho.to_numpy())
    same_count = max(np.count_nonzero(signs > 0), np.count_nonzero(signs < 0))
    relation_support = bool(same_count >= 3 and primary_meta.fixed_p < 0.05 and primary_meta.random_p < 0.05 and np.count_nonzero(primary.q_within_endpoint < 0.05) >= 2)
    cohort_indirect = pd.DataFrame(indirect_rows)
    short_pool = pooled[pooled.effect_type == "SHORT_TYMP_MTAP_MTA"].iloc[0]
    short_cohort = cohort_indirect[cohort_indirect.effect_type == "SHORT_TYMP_MTAP_MTA"]
    short_sig = short_cohort[(short_cohort.bootstrap_ci_low > 0) | (short_cohort.bootstrap_ci_high < 0)]
    pool_short_nonzero = bool(short_pool.bootstrap_ci_low > 0 or short_pool.bootstrap_ci_high < 0)
    short_direction_count = 0 if short_sig.empty else max(np.count_nonzero(short_sig.indirect_effect > 0), np.count_nonzero(short_sig.indirect_effect < 0))
    short_support = bool(pool_short_nonzero and short_direction_count >= 2)
    serial_support = {}
    for effect_type in ["SERIAL_TO_SPP1", "SERIAL_TO_PROGRAM"]:
        prow = pooled[pooled.effect_type == effect_type].iloc[0]
        crow = cohort_indirect[cohort_indirect.effect_type == effect_type]
        sig = crow[(crow.bootstrap_ci_low > 0) | (crow.bootstrap_ci_high < 0)]
        direction_count = 0 if sig.empty else max(np.count_nonzero(sig.indirect_effect > 0), np.count_nonzero(sig.indirect_effect < 0))
        serial_support[effect_type] = bool(short_support and prow.q_serial_outcomes < 0.05 and direction_count >= 2)

    status = {
        "status": "PASS_COMPUTATION_VALIDATION_PENDING", "analysis_version": "phase99_v1_preregistered",
        "datasets": list(DATASETS), "tumor_samples": int(primary.n.sum()),
        "primary_tymp_mtap_same_direction_count": int(same_count),
        "primary_tymp_mtap_fixed_rho": float(primary_meta.fixed_rho),
        "primary_tymp_mtap_fixed_p": float(primary_meta.fixed_p),
        "primary_tymp_mtap_random_rho": float(primary_meta.random_rho),
        "primary_tymp_mtap_random_p": float(primary_meta.random_p),
        "primary_tymp_mtap_cohort_q_lt_0_05": int(np.count_nonzero(primary.q_within_endpoint < 0.05)),
        "reproducible_tymp_mtap_association": relation_support,
        "short_tymp_mtap_mta_support": short_support,
        "serial_support": serial_support,
        "decision": "SERIAL_MODEL_SUPPORTED" if any(serial_support.values()) else "NO_SERIAL_MODEL_SUPPORT",
        "seed": SEED, "bootstrap_replicates": N_BOOT,
    }
    (out / "analysis_status_phase99_v1.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
